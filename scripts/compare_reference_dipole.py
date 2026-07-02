"""Check the milestone map dipoles against a standard HMI synoptic chart.

Computes g10 from a reference full-CR HMI synoptic radial-field map (CEA /
uniform sine-latitude grid, e.g. hmi.synoptic_mr_polfil_720s) and compares
it with the PHI-only / HMI-only / merged milestone products under every
polar-fill assumption. The reference can be a local FITS file
(--reference) or fetched from JSOC by Carrington rotation (--car-rot)
using the same no-registration rs_list/SUMS path as the data downloader.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

sys.path.append(str(Path(__file__).resolve().parents[1]))

from baseline_config import OUT_DIR
from solar_pipeline.dipole import FILL_MODES, axial_dipole_g10, axial_dipole_g10_sinlat

JSOC_BASE = "http://jsoc.stanford.edu"


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, default=None, help="Local reference synoptic FITS")
    parser.add_argument("--car-rot", type=int, default=None, help="Fetch this Carrington rotation from JSOC")
    parser.add_argument(
        "--series",
        type=str,
        default="hmi.synoptic_mr_polfil_720s",
        help="JSOC synoptic series for --car-rot (polar-filled radial field by default)",
    )
    parser.add_argument("--maps-dir", type=Path, default=OUT_DIR / "milestone")
    parser.add_argument("--reference-dir", type=Path, default=OUT_DIR / "reference")
    return parser.parse_args()


def fetch_reference(series: str, car_rot: int, dest_dir: Path) -> Path:
    try:
        import drms
        import requests
    except ImportError as exc:
        raise SystemExit(
            f"Missing download dependency ({exc}). Install with: pip install -e '.[download]'"
        )

    dest_dir.mkdir(exist_ok=True, parents=True)
    dest = dest_dir / f"{series}.{car_rot}.fits"
    if dest.exists():
        print(f"already present: {dest}")
        return dest

    client = drms.Client()
    seg_names = list(client.info(series).segments.index)
    if not seg_names:
        raise SystemExit(f"Series {series} has no segments")
    seg = seg_names[0]

    ds = f"{series}[{car_rot}]"
    print(f"fetching {ds} (segment {seg}) ...")
    _, segs = client.query(ds, key="CAR_ROT", seg=seg)
    if not len(segs):
        raise SystemExit(f"No record for {ds}")
    seg_path = str(segs[seg].iloc[0])
    if not seg_path.startswith("/"):
        raise SystemExit(f"Record {ds} is offline in SUMS; export it manually from JSOC.")

    tmp = dest.with_suffix(".part")
    with requests.get(JSOC_BASE + seg_path, stream=True, timeout=300) as resp:
        resp.raise_for_status()
        with open(tmp, "wb") as f:
            for chunk in resp.iter_content(1 << 20):
                f.write(chunk)
    tmp.rename(dest)
    print(f"saved: {dest}")
    return dest


def load_reference_grid(path: Path) -> np.ndarray:
    with fits.open(path) as hdul:
        for hdu in hdul:
            data = getattr(hdu, "data", None)
            if data is not None and getattr(data, "ndim", 0) == 2:
                return np.asarray(data, dtype=float)
    raise SystemExit(f"No 2D image HDU in {path}")


def main():
    args = parse_args()

    if args.reference is None and args.car_rot is None:
        raise SystemExit("Provide --reference <fits> or --car-rot <N>")

    ref_path = args.reference or fetch_reference(args.series, args.car_rot, args.reference_dir)
    ref_grid = load_reference_grid(ref_path)
    g10_ref = axial_dipole_g10_sinlat(ref_grid)
    finite_frac = float(np.isfinite(ref_grid).mean())
    print(f"\nReference: {ref_path.name}  shape={ref_grid.shape}  finite={finite_frac:.3f}")
    print(f"Reference g10 = {g10_ref:+.4f} G")

    lat_centers = np.load(args.maps_dir / "lat_centers.npy")
    rows = [{"product": "reference", "mode": "sinlat", "g10": g10_ref, "delta_vs_reference": 0.0}]
    for product in ("phi", "hmi", "merged"):
        grid_path = args.maps_dir / f"grid_{product}.npy"
        if not grid_path.exists():
            continue
        grid = np.load(grid_path)
        for mode in FILL_MODES:
            g10 = axial_dipole_g10(grid, lat_centers, mode=mode)["g10"]
            rows.append({
                "product": product,
                "mode": mode,
                "g10": g10,
                "delta_vs_reference": g10 - g10_ref,
            })

    df = pd.DataFrame(rows)
    out_csv = args.maps_dir / "reference_dipole_comparison.csv"
    df.to_csv(out_csv, index=False)
    print("\n=== Dipole vs reference synoptic chart ===")
    print(df.to_string(index=False))
    print(f"\nSaved: {out_csv.resolve()}")


if __name__ == "__main__":
    main()
