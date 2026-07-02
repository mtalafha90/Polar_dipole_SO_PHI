import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
from sunpy.util.exceptions import SunpyMetadataWarning

warnings.filterwarnings("ignore", category=SunpyMetadataWarning)
sys.path.append(str(Path(__file__).resolve().parents[1]))

from baseline_config import (
    PHI_DIR,
    HMI_DIR,
    OUT_DIR,
    ONLY_DATES,
    MAX_TIME_DIFF_SEC,
    R_INNER,
    R_OUTER,
    DISK_FRACTION,
    MU_MIN,
    ALPHA,
    NLAT,
    NLON,
)
from solar_pipeline.io_utils import build_hmi_time_index, expand_date_spec
from solar_pipeline.pipeline import compute_case_fields
from solar_pipeline.carrington import (
    cm_weight,
    bin_br_to_carrington_weighted,
    combine_weighted_grids,
    axial_dipole_from_carrington_grid,
    carrington_fill_fraction,
)
from solar_pipeline.plotting import plot_carrington_map


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Assimilate several visible-disk PHI/HMI cases into one Carrington-style "
            "synoptic map, weighting each case's contribution by closeness to central "
            "meridian at observation time."
        )
    )
    parser.add_argument("--phi-dir", type=Path, default=PHI_DIR)
    parser.add_argument("--hmi-dir", type=Path, default=HMI_DIR)
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    parser.add_argument("--dates", type=str, default=",".join(sorted(ONLY_DATES)))
    parser.add_argument("--max-time-diff-sec", type=float, default=MAX_TIME_DIFF_SEC)
    parser.add_argument("--r-inner", type=float, default=R_INNER)
    parser.add_argument("--r-outer", type=float, default=R_OUTER)
    parser.add_argument("--disk-fraction", type=float, default=DISK_FRACTION)
    parser.add_argument("--mu-min", type=float, default=MU_MIN)
    parser.add_argument("--alpha", type=float, default=ALPHA)
    parser.add_argument("--nlat", type=int, default=NLAT)
    parser.add_argument("--nlon", type=int, default=NLON)
    parser.add_argument(
        "--cm-weight-power",
        type=float,
        default=1.0,
        help="Exponent n in cos(cmd)**n used to weight each case's contribution "
        "by central-meridian distance (higher = sharper preference for near-CM data)",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    only_dates = expand_date_spec(args.dates)

    out_dir = args.out_dir / "synoptic"
    out_dir.mkdir(exist_ok=True, parents=True)

    phi_blos_files = sorted(args.phi_dir.glob("solo_L2_phi-fdt-blos_*.fits"))
    phi_blos_files = [f for f in phi_blos_files if only_dates is None or any(d in f.name for d in only_dates)]
    hmi_files = sorted(args.hmi_dir.glob("hmi.M_720s.*.magnetogram.fits"))

    if not phi_blos_files:
        raise RuntimeError("No PHI blos files found.")
    if not hmi_files:
        raise RuntimeError("No HMI files found.")

    print(f"Found {len(phi_blos_files)} PHI blos files")
    print(f"Found {len(hmi_files)} HMI magnetogram files")

    hmi_index = build_hmi_time_index(hmi_files)

    wsum_list = []
    weight_list = []
    lat_centers = lon_centers = None
    n_ok = 0

    for i, phi_blos_path in enumerate(phi_blos_files, start=1):
        print(f"[{i}/{len(phi_blos_files)}] Processing {phi_blos_path.name}")
        try:
            fields = compute_case_fields(
                phi_blos_path,
                hmi_index,
                max_time_diff_sec=args.max_time_diff_sec,
                r_inner=args.r_inner,
                r_outer=args.r_outer,
                disk_fraction=args.disk_fraction,
                mu_min=args.mu_min,
                alpha=args.alpha,
            )
            weight = cm_weight(fields["cmd"], power=args.cm_weight_power)
            wsum, wt, lat_centers, lon_centers = bin_br_to_carrington_weighted(
                fields["br_merged"],
                fields["lat"],
                fields["lon"],
                fields["valid_merged"],
                weight,
                nlat=args.nlat,
                nlon=args.nlon,
            )
            wsum_list.append(wsum)
            weight_list.append(wt)
            n_ok += 1
        except Exception as exc:
            print(f"  ERROR: {exc}")

    if n_ok == 0:
        raise RuntimeError("No cases were successfully assimilated.")

    grid, total_weight = combine_weighted_grids(wsum_list, weight_list)
    fill = carrington_fill_fraction(total_weight)
    dipole = axial_dipole_from_carrington_grid(grid, lat_centers)

    np.save(out_dir / "synoptic_grid.npy", grid)
    np.save(out_dir / "synoptic_weight.npy", total_weight)
    np.save(out_dir / "lat_centers.npy", lat_centers)
    np.save(out_dir / "lon_centers.npy", lon_centers)

    with open(out_dir / "synoptic_summary.txt", "w") as f:
        f.write("Synoptic assimilation summary\n")
        f.write(f"cases assimilated: {n_ok}/{len(phi_blos_files)}\n")
        f.write(f"cm_weight_power  : {args.cm_weight_power}\n")
        f.write(f"fill_fraction    : {fill:.6f}\n")
        f.write(f"axial_dipole     : {dipole:.6f}\n")

    plots_dir = out_dir / "plots"
    plots_dir.mkdir(exist_ok=True, parents=True)
    plot_carrington_map(
        grid,
        lat_centers,
        lon_centers,
        plots_dir / "synoptic_map.png",
        title=f"CM-weighted synoptic map ({n_ok} cases, dipole={dipole:.3f})",
    )

    print("\nDone.")
    print(f"Cases assimilated: {n_ok}/{len(phi_blos_files)}")
    print(f"Fill fraction    : {fill:.6f}")
    print(f"Axial dipole     : {dipole:.6f}")
    print(f"Saved outputs in : {out_dir.resolve()}")


if __name__ == "__main__":
    main()
