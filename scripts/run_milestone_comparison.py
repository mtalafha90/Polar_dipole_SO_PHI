"""First-milestone comparison (see README section 6b):

Build three Carrington-style synoptic Br maps over the selected subset —
PHI-only (SolO vantage, native geometry), HMI-only (Earth vantage, native
geometry, NOT reprojected through the PHI grid), and merged (PHI+HMI blend
on the PHI grid) — then compute and compare the axial dipole moment g10
from each, under each polar-filling assumption.
"""

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
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
from solar_pipeline.io_utils import build_hmi_time_index, save_synoptic_fits
from solar_pipeline.pipeline import compute_case_fields, compute_native_disk_fields
from solar_pipeline.calibration import calibration_stats
from solar_pipeline.carrington import (
    cm_weight,
    bin_br_to_carrington_weighted,
    bin_max_to_carrington,
    combine_weighted_grids,
    combine_max_grids,
    carrington_fill_fraction,
)
from solar_pipeline.dipole import FILL_MODES, axial_dipole_g10, polar_fill_fractions
from solar_pipeline.plotting import plot_carrington_map


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
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
    parser.add_argument("--cm-weight-power", type=float, default=1.0)
    parser.add_argument(
        "--polar-lat-deg",
        type=float,
        default=60.0,
        help="Latitude bound for polar-cap fill-fraction diagnostics",
    )
    parser.add_argument(
        "--calibrate-phi",
        action="store_true",
        help="Divide PHI Br by the per-case PHI-vs-HMI regression slope before binning",
    )
    parser.add_argument(
        "--calib-min-abs-g",
        type=float,
        default=10.0,
        help="|B| threshold (Gauss) for pixels entering the calibration regression",
    )
    return parser.parse_args()


class Accumulator:
    """Accumulates weighted synoptic bins and max-mu confidence per product."""

    def __init__(self, nlat: int, nlon: int):
        self.nlat, self.nlon = nlat, nlon
        self.wsums, self.weights, self.quals = [], [], []
        self.lat_centers = self.lon_centers = None

    def add(self, br, lat, lon, valid, weight, mu):
        wsum, wt, self.lat_centers, self.lon_centers = bin_br_to_carrington_weighted(
            br, lat, lon, valid, weight, nlat=self.nlat, nlon=self.nlon
        )
        self.wsums.append(wsum)
        self.weights.append(wt)
        self.quals.append(bin_max_to_carrington(mu, lat, lon, valid, nlat=self.nlat, nlon=self.nlon))

    def combine(self):
        grid, total_weight = combine_weighted_grids(self.wsums, self.weights)
        quality = combine_max_grids(self.quals)
        return grid, total_weight, quality


def main():
    args = parse_args()
    only_dates = {d.strip() for d in args.dates.split(",") if d.strip()}

    out_dir = args.out_dir / "milestone"
    out_dir.mkdir(exist_ok=True, parents=True)
    plots_dir = out_dir / "plots"
    plots_dir.mkdir(exist_ok=True, parents=True)

    phi_blos_files = sorted(args.phi_dir.glob("solo_L2_phi-fdt-blos_*.fits"))
    phi_blos_files = [f for f in phi_blos_files if any(d in f.name for d in only_dates)]
    hmi_files = sorted(args.hmi_dir.glob("hmi.M_720s.*.magnetogram.fits"))

    if not phi_blos_files:
        raise RuntimeError("No PHI blos files found.")
    if not hmi_files:
        raise RuntimeError("No HMI files found.")

    hmi_index = build_hmi_time_index(hmi_files)

    acc = {name: Accumulator(args.nlat, args.nlon) for name in ("phi", "hmi", "merged")}
    calib_rows = []
    seen_hmi = set()
    n_ok = 0

    for i, phi_blos_path in enumerate(phi_blos_files, start=1):
        print(f"[{i}/{len(phi_blos_files)}] {phi_blos_path.name}")
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

            calib = calibration_stats(
                fields["phi_blos"].data, fields["hmi_on_phi"].data, fields["mu"],
                mu_min=args.mu_min, min_abs_ref=args.calib_min_abs_g,
            )
            calib_rows.append({"phi_blos_file": phi_blos_path.name, **calib})
            print(
                f"  calib: PHI ~= {calib['slope']:.3f} x HMI "
                f"(r={calib['pearson_r']:.3f}, n={calib['n_pixels']})"
            )

            phi_scale = 1.0
            if args.calibrate_phi and np.isfinite(calib["slope"]) and calib["slope"] != 0.0:
                phi_scale = 1.0 / calib["slope"]

            w_phi = cm_weight(fields["cmd"], power=args.cm_weight_power)

            # PHI-only, SolO vantage (PHI is already in its native geometry)
            acc["phi"].add(
                fields["br_phi"] * phi_scale, fields["lat"], fields["lon"],
                fields["valid_phi"], w_phi, fields["mu"],
            )
            # merged product lives on the PHI grid
            br_merged = fields["br_merged"] * phi_scale if args.calibrate_phi else fields["br_merged"]
            acc["merged"].add(
                br_merged, fields["lat"], fields["lon"],
                fields["valid_merged"], w_phi, fields["mu"],
            )

            # HMI-only, Earth vantage: bin the matched HMI magnetogram from
            # its OWN disk geometry (deduplicated if two PHI files matched
            # the same HMI record)
            if fields["hmi_path"] not in seen_hmi:
                seen_hmi.add(fields["hmi_path"])
                native = compute_native_disk_fields(
                    fields["hmi"],
                    disk_fraction=args.disk_fraction,
                    mu_min=args.mu_min,
                    alpha=args.alpha,
                )
                acc["hmi"].add(
                    native["br"], native["lat"], native["lon"],
                    native["valid"], cm_weight(native["cmd"], power=args.cm_weight_power),
                    native["mu"],
                )

            n_ok += 1
        except Exception as exc:
            print(f"  ERROR: {exc}")

    if n_ok == 0:
        raise RuntimeError("No cases processed successfully.")

    rows = []
    for name in ("phi", "hmi", "merged"):
        grid, weight, quality = acc[name].combine()
        lat_c, lon_c = acc[name].lat_centers, acc[name].lon_centers

        np.save(out_dir / f"grid_{name}.npy", grid)
        np.save(out_dir / f"weight_{name}.npy", weight)
        np.save(out_dir / f"quality_{name}.npy", quality)
        save_synoptic_fits(
            out_dir / f"synoptic_{name}.fits", grid, lat_c, lon_c,
            extra={"product": name, "history": "SO-PHI-SDO-HMI milestone comparison"},
        )
        plot_carrington_map(
            grid, lat_c, lon_c, plots_dir / f"map_{name}.png",
            title=f"{name} synoptic map ({n_ok} cases)",
        )

        polar = polar_fill_fractions(weight, lat_c, polar_lat_deg=args.polar_lat_deg)
        row = {
            "product": name,
            "fill_fraction": carrington_fill_fraction(weight),
            f"polar_fill_north_{args.polar_lat_deg:.0f}": polar["north"],
            f"polar_fill_south_{args.polar_lat_deg:.0f}": polar["south"],
        }
        for mode in FILL_MODES:
            dip = axial_dipole_g10(grid, lat_c, mode=mode)
            row[f"g10_{mode}"] = dip["g10"]
            row[f"g10_north_{mode}"] = dip["g10_north"]
            row[f"g10_south_{mode}"] = dip["g10_south"]
        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "milestone_dipole_comparison.csv", index=False)
    calib_df = pd.DataFrame(calib_rows)
    calib_df.to_csv(out_dir / "calibration_stats.csv", index=False)

    np.save(out_dir / "lat_centers.npy", acc["phi"].lat_centers)
    np.save(out_dir / "lon_centers.npy", acc["phi"].lon_centers)

    print("\n=== Milestone dipole comparison ===")
    print(df.to_string(index=False))
    print("\n=== PHI-vs-HMI calibration ===")
    print(calib_df.to_string(index=False))
    print(f"\nOutputs in: {out_dir.resolve()}")


if __name__ == "__main__":
    main()
