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
from solar_pipeline.io_utils import build_hmi_time_index, expand_date_spec, save_fits
from solar_pipeline.pipeline import run_case, summarize_dataframe
from solar_pipeline.plotting import make_baseline_plots


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run the PHI-HMI baseline pipeline. Defaults reproduce baseline v1 (see README)."
    )
    parser.add_argument("--phi-dir", type=Path, default=PHI_DIR, help="Directory containing PHI blos FITS files")
    parser.add_argument("--hmi-dir", type=Path, default=HMI_DIR, help="Directory containing HMI magnetogram FITS files")
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR, help="Directory to write pipeline outputs")
    parser.add_argument(
        "--dates",
        type=str,
        default=",".join(sorted(ONLY_DATES)),
        help="Comma-separated list of YYYYMMDD dates to include",
    )
    parser.add_argument("--max-time-diff-sec", type=float, default=MAX_TIME_DIFF_SEC)
    parser.add_argument("--r-inner", type=float, default=R_INNER)
    parser.add_argument("--r-outer", type=float, default=R_OUTER)
    parser.add_argument("--disk-fraction", type=float, default=DISK_FRACTION)
    parser.add_argument("--mu-min", type=float, default=MU_MIN)
    parser.add_argument("--alpha", type=float, default=ALPHA)
    parser.add_argument("--nlat", type=int, default=NLAT)
    parser.add_argument("--nlon", type=int, default=NLON)
    return parser.parse_args()


def main():
    args = parse_args()
    only_dates = expand_date_spec(args.dates)

    args.out_dir.mkdir(exist_ok=True, parents=True)

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

    rows = []

    for i, phi_blos_path in enumerate(phi_blos_files, start=1):
        print(f"\n[{i}/{len(phi_blos_files)}] Processing {phi_blos_path.name}")
        try:
            row, arrays = run_case(
                phi_blos_path,
                hmi_index,
                max_time_diff_sec=args.max_time_diff_sec,
                r_inner=args.r_inner,
                r_outer=args.r_outer,
                disk_fraction=args.disk_fraction,
                mu_min=args.mu_min,
                alpha=args.alpha,
                nlat=args.nlat,
                nlon=args.nlon,
            )
            rows.append(row)

            print(
                f"  matched HMI: {row['hmi_file']}\n"
                f"  Δt = {row['time_diff_sec']:.1f} s\n"
                f"  dip_phi = {row['dip_phi']:.6f}, "
                f"dip_hmi = {row['dip_hmi']:.6f}, "
                f"dip_merged = {row['dip_merged']:.6f}"
            )

            case_dir = args.out_dir / Path(phi_blos_path).stem
            case_dir.mkdir(exist_ok=True, parents=True)

            save_fits(case_dir / "merged_smooth_los.fits", arrays["merged"].astype("float32"), arrays["phi_header"])
            np.save(case_dir / "grid_phi.npy", arrays["grid_phi"])
            np.save(case_dir / "grid_hmi.npy", arrays["grid_hmi"])
            np.save(case_dir / "grid_merged.npy", arrays["grid_merged"])
            np.save(case_dir / "lat_centers.npy", arrays["lat_centers"])
            np.save(case_dir / "lon_centers.npy", arrays["lon_centers"])

        except Exception as exc:
            print(f"  ERROR: {exc}")
            rows.append({"phi_blos_file": phi_blos_path.name, "error": str(exc)})

    df = pd.DataFrame(rows)
    full_csv = args.out_dir / "baseline_all_cases.csv"
    summary_csv = args.out_dir / "baseline_summary.csv"
    notes_txt = args.out_dir / "baseline_summary_notes.txt"

    df.to_csv(full_csv, index=False)

    good = df[df.get("error").isna()] if "error" in df.columns else df

    summary_cols = [
        "phi_blos_file",
        "phi_time",
        "hmi_file",
        "hmi_time",
        "time_diff_sec",
        "dip_phi",
        "dip_hmi",
        "dip_merged",
        "merged_minus_phi",
        "merged_minus_hmi",
        "phi_crln_obs",
        "phi_crlt_obs",
        "hmi_crln_obs",
        "hmi_crlt_obs",
        "fill_phi",
        "fill_hmi",
        "fill_merged",
    ]
    summary_cols = [c for c in summary_cols if c in good.columns]
    summary_df = good[summary_cols].copy()
    summary_df.to_csv(summary_csv, index=False)

    stats = summarize_dataframe(good)

    with open(notes_txt, "w") as f:
        f.write("Baseline pipeline summary\n")
        f.write(f"R_INNER={args.r_inner}, R_OUTER={args.r_outer}, MU_MIN={args.mu_min}, ALPHA={args.alpha}\n\n")
        for col, s in stats.items():
            f.write(
                f"{col}: mean={s['mean']:.6f}, std={s['std']:.6f}, "
                f"min={s['min']:.6f}, max={s['max']:.6f}\n"
            )

    print("\nDone.")
    print(f"Saved full table   : {full_csv.resolve()}")
    print(f"Saved summary table: {summary_csv.resolve()}")
    print(f"Saved notes        : {notes_txt.resolve()}")

    print("\nBaseline summary")
    print(f"R_INNER={args.r_inner}, R_OUTER={args.r_outer}, MU_MIN={args.mu_min}, ALPHA={args.alpha}")
    for col, s in stats.items():
        print(
            f"{col}: mean={s['mean']:.6f}, std={s['std']:.6f}, "
            f"min={s['min']:.6f}, max={s['max']:.6f}"
        )

    if len(good) > 0:
        make_baseline_plots(
            good, args.out_dir / "plots", title_suffix=f" (MU_MIN={args.mu_min}, alpha={args.alpha})"
        )


if __name__ == "__main__":
    main()
