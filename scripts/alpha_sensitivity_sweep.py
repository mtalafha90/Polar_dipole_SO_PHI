import argparse
import sys
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
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
    NLAT,
    NLON,
)
from solar_pipeline.io_utils import build_hmi_time_index, expand_date_spec
from solar_pipeline.pipeline import run_case


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Sweep ALPHA in Br ~= Blos / mu**alpha and report how the PHI/HMI/merged "
            "dipole estimates change, to help pick a new baseline default. Does not "
            "change baseline_config.py."
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
    parser.add_argument("--nlat", type=int, default=NLAT)
    parser.add_argument("--nlon", type=int, default=NLON)
    parser.add_argument(
        "--alphas",
        type=str,
        default="0.6,0.7,0.8,0.9,1.0",
        help="Comma-separated ALPHA values to sweep",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    only_dates = expand_date_spec(args.dates)
    alphas = [float(a.strip()) for a in args.alphas.split(",") if a.strip()]

    out_dir = args.out_dir / "alpha_sweep"
    out_dir.mkdir(exist_ok=True, parents=True)

    phi_blos_files = sorted(args.phi_dir.glob("solo_L2_phi-fdt-blos_*.fits"))
    phi_blos_files = [f for f in phi_blos_files if only_dates is None or any(d in f.name for d in only_dates)]
    hmi_files = sorted(args.hmi_dir.glob("hmi.M_720s.*.magnetogram.fits"))

    if not phi_blos_files:
        raise RuntimeError("No PHI blos files found.")
    if not hmi_files:
        raise RuntimeError("No HMI files found.")

    hmi_index = build_hmi_time_index(hmi_files)

    rows = []
    for alpha in alphas:
        print(f"\n=== ALPHA = {alpha} ===")
        for i, phi_blos_path in enumerate(phi_blos_files, start=1):
            print(f"[{i}/{len(phi_blos_files)}] {phi_blos_path.name}")
            try:
                row, _ = run_case(
                    phi_blos_path,
                    hmi_index,
                    max_time_diff_sec=args.max_time_diff_sec,
                    r_inner=args.r_inner,
                    r_outer=args.r_outer,
                    disk_fraction=args.disk_fraction,
                    mu_min=args.mu_min,
                    alpha=alpha,
                    nlat=args.nlat,
                    nlon=args.nlon,
                )
                row["alpha"] = alpha
                rows.append(row)
            except Exception as exc:
                print(f"  ERROR: {exc}")
                rows.append({"phi_blos_file": phi_blos_path.name, "alpha": alpha, "error": str(exc)})

    df = pd.DataFrame(rows)
    all_csv = out_dir / "alpha_sweep_all_cases.csv"
    summary_csv = out_dir / "alpha_sweep_summary.csv"
    df.to_csv(all_csv, index=False)

    good = df[df.get("error").isna()] if "error" in df.columns else df

    summary = good.groupby("alpha")[["dip_phi", "dip_hmi", "dip_merged"]].agg(["mean", "std"])
    summary.columns = ["_".join(c) for c in summary.columns]
    summary = summary.reset_index()
    summary.to_csv(summary_csv, index=False)

    print("\nAlpha sensitivity summary")
    print(summary.to_string(index=False))
    print(f"\nSaved: {all_csv.resolve()}")
    print(f"Saved: {summary_csv.resolve()}")

    plots_dir = out_dir / "plots"
    plots_dir.mkdir(exist_ok=True, parents=True)
    plt.figure(figsize=(7, 5))
    plt.errorbar(
        summary["alpha"], summary["dip_merged_mean"], yerr=summary["dip_merged_std"],
        marker="o", capsize=3, label="Merged",
    )
    plt.errorbar(
        summary["alpha"], summary["dip_phi_mean"], yerr=summary["dip_phi_std"],
        marker="o", capsize=3, label="PHI-only",
    )
    plt.errorbar(
        summary["alpha"], summary["dip_hmi_mean"], yerr=summary["dip_hmi_std"],
        marker="o", capsize=3, label="HMI-on-PHI",
    )
    plt.xlabel("alpha (Br = Blos / mu^alpha)")
    plt.ylabel("Carrington-style dipole (mean ± std across cases)")
    plt.title("Alpha sensitivity")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(plots_dir / "alpha_sensitivity.png", dpi=150)
    plt.close()
    print(f"Saved: {(plots_dir / 'alpha_sensitivity.png').resolve()}")


if __name__ == "__main__":
    main()
