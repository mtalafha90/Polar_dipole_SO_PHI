"""Resolution-degradation test for the PHI-vs-HMI calibration slope.

The per-case calibration slope PHI ~= slope * HMI sits well below 1
(~0.55 in the 2025 campaign) and declines with SolO-Sun distance. One
candidate cause is resolution: PHI's plate scale is coarser than HMI's, so
mixed-polarity fine structure cancels within a PHI pixel and PHI reads low.
This script tests that directly and parameter-free: for a range of Gaussian
smoothing widths it degrades the co-observed HMI map (reprojected onto the
PHI grid) BEFORE the regression and reports how the slope responds. If the
slope climbs toward 1 as HMI is smoothed to PHI-like resolution, the deficit
is resolution, not instrument scale or vantage.

Outputs a per-(case, FWHM) table and an aggregate slope-vs-FWHM curve.

Usage:
    python scripts/calibration_resolution_test.py --dates 20250214-20250301 \
        --fwhms 0,1,2,3,5,8
"""

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
from sunpy.util.exceptions import SunpyMetadataWarning

warnings.filterwarnings("ignore", category=SunpyMetadataWarning)
sys.path.append(str(Path(__file__).resolve().parents[1]))

from baseline_config import (
    PHI_DIR, HMI_DIR, OUT_DIR, MAX_TIME_DIFF_SEC,
    R_INNER, R_OUTER, DISK_FRACTION, MU_MIN, ALPHA,
)
from solar_pipeline.io_utils import build_hmi_time_index, expand_date_spec
from solar_pipeline.pipeline import compute_case_fields
from solar_pipeline.calibration import calibration_stats


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dates", type=str, required=True)
    p.add_argument("--phi-dir", type=Path, default=PHI_DIR)
    p.add_argument("--hmi-dir", type=Path, default=HMI_DIR)
    p.add_argument("--out-dir", type=Path, default=OUT_DIR / "milestone")
    p.add_argument("--fwhms", type=str, default="0,1,2,3,5,8",
                   help="Comma list of Gaussian FWHM (PHI-grid pixels) to smooth HMI by")
    p.add_argument("--max-time-diff-sec", type=float, default=MAX_TIME_DIFF_SEC)
    p.add_argument("--mu-min", type=float, default=MU_MIN)
    p.add_argument("--alpha", type=float, default=ALPHA)
    p.add_argument("--calib-min-abs-g", type=float, default=10.0)
    return p.parse_args()


def nan_gaussian(data, fwhm_pix):
    """Edge-aware Gaussian smooth that ignores (and preserves) NaNs."""
    if fwhm_pix <= 0:
        return data
    sigma = fwhm_pix / 2.3548
    mask = np.isfinite(data).astype(float)
    filled = np.where(np.isfinite(data), data, 0.0)
    num = gaussian_filter(filled, sigma, mode="nearest")
    den = gaussian_filter(mask, sigma, mode="nearest")
    out = np.divide(num, den, out=np.full_like(num, np.nan), where=den > 0)
    out[mask == 0] = np.nan
    return out


def main():
    args = parse_args()
    fwhms = [float(x) for x in args.fwhms.split(",")]
    only_dates = expand_date_spec(args.dates)

    phi_files = sorted(args.phi_dir.glob("solo_L2_phi-fdt-blos_*.fits"))
    phi_files = [f for f in phi_files if only_dates is None or any(d in f.name for d in only_dates)]
    hmi_files = sorted(args.hmi_dir.glob("hmi.M_720s.*.magnetogram.fits"))
    if not phi_files or not hmi_files:
        raise SystemExit("No PHI or HMI files found for the requested dates.")
    hmi_index = build_hmi_time_index(hmi_files)

    rows = []
    for i, phi_path in enumerate(phi_files, start=1):
        try:
            fields = compute_case_fields(
                phi_path, hmi_index, max_time_diff_sec=args.max_time_diff_sec,
                r_inner=R_INNER, r_outer=R_OUTER, disk_fraction=DISK_FRACTION,
                mu_min=args.mu_min, alpha=args.alpha,
            )
        except Exception as exc:
            print(f"[{i}/{len(phi_files)}] {phi_path.name}: skip ({exc})")
            continue
        phi = fields["phi_blos"].data
        hmi = fields["hmi_on_phi"].data
        mu = fields["mu"]
        dsun = float(fields["phi_blos"].meta.get("dsun_obs", np.nan)) / 1.495978707e11
        for fwhm in fwhms:
            calib = calibration_stats(phi, nan_gaussian(hmi, fwhm), mu,
                                      mu_min=args.mu_min, min_abs_ref=args.calib_min_abs_g)
            rows.append({"phi_file": phi_path.name, "dsun_au": dsun, "fwhm_pix": fwhm,
                         "slope": calib["slope"], "pearson_r": calib["pearson_r"],
                         "n_pixels": calib["n_pixels"]})
        s0 = rows[-len(fwhms)]["slope"]
        print(f"[{i}/{len(phi_files)}] {phi_path.name}: slope(fwhm=0)={s0:.3f}")

    if not rows:
        raise SystemExit("No cases processed.")
    df = pd.DataFrame(rows)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out_dir / "calibration_resolution_test.csv", index=False)

    agg = df.groupby("fwhm_pix").agg(
        slope_mean=("slope", "mean"), slope_std=("slope", "std"),
        r_mean=("pearson_r", "mean"), n_cases=("slope", "count")).reset_index()
    print("\n=== slope vs HMI smoothing FWHM (mean over cases) ===")
    print(agg.to_string(index=False))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 4.6))
    ax.errorbar(agg["fwhm_pix"], agg["slope_mean"], yerr=agg["slope_std"],
                marker="o", color="#0072B2", lw=2, capsize=3, label="mean slope +/- std")
    ax.axhline(1.0, color="#7a7a7a", ls="--", lw=1, label="slope = 1 (no flux deficit)")
    ax.set_xlabel("HMI smoothing FWHM [PHI-grid pixels]")
    ax.set_ylabel("PHI-vs-HMI calibration slope")
    ax.set_title("Does degrading HMI to PHI resolution close the slope deficit?")
    ax.grid(alpha=0.3); ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(args.out_dir / "calibration_resolution_test.png", dpi=150)
    plt.close(fig)
    print(f"\nOutputs in: {args.out_dir.resolve()}")


if __name__ == "__main__":
    main()
