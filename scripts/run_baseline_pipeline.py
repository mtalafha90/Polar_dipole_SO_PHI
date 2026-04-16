from pathlib import Path
import sys
from pathlib import Path
import warnings
from sunpy.util.exceptions import SunpyMetadataWarning

warnings.filterwarnings("ignore", category=SunpyMetadataWarning)
sys.path.append(str(Path(__file__).resolve().parents[1]))
import pandas as pd

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
from solar_pipeline.io_utils import build_hmi_time_index, save_fits
from solar_pipeline.pipeline import run_case, summarize_dataframe
from solar_pipeline.plotting import make_baseline_plots


def main():
    OUT_DIR.mkdir(exist_ok=True, parents=True)

    phi_blos_files = sorted(PHI_DIR.glob("solo_L2_phi-fdt-blos_*.fits"))
    phi_blos_files = [f for f in phi_blos_files if any(d in f.name for d in ONLY_DATES)]
    hmi_files = sorted(HMI_DIR.glob("hmi.M_720s.*.magnetogram.fits"))

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
                max_time_diff_sec=MAX_TIME_DIFF_SEC,
                r_inner=R_INNER,
                r_outer=R_OUTER,
                disk_fraction=DISK_FRACTION,
                mu_min=MU_MIN,
                alpha=ALPHA,
                nlat=NLAT,
                nlon=NLON,
            )
            rows.append(row)

            print(
                f"  matched HMI: {row['hmi_file']}\n"
                f"  Δt = {row['time_diff_sec']:.1f} s\n"
                f"  dip_phi = {row['dip_phi']:.6f}, "
                f"dip_hmi = {row['dip_hmi']:.6f}, "
                f"dip_merged = {row['dip_merged']:.6f}"
            )

            case_dir = OUT_DIR / Path(phi_blos_path).stem
            case_dir.mkdir(exist_ok=True, parents=True)

            save_fits(case_dir / "merged_smooth_los.fits", arrays["merged"].astype("float32"), arrays["phi_header"])
            import numpy as np
            np.save(case_dir / "grid_phi.npy", arrays["grid_phi"])
            np.save(case_dir / "grid_hmi.npy", arrays["grid_hmi"])
            np.save(case_dir / "grid_merged.npy", arrays["grid_merged"])
            np.save(case_dir / "lat_centers.npy", arrays["lat_centers"])
            np.save(case_dir / "lon_centers.npy", arrays["lon_centers"])

        except Exception as exc:
            print(f"  ERROR: {exc}")
            rows.append({"phi_blos_file": phi_blos_path.name, "error": str(exc)})

    df = pd.DataFrame(rows)
    full_csv = OUT_DIR / "baseline_all_cases.csv"
    summary_csv = OUT_DIR / "baseline_summary.csv"
    notes_txt = OUT_DIR / "baseline_summary_notes.txt"

    df.to_csv(full_csv, index=False)

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
    summary_cols = [c for c in summary_cols if c in df.columns]
    summary_df = df[summary_cols].copy()
    summary_df.to_csv(summary_csv, index=False)

    good = df[df.get("error").isna()] if "error" in df.columns else df
    stats = summarize_dataframe(good)

    with open(notes_txt, "w") as f:
        f.write("Baseline pipeline summary\n")
        f.write(f"R_INNER={R_INNER}, R_OUTER={R_OUTER}, MU_MIN={MU_MIN}, ALPHA={ALPHA}\n\n")
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
    print(f"R_INNER={R_INNER}, R_OUTER={R_OUTER}, MU_MIN={MU_MIN}, ALPHA={ALPHA}")
    for col, s in stats.items():
        print(
            f"{col}: mean={s['mean']:.6f}, std={s['std']:.6f}, "
            f"min={s['min']:.6f}, max={s['max']:.6f}"
        )

    if len(good) > 0:
        make_baseline_plots(good, OUT_DIR / "plots", title_suffix=f" (MU_MIN={MU_MIN}, alpha={ALPHA})")


if __name__ == "__main__":
    main()