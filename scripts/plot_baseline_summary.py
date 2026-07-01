from pathlib import Path
import pandas as pd

from baseline_config import OUT_DIR, MU_MIN, ALPHA
from solar_pipeline.plotting import make_baseline_plots


def main():
    summary_csv = OUT_DIR / "baseline_summary.csv"
    if not summary_csv.exists():
        raise FileNotFoundError(f"Missing summary CSV: {summary_csv}")

    df = pd.read_csv(summary_csv)
    if "error" in df.columns:
        df = df[df["error"].isna()]
    make_baseline_plots(df, OUT_DIR / "plots", title_suffix=f" (MU_MIN={MU_MIN}, alpha={ALPHA})")
    print(f"Plots saved in: {(OUT_DIR / 'plots').resolve()}")


if __name__ == "__main__":
    main()