import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.append(str(Path(__file__).resolve().parents[1]))

from baseline_config import OUT_DIR, MU_MIN, ALPHA
from solar_pipeline.plotting import make_baseline_plots


def parse_args():
    parser = argparse.ArgumentParser(
        description="Regenerate baseline summary plots from an existing baseline_summary.csv."
    )
    parser.add_argument(
        "--out-dir", type=Path, default=OUT_DIR, help="Directory containing baseline_summary.csv and plots/"
    )
    parser.add_argument("--mu-min", type=float, default=MU_MIN, help="MU_MIN value shown in plot titles")
    parser.add_argument("--alpha", type=float, default=ALPHA, help="ALPHA value shown in plot titles")
    return parser.parse_args()


def main():
    args = parse_args()
    summary_csv = args.out_dir / "baseline_summary.csv"
    if not summary_csv.exists():
        raise FileNotFoundError(f"Missing summary CSV: {summary_csv}")

    df = pd.read_csv(summary_csv)
    if "error" in df.columns:
        df = df[df["error"].isna()]
    make_baseline_plots(df, args.out_dir / "plots", title_suffix=f" (MU_MIN={args.mu_min}, alpha={args.alpha})")
    print(f"Plots saved in: {(args.out_dir / 'plots').resolve()}")


if __name__ == "__main__":
    main()
