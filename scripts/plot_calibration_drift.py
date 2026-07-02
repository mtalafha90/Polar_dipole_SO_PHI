"""Diagnose the PHI-vs-HMI calibration drift.

The CR 2264 run showed the through-origin slope declining monotonically
(0.75 -> 0.50 over eight days) with a distinct low cluster (~0.44-0.46,
r ~ 0.60) for the late-UT observations. This plots slope and r against
time (colored by hour of day to expose observing-program clusters), and —
when the columns are present — against SolO-Sun distance and the
SolO-Earth longitude separation, to identify what the drift tracks.
"""

import argparse
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

sys.path.append(str(Path(__file__).resolve().parents[1]))

from baseline_config import OUT_DIR


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--calib-csv",
        type=Path,
        default=OUT_DIR / "milestone" / "calibration_stats.csv",
        help="calibration_stats.csv from run_milestone_comparison.py",
    )
    parser.add_argument("--out", type=Path, default=None, help="Output figure path")
    return parser.parse_args()


def _time_from_filename(name: str):
    m = re.search(r"_(\d{8}T\d{6})_", name)
    return pd.to_datetime(m.group(1), format="%Y%m%dT%H%M%S", utc=True) if m else pd.NaT


def main():
    args = parse_args()
    df = pd.read_csv(args.calib_csv)

    if "phi_time" in df.columns:
        df["time"] = pd.to_datetime(df["phi_time"], utc=True, format="ISO8601")
    else:  # older CSVs: recover the time from the PHI filename
        df["time"] = df["phi_blos_file"].map(_time_from_filename)
    df = df.sort_values("time").reset_index(drop=True)
    df["hour"] = df["time"].dt.hour + df["time"].dt.minute / 60.0

    have_dsun = "phi_dsun_au" in df.columns and df["phi_dsun_au"].notna().any()
    have_sep = "lon_separation_deg" in df.columns and df["lon_separation_deg"].notna().any()

    ncols = 2 + int(have_dsun) + int(have_sep)
    fig, axes = plt.subplots(1, ncols, figsize=(4.8 * ncols, 4.4))
    axes = list(axes) if ncols > 1 else [axes]

    ax = axes[0]
    sc = ax.scatter(df["time"], df["slope"], c=df["hour"], cmap="twilight", s=45, zorder=3)
    ax.plot(df["time"], df["slope"], color="gray", alpha=0.4, zorder=2)
    fig.colorbar(sc, ax=ax, label="hour of day [UT]")
    ax.set_ylabel("slope (PHI / HMI)")
    ax.set_title("Calibration slope vs time")
    ax.tick_params(axis="x", rotation=30)
    ax.grid(alpha=0.3)

    ax = axes[1]
    ax.scatter(df["time"], df["pearson_r"], c=df["hour"], cmap="twilight", s=45, zorder=3)
    ax.plot(df["time"], df["pearson_r"], color="gray", alpha=0.4, zorder=2)
    ax.set_ylabel("Pearson r")
    ax.set_title("Correlation vs time")
    ax.tick_params(axis="x", rotation=30)
    ax.grid(alpha=0.3)

    k = 2
    if have_dsun:
        ax = axes[k]
        ax.scatter(df["phi_dsun_au"], df["slope"], c=df["hour"], cmap="twilight", s=45)
        ax.set_xlabel("SolO-Sun distance [AU]")
        ax.set_ylabel("slope")
        ax.set_title("Slope vs SolO distance")
        ax.grid(alpha=0.3)
        k += 1
    if have_sep:
        ax = axes[k]
        ax.scatter(df["lon_separation_deg"], df["slope"], c=df["hour"], cmap="twilight", s=45)
        ax.set_xlabel("SolO-Earth longitude separation [deg]")
        ax.set_ylabel("slope")
        ax.set_title("Slope vs vantage separation")
        ax.grid(alpha=0.3)

    fig.suptitle("PHI-vs-HMI calibration drift diagnostics")
    fig.tight_layout()

    out = args.out or (args.calib_csv.parent / "plots" / "calibration_drift.png")
    out.parent.mkdir(exist_ok=True, parents=True)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved: {out.resolve()}")

    # console summary: overall trend + per-hour clusters
    span = (df["time"].iloc[-1] - df["time"].iloc[0]).total_seconds() / 86400.0
    print(f"slope: {df['slope'].iloc[0]:.3f} -> {df['slope'].iloc[-1]:.3f} over {span:.1f} days")
    late = df[df["hour"] >= 20]
    if len(late) and len(late) < len(df):
        rest = df[df["hour"] < 20]
        print(
            f"late-UT (>=20h) cluster: slope {late['slope'].mean():.3f} +/- {late['slope'].std():.3f} "
            f"vs rest {rest['slope'].mean():.3f} +/- {rest['slope'].std():.3f}"
        )


if __name__ == "__main__":
    main()
