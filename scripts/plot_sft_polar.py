"""Figure for the SFT polar-constraint experiment.

Reads the CSV/npy outputs of run_sft_polar_experiment.py and draws the
companion figure: the injected initial conditions (with the PHI-constrained
polar cap shaded) and the axial-dipole evolution g10(t), with reversal
markers and the headline Delta-first-reversal annotated. Reading the saved
outputs means the figure can be regenerated/restyled without re-running the
SFT.

Usage:
    python scripts/plot_sft_polar.py --sft-dir baseline_outputs/cr_2296/milestone/sft_polar
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Okabe-Ito, colour-vision-safe; consistent with run_sft_polar_experiment.py
COLORS = {"hmi": "#E69F00", "phi_constrained": "#0072B2", "phi": "#009E73"}
LABELS = {"hmi": "HMI only", "phi_constrained": "HMI + PHI cap", "phi": "PHI only"}
GRID = "#d8d8d8"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sft-dir", type=Path, required=True,
                   help="run_sft_polar_experiment.py output dir (holds sft_polar_comparison.csv)")
    p.add_argument("--out", type=Path, default=None, help="Default: <sft-dir>/sft_polar_figure.png")
    p.add_argument("--polar-lat-deg", type=float, default=60.0)
    p.add_argument("--title", type=str, default="SFT polar-constraint experiment")
    return p.parse_args()


def _reversals(rev_df, ic):
    row = rev_df[rev_df["ic"] == ic]
    if row.empty or not isinstance(row["reversal_times_yr"].iloc[0], str):
        return []
    s = row["reversal_times_yr"].iloc[0]
    return [float(x) for x in s.split(";") if x]


def main():
    args = parse_args()
    comp = pd.read_csv(args.sft_dir / "sft_polar_comparison.csv")
    rev = pd.read_csv(args.sft_dir / "sft_polar_reversals.csv")
    lat = np.load(args.sft_dir / "sft_latitude.npy")
    ics = [c for c in ("hmi", "phi_constrained", "phi") if (args.sft_dir / f"b_init_{c}.npy").exists()]

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10.5,
                         "axes.edgecolor": "#5a5a5a", "axes.linewidth": 0.8,
                         "axes.titlesize": 11.5, "axes.titleweight": "bold"})
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(12, 4.8))

    # panel A: initial conditions, polar caps shaded
    for c in ics:
        b = np.load(args.sft_dir / f"b_init_{c}.npy")
        ax0.plot(lat, b, color=COLORS[c], lw=2, label=LABELS[c])
    ax0.axvspan(args.polar_lat_deg, 90, color="0.92", zorder=0)
    ax0.axvspan(-90, -args.polar_lat_deg, color="0.92", zorder=0)
    ax0.set_xlabel("latitude [deg]"); ax0.set_ylabel("B init [G]")
    ax0.set_title("A.  Injected initial conditions")
    ax0.set_xlim(-90, 90); ax0.grid(color=GRID, lw=0.6); ax0.set_axisbelow(True)
    ax0.legend(frameon=False, fontsize=9.5)

    # panel B: g10(t) with reversal markers
    for c in ics:
        col = f"g10_{c}"
        if col in comp:
            ax1.plot(comp["time_yr"], comp[col], color=COLORS[c], lw=2, label=LABELS[c])
            for r in _reversals(rev, c):
                ax1.plot(r, 0.0, "o", color=COLORS[c], ms=7, zorder=5)
    ax1.axhline(0, color="#7a7a7a", lw=1)
    ax1.set_xlabel("time [yr]"); ax1.set_ylabel("g10 [G]")
    ax1.set_title("B.  Axial dipole evolution")
    ax1.grid(color=GRID, lw=0.6); ax1.set_axisbelow(True)
    ax1.legend(frameon=False, fontsize=9.5, loc="best")

    # headline: Delta first reversal (phi_constrained vs hmi)
    r_hmi, r_con = _reversals(rev, "hmi"), _reversals(rev, "phi_constrained")
    if r_hmi and r_con:
        d = r_con[0] - r_hmi[0]
        ax1.annotate(f"Δ first reversal = {d:+.2f} yr",
                     xy=(max(r_hmi[0], r_con[0]), 0.0), xytext=(0.5, 0.9),
                     textcoords="axes fraction", fontsize=10, fontweight="bold",
                     color="#1a1a1a", ha="center",
                     bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#bbbbbb"))

    fig.suptitle(args.title, fontsize=13, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = args.out or (args.sft_dir / "sft_polar_figure.png")
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
