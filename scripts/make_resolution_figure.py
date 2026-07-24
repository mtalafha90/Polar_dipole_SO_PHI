"""Resolution-degradation figure (paper.tex Fig. `resolution`).

Aggregates the per-case calibration_resolution_test.py output
(`calibration_stats` re-run with the HMI map smoothed to a range of
Gaussian widths) into the mean PHI-vs-HMI slope and Pearson r versus
FWHM, showing that the sub-unity calibration slope closes to unity as
HMI is degraded to PHI's plate scale -- i.e. the deficit is resolution,
not vantage.

Usage:
    python scripts/make_resolution_figure.py \
        --csv data/calibration_resolution_test.csv \
        --out figures/calibration_resolution.png
"""

import argparse
import csv
import statistics as st
from collections import defaultdict
from pathlib import Path

# Okabe-Ito, colour-vision-safe; consistent with the other paper figures
PHI_C, HMI_C, GRID, MUTED = "#0072B2", "#E69F00", "#d8d8d8", "#5a5a5a"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--csv", type=Path, default=Path("data/calibration_resolution_test.csv"),
                   help="calibration_resolution_test.py output (per case, per FWHM)")
    p.add_argument("--out", type=Path, default=Path("figures/calibration_resolution.png"))
    return p.parse_args()


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    args = parse_args()
    rows = list(csv.DictReader(open(args.csv)))
    by_slope, by_r = defaultdict(list), defaultdict(list)
    for r in rows:
        by_slope[float(r["fwhm_pix"])].append(float(r["slope"]))
        by_r[float(r["fwhm_pix"])].append(float(r["pearson_r"]))
    fw = sorted(by_slope)
    n_cases = len(by_slope[fw[0]])
    ms = [st.mean(by_slope[f]) for f in fw]
    ss = [st.pstdev(by_slope[f]) for f in fw]
    mr = [st.mean(by_r[f]) for f in fw]

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 11,
        "axes.edgecolor": MUTED, "axes.linewidth": 0.8,
        "axes.titlesize": 12, "axes.titleweight": "bold",
    })
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.axhline(1.0, color=GRID, lw=1.4, zorder=0)
    ax.errorbar(fw, ms, yerr=ss, color=PHI_C, lw=2, marker="o", ms=7, capsize=3,
                zorder=3, label="mean slope (PHI/HMI)")
    ax.set_xlabel("HMI Gaussian smoothing FWHM  [PHI-grid px]")
    ax.set_ylabel("calibration slope", color=PHI_C)
    ax.tick_params(axis="y", labelcolor=PHI_C)
    ax.set_ylim(0.5, 1.15)
    ax.grid(color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.annotate("slope crosses unity\nnear 5 px (PHI plate scale)", xy=(5, ms[fw.index(5.0)] if 5.0 in fw else 0.98),
                xytext=(5.4, 0.72), fontsize=9.5, color=MUTED,
                arrowprops=dict(arrowstyle="->", color=MUTED))

    ax2 = ax.twinx()
    ax2.plot(fw, mr, color=HMI_C, lw=2, marker="s", ms=6, ls="--", zorder=2,
             label="mean Pearson $r$")
    ax2.set_ylabel("Pearson $r$", color=HMI_C)
    ax2.tick_params(axis="y", labelcolor=HMI_C)
    ax2.set_ylim(0.80, 0.95)

    l1, la1 = ax.get_legend_handles_labels()
    l2, la2 = ax2.get_legend_handles_labels()
    ax.legend(l1 + l2, la1 + la2, loc="lower right", frameon=False, fontsize=9.5)
    ax.set_title(f"Calibration slope deficit is a resolution effect (CR 2294, {n_cases} cases)")
    fig.tight_layout()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=160)
    print(f"Wrote {args.out}  (FWHM {fw}, slope {[round(x,3) for x in ms]})")


if __name__ == "__main__":
    main()
