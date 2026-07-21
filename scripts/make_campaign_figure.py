"""Regenerate the campaign coverage figure (paper.tex Fig. 1) from the
published per-rotation numbers in Table 2 of paper.tex / Sect. 4.2 of
docs/manuscript.md, without requiring the underlying FITS-derived
Carrington maps (which need SOAR/JSOC downloads; see
scripts/plot_campaign_summary.py for the version driven directly off a
real `cr_<N>/milestone/` pipeline run).

Usage:
    python scripts/make_campaign_figure.py --out figures/campaign_polar_advantage.png
"""

import argparse
from pathlib import Path

import numpy as np

# Okabe-Ito, colour-vision-safe; matches plot_campaign_summary.py
PHI_C, HMI_C = "#0072B2", "#E69F00"
INK, MUTED, GRID = "#1a1a1a", "#5a5a5a", "#d8d8d8"

# Table 2 (paper.tex): rotation, date range, SolO B0 range (deg),
# separation range (deg), N-cap PHI/HMI (%), S-cap PHI/HMI (%).
# CR2298-2299 B0/separation ranges are the paper's own approximate figures.
ROTATIONS = [
    dict(cr=2294, dates="Feb 14 - Mar 1",  b0=(-2, -8),   sep=(15, 21),   n=(5, 0),   s=(30, 38)),
    dict(cr=2295, dates="Mar 2 - Mar 28",  b0=(-8, -17),  sep=(0.3, 60),  n=(0, 0),   s=(64, 46)),
    dict(cr=2296, dates="Mar 31 - Apr 24", b0=(-8, 17),   sep=(80, 165),  n=(51, 3),  s=(12, 41)),
    dict(cr=2297, dates="Apr 26 - May 21", b0=(14, 17),   sep=(168, 180), n=(77, 12), s=(0, 34)),
    dict(cr=2298, dates="May 23 - Jun 19", b0=(14, 6),    sep=(155, 179), n=(64, 23), s=(0, 24)),
    dict(cr=2299, dates="Jun 20 - Jul 16", b0=(6, 0),     sep=(90, 160),  n=(45, 34), s=(4, 13)),
]
# Earth's B0 is not tabulated per rotation in the paper -- only the
# campaign-level trend ("-7 deg toward -2 deg, northern summer"). Shown as
# a linear reference between those two stated endpoints, not a per-rotation
# measurement.
EARTH_B0_START, EARTH_B0_END = -7.0, -2.0


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, default=Path("figures/campaign_polar_advantage.png"))
    return p.parse_args()


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    args = parse_args()
    n = len(ROTATIONS)
    x = np.arange(n)

    solo_lo = np.array([min(r["b0"]) for r in ROTATIONS])
    solo_hi = np.array([max(r["b0"]) for r in ROTATIONS])
    earth_b0 = np.linspace(EARTH_B0_START, EARTH_B0_END, n)
    phi_N = np.array([r["n"][0] for r in ROTATIONS]); hmi_N = np.array([r["n"][1] for r in ROTATIONS])
    phi_S = np.array([r["s"][0] for r in ROTATIONS]); hmi_S = np.array([r["s"][1] for r in ROTATIONS])
    rot_lbl = [f"CR {r['cr']}\n{r['dates']}" for r in ROTATIONS]
    sep_lbl = [f"sep {r['sep'][0]:g}-{r['sep'][1]:g}°" for r in ROTATIONS]

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10.5,
        "axes.edgecolor": MUTED, "axes.linewidth": 0.8,
        "axes.titlesize": 11.5, "axes.titleweight": "bold",
    })
    fig, (ax0, ax1, ax2) = plt.subplots(
        3, 1, figsize=(8.4, 10.2), sharex=True,
        gridspec_kw={"height_ratios": [1.15, 1, 1], "hspace": 0.16})

    # panel A: the orbital driver
    b0lim = max(20.0, float(np.abs(np.r_[solo_lo, solo_hi, earth_b0]).max()) + 4)
    ax0.axhline(0, color=GRID, lw=1.2, zorder=0)
    for i in x:
        ax0.plot([i, i], [solo_lo[i], solo_hi[i]], color=PHI_C, lw=7,
                 solid_capstyle="round", alpha=0.85, zorder=2,
                 label="SolO/PHI B0 range" if i == 0 else None)
        ax0.annotate(f"{solo_hi[i]:+.0f}°", (i, solo_hi[i]), textcoords="offset points",
                     xytext=(9, -1), color=PHI_C, fontsize=9.5, fontweight="bold", va="center")
        ax0.annotate(f"{solo_lo[i]:+.0f}°", (i, solo_lo[i]), textcoords="offset points",
                     xytext=(9, 1), color=PHI_C, fontsize=9.5, va="center")
    ax0.plot(x, earth_b0, "o--", color=HMI_C, lw=2, ms=8, zorder=3,
             label="Earth/HMI B0 (linear ref., campaign endpoints)")
    ax0.set_ylabel("B0  (deg)"); ax0.set_ylim(-b0lim, b0lim)
    ax0.set_title("A.  The orbital driver: SolO's heliolatitude vs Earth's")
    ax0.legend(loc="lower left", frameon=False, fontsize=8.8)
    ax0.grid(axis="y", color=GRID, lw=0.6); ax0.set_axisbelow(True)

    # panels B/C: polar-cap fill
    w = 0.34
    ytop = max(72.0, float(np.r_[phi_N, hmi_N, phi_S, hmi_S].max()) + 14)

    def fill_panel(ax, phi, hmi, title):
        b1 = ax.bar(x - w / 2, phi, w, color=PHI_C, label="PHI (SolO)", zorder=3)
        b2 = ax.bar(x + w / 2, hmi, w, color=HMI_C, label="HMI (SDO)", zorder=3)
        for b in list(b1) + list(b2):
            h = b.get_height()
            ax.annotate(f"{h:.0f}%", (b.get_x() + b.get_width() / 2, h),
                        textcoords="offset points", xytext=(0, 3), ha="center",
                        fontsize=9.5, color=INK, fontweight="bold" if h >= 45 else "normal")
        ax.set_ylabel("cap fill  (%)"); ax.set_title(title); ax.set_ylim(0, ytop)
        ax.grid(axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True)

    fill_panel(ax1, phi_N, hmi_N, "B.  North polar cap (≥ 60°):  PHI vs HMI")
    ax1.legend(loc="upper left", frameon=False, fontsize=9.5, ncol=2)
    fill_panel(ax2, phi_S, hmi_S, "C.  South polar cap (≥ 60°):  PHI vs HMI")

    ax2.set_xticks(x)
    ax2.set_xticklabels([f"{r}\n{s}" for r, s in zip(rot_lbl, sep_lbl)], fontsize=9.2)
    ax2.tick_params(axis="x", pad=6)

    fig.suptitle("Solar Orbiter/PHI polar coverage vs SDO/HMI", fontsize=13, fontweight="bold", y=0.985)
    fig.text(0.5, 0.01,
              "Coverage is calibration-independent. Where the SolO-Earth separation is large, a per-pixel "
              "PHI+HMI merge is\nnot meaningful (opposite hemispheres) — there PHI is a standalone polar constraint.",
              ha="center", va="bottom", fontsize=8.6, color=MUTED)
    fig.subplots_adjust(top=0.93, bottom=0.145, left=0.1, right=0.965)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=160)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
