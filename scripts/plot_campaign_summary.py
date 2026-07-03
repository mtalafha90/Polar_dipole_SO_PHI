"""Summarise a per-rotation campaign in one figure.

Reads the ``cr_<N>/milestone/`` outputs produced by
``run_milestone_by_rotation.py`` and builds the headline figure of the
SolO/PHI polar campaign:

  A. the orbital driver  -- SolO's B0 excursion vs Earth's B0 per rotation
  B. north polar-cap fill -- PHI vs HMI
  C. south polar-cap fill -- PHI vs HMI

Coverage (panels B/C) is the robust, calibration-independent result; the
SolO-Earth separation range is annotated per rotation because it decides
whether a per-pixel PHI+HMI merge is even meaningful (see README 6c).

Usage:
    python scripts/plot_campaign_summary.py --campaign-dir baseline_outputs \
        --out baseline_outputs/campaign_polar_advantage.png
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.append(str(Path(__file__).resolve().parents[1]))

PHI_C, HMI_C = "#0072B2", "#E69F00"          # Okabe-Ito, colour-vision-safe
INK, MUTED, GRID = "#1a1a1a", "#5a5a5a", "#d8d8d8"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--campaign-dir", type=Path, required=True,
                   help="Directory containing cr_<N>/ subdirectories.")
    p.add_argument("--out", type=Path, default=None,
                   help="Output PNG (default: <campaign-dir>/campaign_polar_advantage.png).")
    p.add_argument("--polar-lat-deg", type=float, default=60.0,
                   help="Polar-cap latitude threshold used in the milestone run (default 60).")
    return p.parse_args()


def _month_day(ts: str) -> str:
    t = pd.Timestamp(ts)
    return t.strftime("%b %-d")


def read_rotation(cr_dir: Path, cap: float):
    """Extract polar fill, B0 range and separation range for one rotation."""
    mile = pd.read_csv(cr_dir / "milestone" / "milestone_dipole_comparison.csv").set_index("product")
    calib = pd.read_csv(cr_dir / "milestone" / "calibration_stats.csv")
    ncol, scol = f"polar_fill_north_{cap:.0f}", f"polar_fill_south_{cap:.0f}"

    def fill(product, col):
        v = mile.loc[product, col] if product in mile.index else np.nan
        return float(v) * 100.0 if np.isfinite(v) else 0.0

    times = pd.to_datetime(calib["phi_time"])
    return {
        "phi_N": fill("phi", ncol), "hmi_N": fill("hmi", ncol),
        "phi_S": fill("phi", scol), "hmi_S": fill("hmi", scol),
        "solo_b0_min": float(calib["phi_crlt_obs"].min()),
        "solo_b0_max": float(calib["phi_crlt_obs"].max()),
        "earth_b0": float(calib["hmi_crlt_obs"].mean()),
        "sep_min": float(calib["lon_separation_deg"].min()),
        "sep_max": float(calib["lon_separation_deg"].max()),
        "date_lo": _month_day(times.min()), "date_hi": _month_day(times.max()),
    }


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    args = parse_args()
    cap = args.polar_lat_deg
    cr_dirs = sorted(
        (d for d in args.campaign_dir.glob("cr_*") if (d / "milestone" / "milestone_dipole_comparison.csv").exists()),
        key=lambda d: int(re.search(r"cr_(\d+)", d.name).group(1)),
    )
    if not cr_dirs:
        raise SystemExit(f"No cr_<N>/milestone/ outputs found under {args.campaign_dir}")

    crs = [int(re.search(r"cr_(\d+)", d.name).group(1)) for d in cr_dirs]
    R = [read_rotation(d, cap) for d in cr_dirs]
    x = np.arange(len(R))

    rot_lbl = [f"CR {cr}\n{r['date_lo']} - {r['date_hi']}" for cr, r in zip(crs, R)]
    sep_lbl = [f"sep {r['sep_min']:.0f}-{r['sep_max']:.0f}°" for r in R]
    solo_min = np.array([r["solo_b0_min"] for r in R])
    solo_max = np.array([r["solo_b0_max"] for r in R])
    earth = np.array([r["earth_b0"] for r in R])
    phi_N = np.array([r["phi_N"] for r in R]); hmi_N = np.array([r["hmi_N"] for r in R])
    phi_S = np.array([r["phi_S"] for r in R]); hmi_S = np.array([r["hmi_S"] for r in R])

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10.5,
        "axes.edgecolor": MUTED, "axes.linewidth": 0.8,
        "axes.titlesize": 11.5, "axes.titleweight": "bold",
    })
    fig, (ax0, ax1, ax2) = plt.subplots(
        3, 1, figsize=(8.4, 10.2), sharex=True,
        gridspec_kw={"height_ratios": [1.15, 1, 1], "hspace": 0.16})

    # panel A: the orbital driver
    b0lim = max(20.0, np.abs(np.r_[solo_min, solo_max, earth]).max() + 4)
    ax0.axhline(0, color=GRID, lw=1.2, zorder=0)
    for i in x:
        ax0.plot([i, i], [solo_min[i], solo_max[i]], color=PHI_C, lw=7,
                 solid_capstyle="round", alpha=0.85, zorder=2,
                 label="SolO/PHI B0 range" if i == 0 else None)
        ax0.annotate(f"{solo_max[i]:+.0f}°", (i, solo_max[i]), textcoords="offset points",
                     xytext=(9, -1), color=PHI_C, fontsize=9.5, fontweight="bold", va="center")
        ax0.annotate(f"{solo_min[i]:+.0f}°", (i, solo_min[i]), textcoords="offset points",
                     xytext=(9, 1), color=PHI_C, fontsize=9.5, va="center")
    ax0.plot(x, earth, "o-", color=HMI_C, lw=2, ms=9, zorder=3, label="Earth/HMI B0")
    ax0.set_ylabel("B0  (deg)"); ax0.set_ylim(-b0lim, b0lim)
    ax0.set_title("A.  The orbital driver: SolO's heliolatitude vs Earth's")
    ax0.legend(loc="lower left", frameon=False, fontsize=9.5)
    ax0.grid(axis="y", color=GRID, lw=0.6); ax0.set_axisbelow(True)

    # panels B/C: polar-cap fill
    w = 0.34
    ytop = max(72.0, float(np.r_[phi_N, hmi_N, phi_S, hmi_S].max()) + 14)

    def fill_panel(ax, phi, hmi, title):
        b1 = ax.bar(x - w/2, phi, w, color=PHI_C, label="PHI (SolO)", zorder=3)
        b2 = ax.bar(x + w/2, hmi, w, color=HMI_C, label="HMI (SDO)", zorder=3)
        for b in list(b1) + list(b2):
            h = b.get_height()
            ax.annotate(f"{h:.0f}%", (b.get_x() + b.get_width()/2, h),
                        textcoords="offset points", xytext=(0, 3), ha="center",
                        fontsize=9.5, color=INK, fontweight="bold" if h >= 45 else "normal")
        ax.set_ylabel("cap fill  (%)"); ax.set_title(title); ax.set_ylim(0, ytop)
        ax.grid(axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True)

    fill_panel(ax1, phi_N, hmi_N, f"B.  North polar cap (≥{cap:.0f}°):  PHI vs HMI")
    ax1.legend(loc="upper left", frameon=False, fontsize=9.5, ncol=2)
    fill_panel(ax2, phi_S, hmi_S, f"C.  South polar cap (≥{cap:.0f}°):  PHI vs HMI")

    ax2.set_xticks(x)
    ax2.set_xticklabels([f"{r}\n{s}" for r, s in zip(rot_lbl, sep_lbl)], fontsize=9.5)
    ax2.tick_params(axis="x", pad=6)

    fig.suptitle("Solar Orbiter/PHI polar coverage vs SDO/HMI", fontsize=13, fontweight="bold", y=0.985)
    fig.text(0.5, 0.01,
             "Coverage is calibration-independent. Where the SolO-Earth separation is large, a per-pixel "
             "PHI+HMI merge is\nnot meaningful (opposite hemispheres) — there PHI is a standalone polar constraint.",
             ha="center", va="bottom", fontsize=8.6, color=MUTED)
    fig.subplots_adjust(top=0.93, bottom=0.145, left=0.085, right=0.965)

    out = args.out or (args.campaign_dir / "campaign_polar_advantage.png")
    fig.savefig(out, dpi=160)
    print(f"Wrote {out}  ({len(R)} rotations: {crs})")


if __name__ == "__main__":
    main()
