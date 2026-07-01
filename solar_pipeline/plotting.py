from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def make_baseline_plots(df: pd.DataFrame, plots_dir: Path, title_suffix: str = ""):
    plots_dir.mkdir(exist_ok=True, parents=True)

    data = df.copy()
    data["phi_time"] = pd.to_datetime(data["phi_time"], utc=True)
    data = data.sort_values("phi_time").reset_index(drop=True)
    data["case_label"] = data["phi_time"].dt.strftime("%m-%d %H:%M")

    x = np.arange(len(data))

    plt.figure(figsize=(9, 5))
    plt.plot(x, data["dip_phi"], marker="o", label="PHI-only")
    plt.plot(x, data["dip_hmi"], marker="o", label="HMI-on-PHI")
    plt.plot(x, data["dip_merged"], marker="o", label="Merged")
    plt.xticks(x, data["case_label"], rotation=30, ha="right")
    plt.ylabel("Carrington-style dipole")
    plt.title(f"Baseline dipole series{title_suffix}")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(plots_dir / "dipole_series.png", dpi=150)
    plt.close()

    plt.figure(figsize=(9, 5))
    plt.plot(x, data["merged_minus_phi"], marker="o", label="Merged - PHI")
    plt.plot(x, data["merged_minus_hmi"], marker="o", label="Merged - HMI")
    plt.axhline(0.0, color="k", linestyle="--", linewidth=1)
    plt.xticks(x, data["case_label"], rotation=30, ha="right")
    plt.ylabel("Dipole offset")
    plt.title(f"Baseline merged dipole offsets{title_suffix}")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(plots_dir / "dipole_offsets.png", dpi=150)
    plt.close()

    plt.figure(figsize=(8, 4))
    plt.plot(x, data["time_diff_sec"], marker="o")
    plt.xticks(x, data["case_label"], rotation=30, ha="right")
    plt.ylabel("|Δt| [s]")
    plt.title("PHI-HMI time differences")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(plots_dir / "time_differences.png", dpi=150)
    plt.close()


def plot_carrington_map(grid, lat_centers, lon_centers, path: Path, title: str = ""):
    plt.figure(figsize=(10, 5))
    extent = [
        np.rad2deg(lon_centers[0]),
        np.rad2deg(lon_centers[-1]),
        np.rad2deg(lat_centers[0]),
        np.rad2deg(lat_centers[-1]),
    ]
    vmax = np.nanmax(np.abs(grid)) if np.any(np.isfinite(grid)) else 1.0
    plt.imshow(grid, origin="lower", extent=extent, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    plt.colorbar(label="Br [G]")
    plt.xlabel("Carrington longitude [deg]")
    plt.ylabel("Latitude [deg]")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close()