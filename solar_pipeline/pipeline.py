from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import sunpy.map
from reproject import reproject_interp

from .io_utils import (
    parse_phi_time,
    load_map,
    find_nearest_hmi,
)
from .geometry import build_radius_arrays, estimate_mu_lat_lon
from .blending import smooth_merge
from .radial import los_to_br
from .carrington import (
    bin_br_to_carrington,
    axial_dipole_from_carrington_grid,
    carrington_fill_fraction,
)


def compute_case_fields(
    phi_blos_path: Path,
    hmi_index,
    *,
    max_time_diff_sec: float,
    r_inner: float,
    r_outer: float,
    disk_fraction: float,
    mu_min: float,
    alpha: float,
):
    """Match, reproject, blend, and convert a single PHI/HMI case to Br.

    Returns the per-pixel fields (mu, lat, lon, central-meridian distance,
    Br for PHI/HMI/merged and their validity masks) shared by both the
    per-case baseline pipeline and multi-case synoptic assimilation.
    """
    phi_time = parse_phi_time(phi_blos_path)
    hmi_path, hmi_time, time_diff_sec = find_nearest_hmi(phi_time, hmi_index)

    if time_diff_sec > max_time_diff_sec:
        raise RuntimeError(
            f"No suitable HMI match within {max_time_diff_sec} s. "
            f"Nearest: {hmi_path.name}, Δt = {time_diff_sec:.1f} s"
        )

    phi_blos = load_map(phi_blos_path)
    hmi = load_map(hmi_path)

    reprojected_hmi_data, footprint = reproject_interp(
        hmi, phi_blos.wcs, shape_out=phi_blos.data.shape
    )
    hmi_on_phi = sunpy.map.Map(reprojected_hmi_data, phi_blos.meta.copy())

    crpix1 = float(phi_blos.meta["crpix1"])
    crpix2 = float(phi_blos.meta["crpix2"])
    cdelt1 = float(phi_blos.meta["cdelt1"])
    cdelt2 = float(phi_blos.meta["cdelt2"])
    rsun_arcsec = float(phi_blos.meta.get("rsun_obs", phi_blos.meta.get("rsun_arc")))
    b0_deg = float(phi_blos.meta.get("crlt_obs", 0.0))
    l0_deg = float(phi_blos.meta.get("crln_obs", 0.0))

    _, _, dx, dy, _, rr_norm, rsun_pix = build_radius_arrays(
        shape=phi_blos.data.shape,
        crpix1=crpix1,
        crpix2=crpix2,
        cdelt1=cdelt1,
        cdelt2=cdelt2,
        rsun_arcsec=rsun_arcsec,
    )

    merged, disk_mask = smooth_merge(
        phi_blos.data,
        hmi_on_phi.data,
        rr_norm,
        disk_fraction=disk_fraction,
        r_inner=r_inner,
        r_outer=r_outer,
    )

    mu, lat, lon, cmd = estimate_mu_lat_lon(
        dx=dx, dy=dy, rsun_pix=rsun_pix, b0_deg=b0_deg, l0_deg=l0_deg
    )

    br_phi, valid_phi = los_to_br(phi_blos.data, mu, mu_min=mu_min, alpha=alpha)
    br_hmi, valid_hmi = los_to_br(hmi_on_phi.data, mu, mu_min=mu_min, alpha=alpha)
    br_merged, valid_merged = los_to_br(merged, mu, mu_min=mu_min, alpha=alpha)

    return {
        "phi_blos": phi_blos,
        "hmi": hmi,
        "hmi_path": hmi_path,
        "phi_time": phi_time,
        "hmi_time": hmi_time,
        "time_diff_sec": time_diff_sec,
        "merged": merged,
        "mu": mu,
        "lat": lat,
        "lon": lon,
        "cmd": cmd,
        "br_phi": br_phi,
        "br_hmi": br_hmi,
        "br_merged": br_merged,
        "valid_phi": valid_phi,
        "valid_hmi": valid_hmi,
        "valid_merged": valid_merged,
    }


def run_case(
    phi_blos_path: Path,
    hmi_index,
    *,
    max_time_diff_sec: float,
    r_inner: float,
    r_outer: float,
    disk_fraction: float,
    mu_min: float,
    alpha: float,
    nlat: int,
    nlon: int,
):
    fields = compute_case_fields(
        phi_blos_path,
        hmi_index,
        max_time_diff_sec=max_time_diff_sec,
        r_inner=r_inner,
        r_outer=r_outer,
        disk_fraction=disk_fraction,
        mu_min=mu_min,
        alpha=alpha,
    )

    grid_phi, count_phi, lat_centers, lon_centers = bin_br_to_carrington(
        fields["br_phi"], fields["lat"], fields["lon"], fields["valid_phi"], nlat=nlat, nlon=nlon
    )
    grid_hmi, count_hmi, _, _ = bin_br_to_carrington(
        fields["br_hmi"], fields["lat"], fields["lon"], fields["valid_hmi"], nlat=nlat, nlon=nlon
    )
    grid_merged, count_merged, _, _ = bin_br_to_carrington(
        fields["br_merged"], fields["lat"], fields["lon"], fields["valid_merged"], nlat=nlat, nlon=nlon
    )

    dip_phi = axial_dipole_from_carrington_grid(grid_phi, lat_centers)
    dip_hmi = axial_dipole_from_carrington_grid(grid_hmi, lat_centers)
    dip_merged = axial_dipole_from_carrington_grid(grid_merged, lat_centers)

    phi_blos = fields["phi_blos"]
    hmi = fields["hmi"]

    row = {
        "phi_blos_file": phi_blos_path.name,
        "phi_time": fields["phi_time"].isoformat(),
        "hmi_file": fields["hmi_path"].name,
        "hmi_time": fields["hmi_time"].isoformat(),
        "time_diff_sec": float(fields["time_diff_sec"]),
        "dip_phi": dip_phi,
        "dip_hmi": dip_hmi,
        "dip_merged": dip_merged,
        "merged_minus_phi": dip_merged - dip_phi,
        "merged_minus_hmi": dip_merged - dip_hmi,
        "phi_crln_obs": float(phi_blos.meta.get("crln_obs", np.nan)),
        "phi_crlt_obs": float(phi_blos.meta.get("crlt_obs", np.nan)),
        "hmi_crln_obs": float(hmi.meta.get("crln_obs", np.nan)),
        "hmi_crlt_obs": float(hmi.meta.get("crlt_obs", np.nan)),
        "fill_phi": carrington_fill_fraction(count_phi),
        "fill_hmi": carrington_fill_fraction(count_hmi),
        "fill_merged": carrington_fill_fraction(count_merged),
    }

    arrays = {
        "merged": fields["merged"],
        "phi_header": phi_blos.fits_header,
        "grid_phi": grid_phi,
        "grid_hmi": grid_hmi,
        "grid_merged": grid_merged,
        "lat_centers": lat_centers,
        "lon_centers": lon_centers,
    }

    return row, arrays


def summarize_dataframe(df: pd.DataFrame) -> dict[str, dict[str, float]]:
    out = {}
    for col in [
        "dip_phi",
        "dip_hmi",
        "dip_merged",
        "merged_minus_phi",
        "merged_minus_hmi",
        "fill_phi",
        "fill_hmi",
        "fill_merged",
    ]:
        if col in df.columns:
            vals = pd.to_numeric(df[col], errors="coerce").dropna()
            if len(vals) > 0:
                out[col] = {
                    "mean": float(vals.mean()),
                    "std": float(vals.std()),
                    "min": float(vals.min()),
                    "max": float(vals.max()),
                }
    return out
