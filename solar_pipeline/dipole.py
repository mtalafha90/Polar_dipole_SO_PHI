"""Axial dipole moment (g10) from Carrington-style Br grids.

The baseline pipeline's `axial_dipole_from_carrington_grid` is a weighted
mean over observed bins and is kept unchanged as the baseline-v1 proxy.
This module provides the standard spherical-harmonic coefficient

    g10 = (3 / 4pi) * integral( Br(lat, lon) * sin(lat) * dOmega )

with explicit, selectable assumptions about unobserved bins, which is the
quantity needed for SFT comparisons and polar-filling sensitivity studies.
"""

from __future__ import annotations

import numpy as np

FILL_MODES = ("zero", "project", "polar_extend")


def _bin_areas(lat_centers, nlon):
    """Solid-angle element per bin: cos(lat) * dlat * dlon."""
    nlat = len(lat_centers)
    dlat = np.pi / nlat
    dlon = 2.0 * np.pi / nlon
    return np.cos(lat_centers)[:, None] * dlat * dlon * np.ones((nlat, nlon))


def _polar_extend(grid, lat_centers):
    """Zero-fill everywhere except the polar caps, which are filled with the
    zonal-mean Br of the highest observed latitude band in that hemisphere.

    Only bins poleward of the last observed band are treated as "cap"; gaps
    in longitude at observed latitudes are zero-filled like mode "zero".
    """
    filled = np.where(np.isfinite(grid), grid, 0.0)
    observed_rows = np.where(np.any(np.isfinite(grid), axis=1))[0]
    if len(observed_rows) == 0:
        return filled

    nlat = grid.shape[0]
    # southern cap: below the lowest observed band
    lo = observed_rows[0]
    if lo > 0:
        filled[:lo, :] = np.nanmean(grid[lo, :])
    # northern cap: above the highest observed band
    hi = observed_rows[-1]
    if hi < nlat - 1:
        filled[hi + 1:, :] = np.nanmean(grid[hi, :])
    return filled


def axial_dipole_g10(br_grid, lat_centers, mode: str = "zero") -> dict[str, float]:
    """Compute g10 (units of `br_grid`, typically Gauss) plus its
    hemispheric decomposition (g10 = g10_north + g10_south).

    Modes for unobserved (NaN) bins:
    - "zero": they contribute nothing to the integral. Biases |g10| low
      when coverage is partial, in a way that mimics assuming no polar flux.
    - "project": least-squares projection of the observed bins onto the
      dipole profile Br = g10*sin(lat), area-weighted. Unbiased if the
      underlying field is dipole-dominated; equals the full integral at
      complete coverage.
    - "polar_extend": like "zero" but polar caps poleward of the last
      observed latitude band are filled with that band's zonal mean.
    """
    if mode not in FILL_MODES:
        raise ValueError(f"mode must be one of {FILL_MODES}, got {mode!r}")

    lat2d = lat_centers[:, None]
    area = _bin_areas(lat_centers, br_grid.shape[1])
    north = (lat2d > 0)

    if mode == "project":
        obs = np.isfinite(br_grid)
        if not np.any(obs):
            raise RuntimeError("No observed bins in grid.")
        w = np.where(obs, area, 0.0)
        num = np.where(obs, br_grid * np.sin(lat2d) * w, 0.0)
        den = float(np.sum(np.sin(lat2d) ** 2 * w))
        g10_n = float(np.sum(num[np.broadcast_to(north, br_grid.shape)]) / den)
        g10_s = float(np.sum(num[~np.broadcast_to(north, br_grid.shape)]) / den)
        return {"g10": g10_n + g10_s, "g10_north": g10_n, "g10_south": g10_s}

    filled = np.where(np.isfinite(br_grid), br_grid, 0.0) if mode == "zero" else _polar_extend(br_grid, lat_centers)

    integrand = (3.0 / (4.0 * np.pi)) * filled * np.sin(lat2d) * area
    g10_n = float(np.sum(integrand[np.broadcast_to(north, br_grid.shape)]))
    g10_s = float(np.sum(integrand[~np.broadcast_to(north, br_grid.shape)]))
    return {"g10": g10_n + g10_s, "g10_north": g10_n, "g10_south": g10_s}


def axial_dipole_g10_sinlat(br_grid) -> float:
    """g10 for a full-sphere map on a uniform sine-latitude (CEA) grid,
    e.g. HMI synoptic charts (hmi.Synoptic_Mr*). With s = sin(lat) uniform,
    dOmega = ds dlon, so g10 = 3 * mean(Br * s). NaN bins contribute zero.
    """
    nlat = br_grid.shape[0]
    s_edges = np.linspace(-1.0, 1.0, nlat + 1)
    s = 0.5 * (s_edges[:-1] + s_edges[1:])
    filled = np.where(np.isfinite(br_grid), br_grid, 0.0)
    return float(3.0 * np.mean(filled * s[:, None]))


def polar_fill_fractions(grid_count, lat_centers, polar_lat_deg: float = 60.0) -> dict[str, float]:
    """Fraction of polar-cap bins (|lat| > polar_lat_deg) that are observed,
    per hemisphere. `grid_count` is a per-bin count or weight grid."""
    lat_deg = np.rad2deg(lat_centers)
    out = {}
    for name, sel in (("north", lat_deg > polar_lat_deg), ("south", lat_deg < -polar_lat_deg)):
        cap = grid_count[sel, :]
        out[name] = float(np.count_nonzero(cap > 0) / cap.size) if cap.size else np.nan
    return out
