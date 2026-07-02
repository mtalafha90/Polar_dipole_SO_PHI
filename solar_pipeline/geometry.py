from __future__ import annotations

import numpy as np


def build_radius_arrays(shape, crpix1: float, crpix2: float, cdelt1: float, cdelt2: float, rsun_arcsec: float):
    ny, nx = shape
    y, x = np.indices((ny, nx), dtype=float)

    x0 = crpix1 - 1.0
    y0 = crpix2 - 1.0

    pixscale = 0.5 * (abs(cdelt1) + abs(cdelt2))
    rsun_pix = rsun_arcsec / pixscale

    dx = x - x0
    dy = y - y0
    rr_pix = np.sqrt(dx**2 + dy**2)
    rr_norm = rr_pix / rsun_pix
    return x, y, dx, dy, rr_pix, rr_norm, rsun_pix


def rotate_offsets(dx, dy, crota2_deg: float):
    """Rotate pixel offsets into the solar north-up frame for a map whose
    image axes are rotated by CROTA2 (e.g. HMI's ~180 deg camera rotation).
    No-op for crota2 = 0, which is why the baseline PHI path is unchanged."""
    if abs(crota2_deg) < 1e-6:
        return dx, dy
    rho = np.deg2rad(crota2_deg)
    return dx * np.cos(rho) - dy * np.sin(rho), dx * np.sin(rho) + dy * np.cos(rho)


def estimate_mu_lat_lon(dx, dy, rsun_pix: float, b0_deg: float, l0_deg: float):
    x = dx / rsun_pix
    y = dy / rsun_pix
    r2 = x**2 + y**2

    mu = np.full_like(x, np.nan, dtype=float)
    on_disk = r2 <= 1.0
    mu[on_disk] = np.sqrt(np.clip(1.0 - r2[on_disk], 0.0, 1.0))

    b0 = np.deg2rad(b0_deg)
    l0 = np.deg2rad(l0_deg)
    z = np.where(np.isfinite(mu), mu, np.nan)

    lat = np.full_like(x, np.nan, dtype=float)
    lon = np.full_like(x, np.nan, dtype=float)

    sin_lat = y * np.cos(b0) + z * np.sin(b0)
    sin_lat = np.clip(sin_lat, -1.0, 1.0)
    lat[on_disk] = np.arcsin(sin_lat[on_disk])

    cos_lat = np.cos(lat)

    cmd = np.full_like(x, np.nan, dtype=float)
    valid = on_disk & np.isfinite(cos_lat) & (np.abs(cos_lat) > 1e-6)

    sin_cmd = np.zeros_like(x)
    sin_cmd[valid] = x[valid] / cos_lat[valid]
    sin_cmd = np.clip(sin_cmd, -1.0, 1.0)

    cos_cmd = np.zeros_like(x)
    cos_cmd[valid] = (z[valid] * np.cos(b0) - y[valid] * np.sin(b0)) / cos_lat[valid]
    cos_cmd = np.clip(cos_cmd, -1.0, 1.0)

    cmd[valid] = np.arctan2(sin_cmd[valid], cos_cmd[valid])
    lon[valid] = np.mod(l0 + cmd[valid], 2.0 * np.pi)

    return mu, lat, lon, cmd