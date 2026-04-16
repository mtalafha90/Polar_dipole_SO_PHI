from __future__ import annotations

import numpy as np


def cosine_taper_weight(rr_norm, r_inner: float, r_outer: float):
    w = np.zeros_like(rr_norm, dtype=float)
    w[rr_norm <= r_inner] = 1.0
    w[rr_norm >= r_outer] = 0.0

    mid = (rr_norm > r_inner) & (rr_norm < r_outer)
    x = (rr_norm[mid] - r_inner) / (r_outer - r_inner)
    w[mid] = 0.5 * (1.0 + np.cos(np.pi * x))
    return w


def smooth_merge(phi_data, hmi_data, rr_norm, disk_fraction: float, r_inner: float, r_outer: float):
    disk_mask = rr_norm <= disk_fraction

    w_phi = cosine_taper_weight(rr_norm, r_inner, r_outer)
    w_hmi = 1.0 - w_phi

    w_phi = np.where(disk_mask, w_phi, 0.0)
    w_hmi = np.where(disk_mask, w_hmi, 0.0)

    merged = np.full_like(phi_data, np.nan, dtype=float)

    valid_both = disk_mask & np.isfinite(phi_data) & np.isfinite(hmi_data)
    merged[valid_both] = w_phi[valid_both] * phi_data[valid_both] + w_hmi[valid_both] * hmi_data[valid_both]

    phi_only = disk_mask & np.isnan(merged) & np.isfinite(phi_data)
    merged[phi_only] = phi_data[phi_only]

    hmi_only = disk_mask & np.isnan(merged) & np.isfinite(hmi_data)
    merged[hmi_only] = hmi_data[hmi_only]

    return merged, disk_mask