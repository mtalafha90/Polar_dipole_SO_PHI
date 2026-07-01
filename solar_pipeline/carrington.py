from __future__ import annotations

import numpy as np


def bin_br_to_carrington(br, lat, lon, valid_mask, nlat: int, nlon: int):
    lat_edges = np.linspace(-np.pi / 2, np.pi / 2, nlat + 1)
    lon_edges = np.linspace(0.0, 2.0 * np.pi, nlon + 1)

    lat_centers = 0.5 * (lat_edges[:-1] + lat_edges[1:])
    lon_centers = 0.5 * (lon_edges[:-1] + lon_edges[1:])

    grid_sum = np.zeros((nlat, nlon), dtype=float)
    grid_count = np.zeros((nlat, nlon), dtype=float)

    lat_v = lat[valid_mask]
    lon_v = lon[valid_mask]
    br_v = br[valid_mask]

    ilat = np.digitize(lat_v, lat_edges) - 1
    ilon = np.digitize(lon_v, lon_edges) - 1

    ok = (ilat >= 0) & (ilat < nlat) & (ilon >= 0) & (ilon < nlon)
    ilat = ilat[ok]
    ilon = ilon[ok]
    br_v = br_v[ok]

    np.add.at(grid_sum, (ilat, ilon), br_v)
    np.add.at(grid_count, (ilat, ilon), 1.0)

    grid = np.full((nlat, nlon), np.nan, dtype=float)
    good = grid_count > 0
    grid[good] = grid_sum[good] / grid_count[good]

    return grid, grid_count, lat_centers, lon_centers


def axial_dipole_from_carrington_grid(br_grid, lat_centers):
    lat2d = lat_centers[:, None]
    weights = np.cos(lat2d)

    valid = np.isfinite(br_grid)
    if not np.any(valid):
        raise RuntimeError("No valid bins in Carrington grid.")

    signal = br_grid * np.sin(lat2d)
    w = np.where(valid, weights, 0.0)
    s = np.where(valid, signal * weights, 0.0)

    return float(np.sum(s) / np.sum(w))


def carrington_fill_fraction(grid_count):
    return float(np.count_nonzero(grid_count > 0) / grid_count.size)


def cm_weight(cmd, power: float = 1.0):
    """Per-pixel weight from central-meridian distance `cmd` (radians).

    cos(cmd)**power, clipped to zero beyond +/-90 deg so far-limb pixels
    (already excluded from most cases by MU_MIN/DISK_FRACTION) get no
    influence when combining several cases into one synoptic map.
    """
    w = np.zeros_like(cmd, dtype=float)
    finite = np.isfinite(cmd)
    w[finite] = np.clip(np.cos(cmd[finite]), 0.0, None) ** power
    return w


def bin_br_to_carrington_weighted(br, lat, lon, valid_mask, weight, nlat: int, nlon: int):
    """Like bin_br_to_carrington, but accumulates weight*br and weight per
    bin instead of an unweighted mean, so several cases can later be summed
    together (see combine_weighted_grids) into one assimilated map.
    """
    lat_edges = np.linspace(-np.pi / 2, np.pi / 2, nlat + 1)
    lon_edges = np.linspace(0.0, 2.0 * np.pi, nlon + 1)

    lat_centers = 0.5 * (lat_edges[:-1] + lat_edges[1:])
    lon_centers = 0.5 * (lon_edges[:-1] + lon_edges[1:])

    grid_wsum = np.zeros((nlat, nlon), dtype=float)
    grid_weight = np.zeros((nlat, nlon), dtype=float)

    lat_v = lat[valid_mask]
    lon_v = lon[valid_mask]
    br_v = br[valid_mask]
    w_v = weight[valid_mask]

    ilat = np.digitize(lat_v, lat_edges) - 1
    ilon = np.digitize(lon_v, lon_edges) - 1

    ok = (ilat >= 0) & (ilat < nlat) & (ilon >= 0) & (ilon < nlon) & (w_v > 0)
    ilat = ilat[ok]
    ilon = ilon[ok]
    br_v = br_v[ok]
    w_v = w_v[ok]

    np.add.at(grid_wsum, (ilat, ilon), w_v * br_v)
    np.add.at(grid_weight, (ilat, ilon), w_v)

    return grid_wsum, grid_weight, lat_centers, lon_centers


def combine_weighted_grids(wsum_list, weight_list):
    """Combine per-case (weighted-sum, weight) grids into one assimilated map."""
    total_wsum = np.sum(np.stack(wsum_list), axis=0)
    total_weight = np.sum(np.stack(weight_list), axis=0)

    grid = np.full(total_wsum.shape, np.nan, dtype=float)
    good = total_weight > 0
    grid[good] = total_wsum[good] / total_weight[good]
    return grid, total_weight