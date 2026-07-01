"""Synthetic validation of the milestone comparison machinery.

Ground truth is a pure axial dipole Br = G10_TRUE * sin(lat). Synthetic
observers at B0 = 0 deg (Earth/HMI-like) and B0 = 30 deg (high-latitude
SolO-like vantage) image it as Blos = Br * mu; with alpha = 1 the pipeline's
LOS->Br inversion is exact, so any error is due to geometry/binning — which
is exactly what these tests validate.

Run with: pytest tests/  (or: python tests/test_milestone.py)
"""

import sys
from pathlib import Path

import numpy as np

sys.path.append(str(Path(__file__).resolve().parents[1]))

from solar_pipeline.geometry import build_radius_arrays, estimate_mu_lat_lon
from solar_pipeline.pipeline import compute_native_disk_fields
from solar_pipeline.carrington import (
    cm_weight,
    bin_br_to_carrington_weighted,
    bin_max_to_carrington,
    combine_weighted_grids,
)
from solar_pipeline.dipole import axial_dipole_g10, polar_fill_fractions
from solar_pipeline.calibration import calibration_stats

G10_TRUE = 5.0
NLAT, NLON = 90, 180
MU_MIN = 0.4


class FakeMap:
    """Duck-typed stand-in for sunpy.map.Map (only .data/.meta are used)."""

    def __init__(self, data, meta):
        self.data = data
        self.meta = meta


def make_synthetic_disk(b0_deg: float, l0_deg: float, nx: int = 400, g10: float = G10_TRUE) -> FakeMap:
    crpix = (nx + 1) / 2.0
    rsun_pix = 0.45 * nx
    meta = {
        "crpix1": crpix,
        "crpix2": crpix,
        "cdelt1": 1.0,
        "cdelt2": 1.0,
        "rsun_obs": rsun_pix,  # with cdelt=1, rsun in arcsec == rsun in pixels
        "crlt_obs": b0_deg,
        "crln_obs": l0_deg,
    }
    _, _, dx, dy, _, rr_norm, rp = build_radius_arrays(
        shape=(nx, nx), crpix1=crpix, crpix2=crpix, cdelt1=1.0, cdelt2=1.0, rsun_arcsec=rsun_pix
    )
    mu, lat, lon, cmd = estimate_mu_lat_lon(dx=dx, dy=dy, rsun_pix=rp, b0_deg=b0_deg, l0_deg=l0_deg)
    blos = g10 * np.sin(lat) * mu  # radial dipole seen in the line of sight
    blos[~np.isfinite(mu)] = np.nan
    return FakeMap(blos, meta)


def full_dipole_grid() -> tuple[np.ndarray, np.ndarray]:
    lat_edges = np.linspace(-np.pi / 2, np.pi / 2, NLAT + 1)
    lat_centers = 0.5 * (lat_edges[:-1] + lat_edges[1:])
    grid = G10_TRUE * np.sin(lat_centers)[:, None] * np.ones((NLAT, NLON))
    return grid, lat_centers


def bin_native(fake_map: FakeMap):
    fields = compute_native_disk_fields(fake_map, disk_fraction=0.98, mu_min=MU_MIN, alpha=1.0)
    w = cm_weight(fields["cmd"])
    wsum, weight, lat_c, lon_c = bin_br_to_carrington_weighted(
        fields["br"], fields["lat"], fields["lon"], fields["valid"], w, nlat=NLAT, nlon=NLON
    )
    grid, total_w = combine_weighted_grids([wsum], [weight])
    quality = bin_max_to_carrington(
        fields["mu"], fields["lat"], fields["lon"], fields["valid"], nlat=NLAT, nlon=NLON
    )
    return grid, total_w, quality, lat_c


def test_g10_exact_on_full_grid():
    grid, lat_c = full_dipole_grid()
    for mode in ("zero", "project", "polar_extend"):
        dip = axial_dipole_g10(grid, lat_c, mode=mode)
        assert abs(dip["g10"] - G10_TRUE) < 0.01 * G10_TRUE, (mode, dip)
        # pure axial dipole: hemispheres contribute equally
        assert abs(dip["g10_north"] - dip["g10_south"]) < 0.01 * G10_TRUE


def test_g10_partial_coverage_mode_ordering():
    grid, lat_c = full_dipole_grid()
    grid = grid.copy()
    grid[np.abs(np.rad2deg(lat_c)) > 50.0, :] = np.nan  # no polar data

    g_zero = axial_dipole_g10(grid, lat_c, mode="zero")["g10"]
    g_proj = axial_dipole_g10(grid, lat_c, mode="project")["g10"]
    g_ext = axial_dipole_g10(grid, lat_c, mode="polar_extend")["g10"]

    # projection is unbiased for a pure dipole; zero-fill loses the
    # high-latitude (high sin*cos weight) contribution
    assert abs(g_proj - G10_TRUE) < 0.02 * G10_TRUE
    assert g_zero < 0.75 * G10_TRUE
    assert g_zero < g_ext <= G10_TRUE * 1.02


def test_polar_visibility_earth_vs_solo():
    hmi_like = make_synthetic_disk(b0_deg=0.0, l0_deg=180.0)
    solo_like = make_synthetic_disk(b0_deg=30.0, l0_deg=180.0)

    _, w_hmi, _, lat_c = bin_native(hmi_like)
    _, w_solo, _, _ = bin_native(solo_like)

    hmi_polar = polar_fill_fractions(w_hmi, lat_c, polar_lat_deg=70.0)
    solo_polar = polar_fill_fractions(w_solo, lat_c, polar_lat_deg=70.0)

    # from the ecliptic, mu_min=0.4 cuts everything poleward of ~66 deg
    assert hmi_polar["north"] == 0.0
    assert hmi_polar["south"] == 0.0
    # the B0=30 vantage sees the north cap (pole itself has mu=sin30=0.5)
    assert solo_polar["north"] > 0.3
    # ... at the cost of the south cap
    assert solo_polar["south"] == 0.0


def test_g10_recovery_from_single_vantage():
    for b0 in (0.0, 30.0):
        grid, _, _, lat_c = bin_native(make_synthetic_disk(b0_deg=b0, l0_deg=180.0))
        g_proj = axial_dipole_g10(grid, lat_c, mode="project")["g10"]
        # noise-free dipole: projection recovers the truth from either vantage
        assert abs(g_proj - G10_TRUE) < 0.05 * G10_TRUE, (b0, g_proj)
        # zero-fill on a half-sphere map loses roughly half the signal
        g_zero = axial_dipole_g10(grid, lat_c, mode="zero")["g10"]
        assert g_zero < 0.75 * G10_TRUE


def test_polar_quality_mask():
    _, _, quality, lat_c = bin_native(make_synthetic_disk(b0_deg=30.0, l0_deg=180.0))
    lat_deg = np.rad2deg(lat_c)
    north_cap = quality[lat_deg > 70.0, :]
    # observed north-cap bins are foreshortened: from vantage B0 the best
    # possible mu at latitude lat is cos(lat - B0), i.e. cos(70-30) here —
    # present but well below disk-centre quality, usable as confidence mask
    observed = north_cap[np.isfinite(north_cap)]
    assert observed.size > 0
    assert np.nanmax(quality) > 0.95
    assert observed.max() <= np.cos(np.deg2rad(70.0 - 30.0)) + 0.01


def test_calibration_slope_recovery():
    rng = np.random.default_rng(1)
    solo_like = make_synthetic_disk(b0_deg=10.0, l0_deg=90.0)
    fields = compute_native_disk_fields(solo_like, disk_fraction=0.98, mu_min=MU_MIN, alpha=1.0)
    ref = solo_like.data
    phi = 1.15 * ref + rng.normal(0.0, 0.05, ref.shape)
    stats = calibration_stats(phi, ref, fields["mu"], mu_min=MU_MIN, min_abs_ref=0.5)
    assert abs(stats["slope"] - 1.15) < 0.01
    assert stats["pearson_r"] > 0.99
    assert stats["n_pixels"] > 1000


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"{name}: OK")
    print("All tests passed.")
