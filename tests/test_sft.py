"""Physics validation of the ported SFT model (solar_pipeline/sft.py).

The key check is analytic: with no flow, no source, and no decay term, an
l=1 dipole profile B = B0*sin(lat) is an eigenmode of the surface diffusion
operator and must decay at exactly rate 2*eta/Rsun^2. Any porting error in
the diffusive term, grid, or boundary conditions breaks this.

Run with: pytest tests/  (or: python tests/test_sft.py)
"""

import sys
from pathlib import Path

import numpy as np

sys.path.append(str(Path(__file__).resolve().parents[1]))

from solar_pipeline.sft import (
    RSUN_MM,
    KM2S_TO_MM2_PER_DAY,
    SFTModel,
    HathawayJiangSource,
    apply_polar_constraint,
    axial_dipole_moment,
    polar_cap_mean,
    reversal_times,
    zonal_profile_from_map,
)


def _dipole_grid(nlat=90, nlon=180, observed_max_deg=None):
    lat_edges = np.linspace(-np.pi / 2, np.pi / 2, nlat + 1)
    lat_centers = 0.5 * (lat_edges[:-1] + lat_edges[1:])
    grid = 5.0 * np.sin(lat_centers)[:, None] * np.ones((nlat, nlon))
    if observed_max_deg is not None:
        grid[np.abs(np.rad2deg(lat_centers)) > observed_max_deg, :] = np.nan
    return grid, lat_centers


def test_zonal_nan_mode_masks_unobserved():
    # observed only equatorward of 55 deg
    grid, lat_centers = _dipole_grid(observed_max_deg=55.0)
    model = SFTModel(flowtype=2, u0=0.0, eta=250.0)
    b_nan = zonal_profile_from_map(grid, lat_centers, model.latitude, unobserved="nan")
    assert np.all(np.isnan(b_nan[np.abs(model.latitude) > 60]))       # caps masked
    assert np.all(np.isfinite(b_nan[np.abs(model.latitude) < 45]))    # mid observed


def test_apply_polar_constraint_splices_only_selected_cap():
    model = SFTModel(flowtype=2, u0=0.0, eta=250.0)
    lat = model.latitude
    b_base = np.zeros_like(lat)                       # HMI: nothing in the caps
    b_polar = np.where(np.abs(lat) >= 55, 4.0, np.nan)  # PHI: +4 G in both caps
    out = apply_polar_constraint(b_base, b_polar, lat, polar_lat_deg=60.0,
                                 hemisphere="north", blend_deg=10.0)
    assert np.allclose(out[lat > 65], 4.0)            # north cap took PHI
    assert np.allclose(out[lat < -65], 0.0)           # south cap untouched
    assert np.allclose(out[np.abs(lat) < 40], 0.0)    # mid-latitudes untouched


def test_apply_polar_constraint_keeps_base_where_polar_missing():
    lat = np.linspace(90, -90, 181)
    b_base = np.full_like(lat, 1.0)
    b_polar = np.full_like(lat, np.nan)               # PHI observed nothing
    out = apply_polar_constraint(b_base, b_polar, lat, hemisphere="both")
    assert np.allclose(out, b_base)                   # base is preserved everywhere


def test_dipole_diffusive_decay_rate():
    eta = 250.0  # km^2/s
    model = SFTModel(flowtype=2, u0=0.0, eta=eta, tau_years=None)
    b0 = 5.0 * np.sin(np.deg2rad(model.latitude))
    times, history = model.run(b0, years=4.0, record_every_days=27.0)
    dips = np.array([axial_dipole_moment(b, model.theta) for b in history])

    # fit exponential decay rate over the run
    rate_fit = -np.polyfit(times * 365.25, np.log(dips / dips[0]), 1)[0]  # per day
    rate_true = 2.0 * eta * KM2S_TO_MM2_PER_DAY / RSUN_MM**2
    assert abs(rate_fit - rate_true) / rate_true < 0.02, (rate_fit, rate_true)


def test_flux_conservation_without_decay():
    model = SFTModel(flowtype=2, u0=11.0, eta=250.0, tau_years=None)
    rng = np.random.default_rng(3)
    b0 = rng.normal(0.0, 1.0, model.n)
    times, history = model.run(b0, years=2.0)
    flux = [np.trapezoid(b * np.sin(model.theta), model.theta) for b in history]
    unsigned = np.trapezoid(np.abs(history[0]) * np.sin(model.theta), model.theta)
    # transport terms conserve signed flux up to the scheme's known small
    # near-pole leakage (the original notes latitude0 "tricks" for this);
    # measured drift is ~4e-5 relative — a porting error would be O(1)
    assert abs(flux[-1] - flux[0]) < 1e-3 * unsigned


def test_poleward_transport_builds_polar_field():
    model = SFTModel(flowtype=2, u0=11.0, eta=250.0, tau_years=None)
    # mid-latitude flux band (positive, northern hemisphere)
    b0 = np.exp(-((model.latitude - 30.0) ** 2) / (2 * 8.0**2))
    times, history = model.run(b0, years=6.0)
    cap0 = polar_cap_mean(history[0], model.latitude, 70.0)["north"]
    cap_end = polar_cap_mean(history[-1], model.latitude, 70.0)["north"]
    assert cap0 < 0.05
    assert cap_end > 5 * max(cap0, 1e-6)


def test_source_is_flux_balanced_and_cyclic():
    src = HathawayJiangSource(tau_days=10.0 * 365.25, sigma=0.0)
    lat = np.linspace(90, -90, 181)
    theta = np.deg2rad(90.0 - lat)
    s_mid = src(lat, 4.0 * 365.25)  # mid-cycle: active source
    assert np.max(np.abs(s_mid)) > 0
    net = np.trapezoid(s_mid * np.sin(theta), theta)
    assert abs(net) < 1e-8 * np.max(np.abs(s_mid))


def test_sft_reversal_with_source():
    # the cycle-1 source drives g10 negative (measured), so an initial
    # positive dipole must reverse within the cycle
    model = SFTModel(flowtype=2, u0=11.0, eta=250.0, tau_years=10.0)
    src = HathawayJiangSource(tau_days=10.0 * 365.25, sigma=0.0)
    b0 = 0.5 * np.sin(np.deg2rad(model.latitude))
    times, history = model.run(b0, years=11.0, source_fn=src)
    dips = [axial_dipole_moment(b, model.theta) for b in history]
    revs = reversal_times(times, dips)
    assert len(revs) >= 1
    assert 0.0 < revs[0] < 11.0


def test_zonal_injection_polar_handicap():
    # synthetic map: dipole observed only equatorward of 55 deg (Earth-like)
    nlat, nlon = 90, 180
    lat_edges = np.linspace(-np.pi / 2, np.pi / 2, nlat + 1)
    lat_centers = 0.5 * (lat_edges[:-1] + lat_edges[1:])
    grid = 5.0 * np.sin(lat_centers)[:, None] * np.ones((nlat, nlon))
    grid[np.abs(np.rad2deg(lat_centers)) > 55.0, :] = np.nan

    model = SFTModel(flowtype=2, u0=0.0, eta=250.0)
    b_zero = zonal_profile_from_map(grid, lat_centers, model.latitude, unobserved="zero")
    b_ext = zonal_profile_from_map(grid, lat_centers, model.latitude, unobserved="extend")

    # observed latitudes reproduce the dipole profile
    mid = np.abs(model.latitude) < 50
    assert np.allclose(b_zero[mid], 5.0 * np.sin(np.deg2rad(model.latitude[mid])), atol=0.05)
    # 'zero' zeroes the caps; 'extend' holds the last observed value
    polar = model.latitude > 60
    assert np.all(b_zero[polar] == 0.0)
    assert np.all(b_ext[polar] > 3.5)
    # and the resulting initial dipole differs accordingly
    d_zero = axial_dipole_moment(b_zero, model.theta)
    d_ext = axial_dipole_moment(b_ext, model.theta)
    assert d_zero < 0.8 * d_ext


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"{name}: OK")
    print("All tests passed.")
