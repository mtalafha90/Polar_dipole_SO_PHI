"""The resolution-degradation calibration test (scripts/calibration_resolution_test.py).

Validates the core mechanism: if PHI is a coarser-resolution view of the
same field HMI sees, the through-origin regression slope PHI ~= slope*HMI
sits below 1 (fine mixed-polarity structure cancels within a PHI resolution
element), and degrading HMI to PHI's resolution raises the slope back toward
1. This is exactly the attribution the script is meant to make.
"""

import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter

sys.path.append(str(Path(__file__).resolve().parents[1]))

from scripts.calibration_resolution_test import nan_gaussian
from solar_pipeline.calibration import calibration_stats


def test_nan_gaussian_preserves_mask_and_ignores_nans():
    a = np.ones((20, 20))
    a[5:8, 5:8] = np.nan
    out = nan_gaussian(a, fwhm_pix=3.0)
    assert np.all(np.isnan(out[5:8, 5:8]))          # mask preserved
    assert np.isfinite(out[0, 0])
    assert abs(out[0, 0] - 1.0) < 1e-6              # smoothing a constant is a no-op


def test_degrading_hmi_recovers_the_slope_deficit():
    rng = np.random.default_rng(0)
    truth = rng.normal(0.0, 50.0, (200, 200))       # fine mixed-polarity field
    mu = np.ones_like(truth)
    sigma = 3.0
    phi = gaussian_filter(truth, sigma)             # PHI: coarser view
    hmi = truth                                     # HMI: sharp view

    kw = dict(mu_min=0.0, min_abs_ref=0.0)
    slope0 = calibration_stats(phi, hmi, mu, **kw)["slope"]
    # matched smoothing: nan_gaussian takes FWHM = sigma * 2.3548
    slope_matched = calibration_stats(phi, nan_gaussian(hmi, sigma * 2.3548), mu, **kw)["slope"]

    assert slope0 < 0.9                             # coarse PHI reads low vs sharp HMI
    assert slope_matched > slope0 + 0.1             # degrading HMI closes the deficit
    assert abs(slope_matched - 1.0) < 0.05          # ... essentially to 1


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"{name}: OK")
    print("All tests passed.")
