"""PHI-vs-HMI cross-calibration diagnostics.

PHI-FDT and HMI LOS magnetograms have known scale differences. For each
matched case the two maps share the PHI pixel grid after reprojection, so a
simple regression in the overlap region measures the relative calibration.
"""

from __future__ import annotations

import numpy as np


def calibration_stats(
    phi_data,
    ref_data,
    mu,
    mu_min: float = 0.4,
    min_abs_ref: float = 10.0,
) -> dict[str, float]:
    """Through-origin regression slope of PHI against a reference map (both
    on the same grid), restricted to well-observed pixels: finite in both
    maps, mu >= mu_min, and |ref| >= min_abs_ref (to avoid the noise floor
    dominating the fit).

    Returns slope (PHI ~= slope * ref), Pearson r, and the pixel count.
    A slope near 1 means the two instruments agree in scale.
    """
    mask = (
        np.isfinite(phi_data)
        & np.isfinite(ref_data)
        & np.isfinite(mu)
        & (mu >= mu_min)
        & (np.abs(ref_data) >= min_abs_ref)
    )
    n = int(np.count_nonzero(mask))
    if n < 10:
        return {"slope": np.nan, "pearson_r": np.nan, "n_pixels": n}

    phi = phi_data[mask].astype(float)
    ref = ref_data[mask].astype(float)

    slope = float(np.sum(phi * ref) / np.sum(ref**2))
    r = float(np.corrcoef(phi, ref)[0, 1])
    return {"slope": slope, "pearson_r": r, "n_pixels": n}
