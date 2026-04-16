from __future__ import annotations

import numpy as np


def los_to_br(blos, mu, mu_min: float, alpha: float):
    br = np.full_like(blos, np.nan, dtype=float)
    valid = np.isfinite(blos) & np.isfinite(mu) & (mu >= mu_min)
    br[valid] = blos[valid] / (mu[valid] ** alpha)
    return br, valid