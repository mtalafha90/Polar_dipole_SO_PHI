"""1D surface flux transport (SFT) model.

Python 3 port of the author's `sft/original_transp.py` (see sft/README.md
for provenance and the list of port changes). The physics and
discretization are preserved exactly:

- annular flux density formulation W = Rsun * sin(theta) * B
- dW/dt = d(Wflux)/dx - W/tau + S * Rsun * sin(theta)
- centered differencing for the advective term, conservative centered
  differencing for the diffusive term, forward-Euler time stepping
- polar boundary conditions assuming vanishing third derivative
- the five meridional-flow profiles and the Hathaway (1994) /
  Jiang et al. (2011) source with the author's Oct 2019 tilt treatment

Units: Megameter, day, Gauss (as in the original).
"""

from __future__ import annotations

import numpy as np
import scipy.special as special

RSUN_MM = 695.7  # solar radius in Mm
CYCLE_DAYS = 11.0 * 365.25

# unit conversions from the original CLI conventions
MS_TO_MM_PER_DAY = 8.64e-2      # m/s -> Mm/day
KM2S_TO_MM2_PER_DAY = 8.64e-2   # km^2/s -> Mm^2/day


def meridional_flow(flowtype: int, u0_ms: float, latitude: np.ndarray, theta: np.ndarray) -> np.ndarray:
    """The five flow profiles of the original code. u0 in m/s; returns Mm/day."""
    u0 = u0_ms * MS_TO_MM_PER_DAY
    if flowtype == 1:  # Dikpati et al. 2006 / Cameron et al. 2007
        latitude0 = 90.0
        uc = u0 * np.sin(np.pi * latitude / latitude0)
        uc[np.abs(latitude) > latitude0] = 0.0
    elif flowtype == 2:  # van Ballegooijen 1998 / Jiang et al. 2014
        latitude0 = 75.0
        uc = u0 * np.sin(np.pi * latitude / latitude0)
        uc[np.abs(latitude) > latitude0] = 0.0
    elif flowtype == 3:  # Lemerle et al. 2017
        latitude0 = 89.0
        V, W_, q = 7.0, 1.0, 1
        uc = (
            u0
            * (special.erf(V * np.cos(np.pi / 2 * latitude / latitude0))) ** q
            * special.erf(W_ * np.sin(np.pi / 2 * latitude / latitude0))
        )
        uc[np.abs(latitude) > latitude0] = 0.0
    elif flowtype == 4:  # Whitbread et al. 2017
        p = 3.24
        uc = u0 * (np.sin(theta) ** p) * np.cos(theta)
    elif flowtype == 5:  # Wang 2017
        u0 = 13.0 * MS_TO_MM_PER_DAY
        uc = u0 * np.tanh(np.pi / 2 * latitude / 6.0) * (np.cos(np.pi / 2 * latitude)) ** 2
    else:
        raise ValueError(f"flowtype must be 1..5, got {flowtype}")
    return uc


class HathawayJiangSource:
    """Idealized cycle source: Hathaway (1994) time profile, Jiang et al.
    (2011) latitudinal profile, flux-balanced rings with the author's tilt
    treatment. Port notes: the original drew a random cycle amplitude each
    cycle (gauss(0, 0.13) in log10); here `sigma=0` makes it deterministic,
    and the original's `bjoy=loat(...)` typo is fixed to `float`.
    """

    def __init__(self, tau_days: float, blat: float = 0.0, bjoy: float = 0.0,
                 sigma: float = 0.0, seed: int | None = None,
                 sourcescale2: float = 0.015):
        self.tau_days = tau_days
        self.blat = blat
        self.bjoy = bjoy
        self.sigma = sigma
        self.rng = np.random.default_rng(seed)
        self.sourcescale2 = sourcescale2
        self.sourcescale1 = 0.0015 * np.exp(7.0 / tau_days * 365.25)
        self.sourcescale = self.sourcescale1

    def __call__(self, latitude: np.ndarray, t: float) -> np.ndarray:
        tc = 12.0 * (((t / CYCLE_DAYS) % 1) * CYCLE_DAYS / 365.25)
        ahat, bhat, chat = 0.00185, 48.7, 0.71

        if tc < 0.032:  # start of a cycle: draw (or fix) this cycle's amplitude
            gaussian = self.rng.normal(0.0, self.sigma) if self.sigma > 0 else 0.0
            self.sourcescale = self.sourcescale1 * 10**gaussian

        ampli = self.sourcescale * (ahat * tc**3 / (np.exp(tc**2 / bhat**2) - chat))

        cycleno = int(t // CYCLE_DAYS) + 1
        evenodd = 1 - 2 * (cycleno % 2)

        phase = (t / CYCLE_DAYS) % 1
        lambdan = 14.6 + self.blat * (self.sourcescale - self.sourcescale1) / self.sourcescale1
        lambdai = 26.4 - 34.2 * phase + 16.1 * phase**2
        lambda0 = lambdai * (lambdan / 14.6)
        fwhm = (0.14 + 1.05 * phase - 0.78 * phase**2) * lambdai

        joynorm0 = 1.5
        joynorm = joynorm0
        if tc > 0:
            ampli0 = self.sourcescale1 * ahat * tc**3 / (np.exp(tc**2 / bhat**2) - chat)
            joynorm = joynorm0 * (1 - self.bjoy * ((ampli - ampli0) / ampli0))

        shift = joynorm * np.sin(lambda0 / 180 * np.pi)
        bandn1 = evenodd * ampli * np.exp(-((latitude - lambda0 - shift) ** 2) / 2 / fwhm**2)
        bandn2a = -evenodd * ampli * np.exp(-((latitude - lambda0 + shift) ** 2) / 2 / fwhm**2)
        bands2a = evenodd * ampli * np.exp(-((latitude + lambda0 - shift) ** 2) / 2 / fwhm**2)
        bands1 = -evenodd * ampli * np.exp(-((latitude + lambda0 + shift) ** 2) / 2 / fwhm**2)

        # flux correction on a fine grid so the net flux is zero
        nfine = 181
        thetaf = np.linspace(0, np.pi, nfine)
        latitudef = 90.0 - thetaf * 180 / np.pi
        bandn1f = evenodd * ampli * np.exp(-((latitudef - lambda0 - shift) ** 2) / 2 / fwhm**2)
        bandn2af = -evenodd * ampli * np.exp(-((latitudef - lambda0 + shift) ** 2) / 2 / fwhm**2)
        fluxband1 = np.trapezoid(-np.sin(thetaf) * bandn1f, thetaf)
        fluxband2 = np.trapezoid(-np.sin(thetaf) * bandn2af, thetaf)
        fluxcorrection = -fluxband1 / fluxband2 if ampli != 0 else 1.0

        return bandn1 + fluxcorrection * bandn2a + bands1 + fluxcorrection * bands2a


class SFTModel:
    """The original model's grid and stepper wrapped as a class.

    Parameters follow the original CLI: u0 in m/s, eta in km^2/s,
    tau in years (None or np.inf disables the decay term).
    """

    def __init__(
        self,
        flowtype: int = 2,
        u0: float = 11.0,
        eta: float = 250.0,
        tau_years: float | None = None,
        n: int = 181,
        dt_days: float = 1.0,
    ):
        self.n = n
        self.theta = np.linspace(0, np.pi, n)
        self.latitude = 90.0 - self.theta * 180 / np.pi
        self.dx = np.pi / (n - 1)
        self.dt = dt_days

        self.uc = meridional_flow(flowtype, u0, self.latitude, self.theta)
        self.eta = eta * KM2S_TO_MM2_PER_DAY
        self.tau_days = np.inf if tau_years is None or not np.isfinite(tau_years) else tau_years * 365.25

        # stability guards for the explicit scheme
        diff_limit = self.dx**2 * RSUN_MM**2 / (2.0 * self.eta) if self.eta > 0 else np.inf
        adv_limit = self.dx * RSUN_MM / np.max(np.abs(self.uc)) if np.max(np.abs(self.uc)) > 0 else np.inf
        if self.dt > min(diff_limit, adv_limit):
            raise ValueError(
                f"dt={self.dt} d unstable: diffusion limit {diff_limit:.2f} d, advection limit {adv_limit:.2f} d"
            )

    def step(self, b, w, t, source_fn=None):
        x, dx = self.theta, self.dx
        wfladv = w / RSUN_MM * self.uc
        dwflux = (np.roll(wfladv, -1) - np.roll(wfladv, 1)) / 2
        wfldifr = np.sin(x + dx / 2) * (np.roll(b, -1) - b) / dx
        wfldifl = np.sin(x - dx / 2) * (b - np.roll(b, 1)) / dx
        dwflux += self.eta / RSUN_MM * (wfldifr - wfldifl)

        dw = dwflux / dx
        if np.isfinite(self.tau_days):
            dw -= w / self.tau_days
        if source_fn is not None:
            dw += source_fn(self.latitude, t) * RSUN_MM * np.sin(self.theta)

        w = w + dw * self.dt
        b = b.copy()
        b[1:-1] = w[1:-1] / RSUN_MM / np.sin(x[1:-1])
        # original polar BCs: vanishing third derivative
        b[0] = b[2] + 0.5 * (b[1] - b[3])
        b[-1] = b[-3] + 0.5 * (b[-2] - b[-4])
        return b, w

    def run(self, b_init, years: float, source_fn=None, record_every_days: float = 27.0):
        """Evolve an initial B(latitude) profile; returns (times_yr, B_history)."""
        b = np.asarray(b_init, dtype=float).copy()
        if b.shape != (self.n,):
            raise ValueError(f"b_init must have shape ({self.n},)")
        w = RSUN_MM * np.sin(self.theta) * b

        nt = int(round(years * 365.25 / self.dt))
        rec = max(1, int(round(record_every_days / self.dt)))
        times, history = [0.0], [b.copy()]
        for i in range(1, nt + 1):
            t = i * self.dt
            b, w = self.step(b, w, t, source_fn=source_fn)
            if i % rec == 0 or i == nt:
                times.append(t / 365.25)
                history.append(b.copy())
        return np.array(times), np.array(history)


# --------------------- diagnostics (original definitions) ---------------------

def axial_dipole_moment(b, theta):
    """0.75 * integral(B sin(2 theta) dtheta) — identical to the axisymmetric
    g10 of solar_pipeline.dipole (1.5 * integral(B sin(lat) cos(lat) dlat))."""
    return float(0.75 * np.trapezoid(b * np.sin(2 * theta), theta))


def wso_polar_field(b, theta, thetadeg_max: int = 35):
    """The original's WSO polar field proxy over the first 35 grid points."""
    thetavar = theta[: thetadeg_max + 1]
    integrand = (b * np.sin(theta) ** 2)[: thetadeg_max + 1]
    return float(5.0 * np.trapezoid(integrand, thetavar) / 1.8)


def polar_cap_mean(b, latitude, cap_deg: float = 70.0):
    """Area-weighted mean field in each polar cap."""
    out = {}
    for name, sel in (("north", latitude > cap_deg), ("south", latitude < -cap_deg)):
        wgt = np.cos(np.deg2rad(latitude[sel]))
        out[name] = float(np.sum(b[sel] * wgt) / np.sum(wgt))
    return out


def reversal_times(times, dipoles):
    """Zero crossings of the dipole series (any direction)."""
    d = np.asarray(dipoles)
    t = np.asarray(times)
    idx = np.where(np.diff(np.sign(d)) != 0)[0]
    return [float(t[i] - d[i] * (t[i + 1] - t[i]) / (d[i + 1] - d[i])) for i in idx]


def balance_flux(b, theta):
    """Remove the net signed flux from a B(latitude) profile by subtracting
    the area-weighted mean. A partial-coverage map generally injects
    hemispherically unbalanced flux, which the SFT then transports into
    unphysically large polar asymmetries; balancing is standard practice
    before injection."""
    mean = np.trapezoid(b * np.sin(theta), theta) / np.trapezoid(np.sin(theta), theta)
    return b - mean


# --------------------- map -> initial condition injection ---------------------

def zonal_profile_from_map(grid, lat_centers, sft_latitude, unobserved: str = "zero"):
    """Zonally average a Carrington Br grid into B(latitude) on the SFT grid.

    `unobserved` controls latitudes with no observed bins at all:
    - "zero": B = 0 there (i.e. no polar constraint — this is exactly the
      handicap an Earth-view-only product carries into the SFT)
    - "extend": hold the value of the nearest observed latitude
    """
    if unobserved not in ("zero", "extend"):
        raise ValueError("unobserved must be 'zero' or 'extend'")

    finite = np.isfinite(grid)
    counts = finite.sum(axis=1)
    sums = np.where(finite, grid, 0.0).sum(axis=1)
    zonal = np.divide(sums, counts, out=np.full(len(counts), np.nan), where=counts > 0)
    observed = np.isfinite(zonal)
    if not np.any(observed):
        raise RuntimeError("Map has no observed bins.")

    lat_deg = np.rad2deg(lat_centers)
    if unobserved == "extend":
        zonal = np.interp(lat_deg, lat_deg[observed], zonal[observed])
    else:
        zonal = np.where(observed, zonal, 0.0)

    # SFT latitude runs +90 -> -90; np.interp needs ascending x
    return np.interp(sft_latitude, lat_deg, zonal)
