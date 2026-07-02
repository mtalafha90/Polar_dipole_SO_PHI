# Paper outline (working draft skeleton)

**Title:** Impact of Solar Orbiter Polar-View Magnetic Constraints on
Surface Flux Transport Estimates of the Sun's Axial Dipole Moment

All numbers below are measured from the CR 2264 run (24 PHI-FDT cases,
2022-10-27 to 2022-11-03, matched to hmi.M_720s within 600 s; pipeline
defaults R_INNER=0.70, R_OUTER=0.90, MU_MIN=0.40, ALPHA=0.80, 1x1 deg
Carrington grid, per-case cross-calibration applied, quiet-Sun threshold
50 G). Placeholders marked [2025] await the high-latitude campaign.

## 1. Introduction

- Polar fields as the strongest precursor of the next cycle; axial dipole
  as the SFT-transportable summary quantity.
- Near-ecliptic magnetographs never see the poles at mu > ~0.45; polar
  filling is an assumption, not a measurement.
- Solar Orbiter/PHI provides the first out-of-ecliptic magnetograph;
  question: does adding its vantage measurably improve g10 estimates and
  SFT predictions?

## 2. Observations and data preparation

- 24 SO/PHI-FDT blos L2 files (two observing programs; the Oct 30+
  program ships headers without CUNIT1/2 — patched on load) + matched
  hmi.M_720s magnetograms; time matching Delta-t <= 3 s for 18/24 cases,
  max 189 s.
- Viewing geometry: SolO B0 = +7.09..+7.44 deg vs Earth B0 = +4.72..+4.88
  deg; SolO trails Earth by 40-45 deg in Carrington longitude.
- Per-case cross-calibration (through-origin regression, |B| > 10 G,
  mu >= 0.4): slope declines monotonically 0.751 -> 0.496 over 8 days at
  stable r ~ 0.75; a distinct late-UT (21-22 h) cluster sits at
  slope 0.44-0.46, r ~ 0.60 (observing-program signature).
  [Figure: calibration_drift.png]
- Orientation handling: HMI CROTA2 ~ 180 deg honored in the
  native-geometry path; ignoring it mirrors the map between hemispheres
  and flips the sign of g10 (quantified: -2.11 -> +2.20 G in project
  mode on the two-day subset) — a cautionary systematic for any
  direct-geometry synoptic construction.

## 3. Construction of merged magnetic maps

- Three products on a common 1x1 deg Carrington grid, CM-weighted
  (cos(cmd)) multi-case assimilation:
  - PHI-only (SolO vantage, native geometry), fill 40.5%, north >60 deg
    cap fill 22.0%;
  - HMI-only (Earth vantage, native geometry — deliberately NOT
    reprojected through the PHI grid), fill 43.3%, north cap 17.6%,
    south cap 2.7%;
  - merged (smooth radial blend on the PHI grid), fill 40.5%.
- Map-space consistency on common support (33.7% of sphere):
  r(phi,hmi) = 0.889 over 21,850 bins; r(phi,merged) = 0.976.
- Per-bin max-mu confidence grids accompany every product.

## 4. Axial dipole moment

- g10 = (3/4pi) integral(Br sin(lat) dOmega) under three polar-fill
  assumptions (zero / least-squares dipole projection / polar cap
  extension), with N/S decomposition.
- Full-support results (G): PHI +0.276 / +0.695 / +0.959
  (zero/project/polar_extend); HMI +0.063 / +0.147 / +0.432; merged
  +0.137 / +0.346 / +0.528.
- Common-support results: PHI -0.178 vs HMI -0.231 (zero mode) — the two
  vantages agree to ~0.05 G where they see the same bins, after
  calibration and orientation correction.
- Key argument: the positive (true-sign) dipole signal is carried by the
  coverage UNIQUE to each vantage — dominantly the high-latitude bins —
  which the common support strips out. Vantage diversity is not
  redundancy; it is signal.
- Reference check: g10 vs hmi.synoptic_mr_polfil_720s[2264]
  [run compare_reference_dipole.py --car-rot 2264; insert table].
- Uncertainties: fill-mode spread per product; quiet-Sun (|Br| <= 50 G)
  variants; ALPHA sensitivity sweep [insert from alpha_sweep on real
  data].

## 5. SFT modeling setup

- 1D flux-transport model (author's code; W = R sin(theta) B
  formulation, van Ballegooijen flow u0 = 11 m/s, eta = 250 km^2/s;
  port validated against the analytic l=1 diffusive decay rate 2*eta/R^2
  to <2%).
- Injection: zonally averaged B(lat) per product; unobserved latitudes
  enter as B = 0 (the polar-constraint handicap is carried into the
  simulation); net flux balanced before injection.

## 6. Results

- Decay experiment (source off, 11 yr): PHI-informed initial condition
  injects capN(>70 deg) = +1.80 G vs +0.22 G for HMI-only; asymptotic
  dipole +3.29 G (PHI) vs +1.06 G (HMI-only) vs +3.13 G (merged) — a
  factor ~3 from the added vantage.
- Cycle experiment (source on, tau = 10 yr, 22 yr, deterministic
  amplitude): reversal times differ by [insert from sft_reversals.csv on
  real data] between products.
- [2025] High-B0 window: repeat both experiments where the polar
  visibility advantage is ~17 deg instead of 2.5 deg.

## 7. Uncertainties and limitations

- Br ~ Blos / mu^alpha on AR fields: the quiet-Sun restriction shrinks
  the cross-vantage project-mode gap substantially [quote quiet vs full
  rows]; vector inversions would remove this class of error.
- Calibration drift (Sec. 2) and its late-UT program cluster: cause
  unresolved (candidates: SolO distance / resolution mismatch, FDT
  program differences); slope applied per case, guarded at |r| >= 0.5.
- Coverage: 8 days ~ 105 deg of longitude drift; full-CR coverage per
  vantage awaits denser PHI synoptic programs.
- Fill-mode dependence remains the largest single uncertainty on
  absolute g10; differential (product-vs-product) statements are robust
  against it.
- 1D (axisymmetric) SFT; longitude-dependent assimilation is future
  work.

## 8. Conclusions

- Merged PHI+HMI synoptic maps are constructible and validated
  (map-space r ~ 0.9 between independent vantages; ~0.05 G common-support
  dipole agreement).
- The dipole information added by a second vantage resides in its unique
  (high-latitude) coverage; SFT propagates a factor ~3 difference in
  long-term dipole from it even at a 2.5 deg B0 advantage.
- The methodology, its systematics (orientation, calibration,
  deprojection, polar filling), and open-source tooling are established;
  the decisive polar test is the 2025 high-inclination window.
