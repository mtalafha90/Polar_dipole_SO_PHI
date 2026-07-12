# Paper outline (working draft skeleton)

**Title:** Impact of Solar Orbiter Polar-View Magnetic Constraints on
Surface Flux Transport Estimates of the Sun's Axial Dipole Moment

Numbers in Secs. 2-4 and 6.1 are measured from the CR 2264 run (24 PHI-FDT
cases, 2022-10-27 to 2022-11-03, matched to hmi.M_720s within 600 s;
pipeline defaults R_INNER=0.70, R_OUTER=0.90, MU_MIN=0.40, ALPHA=0.80, 1x1
deg Carrington grid, per-case cross-calibration applied, quiet-Sun
threshold 50 G). Sec. 6.2 reports the 2025 high-B0 campaign (CR 2294-2296),
which delivers the decisive polar test; items marked [pending run] await
SFT/reference runs on the campaign maps.

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
  stable r ~ 0.75, tracking the SolO-Sun distance nearly linearly
  (0.75 at 0.42 AU -> 0.50 at 0.52 AU) — consistent with
  resolution-dependent flux loss as PHI's plate scale coarsens.
  Time/distance/longitude-separation are collinear within one window;
  confirming test: degrade HMI to PHI's per-case resolution before the
  regression (planned). A distinct late-UT (21-22 h) cluster falls
  BELOW the distance trend at slope 0.44-0.46, r ~ 0.60
  (observing-program signature; flagged/excluded in final calibration).
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
- Reference check vs hmi.synoptic_mr_polfil_720s[2264]
  (g10_ref = +0.654 G): PHI-only/project = +0.695 (Delta = +0.041 G, 6%);
  merged/polar_extend = +0.528 (-0.126); HMI-only best (polar_extend) =
  +0.432 (-0.222). The PHI-informed product reproduces the
  community-standard full-CR dipole to 0.04 G from eight days of data;
  the same-pipeline Earth-view product misses by 5x more, and the
  accuracy ordering follows polar coverage (22.0% / annulus-diluted /
  17.6% north-cap fill). Framing: the PHI product carries HMI's
  calibration scale, so the claim is vantage-added accuracy, not
  instrument-alone; and the reference's own poles are model-filled —
  agreement with the standard SFT input, not absolute truth.
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
- Two injection modes: (i) whole-map injection per product
  (`run_sft_from_maps.py`), used when PHI and HMI co-observe so a merged
  map is meaningful; (ii) polar-constraint splice
  (`run_sft_polar_experiment.py`, `apply_polar_constraint`), used at
  high-separation epochs where a per-pixel merge is invalid — the initial
  condition is HMI everywhere except the cap PHI observes, where PHI's
  zonal field is blended in. Mode (ii) is the appropriate one for the 2025
  campaign's decisive rotations (Sec. 6.2).

## 6. Results

### 6.1 Method demonstration (CR 2264, low-B0)

- Decay experiment (source off, 11 yr): PHI-informed initial condition
  injects capN(>70 deg) = +1.80 G vs +0.22 G for HMI-only; asymptotic
  dipole +3.29 G (PHI) vs +1.06 G (HMI-only) vs +3.13 G (merged) — a
  factor ~3 from the added vantage.
- Cycle experiment (source on, tau = 10 yr, 22 yr, deterministic
  amplitude): first polar reversal at 2.46 yr from the HMI-only initial
  condition vs 3.72 yr PHI-informed (merged 3.64) — a 1.26 yr shift in
  predicted reversal timing from the vantage-added polar constraint.
  The memory decays as the source dominates: second-reversal spread
  shrinks to 0.32 yr (14.87/15.19/14.89) and final dipoles converge
  (3.41/3.17/3.40 G).
- Caveat: at a 2.5 deg B0 advantage this demonstrates the propagation
  mechanism, not a large polar-visibility gain; the decisive test is the
  high-inclination window below.

### 6.2 The decisive polar test (2025 high-B0 campaign)

Data: SO/PHI-FDT blos matched to hmi.M_720s over 2025-02-14 to 2025-07-16
(time matches Delta-t = 2-3 s). The window spans six Carrington
rotations, so it is analysed per rotation (`run_milestone_by_rotation.py`)
rather than as one smeared synoptic map, and the merged product excludes
cases beyond a 60 deg SolO-Earth Carrington-longitude separation
(`--max-separation-deg 60`), since a per-pixel PHI+HMI blend is
meaningless once the two spacecraft view different hemispheres.

The orbital driver: SolO's B0 sweeps from near-ecliptic (-2 deg) to deep
south (-16.7 deg, late March) to a high-north plateau (+16.8 deg, late
April), while Earth's B0 stays at -5 to -7 deg. The polar-cap fill
(>=60 deg, in %) tracks it directly:

| rotation | dates | SolO B0 | sep | N-cap PHI/HMI | S-cap PHI/HMI |
|---|---|---|---|---|---|
| CR 2294 | Feb 14 - Mar 1 | -2 -> -8 | 15-21 deg | 5 / 0 | 30 / 38 |
| CR 2295 | Mar 2 - Mar 28 | -8 -> -17 | 0.3-60 deg | 0 / 0 | **64 / 46** |
| CR 2296 | Mar 31 - Apr 24 | -8 -> +17 | 80-165 deg | **51 / 3** | 12 / 41 |
| CR 2297 | Apr 26 - May 21 | +14 -> +17 | 168-180 deg | **77 / 12** | 0 / 34 |
| CR 2298 | May 23 - Jun 19 | descending | far side | **64 / 23** | 0 / 24 |
| CR 2299 | Jun 20 - Jul 16 | descending | far side | 45 / 34 | 4 / 13 |

- Coverage (calibration-independent, the robust headline): PHI overtakes
  HMI on the south cap in March (64% vs 46%) as SolO dives south, and on
  the north cap in April (51% vs 3%, a 16x advantage) as SolO crosses to a
  northern view; the advantage STRENGTHENS across the full CR 2297 rotation
  (77% vs 12% on 203 cases) as SolO holds the +14-17 deg apex in near-perfect
  opposition (sep 168-180 deg, zero co-observed pixels all month), then
  SWITCHES OFF symmetrically as SolO descends: PHI north fill 77 -> 64 -> 45%
  while HMI's climbs 12 -> 23 -> 34% (Earth's B0 rising into northern
  summer), nearly closing the gap by CR 2299. The advantage is governed by
  the DIFFERENCE of the two heliolatitudes. The merged product stays empty
  through CR 2296-2299 (SolO still far-side even as the advantage fades):
  the coverage (B0-driven) and merge-validity (separation-driven) effects
  decouple. [Figure: campaign_polar_advantage.png.]
- Polarity (CR 2297): with 77% north-cap coverage PHI measures
  g10_north = -0.25 G, the SIGN of the reference chart (-0.39 G) and the
  cycle-25 north reversal; HMI's 12% grazing-angle extrapolation gives +0.19
  (wrong sign) -- the coverage advantage as a corrected polar-field polarity.
- Vantage effect on the dipole: on common support in CR 2295 the two
  vantages disagree on the SIGN of the south-polar g10 contribution
  (project mode: PHI +0.21 vs HMI -0.10 G on the same bins). HMI's
  south-cap Br is a near-limb 1/mu^alpha extrapolation there, unreliable
  near solar maximum where the polar field is weak and reversing; PHI sees
  the same bins near disk centre and gives the better-constrained term.
- Key structural finding: PHI's polar advantage is largest exactly when
  the SolO-Earth separation is largest (both are driven by the same
  orbital motion). When separation is small (Feb, and the Mar 13
  conjunction) a per-pixel merge is valid but PHI's polar advantage is
  small; when the advantage is decisive (April north cap) the separation
  is 80-165 deg and the merge is meaningless (the merged product is
  correctly empty). The two vantages are therefore complementary in TIME,
  not co-registered in space: PHI's value is as an independent polar
  constraint at the epochs Earth cannot see a given pole, not as a
  pixel-level merge partner. This reframes the deliverable and is enforced
  by the separation guard.
- SFT impact via the polar-constraint splice
  (`run_sft_polar_experiment.py`; Sec. 5 mode ii; source on, tau = 10 yr,
  22 yr, flux-balanced): the effect is concentrated exactly where the Earth
  view is blind. Constraining the NORTH cap from CR 2296 (HMI cap fill 3%)
  injects capN = -1.51 G where HMI has ~0, flips the initial dipole sign
  (g10(0) = -0.26 vs +0.10 HMI-only), and shifts the first predicted polar
  reversal by -1.68 yr (0.68 vs 2.36 yr). Constraining the SOUTH cap from
  CR 2295 (HMI cap fill 46%, co-observed) changes essentially nothing
  (Delta first-reversal -0.004 yr, Delta g10_final -0.6 mG), because HMI
  already constrains that pole. The 22-yr dipoles converge to ~+3.1 G in
  all cases as the deterministic source dominates — the memory of the polar
  constraint lives in the reversal TIMING, not the asymptote. The PHI polar
  view adds SFT-relevant information only where the Earth view lacks it, so
  its value is maximal at the high-separation north epoch — precisely the
  standalone-constraint regime.
- Per-rotation reference check against hmi.synoptic_mr_polfil_720s
  (`compare_reference_by_rotation.py`): the standard polar-filled chart's
  dipole grows steadily more negative through the window (g10_ref = -0.21,
  -0.17, -0.18, -0.39, -0.57, -0.52 G for CR 2294-2299) — the cycle-25
  axial dipole rebuilding after reversal. Where the true dipole is near
  zero (early window) the partial-coverage products miss its sign (all
  positive, +0.3 to +1.0 G; a near-zero dipole is swamped by
  longitude-sampling), so early absolute g10 is provisional. As the true
  dipole strengthens the products lock on: by CR 2299 both PHI (-0.19) and
  HMI (-0.03) are negative, and PHI is now CLOSER to the reference than HMI
  (Delta +0.37 vs +0.49) — reversing the early ordering and corroborating
  the CR 2297 polarity result: once the north polar field is strong, the
  vantage that sees it recovers more of it than the extrapolation. The
  campaign's robust results remain the polar COVERAGE and the SFT
  reversal-timing DIFFERENTIAL (both calibration- and fill-independent).
- Robustness (Sec. 7): (a) the sub-unity calibration slope is RESOLUTION —
  smoothing HMI to PHI's plate scale raises the mean slope 0.56 -> ~1.0
  (crossing unity near 5 px FWHM, r peaking at 2-3 px;
  `calibration_resolution_test.py`), not an instrument-scale or vantage
  error; (b) the coverage-advantage MAGNITUDE depends on the limb cut
  (~16x at mu_min=0.4 vs ~2.6x at mu_min=0.25, where HMI admits
  low-reliability grazing-angle polar pixels) but the SIGN is robust;
  (c) alpha is sub-dominant — over alpha=0.6-1.0 the PHI dipole varies ~3%,
  HMI ~9% (`alpha_sensitivity_sweep`), below the case scatter.
- Caveat: for CR 2296/2297 there is no common support (opposite
  hemispheres) and no valid calibration, so the April g10 magnitudes are
  coverage-driven, not a controlled vantage comparison; the 51% / 37% vs
  3% / 2% coverage result is robust, the April g10 is provisional.

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

- Merged PHI+HMI synoptic maps are constructible and validated at small
  SolO-Earth separation (map-space r ~ 0.9 between independent vantages;
  ~0.05 G common-support dipole agreement).
- The dipole information added by a second vantage resides in its unique
  (high-latitude) coverage; SFT propagates a factor ~3 difference in
  long-term dipole from it even at a 2.5 deg B0 advantage.
- The 2025 high-B0 campaign delivers the decisive test: PHI's polar-cap
  coverage overtakes HMI by up to 16x (north cap 51% vs 3% in April), and
  the advantage switches on with SolO's heliolatitude exactly as
  hypothesised — south cap in March, north cap in April.
- Central reframing: PHI's polar advantage is largest exactly when the
  SolO-Earth separation is largest, so at the decisive epochs a per-pixel
  merge is invalid. The scientific value of the polar view is as an
  INDEPENDENT polar constraint (a boundary condition for SFT) at times
  Earth cannot see a pole, not as a co-registered merge partner. The two
  vantages are complementary in time, not space.
- Quantified impact: adding PHI's north-cap constraint at the April epoch
  where HMI is blind (3% fill) shifts the first SFT-predicted polar
  reversal by -1.68 yr, while adding its south-cap constraint where HMI
  co-observes (46% fill) changes nothing (-0.004 yr). The polar view moves
  the prediction only where the Earth view is blind — the concrete payoff
  of the standalone-constraint framing.
- The methodology, its systematics (orientation, calibration,
  deprojection, polar filling, co-observation separation), and open-source
  tooling are established end to end, from download through the
  polar-constrained SFT experiment.
