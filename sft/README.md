# SFT model provenance

`original_transp.py` is the author's 1D surface flux transport code,
committed verbatim as uploaded (Python 2). It is the reference
implementation; the version used by this package is the faithful port in
`solar_pipeline/sft.py`.

## Port changes (physics unchanged)

- Python 3: `print` statement, integer divisions in slice indices
  (`(N-1)/2+1` etc.), `zip()` materialization for plotting, `np.trapz` ->
  `np.trapezoid`.
- Bug fix: `bjoy=loat(sys.argv[6])` -> `float(...)` (line 180 of the
  original; it would raise `NameError` on the first time step with
  `tc > 0`).
- The stochastic cycle-amplitude draw (`gauss(0, 0.13)` in log10) is kept
  but exposed as `sigma`/`seed` parameters with a deterministic default
  (`sigma=0`), so data-driven comparison runs are reproducible.
- The CLI/1000-cycle parameter-study driver, per-cycle plotting, and
  `params*.dat` bookkeeping are not ported; `SFTModel.run()` exposes the
  stepper directly and `scripts/run_sft_from_maps.py` is the data-driven
  driver.

## Correspondence checks

- `axial_dipole_moment` (0.75 * integral B sin(2 theta) dtheta) is identical
  to the original's `dipmom` and to the axisymmetric g10 of
  `solar_pipeline/dipole.py`.
- `wso_polar_field` reproduces the original's WSO proxy.
- The discretization (centered advective differencing, conservative
  diffusive term, forward Euler, vanishing-third-derivative polar BCs, the
  five meridional-flow profiles, the flux-balanced Hathaway/Jiang source)
  matches the original line for line; `tests/test_sft.py` additionally
  verifies the l=1 diffusive decay rate against the analytic 2*eta/R^2.
