# Research plan → code map

**Project:** Impact of Solar Orbiter Polar-View Magnetic Constraints on
Surface Flux Transport Estimates of the Sun's Axial Dipole Moment

**Science question:** does adding Solar Orbiter polar-view magnetic
information improve the estimation of the Sun's axial dipole moment compared
with Earth-view magnetogram products alone?

**Hypothesis:** near-polar magnetic flux is not fully captured by
near-ecliptic observations, especially during reversal and polar-field
rebuilding; Solar Orbiter's higher-latitude vantage should reduce this
uncertainty and yield a more reliable dipole estimate in SFT simulations.

## Stage map

| Stage | Requirement | Where it lives | Status |
|---|---|---|---|
| A | Collect overlapping HMI + SO/PHI intervals | `scripts/download_baseline_data.py` | done (needs archive network access) |
| A | Remap to common Carrington grid | `solar_pipeline/carrington.py` | done |
| A | Resolution/projection handling | `compute_case_fields` (reprojection), `compute_native_disk_fields` | done |
| A | Sign/calibration differences | `solar_pipeline/calibration.py`, per-case stats in all outputs | done (diagnostic + optional `--calibrate-phi`) |
| A | Polar confidence mask | `bin_max_to_carrington` max-mu quality grids | done |
| B | HMI-only map | `run_milestone_comparison.py` — **native Earth-view geometry**, not HMI-on-PHI | done |
| B | PHI-only map | `run_milestone_comparison.py` | done |
| B | Merged HMI+PHI map | blend on PHI grid + CM-weighted assimilation | done |
| C | Axial dipole (standard integral) | `solar_pipeline/dipole.py` (`g10`, three polar-fill modes) | done |
| C | North/south contributions | `axial_dipole_g10` hemispheric split | done |
| C | Polar-filling sensitivity | compare `g10_zero` / `g10_project` / `g10_polar_extend` columns | done |
| C | Time evolution | per-case dipole series (`run_baseline_pipeline.py`) | done (proxy); extend to g10 per case when needed |
| D | SFT baseline vs merged input | author's 1D SFT (`sft/original_transp.py`, ported to `solar_pipeline/sft.py`) + `scripts/run_sft_from_maps.py` | done (initial-condition injection) |
| D | Dipole buildup / reversal timing / polar-cap flux / asymmetry | `sft_comparison.csv`: g10(t), polar-cap means N/S, reversal detection per product | done |
| D | Continuous data assimilation into the SFT (beyond IC injection) | future work | open |

## First milestone (defined in the plan)

> Produce one Carrington map from HMI and one co-temporal reprojected map
> from PHI, then compute and compare the axial dipole moment from both and
> from a simple merged version.

One command once data is in `PHI/` and `HMI/`:

```bash
python scripts/run_milestone_comparison.py
```

Synthetic validation (`tests/test_milestone.py`, ground truth g10 = 5 G,
Earth-like B0=0° vs SolO-like B0=30° observers) confirms:

- the projection estimator recovers the true g10 from either vantage
  (error < 0.1% on synthetic data; < 3% through the full script including
  reprojection and merging),
- zero-fill underestimates g10 by ~75% on single-vantage maps — the size of
  the bias the polar data is supposed to remove,
- the B0=30° vantage fills >30% of the >60° north polar cap where the
  Earth-view map has 0%,
- the calibration regression recovers an imposed scale factor to <1%.

## Design decisions taken

- **HMI-only baseline uses HMI's native disk geometry.** Reprojecting HMI
  onto the PHI grid before binning would give the Earth-view product PHI's
  viewing geometry and make the vantage comparison circular.
- **Polar filling is a first-class parameter** (`zero`, `project`,
  `polar_extend`), because "sensitivity to polar filling assumptions" is one
  of the paper's headline metrics.
- **Calibration is reported per case** (through-origin slope, Pearson r) and
  only applied when explicitly requested (`--calibrate-phi`). Note that when
  the two observers are at genuinely different vantages, LOS projection
  differences contribute to the slope alongside instrument scale.
- **Merged maps carry a confidence grid** (per-bin max mu and accumulated
  CM weight) so SFT experiments can weight or mask low-confidence polar bins.

## Next campaign: high-latitude window (2025+)

The CR 2264 subset had SolO only ~2.5 deg above Earth's B0, so the
polar-visibility advantage was small; the measured gains came mostly from
longitude coverage and the +7 deg vantage. Since early 2025 Solar
Orbiter's orbital inclination reaches ~17 deg, and the high-latitude
windows are where the project's core question gets its real test: at
B0 = 17 deg the pole itself is visible at mu = sin(17) ~ 0.29 and the
>60 deg cap fill should roughly double relative to CR 2264.

The tooling is date-generic; a high-latitude campaign is only a data
selection. Check SOAR for phi-fdt-blos availability inside a high-B0
window (the PHI team publishes observing plans), then:

```bash
python scripts/download_baseline_data.py --start <YYYY-MM-DD> --end <YYYY-MM-DD>
python scripts/run_baseline_pipeline.py    --dates <YYYYMMDD-YYYYMMDD>
python scripts/run_milestone_comparison.py --dates <YYYYMMDD-YYYYMMDD> --calibrate-phi --quiet-sun-max-g 50
python scripts/plot_calibration_drift.py
python scripts/compare_reference_dipole.py --car-rot <CR>
python scripts/run_sft_from_maps.py --balance-flux
python scripts/run_sft_from_maps.py --balance-flux --source on --tau 10 --years 22
```

Numbers to compare directly against CR 2264: polar_fill_north/south,
the full-vs-common g10 split, capN at injection, and the SFT asymptotic
dipole and reversal times per product.

## Paper structure (target)

1. Introduction
2. Observations and data preparation (Stage A tools)
3. Construction of merged magnetic maps (Stage B tools)
4. Axial dipole moment calculation (`dipole.py`, milestone outputs)
5. SFT modeling setup (user's SFT code + injection interface)
6. Results
7. Uncertainties and limitations (alpha sweep, calibration stats, fill-mode spread, confidence masks)
8. Conclusions
