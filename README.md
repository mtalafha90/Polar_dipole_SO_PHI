# SO/PHI + SDO/HMI Magnetogram Merging and Axial Dipole Pipeline

This package builds merged **Solar Orbiter/PHI-FDT** + **SDO/HMI**
line-of-sight magnetogram products, estimates the Sun's **axial dipole
moment (g₁₀)** from them, and drives a **1D surface flux transport (SFT)**
model from the resulting maps. Its purpose is to test whether adding Solar
Orbiter's out-of-ecliptic magnetic vantage improves axial-dipole estimates
and SFT predictions compared with Earth-view (HMI) data alone.

The full research plan is in [`docs/RESEARCH_PLAN.md`](docs/RESEARCH_PLAN.md);
the paper skeleton with measured results is in
[`docs/paper_outline.md`](docs/paper_outline.md).

The pipeline provides:

- PHI/HMI time matching and WCS-based reprojection
- smooth radial blending of PHI and HMI into a merged LOS map
- approximate LOS-to-radial conversion (`Br ≈ Blos / mu^alpha`)
- Carrington-style latitude–longitude binning with central-meridian weighting
- **three synoptic products** (PHI-only, HMI-only native, merged) with the
  standard axial dipole g₁₀ under selectable polar-filling assumptions
- PHI-vs-HMI cross-calibration and drift diagnostics
- comparison against standard HMI synoptic charts
- an SFT model (author's code, Python 3 port) driven by the merged maps
- CSV summaries, diagnostic plots, and a synthetic-data test suite

---

## 1. Baseline method and locked defaults

The per-case baseline parameters are locked as **baseline v1** (see §10):

- **Smooth blend:** `R_INNER = 0.70`, `R_OUTER = 0.90`
- **Radial-field approximation:** `Br ≈ Blos / mu^0.8`
- **Threshold:** `MU_MIN = 0.40`
- **Disk mask:** `DISK_FRACTION = 0.98`
- **Grid:** `NLAT = 180`, `NLON = 360`
- **PHI/HMI matching tolerance:** `MAX_TIME_DIFF_SEC = 600`
- **Default subset:** `20221027`, `20221028`

These settings were chosen after sensitivity tests showed that:

- the dipole estimate is **more sensitive to `MU_MIN` than to `alpha`**,
- `MU_MIN = 0.40` is a practical compromise between limb amplification and over-clipping,
- `alpha = 0.80` is a conservative smoother alternative to the standard `1/mu` correction.

All scripts accept flag overrides, so experiments never require editing
`baseline_config.py`.

---

## 2. Folder structure

```text
package/
  README.md
  pyproject.toml
  baseline_config.py
  solar_pipeline/
    __init__.py
    io_utils.py          # loading, time parsing, date-range spec, synoptic FITS export
    geometry.py          # disk geometry, mu/lat/lon, CROTA2 / PC-matrix rotation
    blending.py          # cosine-taper PHI+HMI radial blend
    radial.py            # LOS -> Br approximation
    calibration.py       # PHI-vs-HMI cross-calibration regression
    carrington.py        # Carrington binning, CM weighting, multi-case assimilation
    dipole.py            # standard g10 (fill modes, N/S split), CEA g10
    pipeline.py          # per-case orchestration; native-geometry path
    plotting.py          # baseline + Carrington-map plots
    sft.py               # 1D surface flux transport model (Python 3 port)
  scripts/
    __init__.py
    download_baseline_data.py     # fetch PHI (SOAR) + HMI (JSOC) inputs
    run_baseline_pipeline.py      # per-case dipole proxy table
    plot_baseline_summary.py      # regenerate baseline plots
    alpha_sensitivity_sweep.py    # sweep ALPHA
    run_synoptic_assimilation.py  # single merged synoptic map
    run_milestone_comparison.py   # PHI / HMI-native / merged + g10 comparison
    compare_reference_dipole.py   # g10 vs a standard HMI synoptic chart
    plot_calibration_drift.py     # calibration slope/r vs time, distance, separation
    run_sft_from_maps.py          # SFT experiment driven by the maps
  sft/
    original_transp.py   # author's original SFT code (verbatim, Python 2)
    README.md            # provenance + port notes
  tests/                 # pytest suite (synthetic-data validation)
  docs/
    RESEARCH_PLAN.md
    paper_outline.md
  PHI/
    solo_L2_phi-fdt-blos_*.fits
  HMI/
    hmi.M_720s.*.magnetogram.fits
  baseline_outputs/
```

---

## 3. Installation

Create and activate a virtual environment, then install the package (this
also provides the console commands listed below, usable from any directory):

```bash
python3 -m venv sopyhi_env
source sopyhi_env/bin/activate
python -m pip install --upgrade pip
pip install -e .            # core
pip install -e ".[download]"  # + SOAR/JSOC download support (sunpy-soar, drms, requests)
pip install -e ".[dev]"       # + pytest
```

Console commands installed by `pip install -e .`:
`run-baseline-pipeline`, `plot-baseline-summary`, `download-baseline-data`,
`run-milestone-comparison`, `run-sft-from-maps`. The remaining scripts
(`alpha_sensitivity_sweep`, `run_synoptic_assimilation`,
`compare_reference_dipole`, `plot_calibration_drift`) are run as
`python scripts/<name>.py`.

Without installing the package:

```bash
pip install "sunpy[map]" reproject astropy matplotlib numpy scipy pandas
```

Requires Python ≥ 3.10. Under sunpy ≥ 8, PHI-FDT files whose L2 headers lack
`CUNIT1/2` are repaired automatically on load.

---

## 3a. Downloading the input data

Input magnetograms come from the mission archives:

- **PHI**: `solo_L2_phi-fdt-blos_*.fits` from the ESA Solar Orbiter Archive
  (SOAR, `soar.esac.esa.int`)
- **HMI**: `hmi.M_720s.*.magnetogram.fits` from JSOC (`jsoc.stanford.edu`)

```bash
pip install -e ".[download]"
python scripts/download_baseline_data.py --start 2022-10-27 --end 2022-10-29
```

It downloads all PHI-FDT blos L2 files in the date range from SOAR, then for
**each** PHI file queries JSOC for the nearest `hmi.M_720s` record and
downloads only those magnetograms (not the full 12-minute-cadence series).
HMI records further than `MAX_TIME_DIFF_SEC` from any PHI time produce a
warning, matching the pipeline's own rejection rule.

Useful flags: `--start/--end`, `--phi-only` / `--hmi-only`,
`--hmi-window-min`. Downloads are idempotent (existing files kept).

JSOC records are fetched directly from SUMS over HTTP with **no
registration** where possible; the DRMS keywords (WCS, observer geometry,
CROTA2) are queried alongside and written into each file's header, since raw
SUMS segments carry no metadata of their own. Not all SUMS partitions are
web-exposed (common for **recent data**), so some records' direct paths
return 404 — the downloader falls back to a real export request when
`--jsoc-email` (or `JSOC_EXPORT_EMAIL`) is set with an address registered at
http://jsoc.stanford.edu/ajax/register_email.html, otherwise it records
those and continues. A single unfetchable record never aborts the run, and
because the download is idempotent, re-running (optionally after registering
an email) resumes from what is missing.

HMI files downloaded before header injection existed can be repaired in
place (no re-download):

```bash
python scripts/download_baseline_data.py --fix-headers
```

---

## 4. Running the per-case baseline

Always run from the **package root** (the directory containing
`baseline_config.py`, `solar_pipeline/`, and `scripts/`).

```bash
python scripts/run_baseline_pipeline.py
```

With no flags this reproduces **baseline v1** exactly. All parameters can be
overridden:

```bash
python scripts/run_baseline_pipeline.py \
  --phi-dir PHI --hmi-dir HMI --out-dir baseline_outputs \
  --dates 20221027,20221028 \
  --max-time-diff-sec 600 \
  --r-inner 0.70 --r-outer 0.90 \
  --disk-fraction 0.98 --mu-min 0.40 --alpha 0.80 \
  --nlat 180 --nlon 360
```

`--dates` (all scripts) accepts single days (`20221027`), comma lists,
inclusive ranges (`20221017-20221113`), or `all` — use ranges for
full-Carrington-rotation runs.

To regenerate plots only: `python scripts/plot_baseline_summary.py`.

Outputs land in `baseline_outputs/`: `baseline_all_cases.csv`,
`baseline_summary.csv`, `baseline_summary_notes.txt`, per-case grids/maps,
and `plots/` (dipole series, offsets, time differences).

---

## 5. The full workflow

The end-to-end analysis, in order (defaults reproduce the CR 2264 study):

```bash
# 1. per-case baseline dipole proxy table
python scripts/run_baseline_pipeline.py    --dates 20221017-20221113

# 2. three synoptic products (PHI / HMI-native / merged) + g10 comparison
python scripts/run_milestone_comparison.py --dates 20221017-20221113 \
       --calibrate-phi --quiet-sun-max-g 50

# 3. diagnostics on the milestone outputs
python scripts/plot_calibration_drift.py
python scripts/compare_reference_dipole.py --car-rot 2264

# 4. SFT experiment driven by the merged maps
python scripts/run_sft_from_maps.py --balance-flux                       # decay
python scripts/run_sft_from_maps.py --balance-flux --source on --tau 10 --years 22  # cycle
```

---

## 6. What each stage does

### 6a. Per-case pipeline (`run_baseline_pipeline.py`)

For each PHI LOS magnetogram in the subset: parse its time, find the nearest
HMI record (reject if `Δt > MAX_TIME_DIFF_SEC`), reproject HMI onto the PHI
grid, build a smooth merged LOS map by radial cosine weighting, estimate
heliographic coordinates, convert LOS to approximate radial field
(`Br ≈ Blos / mu^alpha`), bin onto a Carrington grid, and compute an axial
dipole proxy for PHI/HMI/merged.

### 6b. Optional single-tool scripts

- **Alpha sensitivity sweep** — re-runs the pipeline across `ALPHA` values
  and reports how the dipole estimates change
  (`baseline_outputs/alpha_sweep/`):

  ```bash
  python scripts/alpha_sensitivity_sweep.py --alphas 0.6,0.7,0.8,0.9,1.0
  ```

- **Synoptic assimilation** — combines cases into one CM-weighted merged
  synoptic map, each bin weighted by `cos(cmd)^n`
  (`baseline_outputs/synoptic/`):

  ```bash
  python scripts/run_synoptic_assimilation.py
  ```

### 6c. Milestone comparison (`run_milestone_comparison.py`)

Builds **three** synoptic map products and compares their axial dipole:

- **PHI-only** — SolO vantage, PHI's native disk geometry
- **HMI-only** — Earth vantage, HMI's **native** disk geometry (deliberately
  *not* reprojected through the PHI grid, so the vantage comparison is not
  circular; HMI's ~180° `CROTA2` camera rotation is honored)
- **merged** — the smooth PHI+HMI blend on the PHI grid

For each product it computes the standard axial dipole coefficient
`g10 = (3/4π) ∮ Br sin(lat) dΩ` (`solar_pipeline/dipole.py`) under three
polar-filling assumptions — `zero` (unobserved bins contribute nothing),
`project` (least-squares projection onto the dipole profile), `polar_extend`
(caps filled with the last observed band's zonal mean) — with north/south
decomposition, polar-cap fill fractions, and a per-bin max-mu confidence
grid. Per-case PHI-vs-HMI calibration slopes are reported; `--calibrate-phi`
applies them, guarded by `--calib-min-r` (default 0.5) so an uncorrelated
fit never rescales the maps.

The table also reports:

- `*_common` rows — each product's dipole on the intersection of all
  products' observed bins, isolating vantage/calibration effects from
  longitude-coverage effects;
- `*_quiet` / `*_quiet_common` rows (with `--quiet-sun-max-g <G>`) —
  strong-field bins excluded, since active regions are where the
  radial-field assumption fails and the vantages disagree most;
- a **map-space correlation** table (`map_correlations.csv`) between
  products on the common support — a decisive orientation/consistency check.

`--max-separation-deg <D>` guards the **merged** product: cases whose
SolO–Earth Carrington-longitude separation exceeds `D` are excluded from the
blend (and left uncalibrated), because per-pixel PHI+HMI blending is
meaningless once the two spacecraft look at different hemispheres — as in the
2025 high-B0 campaign, where the separation reached ~171°. The PHI and HMI
native-geometry products keep every case, so polar-coverage statistics are
unaffected. `calibration_stats.csv` records each case's `lon_separation_deg`
and both `crlt_obs` (B0) values.

Outputs (`baseline_outputs/milestone/`): `.npy` grids, SFT-ready
`synoptic_*.fits` maps (plate-carrée WCS), `milestone_dipole_comparison.csv`,
`calibration_stats.csv`, `map_correlations.csv`, and map plots.

**Per-rotation runs.** A multi-month window spans several Carrington
rotations; a single combined map smears the dipole across them and cannot be
checked against a per-rotation reference chart. `run_milestone_by_rotation.py`
splits the window by rotation (assigning each day by its noon) and runs the
milestone comparison once per CR into `<out-dir>/cr_<N>/`, forwarding any
flags after a literal `--`:

```bash
python scripts/run_milestone_by_rotation.py --dates 20250211-20250429 \
       -- --calibrate-phi --quiet-sun-max-g 50 --max-separation-deg 60
```

`plot_campaign_summary.py` then turns those per-rotation outputs into the
campaign's headline figure — SolO's B0 excursion and the PHI-vs-HMI north/
south polar-cap fill per rotation, with the separation range annotated:

```bash
python scripts/plot_campaign_summary.py --campaign-dir baseline_outputs \
       --out baseline_outputs/campaign_polar_advantage.png
```

### 6d. Calibration drift and reference checks

- `plot_calibration_drift.py` plots the PHI-vs-HMI slope and Pearson r vs
  time (colored by hour of day), vs SolO–Sun distance, and vs SolO–Earth
  longitude separation, and prints a trend/cluster summary
  (`baseline_outputs/milestone/plots/calibration_drift.png`).
- `compare_reference_dipole.py` computes g₁₀ from a standard HMI synoptic
  chart (local FITS via `--reference`, or fetched by
  `--car-rot <N>` from `hmi.synoptic_mr_polfil_720s`) and tabulates every
  milestone product/fill-mode against it
  (`reference_dipole_comparison.csv`).
- `compare_reference_by_rotation.py` runs that check for every rotation of a
  per-CR campaign — it discovers each `cr_<N>/milestone` directory, fetches
  the matching `--car-rot N` chart, and collects a combined
  `reference_dipole_by_rotation.csv` at the campaign root:

  ```bash
  python scripts/compare_reference_by_rotation.py --campaign-dir baseline_outputs
  ```

### 6e. SFT experiment (`run_sft_from_maps.py`)

The author's 1D surface flux transport model (`sft/original_transp.py`,
committed verbatim; identical-discretization Python 3 port in
`solar_pipeline/sft.py`, see `sft/README.md`) closes the loop from maps to
flux-transport predictions. Each product is zonally averaged into an initial
`B(latitude)` profile and evolved under the same SFT configuration (flow
profiles 1–5, `--u0`, `--eta`, `--tau`, optional cycle source via
`--source on`, `--balance-flux` to remove net injected flux). With the
default `--unobserved zero`, latitudes a product never observed enter with
zero field — so an Earth-view-only product carries its missing-polar-field
handicap into the simulation while merged PHI+HMI input does not, which is
precisely the experiment.

Outputs (`baseline_outputs/sft/`): dipole and polar-cap time series per
product (`sft_comparison.csv`), reversal timing (`sft_reversals.csv`),
injected profiles, and a comparison figure. The port is physics-validated in
`tests/test_sft.py`, including the analytic l=1 diffusive decay rate
`2·eta/Rsun²` recovered to <2%.

**Polar-constraint experiment (`run_sft_polar_experiment.py`).** At the
epochs where PHI has a real polar advantage the SolO–Earth separation is
large (opposite hemispheres), so a per-pixel merge is invalid — PHI is a
*standalone polar constraint*, not a merge partner. This script builds the
initial condition from HMI everywhere **except** the polar cap PHI observes,
where PHI's zonal field is spliced in (`apply_polar_constraint`, with a
`--blend-deg` ramp), and evolves three cases — `hmi` (the polar handicap),
`phi_constrained` (HMI + PHI cap), and `phi` (PHI-only) — reporting the
Δg₁₀ and Δreversal-timing the PHI constraint induces. Run it per rotation on
the pole PHI covers (`--hemisphere south` for the March rotation,
`--hemisphere north` for April):

```bash
python scripts/run_sft_polar_experiment.py --maps-dir out/cr_2296/milestone \
       --hemisphere north --source on --tau 10 --years 22 --balance-flux
```

Outputs (`<maps-dir>/sft_polar/`): `sft_polar_comparison.csv`,
`sft_polar_reversals.csv`, injected profiles, and a figure.

---

## 7. Results on the CR 2264 subset (Oct 27 – Nov 3, 2022)

Measured from 24 PHI-FDT cases matched to `hmi.M_720s` (per-case calibration
applied, quiet-Sun threshold 50 G):

- **Consistency:** map-space correlation PHI vs HMI = **0.89** on common
  support; the two vantages agree on the axial dipole to **~0.05 G** where
  they see the same bins (zero mode: PHI −0.18 vs HMI −0.23 G).
- **Reference accuracy:** against the standard polar-filled HMI synoptic
  chart (g₁₀ = **+0.65 G**), the PHI-informed product (project mode) agrees
  to **+0.04 G (6%)**, while the same-pipeline Earth-view product misses by
  ~5× more — the accuracy ordering follows polar coverage (PHI 22% north-cap
  fill vs HMI 18%).
- **SFT impact:** PHI-informed input injects a north polar-cap field of
  **+1.80 G** vs **+0.22 G** for HMI-only, and evolves to an asymptotic
  dipole of **+3.29 G vs +1.06 G** (factor ~3). With a cycle source, the
  first polar reversal shifts by **1.26 yr** (3.72 vs 2.46 yr) from the
  added polar constraint.
- **Calibration:** the PHI/HMI slope declines monotonically 0.75 → 0.50 over
  eight days, tracking SolO–Sun distance nearly linearly (consistent with
  resolution-dependent flux loss); a late-UT observing-program cluster sits
  below that trend.

Two systematics were found and fixed during this study and are worth
flagging for anyone building direct-geometry synoptic maps: HMI's ~180°
`CROTA2` (ignoring it mirrors the map between hemispheres and flips the sign
of g₁₀), and the absence of WCS metadata in raw JSOC SUMS segments (injected
on download / `--fix-headers`).

---

## 8. Current limitations

- `Br ≈ Blos / mu^alpha` is only an approximation; the cross-vantage
  disagreement concentrates in active-region bins where it fails (hence the
  quiet-Sun diagnostic). No full vector inversion is used.
- Coverage: eight days ≈ 105° of longitude drift per vantage; full-CR
  coverage awaits denser PHI synoptic programs.
- Fill-mode dependence is the largest single uncertainty on *absolute* g₁₀;
  differential (product-vs-product) statements are robust against it.
- The calibration-distance correlation is consistent-with but not yet
  confirmed (the degrade-HMI-to-PHI-resolution test is planned).
- The SFT model is 1D (axisymmetric); longitude-dependent assimilation is
  future work.

---

## 9. Status and next steps

Done and validated on real data:

- ~~command-line interface~~ (§4) and ~~formal packaging~~ (`pyproject.toml`, §3)
- ~~multi-map Carrington assimilation~~ with central-meridian weighting (§6b–c)
- ~~standard g₁₀ with polar-filling options and N/S split~~ (`dipole.py`)
- ~~non-circular HMI-native baseline, calibration, confidence masks~~ (§6c)
- ~~SFT model port and map-driven experiment~~ (§6e)
- ~~full-Carrington-rotation run + diagnostics~~ (§6d, §7)
- ~~reference comparison, reversal timing, paper outline~~ (§6d–e, `docs/`)

Remaining, data-gated rather than code-gated:

- **Degrade-HMI resolution test** — confirm the calibration-distance
  correlation.
- **2025 high-B₀ window** — the decisive polar test, where Solar Orbiter's
  ~17° inclination makes the pole itself visible (see the "Next campaign"
  section of `docs/RESEARCH_PLAN.md`).

---

## 10. Version note

The per-case **baseline v1** defaults (`R_INNER = 0.70`, `R_OUTER = 0.90`,
`MU_MIN = 0.40`, `ALPHA = 0.80`) are kept fixed as the reference
configuration. Analysis scripts layer on top via flags without changing
them, so the baseline remains reproducible while experiments proceed.
