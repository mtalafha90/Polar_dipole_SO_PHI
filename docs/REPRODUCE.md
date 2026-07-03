# Reproducing the full analysis from scratch

A clean, end-to-end run in a fresh folder. Two data windows are used — the
2022 **CR 2264** demonstration and the 2025 **CR 2294–2297** high-B₀
campaign — and both are downloaded into the same `PHI/` and `HMI/` folders;
each analysis command selects its subset with `--dates`.

Expected wall time is dominated by the 2025 download (a few hundred PHI +
matched HMI files). Everything after that is minutes.

## 0. Prerequisites

- Python ≥ 3.10 and `git`.
- For the 2025 download: a **JSOC-registered email** — recent HMI records are
  not web-exposed in SUMS and need the export fallback. Register once at
  <http://jsoc.stanford.edu/ajax/register_email.html>.

## 1. Clone and install

```bash
git clone https://github.com/mtalafha90/SO-PHI-SDO-HMI.git
cd SO-PHI-SDO-HMI
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -e ".[download,dev]"
pytest -q                     # sanity check: expect "27 passed"
```

Run every command below **from the repo root** (the folder with
`baseline_config.py`).

## 2. Download the input magnetograms

```bash
# Demonstration window (CR 2264, 2022) — small and fast
python scripts/download_baseline_data.py --start 2022-10-27 --end 2022-11-03

# 2025 high-B0 campaign (CR 2294-2297) — larger; recent data needs the email
export JSOC_EXPORT_EMAIL="you@example.com"     # your registered address
python scripts/download_baseline_data.py --start 2025-02-11 --end 2025-04-29 \
       --jsoc-email "$JSOC_EXPORT_EMAIL"
```

Normal, non-fatal messages:
- `SunpyUserWarning: ... on-board averaging of 3 individual observations` —
  informational.
- `No suitable HMI match within 600 s` — that PHI frame has no simultaneous
  HMI record; skipped by design.
- A few April records may still fail to fetch; the run continues and is
  **idempotent** (re-run to resume from what is missing).

## 3. Demonstration analysis (CR 2264, 2022)

```bash
python scripts/run_baseline_pipeline.py    --dates 20221027-20221103
python scripts/run_milestone_comparison.py --dates 20221027-20221103 \
       --calibrate-phi --quiet-sun-max-g 50
python scripts/plot_calibration_drift.py
python scripts/compare_reference_dipole.py --car-rot 2264
python scripts/alpha_sensitivity_sweep.py  --alphas 0.6,0.7,0.8,0.9,1.0
python scripts/run_sft_from_maps.py --balance-flux                                  # decay
python scripts/run_sft_from_maps.py --balance-flux --source on --tau 10 --years 22  # cycle
```

Outputs under `baseline_outputs/`: `baseline_*.csv`, `milestone/`,
`alpha_sweep/`, `sft/`.

## 4. The 2025 campaign (the decisive test)

```bash
# Per-rotation milestone with the separation guard -> baseline_outputs/cr_2294..2297/
python scripts/run_milestone_by_rotation.py --dates 20250211-20250429 \
       -- --calibrate-phi --quiet-sun-max-g 50 --max-separation-deg 60

# Headline coverage figure (auto-discovers the cr_* rotations)
python scripts/plot_campaign_summary.py --campaign-dir baseline_outputs \
       --out baseline_outputs/campaign_polar_advantage.png

# Per-rotation reference-dipole check (needs JSOC network)
python scripts/compare_reference_by_rotation.py --campaign-dir baseline_outputs

# SFT polar-constraint experiment, on the pole PHI covers each rotation
python scripts/run_sft_polar_experiment.py --maps-dir baseline_outputs/cr_2296/milestone \
       --hemisphere north --source on --tau 10 --years 22 --balance-flux
python scripts/run_sft_polar_experiment.py --maps-dir baseline_outputs/cr_2295/milestone \
       --hemisphere south --source on --tau 10 --years 22 --balance-flux

# Publication figure for the north run
python scripts/plot_sft_polar.py --sft-dir baseline_outputs/cr_2296/milestone/sft_polar
```

## 5. Robustness (Sec. 7)

```bash
# Attribute the calibration-slope deficit to resolution (co-observing Feb window)
python scripts/calibration_resolution_test.py --dates 20250214-20250301 --fwhms 0,1,2,3,5,8

# Optional: MU_MIN sensitivity of the coverage advantage (kept in a separate dir)
python scripts/run_milestone_by_rotation.py --dates 20250211-20250429 \
       --out-dir baseline_outputs_mu025 \
       -- --calibrate-phi --quiet-sun-max-g 50 --max-separation-deg 60 --mu-min 0.25
```

(The α sweep was already run in step 3.)

## 6. Where the results land

| Result | Path |
|---|---|
| Coverage figure | `baseline_outputs/campaign_polar_advantage.png` |
| Per-rotation dipole + calibration | `baseline_outputs/cr_<N>/milestone/{milestone_dipole_comparison,calibration_stats}.csv` |
| SFT polar experiment | `baseline_outputs/cr_2296/milestone/sft_polar/{sft_polar_reversals.csv, sft_polar_figure.png}` |
| Reference-dipole check | `baseline_outputs/reference_dipole_by_rotation.csv` |
| Resolution attribution | `baseline_outputs/milestone/calibration_resolution_test.{csv,png}` |
| α sweep | `baseline_outputs/alpha_sweep/alpha_sweep_summary.csv` |

Interpretation of every number is in [`manuscript.md`](manuscript.md); a
one-page summary is in [`results_brief.html`](results_brief.html).

## Expected headline numbers (as a check)

- Coverage (≥60° cap, µ_min = 0.4): south **64% / 46%** (CR 2295, PHI/HMI),
  north **51% / 3%** (CR 2296).
- SFT: Δ first reversal **−1.68 yr** (north constraint) vs **−0.004 yr**
  (south).
- Resolution: mean slope **0.56 → ~1.0** as HMI is smoothed 0 → 5 px FWHM.
- α: PHI dipole varies **~3%** over α = 0.6–1.0.
