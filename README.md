# PHI–HMI Baseline Pipeline

This package provides a baseline workflow for comparing and combining **Solar Orbiter/PHI-FDT** and **SDO/HMI** line-of-sight magnetograms over the **27–28 October 2022** subset, then estimating an **approximate Carrington-style axial dipole proxy**.

The current codebase is intended as a **working baseline version** for development and controlled experiments. It includes:

- PHI/HMI time matching
- WCS-based reprojection of HMI onto the PHI grid
- smooth radial blending of PHI and HMI maps
- approximate LOS-to-radial conversion
- Carrington-style latitude–longitude binning
- approximate axial dipole estimation
- CSV summaries and diagnostic plots

---

## 1. Current baseline method

The baseline parameters are locked to:

- **Smooth blend:** `R_INNER = 0.70`, `R_OUTER = 0.90`
- **Radial-field approximation:** `Br ≈ Blos / mu^0.8`
- **Threshold:** `MU_MIN = 0.40`
- **Disk mask:** `DISK_FRACTION = 0.98`
- **Grid:** `NLAT = 180`, `NLON = 360`
- **PHI/HMI matching tolerance:** `MAX_TIME_DIFF_SEC = 600`
- **Subset:** `20221027`, `20221028`

These settings were chosen after sensitivity tests showed that:

- the dipole estimate is **more sensitive to `MU_MIN` than to `alpha`**,
- `MU_MIN = 0.40` is a practical compromise between limb amplification and over-clipping,
- `alpha = 0.80` is a conservative smoother alternative to the standard `1/mu` correction.

---

## 2. Folder structure

Recommended layout:

```text
package/
  README.md
  baseline_config.py
  solar_pipeline/
    __init__.py
    io_utils.py
    geometry.py
    blending.py
    radial.py
    carrington.py
    pipeline.py
    plotting.py
  scripts/
    run_baseline_pipeline.py
    plot_baseline_summary.py
  PHI/
    solo_L2_phi-fdt-blos_*.fits
    solo_L2_phi-fdt-icnt_*.fits
  HMI/
    hmi.M_720s.*.magnetogram.fits
  baseline_outputs/
```

---

## 3. Required Python packages

Create and activate a virtual environment, then install:

```bash
python3 -m venv sopyhi_env
source sopyhi_env/bin/activate
python -m pip install --upgrade pip
pip install "sunpy[map]" reproject astropy matplotlib numpy scipy pandas
```

If SunPy metadata warnings appear but the pipeline still completes successfully, they can be ignored during development or suppressed in the run script.

---

## 4. How to run

Always run from the **package root**, i.e. the directory containing:

- `baseline_config.py`
- `solar_pipeline/`
- `scripts/`

Example:

```bash
cd /path/to/package
python scripts/run_baseline_pipeline.py
```

To regenerate plots only:

```bash
python scripts/plot_baseline_summary.py
```

---

## 5. Main outputs

The baseline pipeline writes results into:

```text
baseline_outputs/
```

Main files:

- `baseline_all_cases.csv` — full case table
- `baseline_summary.csv` — compact summary table
- `baseline_summary_notes.txt` — text summary of campaign statistics
- `plots/dipole_series.png` — PHI/HMI/merged dipole time series
- `plots/dipole_offsets.png` — merged-minus-PHI and merged-minus-HMI offsets
- `plots/time_differences.png` — PHI/HMI time offsets

Per-case folders also contain:

- merged smooth LOS maps
- binned Carrington-style grids
- latitude/longitude arrays

---

## 6. What the pipeline does

For each PHI LOS magnetogram in the selected subset:

1. parse PHI observation time
2. find nearest HMI magnetogram in time
3. reject the match if `Δt > MAX_TIME_DIFF_SEC`
4. reproject HMI onto the PHI grid
5. build a **smooth merged LOS map** using radial cosine weighting
6. estimate heliographic coordinates on the visible disk
7. convert LOS field to approximate radial field using
   `Br ≈ Blos / mu^alpha`
8. bin the visible-disk radial field onto a Carrington-style grid
9. compute an approximate axial dipole proxy
10. save results and diagnostic plots

---

## 7. Scientific interpretation of the current baseline

For the 27–28 October 2022 subset, the baseline results show:

- PHI and HMI LOS maps agree well after reprojection
- smooth blending is better behaved than a hard radial replacement
- the merged dipole is typically **intermediate between PHI and HMI**
- the merged dipole is often **closer to HMI than to PHI**
- the merged series is **more stable** than either PHI or HMI alone

This baseline should therefore be treated as a **working reference implementation**, not as the final physical model.

---

## 8. Current limitations

The main limitations are physical rather than geometric:

- `Br ≈ Blos / mu^alpha` is only an approximation
- no full vector inversion is used
- only the visible disk is used at a time
- the Carrington map is a binned visible-disk approximation, not a full synoptic assimilation product
- metadata cleanup could still be improved for completely warning-free runs

---

## 9. Recommended next steps

Good next directions include:

- improving the LOS-to-radial conversion
- assimilating multiple visible-disk maps into a fuller Carrington map
- extending the analysis beyond 27–28 Oct 2022
- adding a command-line interface for cleaner execution
- packaging the code formally for reuse

---

## 10. Version note

This README describes the current **baseline v1** workflow:

- `R_INNER = 0.70`
- `R_OUTER = 0.90`
- `MU_MIN = 0.40`
- `ALPHA = 0.80`

This version should be kept unchanged as the reference baseline before introducing further physics updates.
