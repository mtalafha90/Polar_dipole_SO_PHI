"""First-milestone comparison (see README section 6b):

Build three Carrington-style synoptic Br maps over the selected subset —
PHI-only (SolO vantage, native geometry), HMI-only (Earth vantage, native
geometry, NOT reprojected through the PHI grid), and merged (PHI+HMI blend
on the PHI grid) — then compute and compare the axial dipole moment g10
from each, under each polar-filling assumption.
"""

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sunpy.util.exceptions import SunpyMetadataWarning

warnings.filterwarnings("ignore", category=SunpyMetadataWarning)
sys.path.append(str(Path(__file__).resolve().parents[1]))

from baseline_config import (
    PHI_DIR,
    HMI_DIR,
    OUT_DIR,
    ONLY_DATES,
    MAX_TIME_DIFF_SEC,
    R_INNER,
    R_OUTER,
    DISK_FRACTION,
    MU_MIN,
    ALPHA,
    NLAT,
    NLON,
)
from solar_pipeline.io_utils import build_hmi_time_index, expand_date_spec, save_synoptic_fits
from solar_pipeline.pipeline import compute_case_fields, compute_native_disk_fields
from solar_pipeline.calibration import calibration_stats
from solar_pipeline.carrington import (
    cm_weight,
    bin_br_to_carrington_weighted,
    bin_max_to_carrington,
    combine_weighted_grids,
    carrington_fill_fraction,
)
from solar_pipeline.dipole import FILL_MODES, axial_dipole_g10, polar_fill_fractions
from solar_pipeline.plotting import plot_carrington_map


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phi-dir", type=Path, default=PHI_DIR)
    parser.add_argument("--hmi-dir", type=Path, default=HMI_DIR)
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    parser.add_argument("--dates", type=str, default=",".join(sorted(ONLY_DATES)))
    parser.add_argument("--max-time-diff-sec", type=float, default=MAX_TIME_DIFF_SEC)
    parser.add_argument("--r-inner", type=float, default=R_INNER)
    parser.add_argument("--r-outer", type=float, default=R_OUTER)
    parser.add_argument("--disk-fraction", type=float, default=DISK_FRACTION)
    parser.add_argument("--mu-min", type=float, default=MU_MIN)
    parser.add_argument("--alpha", type=float, default=ALPHA)
    parser.add_argument("--nlat", type=int, default=NLAT)
    parser.add_argument("--nlon", type=int, default=NLON)
    parser.add_argument("--cm-weight-power", type=float, default=1.0)
    parser.add_argument(
        "--polar-lat-deg",
        type=float,
        default=60.0,
        help="Latitude bound for polar-cap fill-fraction diagnostics",
    )
    parser.add_argument(
        "--calibrate-phi",
        action="store_true",
        help="Divide PHI Br by the per-case PHI-vs-HMI regression slope before binning",
    )
    parser.add_argument(
        "--calib-min-abs-g",
        type=float,
        default=10.0,
        help="|B| threshold (Gauss) for pixels entering the calibration regression",
    )
    parser.add_argument(
        "--calib-min-r",
        type=float,
        default=0.5,
        help="Minimum Pearson r for a per-case calibration slope to be applied; "
        "an uncorrelated fit (e.g. from misaligned data) would otherwise scale "
        "PHI by a meaningless factor",
    )
    parser.add_argument(
        "--quiet-sun-max-g",
        type=float,
        default=None,
        help="If set, also report *_quiet rows: dipoles from bins with "
        "|Br| below this threshold, excluding active regions (where the "
        "radial-field assumption fails and the two vantages disagree most)",
    )
    return parser.parse_args()


class Accumulator:
    """Accumulates weighted synoptic bins and max-mu confidence per product.

    Uses running totals rather than per-case grids so memory stays flat for
    full-Carrington-rotation runs with ~100+ cases.
    """

    def __init__(self, nlat: int, nlon: int):
        self.nlat, self.nlon = nlat, nlon
        self.wsum_total = np.zeros((nlat, nlon))
        self.weight_total = np.zeros((nlat, nlon))
        self.quality = np.full((nlat, nlon), np.nan)
        self.lat_centers = self.lon_centers = None

    def add(self, br, lat, lon, valid, weight, mu):
        wsum, wt, self.lat_centers, self.lon_centers = bin_br_to_carrington_weighted(
            br, lat, lon, valid, weight, nlat=self.nlat, nlon=self.nlon
        )
        self.wsum_total += wsum
        self.weight_total += wt
        q = bin_max_to_carrington(mu, lat, lon, valid, nlat=self.nlat, nlon=self.nlon)
        self.quality = np.fmax(self.quality, q)

    def combine(self):
        grid, total_weight = combine_weighted_grids([self.wsum_total], [self.weight_total])
        return grid, total_weight, self.quality


def main():
    args = parse_args()
    only_dates = expand_date_spec(args.dates)

    out_dir = args.out_dir / "milestone"
    out_dir.mkdir(exist_ok=True, parents=True)
    plots_dir = out_dir / "plots"
    plots_dir.mkdir(exist_ok=True, parents=True)

    phi_blos_files = sorted(args.phi_dir.glob("solo_L2_phi-fdt-blos_*.fits"))
    phi_blos_files = [f for f in phi_blos_files if only_dates is None or any(d in f.name for d in only_dates)]
    hmi_files = sorted(args.hmi_dir.glob("hmi.M_720s.*.magnetogram.fits"))

    if not phi_blos_files:
        raise RuntimeError("No PHI blos files found.")
    if not hmi_files:
        raise RuntimeError("No HMI files found.")

    hmi_index = build_hmi_time_index(hmi_files)

    acc = {name: Accumulator(args.nlat, args.nlon) for name in ("phi", "hmi", "merged")}
    calib_rows = []
    seen_hmi = set()
    n_ok = 0

    for i, phi_blos_path in enumerate(phi_blos_files, start=1):
        print(f"[{i}/{len(phi_blos_files)}] {phi_blos_path.name}")
        try:
            fields = compute_case_fields(
                phi_blos_path,
                hmi_index,
                max_time_diff_sec=args.max_time_diff_sec,
                r_inner=args.r_inner,
                r_outer=args.r_outer,
                disk_fraction=args.disk_fraction,
                mu_min=args.mu_min,
                alpha=args.alpha,
            )

            calib = calibration_stats(
                fields["phi_blos"].data, fields["hmi_on_phi"].data, fields["mu"],
                mu_min=args.mu_min, min_abs_ref=args.calib_min_abs_g,
            )
            phi_crln = float(fields["phi_blos"].meta.get("crln_obs", np.nan))
            hmi_crln = float(fields["hmi"].meta.get("crln_obs", np.nan))
            calib_rows.append({
                "phi_blos_file": phi_blos_path.name,
                "phi_time": fields["phi_time"].isoformat(),
                "phi_dsun_au": float(fields["phi_blos"].meta.get("dsun_obs", np.nan)) / 1.495978707e11,
                "lon_separation_deg": abs(((phi_crln - hmi_crln) + 180.0) % 360.0 - 180.0),
                **calib,
            })
            print(
                f"  calib: PHI ~= {calib['slope']:.3f} x HMI "
                f"(r={calib['pearson_r']:.3f}, n={calib['n_pixels']})"
            )

            phi_scale = 1.0
            if args.calibrate_phi:
                if (
                    np.isfinite(calib["slope"])
                    and calib["slope"] > 0
                    and abs(calib["pearson_r"]) >= args.calib_min_r
                ):
                    phi_scale = 1.0 / calib["slope"]
                else:
                    print(
                        f"  WARNING: calibration unreliable (slope={calib['slope']:.3f}, "
                        f"r={calib['pearson_r']:.3f}); not applied to this case"
                    )

            # compute everything BEFORE accumulating, so a failure partway
            # through a case cannot leave the products inconsistently filled
            native = None
            if fields["hmi_path"] not in seen_hmi:
                # HMI-only, Earth vantage: the matched HMI magnetogram in its
                # OWN disk geometry (deduplicated if two PHI files matched
                # the same HMI record)
                native = compute_native_disk_fields(
                    fields["hmi"],
                    disk_fraction=args.disk_fraction,
                    mu_min=args.mu_min,
                    alpha=args.alpha,
                )

            w_phi = cm_weight(fields["cmd"], power=args.cm_weight_power)

            # PHI-only, SolO vantage (PHI is already in its native geometry)
            acc["phi"].add(
                fields["br_phi"] * phi_scale, fields["lat"], fields["lon"],
                fields["valid_phi"], w_phi, fields["mu"],
            )
            # merged product lives on the PHI grid
            br_merged = fields["br_merged"] * phi_scale if args.calibrate_phi else fields["br_merged"]
            acc["merged"].add(
                br_merged, fields["lat"], fields["lon"],
                fields["valid_merged"], w_phi, fields["mu"],
            )
            if native is not None:
                seen_hmi.add(fields["hmi_path"])
                acc["hmi"].add(
                    native["br"], native["lat"], native["lon"],
                    native["valid"], cm_weight(native["cmd"], power=args.cm_weight_power),
                    native["mu"],
                )

            n_ok += 1
        except Exception as exc:
            print(f"  ERROR: {exc}")

    if n_ok == 0:
        raise RuntimeError("No cases processed successfully.")

    combined = {name: acc[name].combine() for name in ("phi", "hmi", "merged")}
    # bins observed by every product: comparing dipoles on this common
    # support isolates vantage/calibration effects from coverage effects
    # (the PHI and HMI vantage points can be tens of degrees apart in
    # longitude, which otherwise dominates the product differences)
    common_mask = np.logical_and.reduce([combined[n][1] > 0 for n in combined])

    rows = []
    for name in ("phi", "hmi", "merged"):
        grid, weight, quality = combined[name]
        lat_c, lon_c = acc[name].lat_centers, acc[name].lon_centers

        np.save(out_dir / f"grid_{name}.npy", grid)
        np.save(out_dir / f"weight_{name}.npy", weight)
        np.save(out_dir / f"quality_{name}.npy", quality)
        save_synoptic_fits(
            out_dir / f"synoptic_{name}.fits", grid, lat_c, lon_c,
            extra={"product": name, "history": "SO-PHI-SDO-HMI milestone comparison"},
        )
        plot_carrington_map(
            grid, lat_c, lon_c, plots_dir / f"map_{name}.png",
            title=f"{name} synoptic map ({n_ok} cases)",
        )

        polar = polar_fill_fractions(weight, lat_c, polar_lat_deg=args.polar_lat_deg)
        row = {
            "product": name,
            "fill_fraction": carrington_fill_fraction(weight),
            f"polar_fill_north_{args.polar_lat_deg:.0f}": polar["north"],
            f"polar_fill_south_{args.polar_lat_deg:.0f}": polar["south"],
        }
        for mode in FILL_MODES:
            dip = axial_dipole_g10(grid, lat_c, mode=mode)
            row[f"g10_{mode}"] = dip["g10"]
            row[f"g10_north_{mode}"] = dip["g10_north"]
            row[f"g10_south_{mode}"] = dip["g10_south"]
        rows.append(row)

        common_grid = np.where(common_mask, grid, np.nan)
        common_row = {
            "product": f"{name}_common",
            "fill_fraction": float(np.count_nonzero(common_mask) / common_mask.size),
        }
        for mode in FILL_MODES:
            dip = axial_dipole_g10(common_grid, lat_c, mode=mode)
            common_row[f"g10_{mode}"] = dip["g10"]
            common_row[f"g10_north_{mode}"] = dip["g10_north"]
            common_row[f"g10_south_{mode}"] = dip["g10_south"]
        rows.append(common_row)

        if args.quiet_sun_max_g is not None:
            for label, base in ((f"{name}_quiet", grid), (f"{name}_quiet_common", common_grid)):
                quiet_grid = np.where(np.abs(base) <= args.quiet_sun_max_g, base, np.nan)
                quiet_row = {
                    "product": label,
                    "fill_fraction": float(np.count_nonzero(np.isfinite(quiet_grid)) / quiet_grid.size),
                }
                for mode in FILL_MODES:
                    dip = axial_dipole_g10(quiet_grid, lat_c, mode=mode)
                    quiet_row[f"g10_{mode}"] = dip["g10"]
                    quiet_row[f"g10_north_{mode}"] = dip["g10_north"]
                    quiet_row[f"g10_south_{mode}"] = dip["g10_south"]
                rows.append(quiet_row)

    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "milestone_dipole_comparison.csv", index=False)
    calib_df = pd.DataFrame(calib_rows)
    calib_df.to_csv(out_dir / "calibration_stats.csv", index=False)

    # map-space correlation between products on the common support: a
    # decisive orientation/consistency check. Strongly positive is expected;
    # near-zero or negative means one product's map is misplaced (wrong WCS,
    # unhandled rotation, ...).
    corr_rows = []
    grids = {name: combined[name][0] for name in combined}
    for a, b in (("phi", "hmi"), ("phi", "merged"), ("hmi", "merged")):
        sel = common_mask & np.isfinite(grids[a]) & np.isfinite(grids[b])
        r = float(np.corrcoef(grids[a][sel], grids[b][sel])[0, 1]) if sel.sum() >= 10 else float("nan")
        corr_rows.append({"pair": f"{a}-{b}", "pearson_r": r, "n_bins": int(sel.sum())})
    corr_df = pd.DataFrame(corr_rows)
    corr_df.to_csv(out_dir / "map_correlations.csv", index=False)

    np.save(out_dir / "lat_centers.npy", acc["phi"].lat_centers)
    np.save(out_dir / "lon_centers.npy", acc["phi"].lon_centers)

    print("\n=== Milestone dipole comparison ===")
    print(df.to_string(index=False))
    print("\n=== Map-space correlations on common support ===")
    print("(orientation/consistency check: strongly positive expected)")
    print(corr_df.to_string(index=False))
    print("\n=== PHI-vs-HMI calibration ===")
    print(calib_df.to_string(index=False))
    print(f"\nOutputs in: {out_dir.resolve()}")


if __name__ == "__main__":
    main()
