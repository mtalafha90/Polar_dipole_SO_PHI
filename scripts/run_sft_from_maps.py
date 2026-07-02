"""Stage D driver: evolve the SFT model from observed synoptic map products.

Takes the three milestone map products (PHI-only, HMI-only Earth-view,
merged), zonally averages each into an initial B(latitude) profile, evolves
each with the same SFT configuration, and compares the axial dipole moment
evolution, polar-cap fields, and reversal timing.

The `--unobserved zero` default is the experiment itself: latitudes a
product never observed enter the SFT with B = 0, so an Earth-view-only
product carries its missing-polar-field handicap into the simulation, while
a merged product with SolO polar coverage does not.
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.append(str(Path(__file__).resolve().parents[1]))

from baseline_config import OUT_DIR
from solar_pipeline.sft import (
    SFTModel,
    HathawayJiangSource,
    axial_dipole_moment,
    balance_flux,
    polar_cap_mean,
    reversal_times,
    zonal_profile_from_map,
)

PRODUCTS = ("phi", "hmi", "merged")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--maps-dir",
        type=Path,
        default=OUT_DIR / "milestone",
        help="Directory holding grid_{phi,hmi,merged}.npy + lat_centers.npy (from run_milestone_comparison.py)",
    )
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR / "sft")
    parser.add_argument("--years", type=float, default=11.0)
    parser.add_argument("--flowtype", type=int, default=2, choices=[1, 2, 3, 4, 5])
    parser.add_argument("--u0", type=float, default=11.0, help="Meridional flow amplitude [m/s]")
    parser.add_argument("--eta", type=float, default=250.0, help="Diffusivity [km^2/s]")
    parser.add_argument("--tau", type=float, default=None, help="Decay time [yr]; omit for no decay")
    parser.add_argument("--dt-days", type=float, default=1.0)
    parser.add_argument(
        "--source",
        choices=["off", "on"],
        default="off",
        help="'off': pure decay/transport of the injected map. 'on': add the idealized cycle source",
    )
    parser.add_argument("--blat", type=float, default=0.0)
    parser.add_argument("--bjoy", type=float, default=0.0)
    parser.add_argument("--source-sigma", type=float, default=0.0, help="Cycle-amplitude scatter (0 = deterministic)")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--unobserved", choices=["zero", "extend"], default="zero")
    parser.add_argument(
        "--balance-flux",
        action="store_true",
        help="Remove net signed flux from each injected profile (recommended for "
        "partial-coverage maps, whose unbalanced flux otherwise drives "
        "unphysically large polar asymmetries)",
    )
    parser.add_argument("--cap-deg", type=float, default=70.0)
    return parser.parse_args()


def main():
    args = parse_args()
    args.out_dir.mkdir(exist_ok=True, parents=True)

    lat_centers = np.load(args.maps_dir / "lat_centers.npy")

    model = SFTModel(
        flowtype=args.flowtype, u0=args.u0, eta=args.eta,
        tau_years=args.tau, dt_days=args.dt_days,
    )
    source_fn = None
    if args.source == "on":
        tau_days = (args.tau or 1e9) * 365.25
        source_fn = HathawayJiangSource(
            tau_days=tau_days, blat=args.blat, bjoy=args.bjoy,
            sigma=args.source_sigma, seed=args.seed,
        )

    results = {}
    profiles = {}
    for product in PRODUCTS:
        grid_path = args.maps_dir / f"grid_{product}.npy"
        if not grid_path.exists():
            print(f"skipping {product}: {grid_path} not found")
            continue
        grid = np.load(grid_path)
        b_init = zonal_profile_from_map(grid, lat_centers, model.latitude, unobserved=args.unobserved)
        if args.balance_flux:
            b_init = balance_flux(b_init, model.theta)
        profiles[product] = b_init

        times, history = model.run(b_init, years=args.years, source_fn=source_fn)
        dip = np.array([axial_dipole_moment(b, model.theta) for b in history])
        capN = np.array([polar_cap_mean(b, model.latitude, args.cap_deg)["north"] for b in history])
        capS = np.array([polar_cap_mean(b, model.latitude, args.cap_deg)["south"] for b in history])
        results[product] = {"times": times, "dip": dip, "capN": capN, "capS": capS}
        revs = reversal_times(times, dip)
        print(
            f"{product}: g10(0)={dip[0]:+.4f}  g10({args.years:.0f}yr)={dip[-1]:+.4f}  "
            f"capN(0)={capN[0]:+.3f}  reversals={['%.2f' % r for r in revs]}"
        )

    if not results:
        raise RuntimeError(f"No map products found in {args.maps_dir}")

    rows = []
    ref = next(iter(results.values()))
    for k, t in enumerate(ref["times"]):
        row = {"time_yr": float(t)}
        for product, r in results.items():
            row[f"g10_{product}"] = float(r["dip"][k])
            row[f"polar_north_{product}"] = float(r["capN"][k])
            row[f"polar_south_{product}"] = float(r["capS"][k])
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(args.out_dir / "sft_comparison.csv", index=False)

    np.save(args.out_dir / "sft_latitude.npy", model.latitude)
    for product, b in profiles.items():
        np.save(args.out_dir / f"b_init_{product}.npy", b)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    for product in results:
        axes[0].plot(model.latitude, profiles[product], label=product)
        axes[1].plot(results[product]["times"], results[product]["dip"], label=product)
        axes[2].plot(results[product]["times"], results[product]["capN"], label=f"{product} N")
        axes[2].plot(results[product]["times"], results[product]["capS"], "--", label=f"{product} S")
    axes[0].set_xlabel("latitude [deg]")
    axes[0].set_ylabel("B init [G]")
    axes[0].set_title(f"Injected zonal profiles (unobserved={args.unobserved})")
    axes[1].set_xlabel("time [yr]")
    axes[1].set_ylabel("g10 [G]")
    axes[1].set_title("Axial dipole evolution")
    axes[2].set_xlabel("time [yr]")
    axes[2].set_ylabel(f"mean |lat|>{args.cap_deg:.0f} field [G]")
    axes[2].set_title("Polar-cap fields")
    for ax in axes:
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.suptitle(
        f"SFT from maps: flow {args.flowtype}, u0={args.u0} m/s, eta={args.eta} km2/s, "
        f"tau={args.tau or 'inf'} yr, source={args.source}"
    )
    fig.tight_layout()
    fig.savefig(args.out_dir / "sft_comparison.png", dpi=150)
    plt.close(fig)

    print(f"\nOutputs in: {args.out_dir.resolve()}")


if __name__ == "__main__":
    main()
