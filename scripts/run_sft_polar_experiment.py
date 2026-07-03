"""Polar-constraint SFT experiment: does a Solar Orbiter polar view change
the flux-transport outcome vs an Earth-view-only initial condition?

This is the clean form of the Stage-D test for the high-B0 campaign. At the
epochs where PHI has a real polar advantage the SolO-Earth separation is
large (80-165 deg in April 2025), so a per-pixel PHI+HMI merge is not
meaningful. Instead we build the initial condition from HMI everywhere
*except* the polar cap PHI observes, where PHI's zonal field is spliced in
(solar_pipeline.sft.apply_polar_constraint). We then evolve:

  - hmi              -- Earth-view only (the polar handicap: B=0 in the cap)
  - phi_constrained  -- HMI + PHI polar cap
  - phi              -- PHI-only zonal profile, for reference

with the same SFT configuration, and compare g10(t) and reversal timing.

Run per rotation, constraining the pole PHI covers that rotation, e.g.
    # March (CR 2295): PHI dominates the SOUTH cap
    python scripts/run_sft_polar_experiment.py --maps-dir out/cr_2295/milestone \
        --hemisphere south --source on --tau 10 --years 22 --balance-flux
    # April (CR 2296): PHI dominates the NORTH cap
    python scripts/run_sft_polar_experiment.py --maps-dir out/cr_2296/milestone \
        --hemisphere north --source on --tau 10 --years 22 --balance-flux
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
    apply_polar_constraint,
    axial_dipole_moment,
    balance_flux,
    polar_cap_mean,
    reversal_times,
    zonal_profile_from_map,
)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--maps-dir", type=Path, default=OUT_DIR / "milestone",
                   help="Directory holding grid_{hmi,phi}.npy + lat_centers.npy")
    p.add_argument("--out-dir", type=Path, default=None,
                   help="Default: <maps-dir>/sft_polar")
    p.add_argument("--hemisphere", choices=["north", "south", "both"], default="both",
                   help="Which polar cap PHI constrains (the pole it sees this rotation)")
    p.add_argument("--polar-lat-deg", type=float, default=60.0)
    p.add_argument("--blend-deg", type=float, default=10.0)
    # SFT configuration (mirrors run_sft_from_maps.py)
    p.add_argument("--years", type=float, default=11.0)
    p.add_argument("--flowtype", type=int, default=2, choices=[1, 2, 3, 4, 5])
    p.add_argument("--u0", type=float, default=11.0)
    p.add_argument("--eta", type=float, default=250.0)
    p.add_argument("--tau", type=float, default=None)
    p.add_argument("--dt-days", type=float, default=1.0)
    p.add_argument("--source", choices=["off", "on"], default="off")
    p.add_argument("--blat", type=float, default=0.0)
    p.add_argument("--bjoy", type=float, default=0.0)
    p.add_argument("--source-sigma", type=float, default=0.0)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--balance-flux", action="store_true")
    p.add_argument("--cap-deg", type=float, default=70.0)
    return p.parse_args()


def build_initial_conditions(maps_dir, latitude, lat_centers, hemisphere,
                             polar_lat_deg, blend_deg):
    """Return {label: B(latitude)} for the three initial conditions."""
    def load(name):
        path = maps_dir / f"grid_{name}.npy"
        if not path.exists():
            raise SystemExit(f"missing {path}")
        return np.load(path)

    grid_hmi, grid_phi = load("hmi"), load("phi")
    b_hmi = zonal_profile_from_map(grid_hmi, lat_centers, latitude, unobserved="zero")
    b_phi = zonal_profile_from_map(grid_phi, lat_centers, latitude, unobserved="zero")
    # NaN outside PHI's observed latitudes, so the splice touches only the
    # polar cap PHI actually measured
    b_phi_polar = zonal_profile_from_map(grid_phi, lat_centers, latitude, unobserved="nan")
    b_constrained = apply_polar_constraint(
        b_hmi, b_phi_polar, latitude,
        polar_lat_deg=polar_lat_deg, hemisphere=hemisphere, blend_deg=blend_deg,
    )
    return {"hmi": b_hmi, "phi_constrained": b_constrained, "phi": b_phi}


def main():
    args = parse_args()
    out_dir = args.out_dir or (args.maps_dir / "sft_polar")
    out_dir.mkdir(exist_ok=True, parents=True)

    lat_centers = np.load(args.maps_dir / "lat_centers.npy")
    model = SFTModel(flowtype=args.flowtype, u0=args.u0, eta=args.eta,
                     tau_years=args.tau, dt_days=args.dt_days)

    source_fn = None
    if args.source == "on":
        source_fn = HathawayJiangSource(
            tau_days=(args.tau or 1e9) * 365.25, blat=args.blat, bjoy=args.bjoy,
            sigma=args.source_sigma, seed=args.seed,
        )

    ics = build_initial_conditions(
        args.maps_dir, model.latitude, lat_centers, args.hemisphere,
        args.polar_lat_deg, args.blend_deg,
    )
    if args.balance_flux:
        ics = {k: balance_flux(v, model.theta) for k, v in ics.items()}

    results, reversal_rows = {}, []
    for label, b_init in ics.items():
        times, history = model.run(b_init, years=args.years, source_fn=source_fn)
        dip = np.array([axial_dipole_moment(b, model.theta) for b in history])
        capN = np.array([polar_cap_mean(b, model.latitude, args.cap_deg)["north"] for b in history])
        capS = np.array([polar_cap_mean(b, model.latitude, args.cap_deg)["south"] for b in history])
        results[label] = {"times": times, "dip": dip, "capN": capN, "capS": capS}
        revs = reversal_times(times, dip)
        print(f"{label:16s}: g10(0)={dip[0]:+.4f}  g10({args.years:.0f}yr)={dip[-1]:+.4f}  "
              f"capN(0)={capN[0]:+.3f} capS(0)={capS[0]:+.3f}  "
              f"reversals={['%.2f' % r for r in revs]}")
        reversal_rows.append({
            "ic": label, "g10_initial": float(dip[0]), "g10_final": float(dip[-1]),
            "capN_initial": float(capN[0]), "capS_initial": float(capS[0]),
            "n_reversals": len(revs),
            "first_reversal_yr": revs[0] if revs else float("nan"),
            "reversal_times_yr": ";".join(f"{r:.3f}" for r in revs),
        })

    # the headline number: how much the PHI polar constraint shifts things
    base, con = results["hmi"], results["phi_constrained"]
    d_g10_final = con["dip"][-1] - base["dip"][-1]
    rb = reversal_times(base["times"], base["dip"])
    rc = reversal_times(con["times"], con["dip"])
    d_rev = (rc[0] - rb[0]) if (rb and rc) else float("nan")
    print(f"\nPHI polar constraint ({args.hemisphere}) vs HMI-only:")
    print(f"  Delta g10_final   = {d_g10_final:+.4f} G")
    print(f"  Delta first-reversal = {d_rev:+.3f} yr")

    rows = []
    for k, t in enumerate(base["times"]):
        row = {"time_yr": float(t)}
        for label, r in results.items():
            row[f"g10_{label}"] = float(r["dip"][k])
            row[f"polar_north_{label}"] = float(r["capN"][k])
            row[f"polar_south_{label}"] = float(r["capS"][k])
        rows.append(row)
    pd.DataFrame(rows).to_csv(out_dir / "sft_polar_comparison.csv", index=False)
    pd.DataFrame(reversal_rows).to_csv(out_dir / "sft_polar_reversals.csv", index=False)
    np.save(out_dir / "sft_latitude.npy", model.latitude)
    for label, b in ics.items():
        np.save(out_dir / f"b_init_{label}.npy", b)

    colors = {"hmi": "#E69F00", "phi_constrained": "#0072B2", "phi": "#009E73"}
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    for label in ics:
        axes[0].plot(model.latitude, ics[label], color=colors[label], label=label)
        axes[1].plot(results[label]["times"], results[label]["dip"], color=colors[label], label=label)
    axes[0].axvspan(args.polar_lat_deg, 90, color="0.9", zorder=0)
    axes[0].axvspan(-90, -args.polar_lat_deg, color="0.9", zorder=0)
    axes[0].set_xlabel("latitude [deg]"); axes[0].set_ylabel("B init [G]")
    axes[0].set_title(f"Initial conditions (PHI constrains {args.hemisphere} cap)")
    axes[1].axhline(0, color="0.7", lw=1)
    axes[1].set_xlabel("time [yr]"); axes[1].set_ylabel("g10 [G]")
    axes[1].set_title("Axial dipole evolution")
    for ax in axes:
        ax.grid(alpha=0.3); ax.legend(fontsize=9)
    fig.suptitle(f"SFT polar-constraint experiment: source={args.source}, tau={args.tau or 'inf'} yr")
    fig.tight_layout()
    fig.savefig(out_dir / "sft_polar_comparison.png", dpi=150)
    plt.close(fig)

    print(f"\nOutputs in: {out_dir.resolve()}")


if __name__ == "__main__":
    main()
