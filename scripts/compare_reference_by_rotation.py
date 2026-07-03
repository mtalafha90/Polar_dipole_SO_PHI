"""Reference-dipole check for every rotation of a per-CR campaign.

Wraps compare_reference_dipole.py: for each ``cr_<N>/milestone`` directory
produced by run_milestone_by_rotation.py, it fetches the matching HMI
synoptic chart (``--car-rot N``) and tabulates the PHI / HMI / merged
milestone dipoles against it. A combined ``reference_dipole_by_rotation.csv``
is written at the campaign root.

Requires the download extras (drms + requests) and JSOC network access, as
the reference charts are fetched by Carrington rotation.

Usage:
    python scripts/compare_reference_by_rotation.py --campaign-dir out
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.append(str(Path(__file__).resolve().parents[1]))

from baseline_config import OUT_DIR

COMPARE_SCRIPT = Path(__file__).resolve().parent / "compare_reference_dipole.py"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--campaign-dir", type=Path, default=OUT_DIR,
                   help="Directory containing cr_<N>/ subdirectories.")
    p.add_argument("--series", type=str, default="hmi.synoptic_mr_polfil_720s",
                   help="JSOC synoptic series to compare against.")
    return p.parse_args()


def main():
    args = parse_args()
    cr_dirs = sorted(
        (d for d in args.campaign_dir.glob("cr_*")
         if (d / "milestone" / "lat_centers.npy").exists()),
        key=lambda d: int(re.search(r"cr_(\d+)", d.name).group(1)),
    )
    if not cr_dirs:
        raise SystemExit(f"No cr_<N>/milestone/ outputs found under {args.campaign_dir}")

    combined = []
    for d in cr_dirs:
        cr = int(re.search(r"cr_(\d+)", d.name).group(1))
        maps_dir = d / "milestone"
        print(f"\n=== CR {cr}: reference check ({maps_dir}) ===")
        cmd = [
            sys.executable, str(COMPARE_SCRIPT),
            "--car-rot", str(cr),
            "--series", args.series,
            "--maps-dir", str(maps_dir),
            "--reference-dir", str(d / "reference"),
        ]
        rc = subprocess.run(cmd).returncode
        if rc != 0:
            print(f"  CR {cr}: reference check failed (rc={rc}); continuing")
            continue
        csv = maps_dir / "reference_dipole_comparison.csv"
        if csv.exists():
            df = pd.read_csv(csv)
            df.insert(0, "car_rot", cr)
            combined.append(df)

    if combined:
        out = pd.concat(combined, ignore_index=True)
        out_csv = args.campaign_dir / "reference_dipole_by_rotation.csv"
        out.to_csv(out_csv, index=False)
        print(f"\nWrote {out_csv}  ({len(combined)} rotations)")
    else:
        print("\nNo per-rotation reference comparisons were produced.")


if __name__ == "__main__":
    main()
