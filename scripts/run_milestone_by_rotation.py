"""Run the milestone comparison once per Carrington rotation.

A multi-month window (e.g. the 2025 high-B0 campaign) spans several
Carrington rotations. A single combined synoptic map smears the dipole
across rotations and cannot be checked against a per-rotation reference
chart. This wrapper assigns each day in the requested range to its
Carrington rotation and invokes run_milestone_comparison.py once per
rotation, writing outputs to <out-dir>/cr_<N>/.

Extra flags for the milestone script are passed through after a literal
`--`, e.g.:

    python scripts/run_milestone_by_rotation.py --dates 20250211-20250429 \
        -- --calibrate-phi --quiet-sun-max-g 50 --max-separation-deg 60
"""

import argparse
import subprocess
import sys
from collections import defaultdict
from datetime import datetime, timedelta, timezone
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))

from baseline_config import PHI_DIR, HMI_DIR, OUT_DIR
from solar_pipeline.io_utils import expand_date_spec

MILESTONE_SCRIPT = Path(__file__).resolve().parent / "run_milestone_comparison.py"


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dates", type=str, required=True, help="Date spec (range/list; 'all' not supported here)")
    parser.add_argument("--phi-dir", type=Path, default=PHI_DIR)
    parser.add_argument("--hmi-dir", type=Path, default=HMI_DIR)
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    parser.add_argument(
        "passthrough",
        nargs=argparse.REMAINDER,
        help="Flags after `--` are forwarded to run_milestone_comparison.py",
    )
    return parser.parse_args()


def carrington_rotation_of(date_str: str) -> int:
    from astropy.time import Time
    from sunpy.coordinates.sun import carrington_rotation_number
    import numpy as np

    # assign each day by its noon, so a day never straddles two rotations
    t = Time(datetime.strptime(date_str, "%Y%m%d").replace(hour=12, tzinfo=timezone.utc))
    return int(np.floor(float(carrington_rotation_number(t))))


def main():
    args = parse_args()

    dates = expand_date_spec(args.dates)
    if dates is None:
        raise SystemExit("--dates 'all' is not supported here; give an explicit range or list.")

    groups: dict[int, list[str]] = defaultdict(list)
    for d in sorted(dates):
        groups[carrington_rotation_of(d)].append(d)

    passthrough = args.passthrough
    if passthrough and passthrough[0] == "--":
        passthrough = passthrough[1:]

    print(f"{len(dates)} days span {len(groups)} Carrington rotation(s): {sorted(groups)}")

    results = []
    for cr, ds in sorted(groups.items()):
        cr_out = args.out_dir / f"cr_{cr}"
        print(f"\n=== CR {cr}: {len(ds)} days ({ds[0]}..{ds[-1]}) -> {cr_out} ===")
        cmd = [
            sys.executable, str(MILESTONE_SCRIPT),
            "--dates", ",".join(ds),
            "--phi-dir", str(args.phi_dir),
            "--hmi-dir", str(args.hmi_dir),
            "--out-dir", str(cr_out),
            *passthrough,
        ]
        rc = subprocess.run(cmd).returncode
        results.append((cr, len(ds), rc))
        if rc != 0:
            print(f"  CR {cr} exited with code {rc} (continuing)")

    print("\n=== Summary ===")
    for cr, n, rc in results:
        status = "ok" if rc == 0 else f"FAILED (rc={rc})"
        print(f"  CR {cr}: {n} days -> {status}  [{args.out_dir / f'cr_{cr}' / 'milestone'}]")
    print(
        "\nFor a per-rotation reference check, run e.g.:\n"
        "  python scripts/compare_reference_dipole.py --car-rot <CR> "
        "--maps-dir <out-dir>/cr_<CR>/milestone"
    )


if __name__ == "__main__":
    main()
