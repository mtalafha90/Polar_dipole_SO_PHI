import argparse
import os
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))

from baseline_config import PHI_DIR, HMI_DIR, MAX_TIME_DIFF_SEC
from solar_pipeline.io_utils import parse_phi_time

JSOC_BASE = "http://jsoc.stanford.edu"


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Download the baseline input data described in the README: "
            "Solar Orbiter PHI-FDT L2 blos magnetograms from the ESA SOAR archive, "
            "and for each PHI file the nearest SDO/HMI M_720s magnetogram from JSOC. "
            "Requires network access to soar.esac.esa.int and jsoc.stanford.edu, and "
            "the optional download dependencies: pip install -e '.[download]'"
        )
    )
    parser.add_argument("--start", type=str, default="2022-10-27", help="Start date (YYYY-MM-DD, inclusive)")
    parser.add_argument("--end", type=str, default="2022-10-29", help="End date (YYYY-MM-DD, exclusive)")
    parser.add_argument("--phi-dir", type=Path, default=PHI_DIR)
    parser.add_argument("--hmi-dir", type=Path, default=HMI_DIR)
    parser.add_argument(
        "--max-time-diff-sec",
        type=float,
        default=MAX_TIME_DIFF_SEC,
        help="Warn if the nearest available HMI record is further than this from a PHI time",
    )
    parser.add_argument(
        "--hmi-window-min",
        type=float,
        default=30.0,
        help="Half-width (minutes) of the JSOC query window around each PHI time",
    )
    parser.add_argument("--phi-only", action="store_true", help="Only download PHI files")
    parser.add_argument("--hmi-only", action="store_true", help="Only download HMI files (PHI files must already exist)")
    parser.add_argument(
        "--jsoc-email",
        type=str,
        default=os.environ.get("JSOC_EXPORT_EMAIL"),
        help="JSOC-registered email (register at http://jsoc.stanford.edu/ajax/register_email.html). "
        "Only needed as a fallback when a record is offline at JSOC — online records are fetched "
        "directly from SUMS without any registration. Defaults to $JSOC_EXPORT_EMAIL.",
    )
    return parser.parse_args()


def download_phi(start: str, end: str, phi_dir: Path) -> list[Path]:
    try:
        from sunpy.net import Fido, attrs as a
    except ImportError as exc:
        raise SystemExit(
            f"Missing download dependency ({exc}). Install with: pip install -e '.[download]'"
        )
    if not hasattr(a, "soar"):
        # sunpy >= 8 ships its own SOAR client; older sunpy needs sunpy-soar
        try:
            import sunpy_soar  # noqa: F401  (registers a.soar attrs with Fido)
        except ImportError:
            raise SystemExit(
                "No SOAR client available: upgrade to sunpy>=8 or pip install sunpy-soar"
            )

    phi_dir.mkdir(exist_ok=True, parents=True)

    print(f"Querying SOAR for phi-fdt-blos L2 files in [{start}, {end}) ...")
    result = Fido.search(
        a.Time(start, end),
        a.Instrument("PHI"),
        a.Level(2),
        a.soar.Product("phi-fdt-blos"),
    )
    n_found = sum(len(block) for block in result)
    print(f"Found {n_found} PHI blos files on SOAR")
    if n_found == 0:
        return []

    downloaded = Fido.fetch(result, path=str(phi_dir / "{file}"), progress=True)
    if downloaded.errors:
        for err in downloaded.errors:
            print(f"  DOWNLOAD ERROR: {err}")
    files = sorted(Path(f) for f in downloaded)
    print(f"Downloaded {len(files)} PHI files into {phi_dir.resolve()}")
    return files


def tai_str(t: datetime) -> str:
    return t.strftime("%Y.%m.%d_%H:%M:%S_TAI")


def parse_t_rec(t_rec: str) -> datetime:
    # e.g. "2022.10.27_00:12:00_TAI"; treated as UTC for matching purposes,
    # consistent with parse_hmi_time's handling of filename timestamps
    return datetime.strptime(t_rec, "%Y.%m.%d_%H:%M:%S_TAI").replace(tzinfo=timezone.utc)


def nearest_t_rec(t_recs: list[str], phi_time: datetime) -> tuple[str, float]:
    best, best_dt = None, None
    for t_rec in t_recs:
        dt = abs((parse_t_rec(t_rec) - phi_time).total_seconds())
        if best_dt is None or dt < best_dt:
            best, best_dt = t_rec, dt
    if best is None:
        raise RuntimeError("No HMI records returned for query window.")
    return best, best_dt


def _fetch_url(url: str, dest: Path) -> None:
    import requests

    tmp = dest.with_suffix(dest.suffix + ".part")
    with requests.get(url, stream=True, timeout=300) as resp:
        resp.raise_for_status()
        with open(tmp, "wb") as f:
            for chunk in resp.iter_content(1 << 20):
                f.write(chunk)
    tmp.rename(dest)


def download_hmi(
    phi_files: list[Path],
    hmi_dir: Path,
    window_min: float,
    max_time_diff_sec: float,
    email: str | None = None,
) -> list[Path]:
    try:
        import drms
    except ImportError as exc:
        raise SystemExit(
            f"Missing download dependency ({exc}). Install with: pip install -e '.[download]'"
        )

    hmi_dir.mkdir(exist_ok=True, parents=True)
    client = drms.Client(email=email) if email else drms.Client()

    wanted: dict[str, float] = {}
    for phi_path in phi_files:
        phi_time = parse_phi_time(phi_path)
        t0 = phi_time - timedelta(minutes=window_min)
        ds = f"hmi.M_720s[{tai_str(t0)}/{2 * window_min}m]"
        records = client.query(ds, key="T_REC")
        t_recs = list(records["T_REC"]) if len(records) else []
        if not t_recs:
            print(f"  WARNING: no HMI records near {phi_path.name}")
            continue
        t_rec, dt = nearest_t_rec(t_recs, phi_time)
        if dt > max_time_diff_sec:
            print(f"  WARNING: nearest HMI record for {phi_path.name} is {dt:.0f}s away (> {max_time_diff_sec:.0f}s)")
        if t_rec not in wanted or dt < wanted[t_rec]:
            wanted[t_rec] = dt

    print(f"Need {len(wanted)} unique HMI M_720s records")

    files = []
    n_total = len(wanted)
    for i, t_rec in enumerate(sorted(wanted), start=1):
        expected = hmi_dir / f"hmi.M_720s.{parse_t_rec(t_rec).strftime('%Y%m%d_%H%M%S')}_TAI.magnetogram.fits"
        if expected.exists():
            print(f"  [{i}/{n_total}] already present: {expected.name}")
            files.append(expected)
            continue
        ds = f"hmi.M_720s[{t_rec}]{{magnetogram}}"
        print(f"  [{i}/{n_total}] fetching {ds}")

        # No-registration path: rs_list returns the SUMS path for records
        # that are online at JSOC, downloadable over plain HTTP.
        _, segs = client.query(ds, key="T_REC", seg="magnetogram")
        seg_path = str(segs["magnetogram"].iloc[0]) if len(segs) else ""
        if seg_path.startswith("/"):
            _fetch_url(JSOC_BASE + seg_path, expected)
            files.append(expected)
            continue

        # Offline record: needs a real JSOC export request (registered email).
        if not email:
            raise SystemExit(
                f"JSOC record {t_rec} is not online in SUMS, so it needs an export "
                "request, which requires a JSOC-registered email. Register at "
                "http://jsoc.stanford.edu/ajax/register_email.html and re-run with "
                "--jsoc-email <address> (or set JSOC_EXPORT_EMAIL)."
            )
        export = client.export(ds, method="url_quick", protocol="as-is", email=email)
        result = export.download(str(hmi_dir))
        for local in result["download"]:
            local = Path(local)
            if local.name != expected.name:
                local.rename(expected)
            files.append(expected)
    print(f"Downloaded {len(files)} HMI files into {hmi_dir.resolve()}")
    return files


def main():
    args = parse_args()

    if not args.hmi_only:
        download_phi(args.start, args.end, args.phi_dir)

    if not args.phi_only:
        phi_files = sorted(args.phi_dir.glob("solo_L2_phi-fdt-blos_*.fits"))
        # match HMI only for PHI files inside the requested window, so a
        # data directory holding a longer campaign doesn't trigger extra
        # downloads
        d0 = datetime.strptime(args.start, "%Y-%m-%d").replace(tzinfo=timezone.utc)
        d1 = datetime.strptime(args.end, "%Y-%m-%d").replace(tzinfo=timezone.utc)
        phi_files = [f for f in phi_files if d0 <= parse_phi_time(f) < d1]
        if not phi_files:
            raise SystemExit(
                f"No PHI blos files in {args.phi_dir.resolve()} within [{args.start}, {args.end}); "
                "run without --hmi-only first or adjust --start/--end."
            )
        download_hmi(phi_files, args.hmi_dir, args.hmi_window_min, args.max_time_diff_sec, email=args.jsoc_email)

    print("\nDone. You can now run: python scripts/run_baseline_pipeline.py")


if __name__ == "__main__":
    main()
