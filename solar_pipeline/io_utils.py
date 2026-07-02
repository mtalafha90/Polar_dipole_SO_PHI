from __future__ import annotations

import re
from datetime import datetime, timedelta, timezone
from pathlib import Path

import sunpy.map
from astropy.io import fits


def expand_date_spec(spec: str) -> set[str] | None:
    """Expand a --dates specification into a set of YYYYMMDD strings.

    Accepts comma-separated tokens, each either a single date (20221027) or
    an inclusive range (20221017-20221113). Returns None for "all"/"*",
    meaning no date filtering.
    """
    spec = spec.strip()
    if spec.lower() in ("all", "*"):
        return None
    out: set[str] = set()
    for token in spec.split(","):
        token = token.strip()
        if not token:
            continue
        if "-" in token:
            a, b = token.split("-", 1)
            d0 = datetime.strptime(a.strip(), "%Y%m%d")
            d1 = datetime.strptime(b.strip(), "%Y%m%d")
            if d1 < d0:
                raise ValueError(f"Date range end before start: {token}")
            cur = d0
            while cur <= d1:
                out.add(cur.strftime("%Y%m%d"))
                cur += timedelta(days=1)
        else:
            datetime.strptime(token, "%Y%m%d")
            out.add(token)
    if not out:
        raise ValueError(f"Empty date specification: {spec!r}")
    return out


def parse_phi_time(path: Path) -> datetime:
    m = re.search(r"_(\d{8}T\d{6})_", path.name)
    if not m:
        raise ValueError(f"Could not parse PHI time from {path.name}")
    return datetime.strptime(m.group(1), "%Y%m%dT%H%M%S").replace(tzinfo=timezone.utc)


def parse_hmi_time(path: Path) -> datetime:
    m = re.search(r"\.(\d{8})_(\d{6})_TAI\.", path.name)
    if not m:
        raise ValueError(f"Could not parse HMI time from {path.name}")
    return datetime.strptime(m.group(1) + m.group(2), "%Y%m%d%H%M%S").replace(tzinfo=timezone.utc)


def load_map(path: Path) -> sunpy.map.Map:
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")
    return sunpy.map.Map(path)


def build_hmi_time_index(hmi_files: list[Path]) -> list[tuple[Path, datetime]]:
    out: list[tuple[Path, datetime]] = []
    for f in hmi_files:
        out.append((f, parse_hmi_time(f)))
    return out


def find_nearest_hmi(phi_time: datetime, hmi_index: list[tuple[Path, datetime]]) -> tuple[Path, datetime, float]:
    best_file, best_time, best_dt = None, None, None
    for f, t in hmi_index:
        dt = abs((t - phi_time).total_seconds())
        if best_dt is None or dt < best_dt:
            best_file, best_time, best_dt = f, t, dt
    if best_file is None or best_time is None or best_dt is None:
        raise RuntimeError("No HMI files available.")
    return best_file, best_time, best_dt


def save_fits(path: Path, data, header) -> None:
    fits.writeto(path, data=data, header=header, overwrite=True)


def save_synoptic_fits(path: Path, grid, lat_centers, lon_centers, extra: dict | None = None) -> None:
    """Write a Carrington-style (linear latitude x longitude) Br grid as a
    FITS image with a plate-carree WCS, suitable as SFT model input."""
    import numpy as np

    header = fits.Header()
    header["ctype1"] = "CRLN-CAR"
    header["ctype2"] = "CRLT-CAR"
    header["cunit1"] = "deg"
    header["cunit2"] = "deg"
    header["crpix1"] = 1.0
    header["crpix2"] = 1.0
    header["crval1"] = float(np.rad2deg(lon_centers[0]))
    header["crval2"] = float(np.rad2deg(lat_centers[0]))
    header["cdelt1"] = float(np.rad2deg(lon_centers[1] - lon_centers[0]))
    header["cdelt2"] = float(np.rad2deg(lat_centers[1] - lat_centers[0]))
    header["bunit"] = "Gauss"
    for key, val in (extra or {}).items():
        header[key] = val
    fits.writeto(path, data=grid.astype("float32"), header=header, overwrite=True)