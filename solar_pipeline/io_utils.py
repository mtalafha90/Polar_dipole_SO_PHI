from __future__ import annotations

import re
from datetime import datetime, timezone
from pathlib import Path

import sunpy.map
from astropy.io import fits


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