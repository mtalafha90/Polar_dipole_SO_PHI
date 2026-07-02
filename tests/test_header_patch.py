"""The SUMS header patcher must turn a bare JSOC segment file into one the
pipeline can use: WCS keys written, DATE__OBS renamed, junk values skipped,
existing keywords available afterwards."""

import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

sys.path.append(str(Path(__file__).resolve().parents[1]))

from scripts.download_baseline_data import _patch_fits_header


def _write_bare_segment(path: Path) -> None:
    # like a raw SUMS segment: image data, essentially no header
    fits.writeto(path, np.ones((32, 32), dtype="float32"), overwrite=True)


def test_patch_writes_wcs_and_skips_junk(tmp_path):
    f = tmp_path / "hmi.M_720s.20221030_041200_TAI.magnetogram.fits"
    _write_bare_segment(f)

    _patch_fits_header(
        f,
        {
            "CRPIX1": 2048.5,
            "CRPIX2": 2047.75,
            "CDELT1": 0.504,
            "CDELT2": 0.504,
            "CUNIT1": "arcsec",
            "CUNIT2": "arcsec",
            "CROTA2": 180.07,
            "RSUN_OBS": 967.32,
            "CRLT_OBS": 4.83,
            "CRLN_OBS": 108.9,
            "DATE__OBS": "2022-10-30T04:10:52.90",
            "T_REC": "2022.10.30_04:12:00_TAI",
            "QUALITY": None,          # missing value: skipped
            "OBS_VR": float("nan"),   # NaN: skipped
            "CAMERA": "MISSING",      # JSOC missing marker: skipped
        },
    )

    with fits.open(f) as hdul:
        header = hdul[0].header if hdul[0].data is not None else hdul[1].header
    assert header["CRPIX1"] == 2048.5
    assert header["CROTA2"] == 180.07
    assert header["CRLT_OBS"] == 4.83
    assert header["DATE-OBS"] == "2022-10-30T04:10:52.90"
    assert "DATE__OBS" not in header
    assert "QUALITY" not in header
    assert "OBS_VR" not in header
    assert "CAMERA" not in header


if __name__ == "__main__":
    import tempfile

    with tempfile.TemporaryDirectory() as d:
        test_patch_writes_wcs_and_skips_junk(Path(d))
    print("All tests passed.")
