"""The SUMS header patcher must turn a bare JSOC segment file into one the
pipeline can use: WCS keys written, DATE__OBS renamed, junk values skipped,
existing keywords available afterwards."""

import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

sys.path.append(str(Path(__file__).resolve().parents[1]))

from scripts.download_baseline_data import _patch_fits_header, _fetch_url


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


def test_patch_compressed_segment(tmp_path):
    # real SUMS segments are Rice-compressed integer images; the crashed
    # in-place update is avoided by rewrite-and-replace, which must keep the
    # data bit-identical (do_not_scale + lossless int recompression)
    f = tmp_path / "hmi.M_720s.20221030_041200_TAI.magnetogram.fits"
    data = (np.arange(1024, dtype=np.int32) - 512).reshape(32, 32)
    fits.HDUList([fits.PrimaryHDU(), fits.CompImageHDU(data)]).writeto(f)

    _patch_fits_header(f, {"CRPIX1": 16.5, "CROTA2": 180.07, "CRLT_OBS": 4.83})

    with fits.open(f, do_not_scale_image_data=True) as hdul:
        image = next(
            h for h in hdul
            if getattr(h, "data", None) is not None and getattr(h.data, "ndim", 0) == 2
        )
        assert image.header["CRPIX1"] == 16.5
        assert image.header["CROTA2"] == 180.07
        assert np.array_equal(image.data, data)


def test_fetch_url_cleans_partial_on_failure(tmp_path, monkeypatch):
    # a failed fetch (e.g. the 404 seen for some recent JSOC SUMS paths) must
    # not leave a half-written .part file behind, and must raise so the
    # per-record loop can fall back / record the failure
    import requests

    dest = tmp_path / "hmi.M_720s.20250225_150000_TAI.magnetogram.fits"

    class _Resp:
        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def raise_for_status(self):
            raise requests.HTTPError("404 Client Error: Not Found")

        def iter_content(self, n):
            return iter(())

    monkeypatch.setattr(requests, "get", lambda *a, **k: _Resp())

    raised = False
    try:
        _fetch_url("http://jsoc.stanford.edu/SUM/x", dest)
    except requests.HTTPError:
        raised = True

    assert raised
    assert not dest.exists()
    assert not dest.with_suffix(dest.suffix + ".part").exists()


if __name__ == "__main__":
    import tempfile

    with tempfile.TemporaryDirectory() as d:
        test_patch_writes_wcs_and_skips_junk(Path(d))
        test_patch_compressed_segment(Path(d))
    print("All tests passed.")
