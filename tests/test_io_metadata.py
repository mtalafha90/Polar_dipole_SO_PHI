"""Regression test for PHI-FDT files missing CUNIT1/CUNIT2.

Some PHI L2 observing programs (e.g. the Oct 30+ 2022 files) lack coordinate
units in their headers; sunpy >= 8 refuses to build a Map from them
("Image coordinate units for axis 1 not present in metadata"). load_map
falls back to patching the missing units with arcsec. The fallback is
exercised directly so the test is meaningful on sunpy 7 (warning-only) too.
"""

import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

sys.path.append(str(Path(__file__).resolve().parents[1]))

from solar_pipeline.io_utils import load_map, _load_map_patched_units


def _write_cunitless_fits(path: Path) -> None:
    data = np.ones((64, 64), dtype="float32")
    header = fits.Header()
    header["ctype1"] = "HPLN-TAN"
    header["ctype2"] = "HPLT-TAN"
    header["crpix1"] = 32.5
    header["crpix2"] = 32.5
    header["cdelt1"] = 6.0
    header["cdelt2"] = 6.0
    header["crval1"] = 0.0
    header["crval2"] = 0.0
    header["rsun_obs"] = 160.0
    header["crlt_obs"] = 5.0
    header["crln_obs"] = 100.0
    header["date-obs"] = "2022-10-30T04:15:03"
    # deliberately no CUNIT1/CUNIT2
    fits.writeto(path, data, header=header, overwrite=True)


def test_patched_loader_fills_units(tmp_path):
    f = tmp_path / "solo_L2_phi-fdt-blos_20221030T041503_V02.fits"
    _write_cunitless_fits(f)

    smap = _load_map_patched_units(f)
    assert str(smap.meta["cunit1"]) == "arcsec"
    assert str(smap.meta["cunit2"]) == "arcsec"
    assert smap.data.shape == (64, 64)
    # the WCS must now be constructible for reprojection
    assert smap.wcs is not None


def test_load_map_works_on_cunitless_file(tmp_path):
    f = tmp_path / "solo_L2_phi-fdt-blos_20221030T041503_V02.fits"
    _write_cunitless_fits(f)
    # on sunpy < 8 the primary path succeeds with a warning; on sunpy >= 8
    # the fallback kicks in — either way this must return a usable map
    smap = load_map(f)
    assert smap.data.shape == (64, 64)
    assert smap.wcs is not None


def test_load_map_does_not_override_existing_units(tmp_path):
    f = tmp_path / "with_units.fits"
    data = np.ones((16, 16), dtype="float32")
    header = fits.Header()
    header["cunit1"] = "deg"
    header["cunit2"] = "deg"
    header["ctype1"] = "HPLN-TAN"
    header["ctype2"] = "HPLT-TAN"
    header["crpix1"] = 8.5
    header["crpix2"] = 8.5
    header["cdelt1"] = 0.001
    header["cdelt2"] = 0.001
    header["crval1"] = 0.0
    header["crval2"] = 0.0
    header["date-obs"] = "2022-10-30T04:15:03"
    fits.writeto(f, data, header=header, overwrite=True)

    smap = _load_map_patched_units(f)
    assert str(smap.meta["cunit1"]) == "deg"


if __name__ == "__main__":
    import tempfile

    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        test_patched_loader_fills_units(tmp)
        test_load_map_works_on_cunitless_file(tmp)
        test_load_map_does_not_override_existing_units(tmp)
    print("All tests passed.")
