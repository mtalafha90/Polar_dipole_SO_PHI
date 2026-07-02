"""CROTA2 handling: HMI images are stored with ~180 deg camera rotation, so
ignoring the keyword flips north/south in the native-geometry path. A
180-degree-rotated synthetic instrument must recover the same solar-frame
fields as the unrotated one.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.append(str(Path(__file__).resolve().parents[1]))

from solar_pipeline.geometry import rotate_offsets
from solar_pipeline.pipeline import compute_native_disk_fields
from tests.test_milestone import FakeMap, make_synthetic_disk, MU_MIN


def test_rotate_offsets_identity_and_180():
    dx = np.array([1.0, 0.0, -2.0])
    dy = np.array([0.0, 3.0, 1.0])
    rx, ry = rotate_offsets(dx, dy, 0.0)
    assert np.allclose(rx, dx) and np.allclose(ry, dy)
    rx, ry = rotate_offsets(dx, dy, 180.0)
    assert np.allclose(rx, -dx) and np.allclose(ry, -dy)


def test_rotated_instrument_recovers_same_fields():
    upright = make_synthetic_disk(b0_deg=10.0, l0_deg=180.0)

    # the same scene imaged by a camera rotated 180 deg: the image array is
    # flipped in both axes and the header records crota2 = 180
    rotated_data = upright.data[::-1, ::-1].copy()
    rotated_meta = dict(upright.meta)
    rotated_meta["crota2"] = 180.0
    # crpix of the synthetic disk is at the exact array centre, so the pixel
    # flip i -> n-1-i maps offsets to their negatives with the same crpix
    rotated = FakeMap(rotated_data, rotated_meta)

    f_up = compute_native_disk_fields(upright, disk_fraction=0.98, mu_min=MU_MIN, alpha=1.0)
    f_rot = compute_native_disk_fields(rotated, disk_fraction=0.98, mu_min=MU_MIN, alpha=1.0)

    # solar-frame quantities must agree once the rotated arrays are flipped back
    for key in ("br", "lat", "lon", "mu"):
        a = f_up[key]
        b = f_rot[key][::-1, ::-1]
        both = np.isfinite(a) & np.isfinite(b)
        assert both.sum() > 1000
        assert np.allclose(a[both], b[both], atol=1e-8), key
    # and each map's own valid fraction matches
    assert f_up["valid"].sum() == f_rot["valid"].sum()


def test_missing_wcs_gives_helpful_error():
    bare = FakeMap(np.ones((16, 16)), {"date-obs": "2022-10-30T04:15:03"})
    try:
        compute_native_disk_fields(bare, disk_fraction=0.98, mu_min=MU_MIN, alpha=1.0)
        raise AssertionError("should have raised")
    except RuntimeError as exc:
        assert "fix-headers" in str(exc)


if __name__ == "__main__":
    test_rotate_offsets_identity_and_180()
    test_rotated_instrument_recovers_same_fields()
    test_missing_wcs_gives_helpful_error()
    print("All tests passed.")
