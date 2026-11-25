import numpy as np
import pytest

from .measure import WEIGHTS, controls, flat_volume, marked_volume, native_relay


def test_three_clock_design():
    clocks = np.array([[1, 0], [5 / 3, 4 / 3], [5 / 3, -4 / 3]])
    np.testing.assert_allclose(clocks.T @ np.diag(WEIGHTS) @ clocks, np.ones((2, 2)))
    assert abs(WEIGHTS.sum()) < 1e-15


@pytest.mark.parametrize("boost", [1 / 3, 1.0, 3.0])
def test_flat_subdiamond_off_center(boost):
    for start, stop in ((-0.125, 0.125), (-0.125, 0), (0, 0.125)):
        assert marked_volume(start, stop, boost, 0.0, 0.0) == pytest.approx(
            flat_volume(stop - start)
        )


def test_density_cancellation_and_curvature():
    for panel in controls()["panels"]:
        flat = panel["conformal_quadratic"] == 0
        for read in panel["reads"]:
            curvature = read["exact_center_ricci"]
            if flat:
                assert read["mean_read"] == pytest.approx(0.0, abs=1e-9)
                assert read["full_only_mean"] == pytest.approx(
                    0.75 * (1 + read["null_x"]) ** 2, abs=1e-9
                )
            if read["duration"] == 0.03125:
                assert read["mean_read"] == pytest.approx(curvature, abs=0.002)
            assert (
                abs(read["mean_read"] - read["poisson_read"])
                < 5 * read["poisson_standard_error"]
            )


def test_native_relay_uses_all_ports_and_bounds_the_rest_scaling():
    record = native_relay()
    assert record["ports"] == 960
    for sector in record["sharp_support"]:
        assert sector["nonidentity_blocks"] == 4
    for panel in record["panels"]:
        lo, hi = panel["ratio"]
        assert 2 <= lo <= hi < 11
