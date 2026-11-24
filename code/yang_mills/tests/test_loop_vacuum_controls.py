import json
import sys
from fractions import Fraction as F
from math import exp, log
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from loop_vacuum_controls import arbitrary_coupling, run  # noqa: E402


def test_full_haar_response_and_clock():
    path = Path(__file__).parents[1] / "receipts" / "loop_vacuum_controls.json"
    deposited = json.loads(path.read_bytes())
    actual = run()
    np.testing.assert_allclose(
        actual.pop("haar_quadrature"), deposited.pop("haar_quadrature"), atol=2e-14
    )
    for s in actual["magnetic"]:
        lo, hi = map(F, s["rectangle_mean_interval"])
        b = F(s["coupling"])
        assert 0 < lo < F(7, 216) * b * b < hi
        assert F(s["rectangle_variance_floor"]) > F(1, 5)
    assert actual["magnetic"] == deposited["magnetic"]
    assert actual["cofinal_strong_coupling"][-1]["physical_gap_floor"] > 1e6
    assert arbitrary_coupling(4)["log_physical_gap_floor"] < -1e7
    for row in actual["electric_square_correlation"]:
        t = row["time"]
        # Each of the two one-loop sectors has weight 1/16.
        assert (
            abs(row["value"] - (2 * exp(-8 * t) / 16 + 7 * exp(-16 * t) / 48)) < 1e-14
        )
    assert actual["rectangle_tilings"] == deposited["rectangle_tilings"]


def test_all_coupling_density_comparison():
    for j in range(5):
        row = arbitrary_coupling(j)
        alpha, b = map(F, (row["electric_clock"], row["magnetic_ratio"]))
        # The full two-sided ground-amplitude comparison is squared to compare densities.
        log_ratio = 8 * float(b) * row["faces"] + row["links"] * log(5 / 3)
        np.testing.assert_allclose(
            row["log_physical_gap_floor"], log(3 * float(alpha)) - 2 * log_ratio
        )
        np.testing.assert_allclose(
            row["log_physical_gap_ceiling"], log(4 * float(alpha)) + 2 * log_ratio
        )
