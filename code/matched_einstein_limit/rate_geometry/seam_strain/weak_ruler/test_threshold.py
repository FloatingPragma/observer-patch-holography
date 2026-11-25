import json
from math import factorial
from pathlib import Path

import numpy as np

from . import comparison, threshold

HERE = Path(__file__).parent


def test_two_leg_leading_coefficients_match_independent_integer_hamiltonian_powers():
    saved = json.loads((HERE / "threshold.json").read_text())
    monitored = json.loads((HERE / "flight.json").read_text())
    order = np.array(saved["hop_order"])
    for q in ("1", "-1"):
        row = monitored["onsets"][q]
        np.testing.assert_array_equal(order, row["order"])
        for name, values in saved["direct_coefficients"].items():
            actual = np.asarray(row["states"][name]["return_cumulative"])
            for d in (1, 2, 3):
                c0 = factorial(2 * d) ** 2 / (
                    2 ** (4 * d) * factorial(d) ** 4 * factorial(4 * d + 2)
                )
                pooled = (actual + actual.T) / (2 * c0)
                np.testing.assert_allclose(
                    pooled[order == d],
                    np.asarray(values)[order == d],
                    rtol=1e-10,
                    atol=2e-10,
                )


def test_same_source_threshold_links_and_fixed_metric_obstruction():
    fresh = threshold.run()
    saved = json.loads((HERE / "threshold.json").read_text())
    np.testing.assert_allclose(
        fresh["direct_lengths"], saved["direct_lengths"], atol=1e-11
    )
    assert fresh["certified_triangle"] == saved["certified_triangle"]
    assert fresh["triangle"]["excess"] > 0.254
    panels = fresh["orders"]
    assert [p["links"] for p in panels] == [2520, 3600, 4704, 6672]
    assert panels[0]["link_classes"] == {"blank": 2517, "marked": 2, "strain": 1}
    assert panels[0]["exceptional_links"] == {
        "marked": [[0, 1], [38, 39]],
        "strain": [[1, 38]],
    }
    assert all(p["minimum_coefficient"] > 0 for p in panels)
    assert fresh["fixed_one_hop_metric"]["retained_entropy_slope_lower"] == 0.059
    check = fresh["fixed_one_hop_metric"]["finite_checks"][-1]
    assert abs(check["strain_error"]) < 2e-13
    assert 0.059 < check["entropy_secant"] < 0.060
    four = panels[-1]
    assert min(four["eight_gram_eigenvalues"]) < -1
    for row, stored in zip(
        four["tetrahedron"]["finite_steps"],
        saved["orders"][-1]["tetrahedron"]["finite_steps"],
        strict=True,
    ):
        for arm in ("ordinary", "decorrelated"):
            assert row["arms"][arm]["bracketed"]
            assert abs(row["arms"][arm]["volume_error"]) < 1e-10
        np.testing.assert_allclose(
            row["held_residual"], stored["held_residual"], rtol=0.015, atol=1e-5
        )


def test_common_clock_identities_and_held_raw_readings():
    x = np.arange(8, dtype=float)
    d = np.abs(x[:, None] - x)
    g = np.exp(0.03 * x)
    symmetric = g[:, None] * g * d
    controls = comparison.identities(d, symmetric)
    assert controls["reciprocal_species"]["log_tetrad_error"] < 1e-14
    controls = comparison.identities(d, g[:, None] * d)
    assert controls["oriented"]["plus"]["log_cycle_error"] < 1e-14
    fresh = comparison.run()
    saved = json.loads((HERE / "comparison.json").read_text())
    for a, b in zip(fresh, saved, strict=True):
        assert (
            a["clock_tests"]["capped_range"]["oriented"]["plus"]["log_cycle_error"]
            > 0.020
        )
        assert (
            a["clock_tests"]["capped_range"]["reciprocal_species"]["log_tetrad_error"]
            > 0.005
        )
        assert a["six_station_obstruction"]["maximum_distance_change_lower"] > 0.302
        np.testing.assert_allclose(
            a["held_distance_error"], b["held_distance_error"], rtol=1e-12
        )
