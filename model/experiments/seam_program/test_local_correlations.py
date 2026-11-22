import json
from itertools import product
from pathlib import Path

import numpy as np
import pytest

from . import local_correlations, scalar_vacuum
from .clocked_field import eigenphase

HERE = Path(__file__).resolve().parent


def test_categorical_exclusion_and_independent_four_bits():
    weights = np.array([2, 3, 5, 7, 11, 13, 17, 19])
    exclusive = local_correlations.moments(weights, np.eye(8)[:, :4])
    assert exclusive["connected_fourth"] == pytest.approx(
        -6 * np.prod(exclusive["means"])
    )
    independent = local_correlations.moments(
        np.ones(16), list(product((0, 1), repeat=4))
    )
    assert independent["connected_fourth"] == 0
    with pytest.raises(ValueError):
        local_correlations.moments([1, -1], np.zeros((2, 4)))


def test_fresh_full_sector_and_recorded_moments():
    fresh = local_correlations.run()
    saved = json.loads((HERE / "local_correlations.json").read_text())
    assert fresh["source_sha256"] == saved["source_sha256"]
    assert fresh["scalar_shots"] == saved["scalar_shots"]
    for key, values in fresh["scalar_position"].items():
        np.testing.assert_allclose(
            values, saved["scalar_position"][key], rtol=1e-12, atol=1e-14
        )
    for row, old in zip(fresh["pair_detector"], saved["pair_detector"], strict=True):
        assert (row["cells"], row["arm"], row["quiet"]) == (
            old["cells"],
            old["arm"],
            old["quiet"],
        )
        for source in ("recorded", "population"):
            for key, values in row[source].items():
                np.testing.assert_allclose(
                    values, old[source][key], rtol=1e-10, atol=1e-12
                )
        resolved = abs(row["population"]["connected_fourth"]) > 0.001
        assert resolved == (row["arm"] == "sharp")


def test_vacuum_octant_sum_against_all_signed_modes():
    length, n = 4, 8
    momenta = (
        2 * np.pi * np.array(list(product(range(-n // 2, n // 2), repeat=3))) / length
    )
    momenta = momenta[np.any(momenta != 0, axis=1)]
    phase = eigenphase(np.column_stack((momenta * length / n, np.zeros(len(momenta)))))
    omega = phase * n / length
    source = scalar_vacuum.POINTS
    expected = []
    for a, b in scalar_vacuum.PAIRS:
        dt, *dx = source[a] - source[b]
        expected.append(
            np.sum(
                np.cos(momenta @ dx)
                * np.exp(-abs(dt) * omega)
                / (2 * omega * length**3)
            )
        )
    np.testing.assert_allclose(
        scalar_vacuum.periodic_pairs(source, length, n),
        expected,
        rtol=1e-12,
        atol=1e-16,
    )
    saved = json.loads((HERE / "scalar_vacuum.json").read_text())
    np.testing.assert_allclose(
        scalar_vacuum.flat_pairs(source), saved["flat_pairs"], rtol=1e-13
    )
    fresh = scalar_vacuum.run()
    for row, actual in zip(saved["regulators"], fresh["regulators"], strict=True):
        for key, value in actual.items():
            np.testing.assert_allclose(value, row[key], rtol=1e-12, atol=1e-16)
        assert row["wick_fourth"] == pytest.approx(
            scalar_vacuum.wick(row["pairs"]), rel=1e-13, abs=0
        )
    large = saved["regulators"][-1]
    assert large["max_pair_relative_error"] < 0.05
