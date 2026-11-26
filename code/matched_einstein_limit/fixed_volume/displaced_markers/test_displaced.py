import json
from fractions import Fraction as F
from pathlib import Path

import numpy as np
import pytest

from .measure import decode, measure, shifted_count_mean


def test_displaced_reads_remove_an_arbitrary_local_density_jet():
    radius = F(1, 4)
    # Exact second-order normal-coordinate records, with unequal signed curvatures
    # of the marker density and an independent linear term on each axis.
    m, scalar = F(7), F(-11)
    hessian = (F(-5), F(8), F(3))
    gradient = (F(2), F(-5), F(9))
    q = sum(hessian) / 10 - m * scalar / 30
    central = m + radius**2 * q
    half = m + radius**2 * q / 4
    displaced = [half + sign * g * radius + radius**2 * h / 2
                 for h, g in zip(hessian, gradient, strict=True) for sign in (-1, 1)]
    assert 18 * sum(displaced) / 6 + 22 * half - 40 * central == radius**2 * m * scalar
    assert decode(list(map(float, [central, half, *displaced])), float(radius)) == pytest.approx(scalar)
    flat_h = F(7)
    # Planar stencil: a hidden normal Hessian balances a curvature in every recorded jet.
    assert flat_h / 10 - m * (3 * flat_h / m) / 30 == 0


def test_quadratic_density_leaves_false_center_only_curvature_on_flat_space():
    for radius in (.5, .25, .125):
        gamma = .25
        central = 1 + 3 * gamma * radius**2 / 5
        half = 1 + 3 * gamma * radius**2 / 20
        shifted = shifted_count_mean(0, radius, gamma)
        assert shifted == pytest.approx(half + gamma * radius**2)
        assert -40 * (central - half) / radius**2 == pytest.approx(-4.5)
        assert decode([central, half, *[shifted] * 6], radius) == pytest.approx(0, abs=1e-11)


def test_independently_counted_area_and_blind_center_on_five_space_forms():
    receipt = json.loads(Path(__file__).with_name("results.json").read_text())
    fresh = measure()
    assert len(fresh["rows"]) == len(receipt["rows"]) == 24
    for stored, new in zip(receipt["rows"], fresh["rows"], strict=True):
        assert stored["counts"] == new["counts"]
        np.testing.assert_allclose(stored["expected_normalized_counts"], new["expected_normalized_counts"])
        lower, upper = new["counting_interval"]
        assert lower < new["mean_decoder"] < upper
        if new["density_quadratic"] == 0:
            assert new["expected_normalized_counts"][:2] == [1., 1.]
        if new["radius"] == 1 / 16:
            assert abs(new["matched_area_count_scalar"] - new["mean_decoder"]) < .006
            assert abs(6 * new["sectional_curvature"] - new["mean_decoder"]) < .0051
