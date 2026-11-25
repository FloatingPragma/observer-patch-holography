from fractions import Fraction
from itertools import combinations_with_replacement
import json
import math
from pathlib import Path

import pytest

from .measure import measure, sharp_arc_bound


def test_native_forks_and_positive_half_exchange_panel():
    archive = json.loads(Path(__file__).with_name("results.json").read_text())
    result = measure()
    distances = {row["origin"]: row["distances"] for row in result["available_geodesic_rows"]}
    assert result["geodesic_forks"] == archive["geodesic_forks"]
    assert len(result["geodesic_forks"]) == 5
    for fork in result["geodesic_forks"]:
        i = fork["inner"]
        q, r = fork["rim"]
        assert distances[0][i] == distances[q][i] == distances[r][i] == 1
        assert distances[0][q] == distances[0][r] == 2
        assert distances[q][r] == fork["rim_separation"] == 1
    meter = result["clocked_aperture"]
    assert meter["records"] == archive["clocked_aperture"]["records"]
    assert len(meter["records"]) == 55
    classes = {0: {"0"}, 1: {"5/2"}, 2: {"59/16", "4"}}
    for k, values in classes.items():
        assert {row["mean_lit_panels"] for row in meter["records"]
                if row["rim_occupancy"] == k} == values
    for arm in meter["measured_pair_classes"].values():
        lower, upper = arm["simultaneous_99_interval_nats"]
        assert math.log(10) < lower < arm["predicted_pair_capacity_nats"] < upper


def test_unknown_arc_width_offsets_and_positive_mixture():
    result = sharp_arc_bound()
    assert [(row["near"], row["far"]) for row in result["consecutive_classes"]] == [
        (0, 0), (0, 0), (1, 0), (2, 1), (3, 3), (5, 5)]
    width = result["attaining_uniform_arc_length"]
    target = result["required_overlap_per_single"]
    errors = [math.log(5) * (overlap / width - target)
              for overlap in (width - 1, 2 * width - 5)]
    assert errors[0] == pytest.approx(-errors[1])
    assert errors[0] == pytest.approx(result["sharp_minimum_worst_pair_error_nats"])
    assert errors[0] > 0.19385
    probability = result["fitted_full_rim_lottery_probability"]
    single = probability + (1 - probability) / 5
    double = probability + 2 * (1 - probability) / 5
    assert math.log(5) * double / single == pytest.approx(math.log(10))


def test_radial_overlap_factorization_on_ordered_profiles():
    for c, b, a in combinations_with_replacement([Fraction(i, 8) for i in range(9)], 3):
        f = (a, b, c, c, b)
        near, far = (sum(f[i] * f[(i + offset) % 5] for i in range(5))
                     for offset in (1, 2))
        assert near - far == (b - c) * (2 * a - b - c)
        assert (near == far) == (b == c)
