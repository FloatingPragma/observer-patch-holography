from itertools import combinations
import json
import math
from pathlib import Path

import pytest

from .clock_ball import measure, two_flag_capacity


def test_any_two_load_word_has_an_independent_two_comparison_capacity_reader():
    for m in (2, 3, 5, 6, 12):
        for k in (0, 1, 2):
            both = sum(i < k and j < k for i, j in combinations(range(m), 2)) / math.comb(m, 2)
            assert two_flag_capacity(m, k/m, both) == pytest.approx(math.log(math.comb(m,k)))


def test_a_matched_clock_ball_has_a_correlation_transfer_at_identical_one_site_laws():
    result = json.loads(json.dumps(measure()))
    archive = json.loads(Path(__file__).with_name("clock_ball.json").read_text())
    # Counts and assignments are portable. Avoid pinning BLAS transcendental rounding.
    for key in ("native_gates", "hand_hops_to_stations", "independent_interior_reference_counts"):
        assert result[key] == archive[key]
    assert result["retained_ball_word_records"] == 55
    assert result["total_word_records"] == 66
    assert len(result["interior"]) == 6
    assert len(result["shell"]) == 5
    for name, arm in result["arms"].items():
        stored = archive["arms"][name]
        for key in ("one_site_counts", "shell_pair_counts", "independent_single_shell_occupied",
                    "signed_leg_count"):
            assert arm[key] == stored[key]
        assert arm["volume_vacancy_units"] == 9 and arm["shell_vacancy_units"] == 4
        assert arm["capacity_from_pair_flags"] == pytest.approx(arm["capacity_nats"])
    assert result["arms"]["spread"]["one_site_occupied_probability"] == result["arms"]["split"]["one_site_occupied_probability"]
    contrast = result["spread_minus_split"]
    assert contrast["two_rim_flags_probability"] == pytest.approx(1/20)
    assert contrast["capacity_nats"] == pytest.approx(math.log(2/5)/2)
    assert abs(contrast["flag_reconstruction_nats"] - contrast["capacity_nats"]) < 0.01
    assert abs(contrast["interior_entropy_first_variation_nats"] + contrast["capacity_nats"]) < 0.01
