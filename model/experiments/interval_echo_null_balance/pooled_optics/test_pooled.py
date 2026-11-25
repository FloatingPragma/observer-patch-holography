"""Exercise the complete written-field effects and the exterior current read."""

from fractions import Fraction
import json
from pathlib import Path

import numpy as np
import pytest

from . import circulation, measure


@pytest.fixture(scope="module")
def result():
    return measure.run()


def test_exact_all_source_effects_and_common_current(result):
    assert result["words"] == 458
    assert result["common_current"]["value"] == "53/916"
    assert len(result["all_rays"]) == 30
    assert sum(Fraction(r["difference"]) != 0 for r in result["all_rays"]) == 21
    assert sum(Fraction(r["midpoint_gap"]) > 0 for r in result["all_rays"]) == 11
    assert result["suppressed_ring_all_quartic_contrasts_zero"]


def test_exact_bound_and_mixture_switch(result):
    ray = next(r for r in result["all_rays"] if tuple(r["ray"]) == (2, 4))
    delta, gap = map(Fraction, (ray["difference"], ray["midpoint_gap"]))
    assert delta == Fraction(123793, 1374)
    assert gap == Fraction(16873, 916)
    (switch,) = result["selected_switches"]
    p, drop = map(Fraction, (switch["positive_weight"], switch["downward_slope_jump"]))
    assert (
        p == Fraction(54839, 71712)
        and drop == Fraction(35856, 229)
        and switch["port"] == 11
    )
    assert gap == drop * min(p, 1 - p) / 2


def test_faithful_signed_direction(result):
    row = result["full_rank_lottery"]
    assert row["minimum_eigenvalue"] == "1/7328"
    assert row["midpoint_gap"] == "43907/2748"
    assert row["initial_slope"] == "218015/1832"
    assert row["final_slope"] == "-50905/1832"
    assert row["switches"][0]["positive_weight"] == "1901/2430"


def test_native_unitary_replay(result):
    direct = result["direct_unitary_overlap"][-1]
    assert direct["dwell"] == 1 / 256
    assert direct["vacuum_dwell"] == 1 / 384
    assert abs(direct["calibrated_curvature_difference_over_t4"] - 123793 / 1374) < 0.45
    assert (
        abs(direct["calibrated_curvature_midpoint_gap_over_t4"] - 16873 / 916) < 0.145
    )
    assert abs(direct["difference_over_t4"] - 123793 / 1374) < 0.5
    assert abs(direct["midpoint_gap_over_t4"] - 16873 / 916) < 0.15
    expected = json.loads(Path(__file__).with_name("results.json").read_text())
    assert expected["all_rays"] == result["all_rays"]


def test_exterior_current_for_all_local_density_matrices():
    result = circulation.run()
    for row in result["boards"]:
        assert row["currents_to_optical_gains"] == ["sqrt(3)/96", "-sqrt(3)/96", "0"]
    replay = result["finite_mixed_320"]
    assert abs(replay["independent_current"]) > 0.1
    controls = replay["probabilities_difference_over_t4"]
    gain = np.sqrt(3) / 96 * replay["independent_current"]
    np.testing.assert_allclose(controls["native"][:2], [gain, -gain], atol=5e-4, rtol=0)
    np.testing.assert_allclose(controls["exchange_only"], 0, atol=5e-12, rtol=0)
