import json
from itertools import combinations


from . import station_control, station_region, station_region_counts


def test_radar_region_number_weight_control():
    data = station_control.run()
    assert data["ports"] == (0, 1, 38, 39) and data["fixed_readout_count"] == 26
    assert data["augmented_rank"] == 27 and data["drive_derivative_exact"] == "1/190"
    assert data["region_occupations_derivatives"] == [0, 0]
    assert json.loads(data["lottery_slack_lower"])[0] > 0.03


def test_four_station_untagged_entropy_and_center_preserving_arm():
    data = station_region.run()
    original = data["regions"]["stations"]["fair"]["original_equal_work"]
    assert len(original["spectrum"]) == 9
    assert min(json.loads(v)[0] for v in original["spectrum"]) > 0.001
    assert json.loads(original["slope"])[0] > 0.0002
    assert json.loads(original["finite_change_0_to_0p01"])[0] > 1.5e-6
    assert json.loads(original["extra_conditional_finite_change_0_to_0p01"])[0] > 1e-6
    fixed = data["regions"]["stations"]["fair"]["fixed_all_number_weights"]
    for key in ("slope", "finite_change_0_to_0p01"):
        lo, hi = json.loads(fixed[key])
        assert lo * hi > 0
    witness = data["fixed_number_midpoint_witness"]
    assert (
        json.loads(witness["minimum_midpoint_entropy_deficit_at_symmetric_0p01"])[0] > 0
    )
    lo, hi = json.loads(witness["fair_mode_slope"])
    assert lo * hi > 0
    for q in ("vacancy", "excess"):
        assert (
            len(data["regions"]["stations"][q]["original_equal_work"]["spectrum"]) == 5
        )


def test_four_station_detector_modes_and_sampled_positive_contrast():
    # Check every pair analyzer against the X/Y sign convention.
    effects = station_region.analyzers()
    for (i, j), (x, y) in zip(
        combinations(range(4), 2), zip(effects[4::2], effects[5::2])
    ):
        assert x[i, j] == 1 and complex(y[i, j]) == -1j
    result = station_region_counts.run()
    assert (
        result["source_endpoints"] == [-0.01, 0.01]
        and not result["charge_tag_retained"]
    )
    assert json.loads(result["signed_two_arm_entropy_difference_nats_99pct"])[0] > 4e-6
    extra = result["fixed_all_number_weights"]
    assert (
        extra["source_endpoints"] == [0, 0.01] and extra["trials_per_setting"] == 2**72
    )
    assert json.loads(extra["signed_two_arm_entropy_difference_nats_99pct"])[1] < -3e-8
    for arm in result["arms"].values():
        settings = arm["setting_counts"]
        assert len(settings) == 13 and all(
            sum(row) == 2**60 for row in settings.values()
        )
        assert len(settings["loads"]) == 9
