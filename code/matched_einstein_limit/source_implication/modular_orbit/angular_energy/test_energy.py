import json

from .measure import run


def test_calibrated_native_energies_and_held_triplet_contrast():
    result = run()
    stationary = result["stationary_bounds"]
    assert json.loads(stationary["rim_excess_over_constant_interior"])[0] > 0.345
    assert json.loads(stationary["modular_witness"])[0] > 0.067
    assert result["rounded_stationary_source_commutator"] < 1e-14
    assert (
        result["rounded_stationary_charge_change"]
        < json.loads(stationary["log_change_norm_bound"])[1]
    )
    trial = result["stationary_triplet"]
    assert abs(trial["expectation"] - trial["difference"]) < 4 * trial["standard_error"]
    assert result["stationary_waits"]["two"]["fair_waits"] == 153
    intervals = result["true_bounds"]
    assert json.loads(intervals["modular_witness"])[0] > 0.06786
    assert (
        json.loads(intervals["rim_weight"])[0]
        > json.loads(intervals["spoke_weight"])[1]
        > 0
    )
    assert json.loads(intervals["total_energy_weight"])[1] < 0
    assert result["held_stars_operator_max"] < 2e-10
    assert result["centered_stars_operator_max"] < 2e-10
    assert result["unchanged_weights_untwirled_error"] > 0.001
    assert result["exchange_only_control"]["filter_effect_error"] < 1e-12
    assert result["exchange_only_control"]["uniform_column_norm"] < 1e-14
    arm = result["equal_energy_triplet"]
    assert abs(arm["ordinary_energy_difference"]) < 1e-11
    for label, read in arm["sampled"].items():
        assert (
            abs(read["difference"] - arm["expectation"][label])
            < 4 * read["standard_error"]
        )
        for probe in ("axial", "transverse"):
            assert all(
                sum(item["counts"]) == arm["shots_per_setting"] for item in read[probe]
            )
    assert arm["expectation"]["modular"] < -0.23
    assert arm["expectation"]["spokes"] + arm["expectation"]["rim"] > 3.5
