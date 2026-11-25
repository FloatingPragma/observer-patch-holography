from fractions import Fraction
from functools import lru_cache
from itertools import combinations

import numpy as np

from . import measure as tracer


@lru_cache(maxsize=1)
def run():
    return tracer.measure()


def test_receiver_and_ordinary_calibration():
    receipt = run()
    active = receipt["active"]
    assert active["joint_records"] == 660
    assert active["acquisition"]["verified_routed_transpositions"] == 25
    for hand, exposure in zip(receipt["independent_hand"], [.25, .5, 1.], strict=True):
        assert abs(hand["phase_radians"]-exposure) < .002
    raw = active["exposures"][1]
    np.testing.assert_allclose(raw["held_predicted_nats"], [-.88992533, -.82809624], atol=1e-7)
    assert all(hi < receipt["true_held_completion_nats"] for _, hi in raw["sampled_99_held_nats"])
    assert receipt["exchange_control"]["exposures"][1]["held_predicted_nats"][0] > 0
    for frame in active["exposures"]:
        f, r = frame["inward_and_scaled_outward_at_uniform_reference"]
        assert abs(f-r) < 2e-12


def test_exact_low_orders_and_distinct_equal_capacity_words():
    h, basis, rim, outside = tracer.source(8)
    first, fourth = [tracer.meters(k, basis, outside) for k in tracer.taylor(h)]
    for word in combinations(rim, 2):
        assert first[word][0] == 3
    assert Fraction(fourth[4, 9][0], 192) == -175
    assert Fraction(fourth[2, 9][0], 192) == Fraction(-6559, 96)
    assert run()["active"]["exact_t4_calibration_and_held"][1][0] == "581/24"
    assert run()["active"]["randomized_transmission"]["t6_determinant"] == "42811/192"


def test_exterior_only_prediction_with_unused_outward_exposure():
    receipt = run()
    read = receipt["active"]["randomized_transmission"]
    assert read["class_sizes"] == [15, 30, 10]
    matrix = np.array(read["class_calibration_matrix"])
    np.testing.assert_allclose(np.asarray(read["nats_per_record"]) @ matrix, np.log([1, 5, 10]))
    lo, hi = read["simultaneous_99_charge_interval"]
    assert lo < receipt["true_held_completion_nats"] < hi < 0
    for bank in read["arms"].values():
        assert abs(bank["predicted_outward_rate"]-bank["observed_outward_rate"]) < read[
            "simultaneous_99_outward_discrepancy_bound"]
    counted = read["counted_calibration"]
    lo, hi = counted["simultaneous_99_charge_interval"]
    assert lo < receipt["true_held_completion_nats"] < hi < 0
    assert counted["least_singular_value"] > counted["matrix_error_bound"]
    for bank in counted["arms"].values():
        assert abs(bank["predicted_outward_rate"]-bank["observed_outward_rate"]) < counted[
            "simultaneous_99_outward_discrepancy_bound"]
