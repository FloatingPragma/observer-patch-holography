"""Clock-only inversion, hardware disjointness, whole-walk routing and held signs."""

import json
from math import pi
import pytest
from . import clock_transfer as clock


def archived():
    return json.loads(clock.HERE.joinpath("clock_results.json").read_text())


def test_native_clock_training_and_no_aliased_selection():
    a = archived()
    regenerated = clock.run()
    assert regenerated["material_hands"] == a["material_hands"]
    assert a["joint_clock_failure_probability_upper"] < 0.01 / 33
    assert len(a["material_hands"]) == 2 and all(
        h["degree"] == 5 for h in a["material_hands"]
    )
    for arm in clock.ARMS:
        for name, tau in clock.design().items():
            broad = [0.0, 2 * pi + 1e-12 if name == "full exchange" else pi + 1e-12]
            for new, row in zip(
                regenerated["duration_readings"][arm][name]["stages"],
                a["duration_readings"][arm][name]["stages"],
                strict=True,
            ):
                n = row["repeats"]
                shots = row["right"] + row["left"]
                assert new["right"] == row["right"] and new["left"] == row["left"]
                assert row["probability_tickets"] == clock.canonical_probability(tau, n)
                broad = clock.refine(broad, n, row["right"], shots)
                assert broad[0] <= float(tau.mid) <= broad[1]
            assert broad[1] - broad[0] < (3e-14 if name == "full exchange" else 2.2e-10)
    with pytest.raises(ArithmeticError, match="compatible branches"):
        clock.refine([0.0, pi], 1021, clock.SHOTS // 2, clock.SHOTS)


def test_full_records_and_99_percent_held_predictions():
    x = archived()
    receipt = json.loads(clock.HERE.joinpath("results.json").read_text())
    assert len(x["panels"]) == len(receipt["panels"]) == 32
    for k, v in zip(x["panels"], receipt["panels"], strict=True):
        assert (k["name"], k["probe"], k["face"]) == (
            v["setting"]["name"],
            v["setting"]["probe"],
            v["setting"]["face"],
        )
        assert k["max_source_routing_pulses_per_visit"] in (101, 105)
        assert (
            0
            < k["source_routing_bias_per_release_upper"]
            == 2048 * k["source_record_total_variation_upper"]
            < 5.4e-6
        )
        assert 0 < k["exact_source_timing_bias_per_release_upper"] < 4.1e-6
        assert k["timing_bias_per_release_upper"] < 9.4e-6
        assert k["observed_deviation"] < k["combined99_half_width"] < 0.000129
        if k["name"] in ("matched", "blind"):
            assert abs(v["mean_per_release"]) > k["combined99_half_width"]
    assert x["largest_timing_bias_per_release"] == max(
        p["timing_bias_per_release_upper"] for p in x["panels"]
    )
