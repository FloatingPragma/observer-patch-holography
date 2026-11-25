import json

import numpy as np

from ...action_response.four_ticks import parse_interval
from . import measure, strain


def test_positive_duration_impulse_matches_signed_force_calculation():
    report = measure.run()
    enclosure = json.loads((measure.HERE / "certificate.json").read_text())["charges"]
    for q, surveys in report["force"].items():
        for name, methods in surveys.items():
            for arm, true in methods["commutator"].items():
                approx = methods["pulse_1_1024"][arm]
                assert abs(approx - true) < 0.0021
                assert true in parse_interval(enclosure[q]["observables"][name][arm])
        for arm in ("ordinary", "pairing"):
            assert (
                surveys["right_front"]["commutator"][arm]
                * surveys["right_hops"]["commutator"][arm]
                < 0
            )


def test_equal_strain_and_different_native_work():
    report = strain.run()
    ranges = json.loads((measure.HERE / "certificate.json").read_text())["charges"]
    assert report["leading_flight_comparison"]["compared_integer_entries"] == 6016
    index = report["independent_material_standard"]["records"]
    assert (
        index["full_seam_0.75"]["material_index"][0]
        < 8
        < index["full_seam_0.75"]["material_index"][1]
    )
    assert index["exchange_only_0.75"]["material_index"][0] > 10
    for q, data in report["charges"].items():
        clocks = data["neighbor_clocks"]["settings"]
        for i in (0, 2):
            assert clocks["bridge_stretch"]["frequency_intervals"][i][0] > 0.99
            assert clocks["common_lapse"]["frequency_intervals"][i][1] < 0.76
        work = data["means"]
        for key in ("reference", "ordinary", "pairing"):
            np.testing.assert_allclose(
                work["full_seam"][key] - work["exchange_only"][key],
                8 * work["contact_occupation"][key],
                atol=1e-13,
            )
        ratio = (
            work["exchange_only"]["pairing"] / work["exchange_only"]["ordinary"]
        ) / (work["full_seam"]["pairing"] / work["full_seam"]["ordinary"])
        assert ratio in parse_interval(ranges[q]["ordinary_calibrated_work_ratio"])
        assert 1 not in parse_interval(ranges[q]["ordinary_calibrated_work_ratio"])
