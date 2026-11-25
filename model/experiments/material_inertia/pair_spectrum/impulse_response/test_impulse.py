import copy

import pytest

from . import replay


def test_recoil_of_the_same_rest_pair_and_control():
    record = replay.verify()
    for source in record["preparations"]:
        assert source["controls"]["native"] > 0.99
        assert source["controls"]["contrast_removed"] < 0.32
        assert source["controls"]["uniform_common_response"] > 0.01
        for label, sign in (("common", -1), ("signed_load", 1)):
            case = source["cases"][label]
            low, high = case["derivative_intervals"]["all_ports"]
            assert min(sign * low, sign * high) > 0
            assert case["derivative_intervals"]["within_cluster"][0] > 0
            for key, goal in (("free", sign), ("band_selected", 1)):
                panel = case["finite"][key]
                estimate = panel["signed_position_response"]
                assert goal * estimate > panel["halfwidth_99_simultaneous"]
                assert (
                    abs(estimate - panel["nominal_response"])
                    < panel["halfwidth_99_simultaneous"]
                )
                for shot in panel["records"]:
                    assert shot["enrolled"] == shot["failed_preparation"] + shot[
                        "failed_band_filter"
                    ] + sum(shot["counts"])


def test_custody_rejects_one_reassigned_count_at_fixed_total():
    source = replay.archived("results.json", replay.RESULT_SHA256)
    law = replay.archived("acquisition.json", replay.POPULATIONS_SHA256)[0]
    shot = copy.deepcopy(
        source["preparations"][0]["cases"]["common"]["finite"]["free"]["records"][0]
    )
    shot["counts"][0] += 1
    shot["counts"][1] -= 1
    with pytest.raises(ValueError, match="Seeded impulse counts differ"):
        replay.shot_record(shot, law, 1.0)
