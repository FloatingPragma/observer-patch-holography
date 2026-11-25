import json

import numpy as np
import pytest

from . import certificate, measure


@pytest.fixture(scope="module")
def experiment():
    return measure.experiment()


def test_phase_channels_and_changed_law(experiment):
    np.testing.assert_array_equal(
        experiment["phase_arm_transfer"], [[-1 / 8, 0], [0, -9 / 64]]
    )
    np.testing.assert_array_equal(
        experiment["exchange_only_phase_transfer"], np.zeros((2, 2))
    )
    assert experiment["blank_three_hop_standard"] == [13, 792]


def test_held_transfer(experiment):
    left, right = experiment["channels"]
    assert (
        0.96
        < left["calibrations"]["full"]["held"]["unpairing"]["observed_over_predicted"]
        < 0.98
    )
    assert (
        0.53
        < right["calibrations"]["full"]["held"]["unpairing"]["observed_over_predicted"]
        < 0.54
    )
    for row in (left, right):
        assert (
            row["threshold_slopes"]["ordinary"] > 0
            and row["threshold_slopes"]["phase"] < 0
        )


def test_nonzero_count_responses(experiment):
    replay = experiment["finite_count_replay"]
    assert len(replay["rows"]) == 8
    for row in replay["rows"]:
        assert row["phase_slope_count"] > row["simultaneous_99pct_halfwidth"]
        assert (
            abs(row["phase_slope_count"] - row["phase_slope_population"])
            < row["simultaneous_99pct_halfwidth"]
        )
        assert all(
            n + m == replay["shots_per_setting"]
            for n, m in zip(row["click"], row["miss"], strict=True)
        )


def test_reenclose_true_station():
    import ast

    result = certificate.run()
    for row in result["apertures"]:
        slope = row["slope"]
        assert (
            ast.literal_eval(slope["ordinary"])[0]
            > 0
            > ast.literal_eval(slope["phase"])[1]
        )
        for name in ("full", "collar"):
            assert ast.literal_eval(row["transfer_ratio"][name]["unpairing"])[1] < 1
    expected = json.loads((measure.HERE / "certificate.json").read_text())
    assert result["phase_weights"] == expected["phase_weights"]
