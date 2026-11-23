"""Recover the held forecasts from archived clicks and check the refinement bounds."""

import json
from itertools import combinations

import numpy as np
import pytest

from model.executable_model.refinement import DEPOSIT_SHA256
from .measure import HERE, instruments


def decode(counters, meter, shots):
    readings = np.zeros((20, 12))
    offset = 0
    for pole, (ports, effect) in enumerate(zip(meter.stars, meter.effects, strict=True)):
        weights = [np.diag(effect)[list(ports)]]
        weights += [np.array([effect[a, b], -effect[a, b]]) for a, b in combinations(ports, 2)]
        for weight in weights:
            width = 20 * len(weight) + 1
            setting = counters[offset:offset+width]
            assert setting.sum() == shots
            assert (setting >= 0).all()
            readings[:, pole] += setting[:-1].reshape(20, -1) @ weight / shots
            offset += width
    assert offset == len(counters)
    return readings


@pytest.mark.parametrize("suffix, expected_best", [
    ("", [.936741, .297692, .142175]),
    ("_holdout_3", [.080197]),
    ("_holdout_4", [.140832]),
    ("_holdout_4_precision", [.140832]),
])
def test_archived_coarse_observers_recover_each_held_difference(suffix, expected_best):
    receipt = json.loads((HERE / f"results{suffix}.json").read_text())
    if suffix:
        level = int(suffix.split("_")[2])
        assert receipt["deposited_board_sha256"] == DEPOSIT_SHA256[level]
    meter = instruments()
    assert meter.slow_error < 1e-10
    assert len(receipt["observations"]) == (3 if suffix else 12)
    with np.load(HERE / f"counts{suffix}.npz") as archive:
        for row in receipt["observations"]:
            q = decode(archive[row["count_key"]], meter, receipt["shots_per_setting"])
            diff = q[:, meter.held] - q[:, meter.trained] @ meter.forecast.T
            np.testing.assert_allclose(diff, row["sampled_held_differences"], atol=1e-10)
            if suffix.endswith("_precision") and row["law"] == "native":
                contrast = q - q.mean(axis=1, keepdims=True)
                ratio = np.sqrt(np.mean(diff**2) / np.mean(contrast**2))
                np.testing.assert_allclose(ratio, .76170, atol=1e-5)
                assert row["max_sampled_held_difference_in_standard_errors"] > 14.9
            assert row["sum_worst_held_absolute_error"] < row["angular_bound"] + 1e-10
            if row["source"] == "smooth":
                assert row["native_kinetic_energy"] <= 5 / (16*row["cells_per_parent"])
            if row["law"] == "input":
                assert row["angular_signal_rms"] < 1e-12
                assert row["relative_target_pass"] is None
            else:
                assert row["relative_target_pass"] is False
    smooth = [row for row in receipt["observations"] if row["law"] == "native" and row["source"] == "smooth"]
    np.testing.assert_allclose([row["best_all_direction_relative_error"] for row in smooth],
                               expected_best, atol=1e-6)
