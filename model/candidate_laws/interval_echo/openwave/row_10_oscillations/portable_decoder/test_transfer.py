import copy
import json

import pytest

from .measure import HERE, TRAINING, calibrate, decode, score


def test_single_training_aperture_and_missing_arrivals():
    calibration = json.loads(TRAINING.read_text())
    selected = [
        row
        for row in calibration["banks"]
        if row["key"][:2] == [12, 256] and row["key"][-1] == "training"
    ]
    assert calibrate(calibration, (0, 1, 2)) == calibrate(
        {"banks": selected}, (0, 1, 2)
    )
    empty = {(j, panel): [0, 0, 0, 10, 0, 0, 0] for j in (1, 2) for panel in (-1, 1)}
    result = score(empty, calibrate(calibration, (1, 2)), (1, 2), 10)
    assert result["limits"] is None and result["conditional_accuracy"] is None
    assert result["extra_correct_per_trillion"] == 0


def test_all_launches_and_held_sign_controls():
    data = json.loads((HERE / "counts.json").read_text())
    training = json.loads(TRAINING.read_text())
    result = decode(data, training)
    n = data["starts"]
    for row in result["rows"]:
        for key, choices in [("three", 3), ("two", 2)]:
            sample = row[key]
            correct = sum(c["correct"] for c in sample["cohorts"])
            pairs = sum(c["pairs"] for c in sample["cohorts"])
            gain = (choices * correct - pairs) * 10**12 / (2 * choices**2 * n)
            assert sample["extra_correct_per_trillion"] == pytest.approx(gain)

    def wins(mode, key, cut):
        return sum(
            r[key]["limits"][0] > cut for r in result["rows"] if r["mode"] == mode
        )

    assert wins("native", "three", 1 / 3) == 3
    assert wins("native", "two", 0.5) == 2
    assert wins("which_head", "two", 0.5) == 0
    missing = copy.deepcopy(data)
    missing["banks"][0]["counts"][-1] -= 1
    with pytest.raises(ValueError):
        decode(missing, training)
