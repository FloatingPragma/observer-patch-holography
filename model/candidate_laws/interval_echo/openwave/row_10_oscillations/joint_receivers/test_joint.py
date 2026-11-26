import json
from pathlib import Path

import numpy as np

from .. import eraser
from .apparatus import access, pointer
from .measure import BENCH, MODES, STARTS
from .readback import interval, read


def record():
    return json.loads(Path(__file__).with_name("counts.json").read_text())


def test_native_crossings_select_two_paths_with_the_registered_phase():
    for sense in (-1, 1):
        v = np.array([0, 1, np.exp(sense * 2j * np.pi / 3), 0]) / 2
        np.testing.assert_allclose(pointer(sense)[0], np.outer(v, v.conj()), atol=3e-13)
        np.testing.assert_allclose(
            pointer(sense, True)[0], np.diag([0, 1 / 4, 0, 0]), atol=3e-13
        )
    assert access() == record()["source_apparatus"]


def test_all_material_panel_standards_including_coherent_source_refusals():
    data = record()
    assert len(data["bench"]) == 24
    for row in data["bench"]:
        sense, cut, word = row["key"]
        k, _, _, refused = row["counts"]
        assert sum(row["counts"]) == BENCH
        if word < 4:
            p = 0.25 * (word == 1 or (word == 2 and not cut))
            q = 0
        else:
            same = sense == (-1, 1)[word - 4]
            p = 1 / 4 if same and not cut else 1 / 16
            q = 0.5
        assert interval(k, BENCH)[0] <= p <= interval(k, BENCH)[1]
        assert interval(refused, BENCH)[0] <= q <= interval(refused, BENCH)[1]


def test_physical_click_agrees_with_independent_source_and_guide_forecast():
    data = record()
    original = eraser.forecast()
    held = {tuple(r["key"]): r["counts"] for r in data["held"]}
    for row in original["held"]:
        for sign in (1, 2):
            outcome = row["outcomes"][sign]
            x = outcome["granted_joint"]["X"]
            y = outcome["granted_joint"]["Y"]
            for panel in (-1, 1):
                # Half of a one-hot projector, at angle +/-2pi/3.
                p = (
                    (x[1] + x[2])
                    + np.cos(2 * np.pi / 3) * (x[1] - x[2])
                    + panel * np.sin(2 * np.pi / 3) * (y[1] - y[2])
                ) / 4
                for mode in MODES[:2]:
                    if mode == "which_head":
                        p = (x[1] + x[2]) / 4
                    counts = held[row["owner"], row["time"], sign, panel, mode]
                    low, high = interval(counts[0], STARTS)
                    assert low <= p <= high
                assert abs(sum(x) / row["original_pair_probability"] - 1) < 1e-10


def test_original_start_banks_and_conjugate_sign_controls():
    results = read(record())
    assert len(results) == 48 and not any(
        r["far_marginal"]["resolved"] for r in results
    )
    for mode, n in zip(MODES, (9, 0, 0, 12), strict=True):
        assert (
            sum(r["joint_given_far"]["resolved"] for r in results if r["mode"] == mode)
            == n
        )


def test_blind_labels_are_frozen_by_training_and_the_validation_is_complete():
    from .blind import decode

    data = json.loads(Path(__file__).with_name("blind_counts.json").read_text())
    results = decode(data)
    assert len(results) == 6
    assert all(r["label_count"] == 3 and r["limits"][0] > 1 / 3 for r in results)
    assert sum(len(r["raw_pairwise_separations"]) == 3 for r in results) == 4
    # Rolling the validation source labels must leave every decision unchanged.
    for row in data["banks"]:
        if row["key"][-1] == "validation":
            row["key"][2] = (row["key"][2] + 1) % 3
    reassigned = decode(data)
    assert [r["decision"] for r in results] == [r["decision"] for r in reassigned]
    assert not any(
        x["balanced_success"] == y["balanced_success"]
        for x, y in zip(results, reassigned, strict=True)
    )
