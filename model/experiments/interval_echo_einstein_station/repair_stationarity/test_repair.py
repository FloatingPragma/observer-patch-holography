"""Recover the finite projection, physical cut and entropy variations independently."""

from fractions import Fraction
import json
from math import log

import numpy as np

from .coarse import coarse_patch
from .measure import HERE, SHOTS, STEP, entropy, run


def test_entropies_and_parity_from_archived_word_counts():
    receipt = json.loads((HERE / "results.json").read_text())
    replay = run()
    assert replay["apparatus"] == receipt["apparatus"]
    assert replay["outcome_order"] == receipt["outcome_order"]
    for observed, stored in zip(replay["acquisitions"], receipt["acquisitions"], strict=True):
        for key in ("joint_minus_counts", "joint_plus_counts", "routed_minus_counts", "routed_plus_counts"):
            assert observed[key] == stored[key]
        np.testing.assert_allclose(observed["exact_entropy_slope"], stored["exact_entropy_slope"], atol=1e-13)
    assert receipt["feasible_joint_tangent_dimension"] == 4
    for board in receipt["apparatus"]:
        assert board["occupied_region_loads"] == [6, 6, 6]
        assert sum(Fraction(n) for n in board["common_cut_distribution"].values()) == 1
        assert board["accessible_tables"] == [["3/8", "1/8", "1/8", "3/8"]] * 2
    for row in receipt["acquisitions"]:
        plus = np.array(row["joint_plus_counts"])
        minus = np.array(row["joint_minus_counts"])
        assert plus.sum() == minus.sum() == SHOTS
        measured = (entropy(plus / SHOTS) - entropy(minus / SHOTS)) / (2 * float(STEP))
        np.testing.assert_allclose(measured, row["measured_entropy_slope"], atol=1e-12)
        target = 0. if row["theta"] == "1/16" else -log(5.)
        np.testing.assert_allclose(row["exact_entropy_slope"], target, atol=1e-12)
        assert abs(measured - row["finite_entropy_slope"]) < 3 * row["entropy_slope_standard_error"]
        assert abs(row["measured_routed_mismatch_slope"] + 2.) < 3 * row["routed_slope_standard_error"]


def test_coarse_residual_extensive_on_native_regions():
    rows = json.loads((HERE / "coarse.json").read_text())
    actual = [coarse_patch(level) for level in range(3)]
    for replay, stored in zip(actual, rows, strict=True):
        assert [list(group) for group in replay["groups"]] == stored["groups"]
        assert replay["cut_contacts"] == stored["cut_contacts"]
        np.testing.assert_allclose(replay["entropy_slope_at_one_eighth"], stored["entropy_slope_at_one_eighth"])
    for row, cells in zip(rows, (18, 66, 279), strict=True):
        assert row["region_cells"] == cells
        assert len(set(c for group in row["groups"] for c in group)) == cells
        assert abs(row["entropy_slope_at_one_eighth"] / cells + log(5.)/3) < 1e-12
