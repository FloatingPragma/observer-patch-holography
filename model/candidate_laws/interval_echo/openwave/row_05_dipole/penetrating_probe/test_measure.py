from collections import defaultdict
from copy import deepcopy
import json
import sys

import numpy as np
import pytest

from . import measure


def test_exact_reflection_reverses_the_actual_compact_source(monkeypatch):
    certificate, _, _ = measure.reflection_certificate()
    assert certificate["exact_compact_action"] == [["0", "1"], ["1", "0"]]
    assert certificate["compact_determinant"] == -1 and certificate["source_mapping_error"] < 2e-12
    monkeypatch.setattr(measure, "REFLECTION", tuple(range(12)))
    with pytest.raises(ArithmeticError, match="orientation reversing"):
        measure.reflection_certificate()


def test_principal_restriction_matches_independently_assembled_native_columns():
    board = measure.federation(0)
    x, y, _, kinetic, h, _ = measure.neutral_defect.generators(board)
    xx, yy, zero, native = measure.region(board, x, y, kinetic, h, (0, 4, 6))
    positions = {(int(a), int(b)): i for i, (a, b) in enumerate(zip(xx, yy, strict=True))}
    law = measure.neutral_defect.law
    for strength, matrix in ((0, zero), (8, native)):
        local = law.seam_term(strength)
        for excess, vacancy in ((82, 5), (48, 11), (2, 3), (2, 48)):
            source = positions[excess, vacancy]
            direct = defaultdict(float)
            direct[excess, vacancy] = strength*(3*(6-board.degree[excess])-(6-board.degree[vacancy]))-12*strength
            incident = {tuple(sorted((p, q))) for p in (excess, vacancy) for q in board.neighbours[p]}
            for a, b in incident:
                loads = {excess: 2, vacancy: 0}
                left, right = loads.get(a, 1), loads.get(b, 1)
                word = law.SRC.index((left, right))
                direct[excess, vacancy] += local[word, word]
                if left != right:
                    swap = {a: b, b: a}
                    target = swap.get(excess, excess), swap.get(vacancy, vacancy)
                    if target in positions:
                        direct[target] += local[law.SRC.index((right, left)), word]
            expected = np.zeros(len(xx))
            for pair, value in direct.items():
                expected[positions[pair]] = value
            np.testing.assert_array_equal(matrix[:, source], expected)


def test_one_opening_null_and_two_opening_changed_law_response():
    receipt = json.loads(measure.ARCHIVE.read_text())
    for setting in receipt["studies"][0]["settings"]:
        assert setting["exact_reflection_commutes"]
        for row in setting["readings"]:
            assert row["one_opening_state_symmetry_error"] < 2e-11
            assert abs(row["excess_probability_difference"]) < 2e-11
            np.testing.assert_allclose(row["source_plus_minus"][0]["port_load_zero_host_excess"],
                                       row["source_plus_minus"][1]["port_load_zero_host_excess"],
                                       rtol=0, atol=2e-11)
    native, changed = receipt["studies"][1]["settings"]
    assert native["readings"][-1]["excess_probability_difference"] < -.02
    assert changed["readings"][-1]["excess_probability_difference"] > .007
    for arm in native["readings"][-1]["source_plus_minus"]:
        assert arm["correct_source_and_exterior_then_complement"][0] < .25
        assert arm["probe_exterior_inside"][1] > .15


def test_every_detector_and_source_complement_is_retained():
    receipt = json.loads(measure.ARCHIVE.read_text())
    for study in receipt["studies"]:
        for setting in study["settings"]:
            assert setting["numerical_probability_allowance"] < 1e-9
            for row in setting["readings"]:
                for arm in [*row["source_plus_minus"], row["equal_source_mixture"]]:
                    for values in arm.values():
                        assert min(values) >= 0 and abs(sum(values)-1) < 2e-12


def test_failed_serialization_keeps_the_destination_absent(tmp_path, monkeypatch):
    destination = tmp_path/"receipt.json"
    monkeypatch.setattr(measure, "run", lambda: {"invalid": np.bool_(True)})
    monkeypatch.setattr(sys, "argv", ["measure", "--record", str(destination)])
    with pytest.raises(TypeError):
        measure.main()
    assert not destination.exists()


def test_replay_rejects_a_relabelled_local_detector():
    original = measure.ARCHIVE.read_bytes()
    receipt = json.loads(original)
    changed = deepcopy(receipt)
    changed["studies"][1]["settings"][0]["readings"][-1]["source_plus_minus"][0]["port_load_zero_host_excess"].reverse()
    with pytest.raises(ArithmeticError):
        measure.compare(receipt, changed)
    assert measure.ARCHIVE.read_bytes() == original
