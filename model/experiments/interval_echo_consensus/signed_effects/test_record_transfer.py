"""Native pointer transfer checked against exact gates and reference-entangled maps."""

from itertools import product
import json

import numpy as np
import sympy as s

from model.executable_model import cell
from . import record_transfer as transfer, sequential


def test_exact_binary_refocusing_and_common_signed_phase():
    # Construct the four binary endpoint words independently of the producer.
    swap = s.Matrix([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]])
    contrast = s.diag(0, 1, 1, 0)
    h = (s.eye(4) - swap) / 2 + 8 * contrast
    parity = s.diag(1, 1, -1, -1)
    # Sum of toggled generators is diagonal, and the halves commute exactly.
    assert h * (parity * h * parity) == (parity * h * parity) * h
    assert h + parity * h * parity == 17 * contrast
    phase = s.diag(1, -s.I, -s.I, 1)
    correction = s.diag(1, s.I, s.I, -1)
    assert correction * phase == s.diag(1, 1, 1, -1)
    for q in (-8, 24):
        assert s.simplify(s.exp(-s.I * q * s.pi / 8)) == -1
        assert s.simplify(s.exp(-s.I * q * s.pi / 16)) == s.I


def test_native_contact_geometry_full_sector_and_all_encoded_words():
    required = {(0, 1), (0, 10), (6, 10), (0, 5), (4, 5)}
    assert len(transfer.STATES) == 20
    for a, b in required:
        assert b in cell().neighbours[a]
    for sign, pointer in product((-1, 1), (0, 1)):
        u = transfer.execute(sign, transfer.copy_word(pointer))
        np.testing.assert_allclose(u.conj().T @ u, np.eye(20), atol=4e-13, rtol=0)
        for bits in product((0, 1), repeat=2):
            v = transfer.embedding(bits)
            target = np.zeros_like(v)
            for source in (0, 1):
                changed = list(bits)
                changed[pointer] ^= int(source == 0)
                target[:, source] = -transfer.embedding(tuple(changed))[:, source]
            np.testing.assert_allclose(u @ v, target, atol=4e-13, rtol=0)


def test_deferred_final_read_matches_every_reference_entangled_branch():
    bell = np.array([1, 0, 0, 1], complex) / np.sqrt(2)
    for sign, swapped in product((-1, 1), (False, True)):
        actual = transfer.qubit_maps(sign, swapped)
        direct = sequential.instrument(sign)
        x = np.array([[0, 1], [1, 0]]) if swapped else np.eye(2)
        for (first, second), k in actual.items():
            expected = direct[second][:2, :2] @ x @ direct[first][:2, :2]
            a = np.kron(k, np.eye(2)) @ bell
            b = np.kron(expected, np.eye(2)) @ bell
            np.testing.assert_allclose(np.outer(a, a.conj()), np.outer(b, b.conj()),
                                       atol=4e-13, rtol=0)
        # The complete material record contains invalid-pointer words as well.
        maps = transfer.branch_maps(sign, swapped)
        reference = np.zeros((2, 2), complex)
        for k in maps:
            v = (np.kron(k, np.eye(2)) @ bell).reshape(20, 2)
            reference += v.T @ v.conj()
        np.testing.assert_allclose(reference, np.eye(2) / 2, atol=4e-13, rtol=0)


def test_unread_undo_recovers_phase_and_pointer_dephasing_prevents_recovery():
    v = transfer.embedding()
    psi = v @ (np.array([1, 1j]) / np.sqrt(2))
    for sign in (-1, 1):
        word = transfer.copy_word(0)
        unit = transfer.execute(sign, word)
        undo = transfer.execute(sign, transfer.inverse(word))
        carried = unit @ psi
        np.testing.assert_allclose(undo @ carried, psi, atol=4e-13, rtol=0)
        # Ignore a complete first-pointer occupation read, then reverse the same word.
        restored = np.zeros((20, 20), complex)
        for label in range(4):
            branch = carried.copy()
            for i, marked in enumerate(transfer.STATES):
                a, b = transfer.POINTERS[0]
                if int(a in marked) + 2 * int(b in marked) != label:
                    branch[i] = 0
            branch = undo @ branch
            restored += np.outer(branch, branch.conj())
        np.testing.assert_allclose(v.conj().T @ restored @ v, np.eye(2) / 2,
                                   atol=4e-13, rtol=0)
        assert abs(np.vdot(psi, restored @ psi) - 0.5) < 4e-13


def test_complete_record_banks_and_native_interaction_controls():
    saved = json.loads(transfer.RECEIPT.read_text())
    sequential.compare(saved, transfer.report())
    banks = [row for row in saved["cases"] if "shots" in row]
    assert len(banks) == 16 and sum(row["shots"] for row in banks) == 131072
    for row in banks:
        assert len(row["all_16_pointer_word_counts"]) == 16
        assert sum(row["all_16_pointer_word_counts"]) == row["shots"]
        if row["arm"] == "nondemolition":
            assert sum(row["all_16_pointer_word_counts"][j] for j in (5, 10)) == row["shots"]
    for row in saved["cases"]:
        if "controls" in row:
            c = row["controls"]
            assert c["no_refocus_population_outside_one_mark_per_dimer"] > 0.001
            assert c["nondemolition_population_outside_one_mark_per_dimer"] < 1e-24
            assert c["native_source_record_total_variation"] > 1 - 3e-12
            assert c["cross_contact_omitted_source_record_total_variation"] < 3e-12
