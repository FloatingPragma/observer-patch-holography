"""Independent finite map and observer-history checks for the signed instrument."""

import copy
import json

import numpy as np
import pytest
import sympy as s

from . import sequential


def exact_active(sign):
    # Reconstruct from the two-level splitter and its inverse, independently
    # of the producer's twelve-position exponentials and native block loader.
    eye = s.eye(2)
    x = s.Matrix([[0, 1], [1, 0]])
    splitter = (eye + s.I * x) / s.sqrt(2)
    phase = (-2 * s.pi / 5) if sign == -1 else (6 * s.pi / 5)
    value = s.expand_complex(s.exp(-s.I * phase))
    u = s.simplify(splitter.conjugate().T * s.diag(value, 1) * splitter)
    return [s.simplify(u.conjugate().T * s.diag(1 - j, j) * u) for j in (0, 1)]


def test_exact_phase_effects_and_distinct_poststates():
    x = s.Matrix([[0, 1], [1, 0]])
    for sign in (-1, 1):
        exact = exact_active(sign)
        ops = sequential.instrument(sign)
        changed = sequential.instrument(sign, swapped=True)
        for outcome in (0, 1):
            e = exact[outcome]
            assert s.simplify(e * e - e) == s.zeros(2)
            assert s.simplify(e.trace() - 1) == 0
            assert e[0, 1] != 0
            np.testing.assert_allclose(ops[outcome][:2, :2], np.array(e.evalf(), complex),
                                       atol=3e-13, rtol=0)
            assert s.simplify((x * e).conjugate().T * (x * e) - e) == s.zeros(2)
        disagreement = s.simplify((exact[1] * x * exact[0] * x).trace())
        assert s.simplify(disagreement - (3 + sign * s.sqrt(5)) / 8) == 0


def test_every_matrix_unit_effect_completeness_and_repeat_map():
    for sign in (-1, 1):
        plain = sequential.instrument(sign)
        changed = sequential.instrument(sign, swapped=True)
        for j in range(12):
            for k in range(12):
                matrix = np.zeros((12, 12), complex)
                matrix[j, k] = 1
                for y, (a, b) in enumerate(zip(plain, changed, strict=True)):
                    np.testing.assert_allclose(np.trace(a @ matrix @ a.conj().T),
                                               np.trace(b @ matrix @ b.conj().T),
                                               atol=3e-13, rtol=0)
                    for z in range(12):
                        target = a @ matrix @ a.conj().T if z == y else 0
                        result = plain[z] @ a @ matrix @ a.conj().T @ plain[z].conj().T
                        np.testing.assert_allclose(result, target, atol=3e-13, rtol=0)


def test_untouched_reference_and_same_effects_different_joint_states():
    state = np.zeros(24, complex)
    state[0] = state[3] = 1 / np.sqrt(2)
    initial_reference = np.eye(2) / 2
    for sign in (-1, 1):
        plain = sequential.instrument(sign)
        changed = sequential.instrument(sign, swapped=True)
        reference = np.zeros((2, 2), complex)
        distance = 0.0
        for a, b in zip(plain, changed, strict=True):
            first = np.kron(a, np.eye(2)) @ state
            rival = np.kron(b, np.eye(2)) @ state
            rho = np.outer(first, first.conj())
            other = np.outer(rival, rival.conj())
            np.testing.assert_allclose(np.trace(rho), np.trace(other), atol=3e-13)
            reference += np.trace(rho.reshape(12, 2, 12, 2), axis1=0, axis2=2)
            distance += np.linalg.svd(rho - other, compute_uv=False).sum() / 2
        np.testing.assert_allclose(reference, initial_reference, atol=3e-13)
        assert distance > 0.3


def test_complete_banks_and_replay_tamper_detection():
    saved = json.loads(sequential.RECEIPT.read_text())
    actual = sequential.report()
    sequential.compare(saved, actual)
    assert sum(row["shots"] for row in saved["cases"]) == 262144
    for row in saved["cases"]:
        bins = np.array(row["joint_address_counts"])
        assert bins.shape == (12, 12) and bins.sum() == row["shots"]
        assert bins[2:].sum() == bins[:, 2:].sum() == 0
        if row["arm"] == "repeatable":
            assert bins.sum() == np.trace(bins)
        else:
            assert bins.sum() - np.trace(bins) > 1000
    damaged = copy.deepcopy(saved)
    damaged["cases"][0]["joint_address_counts"][0][0] += 1
    with pytest.raises(ValueError, match="exact source or discrete"):
        sequential.compare(damaged, actual)
    damaged = copy.deepcopy(saved)
    damaged["source_sha256"][next(iter(damaged["source_sha256"]))] = "0" * 64
    with pytest.raises(ValueError, match="exact source or discrete"):
        sequential.compare(damaged, actual)
    damaged = copy.deepcopy(saved)
    damaged["cases"][0]["repeat_disagreement"] += 1e-13
    sequential.compare(damaged, actual)
