"""Check the supplied conditional expectation and distinguish it from the old meter."""

from copy import deepcopy
import json

import numpy as np
import pytest

from openwave.interval_echo.row_19.recovery_identity import measure


def test_conditional_expectation_of_general_effect_and_support():
    ref, product, _, _, score, audit = measure.calibration()
    assert audit["preparations"]["conditional_product"]["realizable"]
    left, right = measure.head.marginals(ref)
    rng = np.random.default_rng(29)
    for b in range(2):
        tau = right[b] / np.trace(right[b])

        def expectation(a):
            reshaped = a.reshape(2, 2, 2, 2)
            local = np.einsum("ijkl,lj->ik", reshaped, tau)
            return np.kron(local, np.eye(2))

        a = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
        a = a + a.conj().T
        np.testing.assert_allclose(
            expectation(expectation(a)), expectation(a), atol=2e-14
        )
        np.testing.assert_allclose(expectation(np.eye(4)), np.eye(4), atol=1e-14)
        np.testing.assert_allclose(
            np.trace(ref[b] @ expectation(a)), np.trace(product[b] @ a), atol=1e-14
        )
        assert np.linalg.matrix_rank(score[b]) == 4
        assert (
            np.linalg.eigvalsh(score[b]).min() < 0 < np.linalg.eigvalsh(score[b]).max()
        )


def test_signed_views_and_full_recovery():
    data = json.loads(measure.RECEIPT.read_text())
    for key, g in data["gates"].items():
        assert (
            measure.decode(g["counts"], g["pauli_coefficients"], float(key))
            == g["measured"]
        )
        for mode, expected in g["expected"].items():
            m = g["measured"][mode]
            assert abs(m["estimate"] - expected) <= m["halfwidth"]
    middle = data["gates"]["0.5"]["measured"]
    assert (
        middle["flagged"]["estimate"] + middle["flagged"]["halfwidth"]
        < middle["unflagged"]["estimate"] - middle["unflagged"]["halfwidth"]
    )
    full = data["gates"]["1.0"]
    assert full["reference_trace_distance"] == 0
    assert full["expected"]["flagged"] == full["expected"]["unflagged"] == 0
    assert (
        full["measured"]["original"]["estimate"]
        > full["measured"]["original"]["halfwidth"]
    )
    broken = deepcopy(full["counts"])
    broken[0][0][1][0] -= 1
    with pytest.raises(ValueError, match="Incomplete"):
        measure.decode(broken, full["pauli_coefficients"], 1.0)
