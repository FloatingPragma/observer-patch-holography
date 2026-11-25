"""Check the signed entropy inequality and the complete-cell shielding measurement."""

from fractions import Fraction

import numpy as np
from scipy.linalg import eigh

from model.executable_model import cell, federation
from .conditional_information import measure, native_front


def test_all_orders_signed_bound_with_a_changing_bath():
    rng = np.random.default_rng(789445)
    x = rng.normal(size=(8, 8)) + 1j*rng.normal(size=(8, 8))
    state = x @ x.conj().T + 3*np.eye(8)
    state /= np.trace(state)
    y = rng.normal(size=(8, 8)) + 1j*rng.normal(size=(8, 8))
    tangent = y + y.conj().T
    tangent -= np.eye(8)*np.trace(tangent)/8
    step = 0.4*eigh(state, eigvals_only=True)[0]/np.linalg.norm(tangent, 2)

    def reductions(matrix):
        blocks = matrix.reshape((2,)*6)
        ab = np.trace(blocks, axis1=2, axis2=5)
        bc = np.trace(blocks, axis1=0, axis2=3)
        b = np.trace(bc, axis1=1, axis2=3)
        return [ab.reshape(4, 4), bc.reshape(4, 4), b, matrix]

    signs = np.array([1, 1, -1, -1])
    def logm(matrix):
        values, vectors = eigh(matrix)
        assert min(values) > 0
        return (vectors*np.log(values)) @ vectors.conj().T

    def trace_pair(x, y):
        return np.trace(x@y).real

    logs = [logm(s) for s in reductions(state)]
    info = -signs @ [trace_pair(s, k) for s, k in zip(reductions(state), logs, strict=True)]
    slopes = [-trace_pair(d, k) for d, k in zip(reductions(tangent), logs, strict=True)]
    slope = signs @ slopes
    assert abs(slopes[1]-slopes[2]) > 0.01
    q = trace_pair(tangent @ tangent, np.linalg.inv(state))
    for sign in (-1, 1):
        changes = reductions(state + sign*step*tangent)
        rel = [trace_pair(x, logm(x)-k) for x, k in zip(changes, logs, strict=True)]
        observed = -signs @ [trace_pair(x, logm(x)) for x in changes]
        np.testing.assert_allclose(observed, info+sign*step*slope-signs@rel, atol=1e-13)
        assert 0 <= observed <= info+sign*step*slope+rel[-1]+1e-13
        assert rel[-1] <= np.log1p(step**2*q)+1e-13
    assert abs(slope) <= info/step + step*q


def test_quantum_screen_and_commuting_control():
    no_buffer, partial, full = measure()['records']
    assert no_buffer['conditional_information'] > partial['conditional_information'] > full['conditional_information']
    np.testing.assert_allclose(full['conditional_information'], 4.56504390e-5, rtol=1e-7)
    np.testing.assert_allclose(full['max_hidden_diagonal_slope'], 2.12365332e-5, rtol=2e-5)
    assert abs(full['contrast_only_information']) < 1e-12
    assert abs(full['contrast_only_hidden_slope']) < 1e-10
    assert partial['contrast_only_hidden_slope'] > 0.004
    for record in (no_buffer, partial, full):
        assert record['reset_signed_error'] <= -record['conditional_information']
        for name, error in record['signed_error'].items():
            assert abs(error) <= record['certified_absolute_error'][name]
    np.testing.assert_allclose(full['signed_error']['swap'], 3.3246669e-5, rtol=1e-7)
    assert 0 < full['signed_error']['parity'] < full['max_hidden_diagonal_slope']

    np.testing.assert_allclose(full['reset_signed_error'], -0.00014718082127, atol=1e-12)


def test_native_propagation_through_buffers():
    for board, distance in ((cell(), 3), (federation(0), 14), (federation(1), 32)):
        fronts = native_front(board, 0)
        assert len(fronts)-1 == distance
        assert sum(len(layer) for layer in fronts) == board.n
        assert all(p > 0 for layer in fronts for p in layer.values())
        if board.n == 12:
            assert fronts[3][3] == Fraction(5, 24)
