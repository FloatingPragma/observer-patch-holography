from fractions import Fraction as F
from unittest.mock import patch

import numpy as np
import pytest
from scipy.linalg import expm

from model.field_protocols import protocol as p
from model.full_fock_field.field import FieldGeometry, FullFock
from model.observer.source_scalar import compiler


@pytest.fixture(scope='module')
def tree():
    return p.adaptive()


def test_all_inputs_and_references_are_complete(tree):
    assert set(tree['success']) == {(0,0), (0,1), (1,0), (1,1)}
    assert tree['trace_preservation_residual'] < 2e-10
    assert tree['numerical_dilation_difference_norm'] < float(tree['ideal_word_half_diamond_bound'])
    rng = np.random.default_rng(741)
    state = rng.normal(size=(16, 3))+1j*rng.normal(size=(16, 3))
    state /= np.linalg.norm(state)
    gram = sum((k@state).conj().T@(k@state) for k in tree['success'].values())
    gram += sum(((k@state).conj().T@(k@state) for branches in tree['aborts'].values()
                 for k in branches), np.zeros((3,3), complex))
    np.testing.assert_allclose(gram, state.conj().T@state, atol=2e-10)


def test_independent_matrix_exponential_and_wrong_feedback_control(tree):
    space = FullFock(4)
    h = space.hamiltonian(FieldGeometry((2,1,1))).toarray()
    u = expm(-.25j*h)
    n0, n1 = (space.number(i).toarray() for i in (0,1))
    first = (np.eye(16)-n0, n0)
    second = (np.eye(16)-n1, n1)
    for bit in (0,1):
        for result in (0,1):
            expected = second[result]@(u if bit == 0 else u.conj().T)@first[bit]@u
            np.testing.assert_allclose(tree['ideal'][bit,result], expected, atol=1e-14)
    wrong = second[0]@u@first[1]@u
    assert np.linalg.norm(wrong-tree['ideal'][1,0]) > .01


def test_injected_failure_survives_all_outcomes_and_reference():
    primitive = compiler.native_pulse
    def leaky(kind, angle):
        unitary = primitive(kind, angle)
        # A complete source unitary with a deliberately nonzero code exit.
        rotation = np.eye(6, dtype=complex)
        c, s = np.cos(.002), np.sin(.002)
        rotation[0,0] = rotation[5,5] = c
        rotation[0,5], rotation[5,0] = -s, s
        return rotation@unitary
    with patch.object(compiler, 'native_pulse', side_effect=leaky):
        tree = p.adaptive(FieldGeometry((1,1,1)), p.Read(F(0), 0, 1),
                          (p.Read(F(0), 1, 1), p.Read(F(0), 1, 1)), diagnostic_only=True)
    failure = tree['failure_gram']
    assert np.trace(failure).real > 1e-5
    success = sum(k.conj().T@k for k in tree['success'].values())
    np.testing.assert_allclose(success+failure, np.eye(4), atol=2e-10)
    assert np.linalg.norm(success-np.eye(4)) > 1e-5
    assert any(prefix == (0,) for prefix, _ in tree['aborts'])
    assert any(prefix == (1,) for prefix, _ in tree['aborts'])
    assert tree['diagnostic_only']
    assert not tree['ideal_source_error_budget_satisfied']


def test_wrong_unitary_source_cannot_pass_on_trace_preservation_alone():
    original = compiler.native_pulse
    def changed(kind, angle):
        return original(kind, -angle)
    with patch.object(compiler, 'native_pulse', side_effect=changed):
        with pytest.raises(ArithmeticError, match='error budget'):
            p.adaptive()


def test_zero_split_error_does_not_accept_leakage_as_ideal_source():
    primitive = compiler.native_pulse
    def leaky(kind, angle):
        rotation = np.eye(6, dtype=complex)
        rotation[0,0] = rotation[5,5] = np.cos(.002)
        rotation[0,5] = -np.sin(.002); rotation[5,0] = np.sin(.002)
        return rotation@primitive(kind, angle)
    with patch.object(compiler, 'native_pulse', side_effect=leaky):
        with pytest.raises(ArithmeticError, match='error budget'):
            p.adaptive(FieldGeometry((1,1,1)), p.Read(F(0), 0, 1),
                       (p.Read(F(0), 1, 1), p.Read(F(0), 1, 1)))


@pytest.mark.parametrize('bad', [True, float('nan'), .5])
def test_exact_time_admission(bad):
    with pytest.raises(ValueError):
        p.Read(bad, 0)


def test_missing_branch_rejects():
    with pytest.raises(ValueError):
        p.adaptive(second=(p.Read(F(1),0),))
