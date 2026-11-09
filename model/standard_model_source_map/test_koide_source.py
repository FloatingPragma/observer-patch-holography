"""Complete coherent-source and unconditioned weak-feedback controls."""
from fractions import Fraction as F
import math

import numpy as np
import pytest
import sympy as s
from scipy.linalg import expm

from . import koide_source as k


def test_full_native_partition_including_zeros():
    rows = k.native_partition()
    assert len(rows) == 128
    assert sum(s.Rational(r['probability']) for r in rows) == 1
    assert sum(r['probability'] == '0' for r in rows) == 104
    assert [sum(s.Rational(r['probability']) for r in rows if r['group'] == j)
            for j in range(5)] == list(k.WEIGHTS)


def test_native_rotation_is_exactly_the_cyclic_frame():
    r, o = k._native_rotation()
    assert k._clean(r*o-o*k.P) == s.zeros(3)
    assert k._clean(o.T*o) == s.eye(3)


@pytest.mark.parametrize('z', [s.S.One, (s.Integer(99)+20*s.I)/101])
def test_complete_analyzer_and_balanced_action(z):
    src = k.exact_source(z)
    assert src['positive_root_chamber_certified']
    assert src['signed_root_Koide_ratio'] == s.Rational(2, 3)
    assert s.trace(src['action_coefficient']) == s.Rational(5, 18)
    a = k.analyzer_block(0, 0, z)[:3, :3]
    assert k._clean(a.H*a-src['action_coefficient']) == s.zeros(3)
    assert [r['probability'] for r in k.exact_analyzer_records(z)] == [
        '5/54', '5/54', '5/54', '1/36', '25/36', '0', '0', '0']
    total = sum((k.analyzer_block(j, 0, z).H*k.analyzer_block(j, 0, z)
                 for j in range(8)), s.zeros(4))
    assert k._clean(total) == s.eye(4)


def test_positive_chamber_is_not_a_signed_identity():
    assert k.exact_source(-1)['signed_root_Koide_ratio'] == s.Rational(2, 3)
    assert not k.exact_source(-1)['physical_action_Koide_ratio_certified']


def test_all_active_flags_give_scalar_and_one_flag_is_input_dependent():
    effects = [k._clean(k.analyzer_block(j, 0)[:3, :3].H*k.analyzer_block(j, 0)[:3, :3])
               for j in range(3)]
    assert k._clean(sum(effects, s.zeros(3))) == s.Rational(5, 18)*s.eye(3)
    assert effects[0] != effects[1]
    assert effects[0] != s.trace(effects[0])*s.eye(3)/3


def test_weak_kick_is_complete_on_all_six_inputs():
    ks = k.weak_kraus()
    assert len(ks) == 8
    assert k._clean(sum((a.H*a for a in ks), s.zeros(6))) == s.eye(6)
    y = k.exact_source()['action_coefficient']; q = s.Rational(1, 64)
    expected = s.kronecker_product(s.eye(2), 4*q*q/(1+q*q)*(y-y*y))
    assert k._clean(sum((a.H*a for a in ks[1:]), s.zeros(6))-expected) == s.zeros(6)
    assert sum(s.sympify(r['probability']) for r in k.weak_records()) == 1


def test_two_complete_rounds_with_reference_obey_original_bound():
    ks = [np.asarray(a, complex) for a in k.weak_kraus()]
    y = np.asarray(k.exact_source()['action_coefficient'], complex)
    h = 2*math.atan(1/64); unitary = expm(-2j*h*np.kron([[0, 1], [1, 0]], y))
    psi = np.eye(6).reshape(36)/math.sqrt(6)
    desired = np.kron(unitary, np.eye(6))@psi
    error = 0.; total = 0.
    for a in range(8):
        for b in range(8):
            out = np.kron(ks[b]@ks[a], np.eye(6))@psi
            rho = np.outer(out, out.conj()); total += np.vdot(out, out).real
            if a == b == 0: rho -= np.outer(desired, desired.conj())
            error += abs(np.linalg.eigvalsh(rho)).sum()
    assert abs(total-1) < 2e-13
    assert error <= 9*2*h*h/4


def test_native_words_execute_all_input_columns_and_keep_failures():
    row = k.native_execution()
    assert row['complete_analyzer_columns'] == 32
    assert row['complete_weak_round_columns'] == 64
    assert row['retained_reset_flags'] == 8
    assert row['indexed_first_service_failure_slots'] > 0
    assert row['valid_code_failure_weights_exact'] == '0'
    assert not row['numerical_admission_is_a_rigorous_roundoff_certificate']


def test_bound_keeps_native_error_and_does_not_renormalize():
    row = k.feedback_bound(1024, F(1), F(1, 10**8))
    assert row['ideal_nonzero_archive_probability_upper'] == F(1, 4096)
    assert row['nonzero_archive_probability_upper'] == F(1, 4096)+F(1024, 10**8)
    assert row['complete_unhalved_trace_distance_upper'] == F(9, 4096)+F(2048, 10**8)
    assert not row['zero_record_postselection_used']


@pytest.mark.parametrize('z', [0, 2, .5, s.Symbol('z')])
def test_phase_admission(z):
    with pytest.raises(ValueError): k.exact_source(z)


@pytest.mark.parametrize('n,t', [(0, F(1)), (True, F(1)), (2, -1), (2, .1)])
def test_clock_admission(n,t):
    with pytest.raises(ValueError): k.feedback_bound(n,t)
