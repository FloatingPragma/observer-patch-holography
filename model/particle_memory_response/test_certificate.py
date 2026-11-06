"""Independent covariance, chronological generator and source-bound controls."""
from fractions import Fraction as F
import itertools

import numpy as np
import pytest
import scipy.linalg as la
import sympy as s

from . import certificate as c
from ..common_source import native_family_record_memory as memory


@pytest.mark.parametrize('policy', c.POLICIES)
def test_group_inverse_replays_full_native_source(policy):
    values = c.response(policy=policy)
    T, pi, R = values['T'], values['pi'], values['R']
    Pi = pi*s.ones(1, 3)
    assert c.zero((s.eye(3)-T)*R-(s.eye(3)-Pi))
    assert c.zero(R*(s.eye(3)-T)-(s.eye(3)-Pi))
    assert c.zero(R*pi)
    for n in (1, 2, 5):
        partial = sum((T**k-Pi for k in range(n)), s.zeros(3))
        assert c.zero(R-partial-T**n*R)


def direct_reward_covariance(T, E, pi, n):
    """Forward first/second reward moments, without lag formulas or resolvents."""
    means, moments = [s.zeros(3, 1) for _ in range(3)], [s.zeros(3) for _ in range(3)]
    distribution = pi.copy()
    for _ in range(n):
        next_mean, next_moment = [s.zeros(3, 1) for _ in range(3)], [s.zeros(3) for _ in range(3)]
        for new in range(3):
            e = s.zeros(3, 1); e[new] = 1
            for old in range(3):
                next_mean[new] += T[new, old]*means[old]+E[new, old]*distribution[old]*e
                next_moment[new] += T[new, old]*moments[old]+E[new, old]*(
                    e*means[old].T+means[old]*e.T+distribution[old]*e*e.T)
        means, moments = next_mean, next_moment
        distribution = T*distribution
    mean = sum(means, s.zeros(3, 1))
    return c.clean((sum(moments, s.zeros(3))-mean*mean.T)/n)


@pytest.mark.parametrize('policy', c.POLICIES)
def test_full_finite_reward_moments_match_lag_sum_with_correct_order(policy):
    values = c.response(policy=policy)
    T, E, pi, w = (values[k] for k in ('T', 'E', 'pi', 'w'))
    C0 = s.diag(*w)-w*w.T
    for n in (1, 2, 5):
        K = C0/2+sum(((1-s.Rational(k, n))*(E*T**(k-1)*s.diag(*w)-w*w.T)
                     for k in range(1, n)), s.zeros(3))
        covariance = direct_reward_covariance(T, E, pi, n)
        assert c.zero(covariance-K-K.T)
        for size in (1, 2, 3):
            for indices in itertools.combinations(range(3), size):
                assert covariance.extract(indices, indices).det() >= 0


def test_chronological_antisymmetric_term_is_not_removed_by_symmetric_covariance():
    values = c.response()
    assert not c.zero(values['K']-values['K'].T)
    assert not c.zero(values['A'])
    assert values['A'][0, 2] == 0
    expected = s.Rational(603, 6666724)*s.I*s.Matrix(
        [[0, 3+s.sqrt(5), 0], [-3-s.sqrt(5), 0, -2], [0, 2, 0]])
    assert c.zero(values['A']-expected)
    assert c.zero(values['A']-values['A'].H)
    assert c.zero(values['A']+values['A'].T)


def test_direct_slow_block_taylor_coefficient_matches_complete_lag_formula():
    # Reconstruct Taylor blocks and the fast inverse independently on memory
    # tensor all36 field matrix entries, without collapsing to the K formula.
    values = c.response(policy='C')
    T, E, pi = c.number(values['T']).real, c.number(values['E']).real, c.number(values['pi']).real[:, 0]
    G0, G1, Li, _ = c.lindblad_generator(values)
    J, L = np.kron(pi[:, None], np.eye(36)), np.kron(np.ones((1, 3)), np.eye(36))
    Pi = J@L
    M0 = np.kron(T, np.eye(36))
    fast = la.inv(np.eye(108)-M0+Pi)-Pi
    M1 = np.block([[E[i, j]*Li[i] for j in range(3)] for i in range(3)])
    M2 = np.block([[E[i, j]*(Li[i]@Li[i])/2 for j in range(3)] for i in range(3)])
    expected0 = L@M1@J
    expected1 = L@M2@J+L@M1@fast@M1@J-expected0@expected0/2
    assert la.norm(expected0-G0, 2) < 3e-13
    assert la.norm(expected1-G1, 2) < 3e-13
    assert la.norm(G1, 2) > 1e-3


@pytest.mark.parametrize('args', [
    (0, F(1), F(1, 12), F(1)),
    (2, F(-1), F(1, 12), F(1)),
    (2, F(1), F(1), F(1)),
    (2, F(1), F(1, 12), F(2)),
    (2, F(1), .1, F(1)),
    (2, F(1), F(1, 12), F(1), F(-1)),
])
def test_invalid_bound_contract_fails_closed(args):
    with pytest.raises(ValueError):
        c.finite_time_bound(*args)


def test_rational_refinement_bound_and_zero_limits():
    errors = [F(c.finite_time_bound(n, F(1), F(1, 12), F(1))['unclipped_bound_exact'])
              for n in (32, 64, 128)]
    assert all(a > b > 0 for a, b in zip(errors, errors[1:]))
    assert all(4 < a/b < 5 for a, b in zip(errors, errors[1:]))
    assert 0.05774 < float(errors[1]) < 0.05775
    assert c.finite_time_bound(8, 0, F(1, 12), F(1))['unhalved_diamond_bound_exact'] == '0'
    assert c.finite_time_bound(8, 1, F(1, 12), 0)['unhalved_diamond_bound_exact'] == '0'


def test_no_C_event_does_not_disable_an_independently_declared_E_duty():
    always = c.response(0, s.Rational(1, 2), 'always')
    e = c.response(0, s.Rational(1, 2), 'E')
    assert c.zero(always['K'])
    assert s.trace(e['D']) == s.Rational(1, 4)
    assert c.zero(e['A'])


def test_literal_source_and_uniform_corrected_convergence():
    assert c.source_join()['actual_source_transfer_residual_diagnostic'] < 3e-12
    rows = c.convergence_cases()
    always = [row for row in rows if row['policy'] == 'always']
    for row in rows:
        assert row['corrected_semigroup_normalized_Choi_trace_error'] < float(F(row['unhalved_diamond_bound_exact']))
    for row in always:
        assert row['corrected_semigroup_normalized_Choi_trace_error'] < row['mean_unitary_normalized_Choi_trace_error']/15
    errors = [row['corrected_semigroup_normalized_Choi_trace_error'] for row in always]
    assert all(3.8 < a/b < 4.2 for a, b in zip(errors, errors[1:]))


def test_asymmetric_lag_orientation_away_from_native_rates():
    values = c.response(s.Rational(2, 5), s.Rational(3, 7), 'E')
    T, E, R, w, expected = (c.number(values[key]).real for key in ('T', 'E', 'R', 'w', 'K'))
    w = w[:, 0]
    summed = (np.diag(w)-np.outer(w, w))/2
    power = np.eye(3)
    for _ in range(159):
        summed += E@power@np.diag(w)-np.outer(w, w)
        power = T@power
    assert la.norm(summed-expected, 2) < 2e-12
    wrong = (np.diag(w)-np.outer(w, w))/2+R@E@np.diag(w)
    assert la.norm(wrong-expected, 2) > 1e-3
