"""Independent spectral controls and mutation-sensitive record comparisons."""
import itertools

import numpy as np
import pytest
import sympy as s

from . import certificate as c


def test_record_intertwining_and_unread_counterexample():
    result = c.observation_certificate()
    assert result['finite_words_checked'] == 127
    assert result['unread_channel_equality_counterexample_tv_at_one_step'] == 1


def test_initialization_is_required_even_for_exact_dynamics():
    instrument = {'a': s.diag(1, 0), 'b': s.diag(0, 1)}
    assert c.observation_intertwining(s.eye(2), instrument, instrument)['exact']
    assert c.word_probability(instrument, [1, 0], ['a']) == 1
    assert c.word_probability(instrument, [0, 1], ['a']) == 0


def test_record_defect_bounds_full_word_law():
    p, q = s.Rational(1, 3), s.Rational(2, 5)
    first = {'a': s.Matrix([[p]]), 'b': s.Matrix([[1-p]])}
    second = {'a': s.Matrix([[q]]), 'b': s.Matrix([[1-q]])}
    defect = c.observation_intertwining(s.eye(1), first, second)['defect_l1']
    assert defect == 2*abs(p-q)
    for n in range(1, 7):
        tv = sum(abs(c.word_probability(first, [1], word)-c.word_probability(second, [1], word))
                 for word in itertools.product(first, repeat=n))/2
        assert tv <= c.total_variation_bound(0, defect, n)


@pytest.mark.parametrize('instrument', [
    {'a': s.Matrix([[s.Rational(1, 2)]])},
    {'a': s.Matrix([[-1]]), 'b': s.Matrix([[2]])},
    {'a': s.Matrix([[1.0]])},
])
def test_incomplete_negative_and_float_instruments_rejected(instrument):
    with pytest.raises(ValueError):
        c.observation_intertwining(s.eye(1), instrument, {'a': s.eye(1)})


def test_nonconserving_coarse_map_rejected():
    with pytest.raises(ValueError):
        c.observation_intertwining(s.Matrix([[2]]), {'a': s.eye(1)}, {'a': s.eye(1)})


@pytest.mark.parametrize('p,h', [(-.1, .1), (.2, .1), (.01, 0), (.01, .51), (float('nan'), .1)])
def test_neutral_bound_domain_rejected(p, h):
    with pytest.raises(ValueError):
        c.neutral_roots(p, h)


def test_neutral_all_roots_against_direct_matrix():
    rng = np.random.default_rng(1612)
    P, Q = np.diag([1., 0.]), np.array([[1., 2.], [2., 4.]])/5
    for _ in range(30):
        p, h = rng.uniform(0, 1/12), rng.uniform(.01, .5)
        U, V = np.eye(2)+(np.exp(-1j*h)-1)*P, np.eye(2)+(np.exp(-1j*h)-1)*Q
        F = np.block([[(1-p)*U, U], [p*V, np.zeros((2, 2))]])
        roots = c.neutral_roots(p, h)
        expected = np.linalg.eigvals(F)
        assert min(max(abs(roots[list(order)]-expected)) for order in itertools.permutations(range(4))) < 3e-12
        omega = min(-np.angle(roots[:2])/h)
        assert 0 <= omega <= c.phase_bound()['omega_upper']


def test_measurement_keeps_first_backaction_and_dark_residue():
    assert c.detector_certificate()['retained_joint_outcomes'] == 9


def test_exact_quantum_link_and_clifford_counterexamples():
    result = c.algebra_certificate()
    assert result['quantum_link_naive_number_eigenvalues']['1/2'] == 2
    assert result['potential_projection_zero_exact']


def test_native_passive_zero_is_not_lifted_by_positive_metric():
    result = c.passive_certificate()
    S = s.Matrix(result['S0']).applyfunc(s.sympify)
    Z = s.Matrix(result['Z']).applyfunc(s.sympify)
    assert S.det() == 0 and Z.det() != 0
    assert c.clean(Z.inv()*S).rank() == 2


def test_native_neutral_symbolic_join_and_branch_certificate():
    result = c.neutral_certificate()
    assert result['native_projector_overlap_exact'] == '1/5'
    assert result['max_eigenvalue_residual'] < c.TOLERANCE
