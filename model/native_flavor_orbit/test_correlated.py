"""Independent full histories, Choi states, antiunitary and alias controls."""
import math

import numpy as np
import pytest
from scipy.linalg import expm, logm
import sympy as s

from ..common_source import native_family_record_memory as source
from . import algebra, correlated as c


def source_histories(p, q, epsilon, steps, duty, initial):
    data = source.architecture_from_factors(p, q)
    pulses = source.feedback_unitaries(epsilon, data)
    result = [(float(probability), old, np.eye(6, dtype=complex)) for old, probability in initial.items()]
    flags = ((False, False, (1-p)*(1-q)), (False, True, (1-p)*q),
             (True, False, p*(1-q)), (True, True, p*q))
    for _ in range(steps):
        next_rows = []
        for probability, old, U in result:
            for flag_c, flag_e, w in flags:
                new = data['gC'][old] if flag_c else 2
                if flag_c and flag_e:
                    new = data['gE'][new]
                if new in (0, 1):
                    new = 2
                active = duty == 'always' or duty == 'C' and flag_c or duty == 'E' and flag_e or duty == 'CE' and flag_c and flag_e
                next_rows.append((probability*float(w), new, pulses[new]@U if active else U))
        result = next_rows
    return result


def choi(histories):
    out = np.zeros((36, 36), complex)
    for probability, _, U in histories:
        vector = U.reshape(36)/math.sqrt(6)
        out += probability*np.outer(vector, vector.conj())
    return out


@pytest.mark.parametrize('policy', algebra.POLICIES)
def test_stationary_one_step_full_Choi_matrix_matches_exact_defect(policy):
    p, q = s.Rational(1, 12), s.Rational(1, 6)
    data = source.architecture_from_factors(p, q)
    initial = dict(enumerate(data['stationary']))
    epsilon = math.atan2(20, 99)
    matrix = choi(source_histories(p, q, epsilon, 1, policy, initial))
    defect = 1-float(np.trace(matrix@matrix).real)
    exact = float(c.stationary_defect(p, q, s.Rational(99, 101), policy))
    assert abs(np.trace(matrix)-1) < 1e-12
    assert abs(defect-exact) < 1e-12
    assert exact > 0 and np.linalg.matrix_rank(matrix, tol=1e-11) > 1
    # The maximally entangled probe marginal is retained, with trace I/6.
    marginal = np.einsum('abad->bd', matrix.reshape(6, 6, 6, 6))
    np.testing.assert_allclose(marginal, np.eye(6)/6, atol=1e-12, rtol=0)


@pytest.mark.parametrize('steps', (2, 3))
@pytest.mark.parametrize('policy', algebra.POLICIES)
def test_complete_histories_from_transient_initial_memory_obey_finite_bounds(steps, policy):
    p, q = s.Rational(1, 12), s.Rational(1, 6)
    epsilon = math.atan2(20, 99)
    initial = {0: s.Rational(1, 5), 3: s.Rational(3, 10), 11: s.Rational(1, 2)}
    histories = source_histories(p, q, epsilon, steps, policy, initial)
    assert len(histories) == 3*4**steps
    matrix = choi(histories)
    defect = 1-float(np.trace(matrix@matrix).real)
    bound = float(c.finite_history_bound(p, q, s.Rational(99, 101), steps, policy))
    assert 0 < bound <= defect+1e-12
    assert abs(sum(row[0] for row in histories)-1) < 1e-12


def test_correlated_joint_channel_symmetry_retains_a_reference():
    p, q = s.Rational(1, 12), s.Rational(1, 6)
    data = source.architecture_from_factors(p, q)
    pulses = [np.kron(U, np.eye(2)) for U in source.feedback_unitaries(.23, data)]
    rng = np.random.default_rng(43139)
    values = rng.normal(size=(12, 12, 4))+1j*rng.normal(size=(12, 12, 4))
    blocks = np.array([v@v.conj().T for v in values])
    blocks /= sum(np.trace(x).real for x in blocks)
    Z = np.kron(np.asarray(c.Z6, complex), np.eye(2))
    X = np.kron(np.asarray(c.X6, complex), np.eye(2))
    T = np.asarray(data['T'], float)

    def act(state):
        return np.array([pulses[new]@sum((T[new, old]*state[old] for old in range(12)),
                                        np.zeros((12, 12), complex))@pulses[new].conj().T
                         for new in range(12)])

    def theta(state):
        return np.array([Z@block.conj()@Z for block in state])

    np.testing.assert_allclose(act(theta(blocks)), theta(act(blocks)), atol=1e-12, rtol=0)
    np.testing.assert_allclose(act(np.array([X@b for b in blocks])),
                               np.array([X@b for b in act(blocks)]), atol=1e-12, rtol=0)
    assert np.linalg.norm(blocks[3]-np.trace(blocks[3]).real*sum(blocks)) > 1e-4


def test_every_actual_projector_and_reachable_word_has_the_correct_log_grading():
    p, q = s.Rational(1, 12), s.Rational(1, 6)
    histories = source_histories(p, q, .2, 3, 'always', {2: 1})
    X, Z = np.asarray(c.X6, complex), np.asarray(c.Z6, complex)
    diagonal_nonzero = False
    for probability, _, U in histories:
        if probability == 0:
            continue
        np.testing.assert_allclose(X@U, U@X, atol=1e-12, rtol=0)
        np.testing.assert_allclose(Z@U.conj()@Z, U, atol=1e-12, rtol=0)
        H = 1j*logm(U)/.6
        np.testing.assert_allclose(H, H.conj().T, atol=1e-11, rtol=0)
        np.testing.assert_allclose(Z@H.conj()@Z, -H, atol=1e-11, rtol=0)
        assert np.linalg.norm(H[:3, 3:].imag) < 1e-11
        diagonal_nonzero |= np.linalg.norm(H[:3, :3]) > 1e-3
    assert diagonal_nonzero


def test_full_memory_has_zero_rows_and_recurrent_negative_eigenvalue():
    p, q = s.Rational(1, 12), s.Rational(1, 6)
    data = source.architecture_from_factors(p, q)
    quantum = source.joint_channel(.23, data)
    assert quantum.shape == (432, 432)
    np.testing.assert_array_equal(quantum[:72], np.zeros((72, 432)))
    ports = [2, 3, 9, 10]
    a = p*(1-q)
    T = data['T'].extract(ports, ports)
    vector = (T+a*s.eye(4)).nullspace()[0]
    blocks = np.zeros((12, 36), complex)
    for port, weight in zip(ports, vector, strict=True):
        blocks[port] = float(weight)*np.eye(6).reshape(36)
    np.testing.assert_allclose(quantum@blocks.reshape(432), -float(a)*blocks.reshape(432),
                               atol=1e-12, rtol=0)
    assert c.memory_log_certificate()['rank_upper_bound'] == 360


def test_two_step_coherent_projection_is_not_a_pure_mass_service():
    result = c.exact_controls()
    assert s.sympify(result['two_step_diagonal_A02_exact']).is_positive
    first = s.diag(1, 2, 3)
    second = s.Matrix([[1, 1, s.I], [1, 2, 1], [-s.I, 1, 3]])
    H = s.kronecker_product(s.eye(2), s.I*s.im(second))+s.kronecker_product(c.X, s.re(second))
    A, S = c.hermitian_grades(H)
    assert not algebra.zero(H-s.kronecker_product(c.X, second))
    assert algebra.cp_invariant(first, second) == -34560*s.I
    assert algebra.cp_invariant(first, S) == 0
    assert not algebra.zero(A)


def test_Hilbert_Schmidt_projection_formula_against_all_hamiltonian_directions():
    projectors = source._projectors(False)
    u = c.pulse(projectors[3], s.Rational(99, 101), s.Rational(20, 101))
    v = c.pulse(projectors[2], s.Rational(99, 101), s.Rational(20, 101))
    words = [(s.Rational(2, 5), u*v), (s.Rational(3, 5), v*u)]
    H = np.asarray(c.coherent_projection(words), complex)
    phi = sum((float(w)*np.kron(np.asarray(U, complex), np.asarray(U.conjugate(), complex))
               for w, U in words), np.zeros((36, 36), complex))
    residual = phi+1j*(np.kron(H, np.eye(6))-np.kron(np.eye(6), H.T))
    directions = []
    for i in range(6):
        for j in range(i):
            E = np.zeros((6, 6), complex); E[i, j] = 1
            directions.extend((E+E.conj().T, 1j*(E-E.conj().T)))
    for i in range(5):
        E = np.zeros((6, 6), complex); E[i, i] = 1; E[5, 5] = -1
        directions.append(E)
    assert len(directions) == 35
    for E in directions:
        direction = -1j*(np.kron(E, np.eye(6))-np.kron(np.eye(6), E.T))
        assert abs(np.vdot(direction, residual).real) < 1e-11


def test_complex_projector_and_nonprincipal_alias_are_countercontrols():
    vector = s.Matrix([1, s.I, 0])/s.sqrt(2)
    U = c.pulse(vector*vector.H, s.Rational(3, 5), s.Rational(4, 5))
    assert not algebra.zero(c.Z6*U.conjugate()*c.Z6-U)
    skew = s.Matrix([[0, 1, 0], [-1, 0, 0], [0, 0, 0]])
    alias = 2*s.pi*s.kronecker_product(c.X, s.I*skew)
    np.testing.assert_allclose(expm(-1j*np.asarray(alias, complex)), np.eye(6), atol=1e-12, rtol=0)
    with pytest.raises(ValueError):
        c.hermitian_grades(alias)


@pytest.mark.parametrize('policy', algebra.POLICIES)
def test_identity_pulse_and_gate_boundaries_do_not_claim_noise(policy):
    assert c.stationary_defect(s.Rational(1, 12), s.Rational(1, 6), 1, policy) == 0
    assert c.finite_history_bound(s.Rational(1, 12), s.Rational(1, 6), 1, 3, policy) == 0
    assert c.stationary_defect(0, 0, s.Rational(3, 5), policy) == 0


@pytest.mark.parametrize('steps', (0, 1, True, 2.0))
def test_finite_history_domain_rejects_unproved_step_counts(steps):
    with pytest.raises(ValueError):
        c.finite_history_bound(s.Rational(1, 12), s.Rational(1, 6), s.Rational(3, 5), steps)


def test_all_source_default_stationary_intervals_are_strictly_positive():
    report = c.certificate()
    assert len(report['stationary_interval_witnesses']) == 12
    for row in report['stationary_interval_witnesses']:
        assert s.Rational(row['normalized_Choi_purity_defect_interval'][0]) > 0
