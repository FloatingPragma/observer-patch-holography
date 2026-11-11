"""Independent full-channel and invariant-graph controls."""
import copy

import numpy as np
import pytest
import scipy.linalg as la
import sympy as sp

from . import response as r


def test_complete_correlation_sum_and_noise_positivity():
    source = r.exact_source()
    T, pi = r.numeric(source['T']).real, r.numeric(source['pi']).real.ravel()
    Pi = np.outer(pi, np.ones(12))
    direct = (np.diag(pi)-np.outer(pi, pi))/2
    power = np.eye(12)
    for _ in range(30):
        power = T@power
        direct += (power-Pi)@np.diag(pi)
    assert np.max(abs(direct-r.numeric(source['C']))) < 2e-14
    assert np.linalg.eigvalsh(r.numeric(source['S'])).min() > -1e-15
    assert source['S'].rank() == 3
    assert any(x != 0 for x in source['J'])


def test_hamiltonian_noise_split_and_pairing_grade():
    G0, G2, noise, H = r.effective_generators()
    source = r.exact_source()
    X = np.kron([[0, 1], [1, 0]], np.eye(3))
    assert np.max(abs(X@H-H@X)) < 1e-14
    assert np.max(abs(H[:3, 3:])) == 0
    assert np.max(abs(r.hamiltonian_projection(noise))) < 1e-14
    assert np.max(abs(r.hamiltonian_projection(G2)-H)) < 1e-14
    assert np.max(abs(r.hamiltonian_projection(G0)-np.kron([[0, 1], [1, 0]], r.numeric(source['Y'])))) < 1e-14


def test_third_order_against_full_invariant_graph_recurrence():
    source = r.exact_source()
    T, pi = r.numeric(source['T']), r.numeric(source['pi']).ravel()
    X = np.array([[0, 1], [1, 0]])
    L = [r.commutator_generator(np.kron(X, r.numeric(P))) for P in source['projectors']]
    M = [np.block([[T[j, i]*np.linalg.matrix_power(L[j], k)/float(sp.factorial(k))
                   for i in range(12)] for j in range(12)]) for k in range(1, 4)]
    E0 = np.vstack([p*np.eye(36) for p in pi])
    marginal = np.hstack([np.eye(36) for _ in range(12)])
    Qsolve = np.kron(r.numeric(source['fundamental']-source['stationary']), np.eye(36))
    A1 = marginal@M[0]@E0
    E1 = Qsolve@M[0]@E0
    A2 = marginal@(M[0]@E1+M[1]@E0)
    E2 = Qsolve@(M[0]@E1+M[1]@E0-E1@A1)
    A3 = marginal@(M[0]@E2+M[1]@E1+M[2]@E0)
    G2 = A2-A1@A1/2
    G3 = A3-(A1@A2+A2@A1)/2+A1@A1@A1/3
    assert np.max(abs(G2-r.effective_generators()[1])) < 1e-14
    Y3, coefficients = r.third_order_pairing()
    assert len(coefficients) == 64
    expected = np.kron(X, r.numeric(Y3))
    assert np.max(abs(r.hamiltonian_projection(G3)-expected)) < 2e-14


def test_complete_reduced_channel_is_trace_preserving_and_positive():
    E = r.complete_reduced_channel(4)
    identity = np.eye(6).reshape(-1)
    assert np.max(abs(E@identity-identity)) < 1e-14
    assert np.max(abs(identity@E-identity)) < 1e-14
    choi = E.reshape(6, 6, 6, 6).transpose(0, 2, 1, 3).reshape(36, 36)
    assert np.max(abs(choi-choi.conj().T)) < 1e-14
    assert la.eigvalsh(choi).min() > -1e-14


def test_finite_clock_correction_converges_and_is_not_an_exact_replacement():
    G0, G2, _, _ = r.effective_generators()
    Y3, _ = r.third_order_pairing()
    errors, odd_errors, product_errors = [], [], []
    for h in (1/8, 1/16):
        G, residual = r.slow_generator(h)
        assert residual < 1e-12
        errors.append(la.norm(G-G0-h*G2, 2))
        H = r.hamiltonian_projection(G)
        odd_errors.append(la.norm((H[:3, 3:]+H[3:, :3])/2-r.numeric(r.exact_source()['Y'])-h*h*r.numeric(Y3), 2))
        exact = r.complete_reduced_channel(round(1/h))
        product_errors.append(la.norm(exact-la.expm(G0+h*G2), 2))
    assert 3.8 < errors[0]/errors[1] < 4.2
    assert 14 < odd_errors[0]/odd_errors[1] < 18
    assert 3 < product_errors[0]/product_errors[1] < 5
    assert min(errors+odd_errors+product_errors) > 0


def test_receipt_mutations_are_rejected():
    from tools import strict_json
    report = strict_json.load(r.OUTPUT)
    altered = copy.deepcopy(report)
    altered['scope']['computed_J_is_a_physical_Yukawa'] = True
    with pytest.raises(ValueError):
        r.compare(altered, report)
    altered = copy.deepcopy(report)
    altered['ordered_correlation_J_exact'][0][1] = '0'
    with pytest.raises(ValueError):
        r.compare(altered, report)
