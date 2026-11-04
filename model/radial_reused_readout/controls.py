"""Independent complete-reference and full-Fock controls of reuse identities."""
from itertools import product

import numpy as np

from ..radial_geometry_join.source import kraus
from .reuse import history_coefficients, unread_coefficients


def quantum_histories():
    generator = np.random.default_rng(65013)
    amplitude = generator.normal(size=(4, 3))+1j*generator.normal(size=(4, 3))
    amplitude /= np.linalg.norm(amplitude)
    vector = amplitude.reshape(-1)
    rho = np.outer(vector, vector.conj())
    D = np.kron(np.diag([1., 0., 0., 1.]), np.eye(3))
    P = np.eye(12)-D
    operators = tuple(np.kron(k, np.eye(3)) for k in kraus(np.diag([1., 0., 0., 1.])))
    maximum, conditional_cross, checked = 0., 0., 0
    for rounds in range(1, 6):
        unread = np.zeros_like(rho)
        probability = 0.
        for history in product((0, 1), repeat=rounds):
            operator = np.eye(12, dtype=complex)
            for bit in history:
                operator = operators[bit]@operator
            actual = operator@rho@operator.conj().T
            row = history_coefficients(rounds, sum(history))
            target = (float(row['parity_zero'])*(P@rho@P)
                      +float(row['parity_one'])*(D@rho@D)
                      +float(row['cross_parity'])*(P@rho@D+D@rho@P))
            maximum = max(maximum, float(np.linalg.norm(actual-target)))
            conditional_cross = max(conditional_cross, float(np.linalg.norm(P@actual@D)))
            unread += actual
            probability += float(np.trace(actual).real)
            checked += 1
        multiplier = float(unread_coefficients(rounds)['parity_cross_multiplier'])
        target = P@rho@P+D@rho@D+multiplier*(P@rho@D+D@rho@P)
        maximum = max(maximum, float(np.linalg.norm(unread-target)), abs(probability-1))
    if maximum > 3e-12 or conditional_cross <= 0.01:
        raise ArithmeticError('complete conditional coherences or reference normalization failed')
    return {'full_field_local_rows': 4, 'untouched_reference_rows': 3,
            'all_ordered_histories_through_rounds': 5,
            'ordered_histories_checked': checked,
            'complete_reference_history_identities_within_3e_minus_12': True,
            'conditional_cross_parity_block_is_nonzero': True,
            'input_dephased_before_readout': False}


def _annihilator(mode):
    result = np.zeros((64, 64))
    for state in range(64):
        if state & (1 << mode):
            result[state ^ (1 << mode), state] = (-1)**(state & ((1 << mode)-1)).bit_count()
    return result


def full_fock_work():
    """All 64 Fock states of the independent three-site radial chain."""
    sigma = np.array([[0., 1.], [1., 0.]])
    h = np.zeros((6, 6), complex)
    h[:2, :2], h[4:, 4:] = 64*sigma, -64*sigma
    for left in (0, 2):
        h[left:left+2, left+2:left+4] = 0.5j*sigma
        h[left+2:left+4, left:left+2] = -0.5j*sigma
    a = [_annihilator(j) for j in range(6)]
    H, star = np.zeros((64, 64), complex), np.zeros((64, 64), complex)
    for i in range(6):
        for j in range(6):
            term = h[i, j]*(a[i].T@a[j])
            H += term
            if (i in (2, 3)) != (j in (2, 3)):
                star += term
    identity = np.eye(64)
    n = [value.T@value for value in a]
    for site in (0, 2, 4):
        H += (n[site]-identity/2)@(n[site+1]-identity/2)/65536
    D = (n[2]+n[3]-identity)@(n[2]+n[3]-identity)
    operators = kraus(D)
    _, vectors = np.linalg.eigh(H)
    ground = vectors[:, :2]@vectors[:, :2].conj().T/2
    after, adjoint = ground.copy(), H.copy()
    maximum = 0.
    for rounds in range(1, 33):
        after = sum(k@after@k.conj().T for k in operators)
        adjoint = sum(k.conj().T@adjoint@k for k in operators)
        coefficient = float(unread_coefficients(rounds)['origin_hopping_energy_injection_coefficient'])
        target = H+coefficient*star
        maximum = max(maximum, float(np.linalg.norm(adjoint-target, 2)))
        actual_work = float(np.trace(H@(after-ground)).real)
        predicted = coefficient*float(np.trace(star@ground).real)
        maximum = max(maximum, abs(actual_work-predicted))
    q_error = abs(float(np.trace(D@(after-ground)).real))
    survival = float(np.trace(2*ground@after).real)
    if maximum > 2e-10 or q_error > 2e-12 or not 0 < survival < 1-1e-8:
        raise ArithmeticError('full-Fock repeated work or spectral-disturbance control failed')
    return {'source_chain_sites': 3, 'complete_Fock_rows': 64,
            'unread_rounds_checked': 32,
            'complete_work_operator_and_state_identities_within_2e_minus_10': True,
            'parity_expectation_preserved': True,
            'actual_postmeasurement_ground_support_probability_below_one': True,
            'ground_sector_protection_assumed': False,
            'chain_control_substituted_for_27_site_theorem': False}


def certificate():
    return {'complete_quantum_histories': quantum_histories(),
            'complete_Fock_measurement_work': full_fock_work()}
