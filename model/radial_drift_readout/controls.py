"""Independent complete-Fock controls; floating results are diagnostics."""
from fractions import Fraction as F
from itertools import product
from math import fsum

import numpy as np
from scipy.linalg import expm

from model.radial_geometry_join.source import kraus


def annihilator(mode):
    result = np.zeros((64, 64))
    for state in range(64):
        if state & (1 << mode):
            result[state ^ (1 << mode), state] = (-1)**(state & ((1 << mode)-1)).bit_count()
    return result


def full_chain(mass=64., coupling=1/65536):
    """Literal three-site bounded radial action, including every number sector."""
    sigma = np.array([[0., 1.], [1., 0.]])
    h = np.zeros((6, 6), complex)
    h[:2, :2], h[4:, 4:] = mass*sigma, -mass*sigma
    for left in (0, 2):
        h[left:left+2, left+2:left+4] = 0.5j*sigma
        h[left+2:left+4, left:left+2] = -0.5j*sigma
    a = [annihilator(j) for j in range(6)]
    n = [x.T@x for x in a]
    identity = np.eye(64)
    H, star = np.zeros((64, 64), complex), np.zeros((64, 64), complex)
    for i, j in product(range(6), repeat=2):
        term = h[i, j]*(a[i].T@a[j])
        H += term
        if (i in (2, 3)) != (j in (2, 3)):
            star += term
    for site in (0, 2, 4):
        H += coupling*(n[site]-identity/2)@(n[site+1]-identity/2)
    D = (n[2]+n[3]-identity)@(n[2]+n[3]-identity)
    W0 = coupling*(D/2-identity/4)
    far = H-star-W0
    if np.linalg.norm(far@D-D@far) > 1e-12 or np.linalg.norm(H-H.conj().T) > 1e-12:
        raise ArithmeticError('the literal source split is not Hermitian and local')
    return H, far, star, D


def squared_norm(value):
    return fsum(float(abs(z)**2) for z in np.asarray(value).reshape(-1))


def reference_input():
    generator = np.random.default_rng(101341)
    state = generator.normal(size=(64, 2))+1j*generator.normal(size=(64, 2))
    return state/np.sqrt(squared_norm(state))


def interleaved_histories():
    H, far, star, D = full_chain()
    source = reference_input()
    operators = kraus(D)
    beta = np.sqrt(2+(1/65536)**2/16)
    if abs(np.linalg.norm(star, 2)-np.sqrt(2)) > 1e-12:
        raise ArithmeticError('independent full-chain star norm is incorrect')
    if abs(np.linalg.norm(H-far, 2)-beta) > 2e-12:
        raise ArithmeticError('independent full-chain local perturbation norm is incorrect')
    count, drift_visible, cross_visible = 0, False, False
    for rounds in (1, 2, 4):
        for elapsed in (F(1, 10000), F(1, 100)):
            U = expm(-1j*float(elapsed/rounds)*H)
            reference_end = expm(-1j*float(elapsed)*far)
            actual, reference = [source], [source]
            for _ in range(rounds):
                actual = [k@U@branch for branch in actual for k in operators]
                reference = [k@branch for branch in reference for k in operators]
            reference = [reference_end@branch for branch in reference]
            a, r = np.concatenate(actual), np.concatenate(reference)
            error = np.sqrt(squared_norm(a-r))
            if error > beta*float(elapsed)+4e-12:
                raise ArithmeticError('complete source/reference/archive violates the drift bound')
            if max(abs(squared_norm(a)-1), abs(squared_norm(r)-1)) > 4e-12:
                raise ArithmeticError('a coherent or low-weight history was discarded')
            probabilities = [squared_norm(x) for x in actual]
            frozen = [squared_norm(x) for x in reference]
            drift_visible |= max(abs(x-y) for x, y in zip(probabilities, frozen, strict=True)) > 1e-7
            cross_visible |= np.linalg.norm(D@actual[0]@actual[0].conj().T@(np.eye(64)-D)) > 1e-5
            count += len(actual)
    if not drift_visible or not cross_visible:
        raise ArithmeticError('the drift or retained cross-parity control is missing')
    return {'full_chain_Fock_rows': 64, 'untouched_reference_columns': 2,
            'complete_ordered_histories_checked': count,
            'round_counts': [1, 2, 4], 'elapsed_times': ['1/10000', '1/100'],
            'coherent_full_archive_bound_within_4e_minus_12': True,
            'normalization_uses_fsum_and_all_histories': True,
            'actual_record_probabilities_differ_from_frozen_mixture': True,
            'conditional_cross_parity_blocks_retained': True,
            'source_chain_substitutes_for_full_cube_spectral_inputs': False}


def simultaneous_control():
    """A complete controlled-parity rotation while the full field H acts."""
    H, far, _, D = full_chain()
    source = reference_input()
    initial = np.zeros((128, 2), complex)
    initial[::2] = source
    X = np.array([[0., 1.], [1., 0.]])
    angle = np.arctan(1/7)*np.eye(64)+(np.arctan(7)-np.arctan(1/7))*D
    integrated_control = -np.kron(angle, X)
    if np.linalg.norm(np.kron(far, np.eye(2))@integrated_control-integrated_control@np.kron(far, np.eye(2))) > 1e-12:
        raise ArithmeticError('the control did not commute with the remote-source reference')
    beta = np.sqrt(2+(1/65536)**2/16)
    for elapsed in (F(1, 10000), F(1, 100)):
        actual = expm(-1j*(float(elapsed)*np.kron(H, np.eye(2))+integrated_control))@initial
        reference = expm(-1j*float(elapsed)*np.kron(far, np.eye(2)))@expm(-1j*integrated_control)@initial
        if np.sqrt(squared_norm(actual-reference)) > float(elapsed)*beta+4e-12:
            raise ArithmeticError('simultaneous control/source drift exceeds the local bound')
        if abs(squared_norm(actual)-1) > 4e-12:
            raise ArithmeticError('a complete pointer output was omitted')
    return {'complete_chain_plus_pointer_rows': 128, 'untouched_reference_columns': 2,
            'elapsed_times': ['1/10000', '1/100'],
            'simultaneous_local_control_and_full_H_tested': True,
            'complete_coherent_pointer_output_retained': True,
            'local_Duhamel_bound_within_4e_minus_12': True,
            'this_control_is_an_executed_40_pulse_native_word': False,
            'purpose': 'independent full-CAR continuous-control theorem check'}


def certificate():
    return {'interleaved_actual_source': interleaved_histories(),
            'simultaneous_continuous_control': simultaneous_control()}
