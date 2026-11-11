"""Derive the source memory correction without an inserted coherent phase.

The twelve-port transition and its logical projector pulses are the declared
native_family_record_memory experiment. The expansion uses its complete
correlated channel, rather than independent samples of its stationary mean.
Its clock and feedback service do not identify physical fermion fields.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import scipy.linalg as la
import sympy as sp

from model.common_source import native_family_record_memory as memory
from model.particle_memory_response import certificate as memory_theorem

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')


def numeric(value):
    return np.asarray(value, dtype=complex)


def commutator_generator(H):
    H = np.asarray(H, complex)
    I = np.eye(len(H))
    return -1j * (np.kron(H, I) - np.kron(I, H.T))


@lru_cache(maxsize=1)
def exact_source():
    data = memory.architecture()
    T, pi = sp.Matrix(data['T']), sp.Matrix(data['stationary'])
    stationary = pi * sp.ones(1, 12)
    fundamental = (sp.eye(12) - T + stationary).inv()
    # sum_{lag>=1}(T**lag-Pi)=fundamental-I; no tail is truncated.
    C = (fundamental - sp.eye(12)/2 - stationary/2) * sp.diag(*pi)
    S, A = (C+C.T)/2, (C-C.T)/2
    P = data['projectors']
    J = sum((A[j, i] * (P[j]*P[i]-P[i]*P[j])
             for j in range(12) for i in range(12)), sp.zeros(3))/(2*sp.I)
    J = J.applyfunc(sp.simplify)
    support = [j for j in range(12) if pi[j] != 0]
    principal_minors = []
    for count in range(1, len(support)+1):
        for subset in itertools.combinations(support, count):
            determinant = sp.factor(S.extract(subset, subset).det())
            if determinant.is_nonnegative is not True:
                raise ArithmeticError('integrated symmetric covariance is not positive')
            principal_minors.append((subset, determinant))
    if (T*pi != pi or C*sp.ones(12, 1) != sp.zeros(12, 1)
            or sp.ones(1, 12)*C != sp.zeros(1, 12)
            or any(sp.simplify(x) != 0 for x in J-J.H)
            or any(sp.simplify(x) != 0 for x in J+J.conjugate())):
        raise ArithmeticError('source stationary or coherent identities failed')
    prior = memory_theorem.response()
    collapse = sp.zeros(3, 12)
    collapse[0, 2] = collapse[1, 3] = collapse[2, 9] = collapse[2, 10] = 1
    if (any(sp.simplify(x) != 0 for x in prior['A']-J)
            or any(sp.simplify(x) != 0 for x in prior['Y']-data['Y'])
            or any(sp.simplify(x) != 0 for x in prior['K']-collapse*C*collapse.T)):
        raise ArithmeticError('the owning all-lag theorem does not match this source')
    return dict(T=T, pi=pi, stationary=stationary, fundamental=fundamental,
                C=C, S=S, A=A, J=J, Y=data['Y'], projectors=P,
                support=support, principal_minors=principal_minors)


def effective_generators():
    source = exact_source()
    X = np.array([[0, 1], [1, 0]], complex)
    V = [np.kron(X, numeric(P)) for P in source['projectors']]
    L = [commutator_generator(H) for H in V]
    G0 = commutator_generator(np.kron(X, numeric(source['Y'])))
    G2 = sum((float(source['C'][j, i]) * L[j]@L[i]
              for j in range(12) for i in range(12)), np.zeros((36, 36), complex))
    dissipator = sum((float(source['S'][j, i]) * L[j]@L[i]
                     for j in range(12) for i in range(12)), np.zeros((36, 36), complex))
    Hcorr = np.kron(np.eye(2), numeric(source['J']))
    if np.max(abs(G2-dissipator-commutator_generator(Hcorr))) > 2e-14:
        raise ArithmeticError('ordered covariance decomposition failed')
    return G0, G2, dissipator, Hcorr


@lru_cache(maxsize=1)
def third_order_pairing():
    """Exact X-coefficient of the order-h-squared slow Hamiltonian.

    All correlation lags enter through the fundamental matrix. The result
    is a source-derived correction to the specified logical feedback service.
    It does not turn that service into an identified fermion-Higgs vertex.
    """
    source = exact_source()
    pi, F, Pi = source['pi'], source['fundamental'], source['stationary']
    D = F-sp.eye(12)
    e = (F-Pi)*sp.diag(*pi)
    Te = D*sp.diag(*pi)
    De = D*e
    B = Te+sp.diag(*pi)/2
    coefficients = {}
    Y3 = sp.zeros(3)
    P = source['projectors']
    for a, b, c in itertools.product(source['support'], repeat=3):
        delta_ab, delta_bc = int(a == b), int(b == c)
        value = (D[a, b]*Te[b, c]-De[a, b]*pi[c]
                 +D[a, b]*pi[b]*delta_bc/2
                 +delta_ab*Te[a, c]/2+pi[a]*delta_ab*delta_bc/6
                 -(pi[a]*B[b, c]+B[a, b]*pi[c])/2
                 +pi[a]*pi[b]*pi[c]/3)
        coefficients[a, b, c] = value
        word = (-(P[a]*P[b]*P[c]+P[c]*P[b]*P[a])/2
                -(P[a]*sp.trace(P[b]*P[c])+P[b]*sp.trace(P[c]*P[a])
                  +P[c]*sp.trace(P[b]*P[a]))/3)
        Y3 += value*word
    Y3 = Y3.applyfunc(sp.simplify)
    if Y3 != Y3.T or any(sp.im(x) != 0 for x in Y3):
        raise ArithmeticError('third-order pairing must be real symmetric')
    return Y3, coefficients


def slow_generator(step):
    """Numerical diagnostic of the full432-dimensional channel's slow sector.

    Marginal coordinates fix the invariant-graph gauge. This extracted slow
    map is not asserted to be a complete positive channel for arbitrary
    marginal inputs. The original memory channel is complete positive.
    """
    if type(step) not in (int, float) or not 0 < step <= 1/4:
        raise ValueError('a positive coherent step at most1/4 is required')
    M = memory.joint_channel(float(step))
    schur, frame, dimension = la.schur(M, output='complex', sort=lambda z: abs(z-1) < 1/2)
    if dimension != 36:
        raise ArithmeticError('the complete slow space must have dimension36')
    frame = frame[:, :dimension]
    marginal_frame = frame.reshape(12, 36, 36).sum(axis=0)
    inverse_frame = la.solve(marginal_frame, np.eye(36))
    reduced = marginal_frame@schur[:36, :36]@inverse_frame
    generator = la.logm(reduced)/step
    embedding = frame@inverse_frame
    error = np.max(abs(M@embedding-embedding@reduced))
    return generator, float(error)


def complete_reduced_channel(steps, time=1):
    if type(steps) is not int or steps < 1 or not 0 < time/steps <= 1/4:
        raise ValueError('positive steps with time/steps at most1/4 required')
    pi = numeric(exact_source()['pi']).real.ravel()
    M = memory.joint_channel(time/steps)
    embedded = np.vstack([p*np.eye(36) for p in pi])
    evolved = np.linalg.matrix_power(M, steps)@embedded
    return evolved.reshape(12, 36, 36).sum(axis=0)


def hamiltonian_projection(generator):
    """Canonical traceless Hamiltonian component in Hilbert-Schmidt coordinates.

    This algebraic projection does not certify that generator is GKSL.
    Traceless-jump dissipators have zero contribution to this projection.
    """
    value = np.asarray(generator, complex)
    if value.shape != (36, 36):
        raise ValueError('the complete six-state Liouville operator is required')
    contracted = np.einsum('ijkj->ik', value.reshape(6, 6, 6, 6))
    return 1j*(contracted-contracted.conj().T)/12


def source_bindings():
    paths = set(memory.feedback.transport.source_bindings())
    paths.update(('model/common_source/native_family_record_memory.py',
                  'model/common_source/test_native_family_record_memory.py',
                  'model/common_source/native_gauge_object.py',
                  'model/common_source/native_repair_mass.py',
                  'model/common_source/native_repair_mass_readout.py',
                  'model/particle_memory_response/certificate.py',
                  'model/particle_memory_response/PROOF.md'))
    paths.update(str(p.relative_to(ROOT)) for p in Path(__file__).parent.glob('*')
                 if p.suffix in ('.py', '.md'))
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            for name in sorted(paths)}


def report():
    source = exact_source()
    G0, G2, dissipator, Hcorr = effective_generators()
    Y3, _ = third_order_pairing()
    rows = []
    for step in (1/4, 1/8, 1/16, 1/32):
        G, invariant_error = slow_generator(step)
        H = hamiltonian_projection(G)
        odd_coefficient = (H[:3, 3:]+H[3:, :3])/2
        odd_eigenvalues = la.eigvalsh(odd_coefficient)
        rows.append({'step': step,
                     'slow_invariant_graph_error': invariant_error,
                     'generator_minus_mean_operator_norm': float(la.norm(G-G0, 2)),
                     'second_order_remainder_operator_norm': float(la.norm(G-G0-step*G2, 2)),
                     'remainder_divided_by_step_squared': float(la.norm(G-G0-step*G2, 2)/step**2),
                     'coherent_X_coefficient_eigenvalues_diagnostic': odd_eigenvalues.tolist(),
                     'coherent_X_coefficient_ratios_diagnostic': (odd_eigenvalues/odd_eigenvalues[0]).tolist(),
                     'coherent_X_coefficient_shift_norm': float(la.norm(odd_coefficient-numeric(source['Y']), 2)),
                     'coherent_X_third_order_remainder_norm': float(la.norm(odd_coefficient-numeric(source['Y'])-step**2*numeric(Y3), 2))})
    finite = []
    for steps in (4, 8, 16):
        E = complete_reduced_channel(steps)
        finite.append({'steps': steps, 'coherent_time': 1,
                       'mean_channel_HS_induced_error': float(la.norm(E-la.expm(G0), 2)),
                       'corrected_channel_HS_induced_error': float(la.norm(E-la.expm(G0+G2/steps), 2))})
    eigenvalues = la.eigvalsh(numeric(source['Y']))
    return {
        'schema': 'oph.yukawa_source_feedback.v1',
        'source': 'native_family_record_memory complete reset2 baseline',
        'native_event_probabilities': {'C': '1/12', 'E': '1/6'},
        'complete_raw_labels_per_input': 128,
        'memory_inputs': 12,
        'logical_work_dimension': 6,
        'complete_CQ_superoperator_dimension': 432,
        'source_mean_Y_exact': [[str(x) for x in row] for row in source['Y'].tolist()],
        'source_mean_Y_eigenvalues_diagnostic': eigenvalues.tolist(),
        'source_mean_Y_ratios_diagnostic': (eigenvalues/eigenvalues[0]).tolist(),
        'ordered_correlation_J_exact': [[str(x) for x in row] for row in source['J'].tolist()],
        'ordered_correlation_J_HS_squared_exact': str(sp.simplify(sp.trace(source['J']**2))),
        'ordered_correlation_J_operator_norm_diagnostic': float(la.norm(numeric(source['J']), 2)),
        'integrated_covariance_rank': int(source['S'].rank()),
        'integrated_covariance_nonnegative_principal_minors': [
            {'indices': list(indices), 'value': str(value)} for indices, value in source['principal_minors']],
        'coherent_correction_operator': 'I_2 tensor J',
        'mean_pairing_operator': 'X tensor Y',
        'first_order_correction_to_Y': 'zero: Tr_Nambu[(X tensor I)(I tensor J)]/2=0',
        'second_order_pairing_correction_Y3_exact': [[str(x) for x in row] for row in Y3.tolist()],
        'dissipator_HS_induced_norm_diagnostic': float(la.norm(dissipator, 2)),
        'slow_generator_controls': rows,
        'complete_reduced_channel_controls': finite,
        'reused_finite_theorem': {
            'owner': 'model/particle_memory_response/PROOF.md',
            'twelve_port_covariance_exactly_intertwines_three_ray_theorem': True,
            'stationary_product_unhalved_diamond_bounds': [
                {'steps': steps, **memory_theorem.finite_time_bound(steps, 1, '1/12', 1)}
                for steps in (128, 512, 2048)],
            'these_bounds_control_G0_plus_hG2_not_third_order_truncation': True},
        'scope': {
            'empirical_mass_or_mixing_input': False,
            'phase_inserted_as_function_of_probability': False,
            'all_correlation_lags_summed_exactly': True,
            'all_original_records_retained_by_source': True,
            'stationary_product_memory_for_reduced_controls': True,
            'feedback_pulse_and_clock_are_declared_services': True,
            'numerical_controls_are_rigorous_error_intervals': False,
            'HS_induced_error_is_diamond_error': False,
            'finite_step_slow_generator_claimed_GKSL': False,
            'computed_J_is_a_physical_Yukawa': False,
            'physical_particle_spectrum_derived': False},
        'source_bindings': source_bindings()}


def compare(actual, expected, path='receipt'):
    if type(actual) is not type(expected):
        raise ValueError(path+' type differs')
    if isinstance(expected, dict):
        if actual.keys() != expected.keys():
            raise ValueError(path+' keys differ')
        for key in expected:
            compare(actual[key], expected[key], path+'.'+key)
    elif isinstance(expected, list):
        if len(actual) != len(expected):
            raise ValueError(path+' length differs')
        for index, (a, b) in enumerate(zip(actual, expected, strict=True)):
            compare(a, b, path+f'[{index}]')
    elif isinstance(expected, float):
        if not np.isfinite(actual) or abs(actual-expected) > 2e-9+2e-7*abs(expected):
            raise ValueError(path+' numerical replay differs')
    elif actual != expected:
        raise ValueError(path+' exact value differs')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = report()
    if args.write:
        OUTPUT.write_text(json.dumps(result, indent=2)+'\n')
    if args.check:
        from tools import strict_json
        compare(result, strict_json.load(OUTPUT))
    print(json.dumps({'ordered_correction': result['ordered_correlation_J_exact'],
                      'slow_generator_controls': result['slow_generator_controls'],
                      'physical_particle_spectrum_derived': False}, indent=2))


if __name__ == '__main__':
    main()
