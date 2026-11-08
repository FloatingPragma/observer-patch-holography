"""Complete native regional instruments, retaining every failure record."""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import math

import numpy as np
from scipy.linalg import expm

from model.centered_native_preparation import tree
from model.certified_reference import certificate as cert
from model.native_spin_observer import algebra, circuits, execute
from . import regional


def physical_label(label):
    if len(label) != 2 or any(axis not in 'IXYZ' for axis in label):
        raise ValueError('complete two-spin Pauli label required')
    return ''.join(axis+('X' if axis in 'XY' else 'I') for axis in label)


def bound_word(template):
    word = list(circuits.emit(template))
    end, angle, variables = cert.consume(word, 0, template)
    if end != len(word) or variables:
        raise ArithmeticError('native instrument has an unbound control')
    port = tree.port_bound(Counter((e['kind'], e['angle']) for e in word))
    return word, angle, port


def complete(template, bits, initial, target, *, prefix_word=(), prefix_error=F(0)):
    word, angle, port = bound_word(template)
    if type(prefix_error) is not F or prefix_error < 0:
        raise ValueError('nonnegative exact preparation error required')
    word = list(prefix_word)+word
    # The contraction kernel accepts32 columns. Both exterior blocks are run
    # without normalization; their squared amplitudes add for the same records.
    blocks = [execute.complete(word, bits, np.ascontiguousarray(initial[:, j:j+32]))
              for j in range(0, initial.shape[1], 32)]
    first = blocks[0]
    success = np.hstack([block.success for block in blocks])
    weights = sum((block.failure_weights for block in blocks),
                  np.zeros_like(first.failure_weights))
    run = execute.Run(bits, initial, success, weights, first.offsets, first.word_sha256,
                      first.native_failure_slots,
                      abs(float(np.vdot(success, success).real)+float(weights.sum())
                          -float(np.vdot(initial, initial).real)),
                      sum(block.rounded_negative_count for block in blocks),
                      sum(block.rounded_negative_mass for block in blocks))
    distance = math.sqrt(float(np.linalg.norm(run.success-target)**2)
                         + run.failure_probability)
    if distance > float(prefix_error+angle+port)+2e-11:
        raise ArithmeticError('literal full native instrument missed its ideal action')
    return run, {'pulses': len(word), 'word_sha256': cert.word_digest(word),
                 'angle_error': str(angle), 'full_generator_error': str(port),
                 'prefix_preparation_error_upper': str(prefix_error),
                 'complete_isometry_error_bound': str(prefix_error+angle+port),
                 'complete_output_distance_diagnostic': distance,
                 'native_first_failure_slots': run.native_failure_slots,
                 'native_first_failure_weights': run.failure_weights.tolist(),
                 'complete_first_failure_probability': run.failure_probability,
                 'normalization_residual': run.normalization_residual,
                 'rounded_negative_count': run.rounded_negative_count,
                 'rounded_negative_mass': run.rounded_negative_mass,
                 'postselected_or_renormalized': False}


def probe_template(label):
    physical_label(label)
    basis = []
    for site, axis in enumerate(label):
        if axis == 'X':
            basis += circuits.spin_rotation('Y', cert._pi(-1, 4), site, 2)
        if axis == 'Y':
            basis += circuits.spin_rotation('X', cert._pi(1, 4), site, 2)
    parity = [event for site, axis in enumerate(label) if axis != 'I'
              for event in cert.cnot(2*site, 4)]
    return cert.inverse(basis)+parity+basis


def regional_purification():
    """Certified integer ground, with all64 exterior columns retained."""
    data = regional.state()
    v = np.asarray(data['vector'], dtype=float)/math.sqrt(data['vector_norm_squared'])
    return v.reshape(64, 4).T


@lru_cache(maxsize=1)
def detector():
    selector_template = [('Y', (bit,), cert._pi(-1, 4)) for bit in range(4)]
    blank = np.zeros((16, 1), complex)
    blank[0, 0] = 1
    coins, coin_report = complete(selector_template, 4, blank, np.ones((16, 1))/4)
    coin_probabilities = np.abs(coins.success[:, 0])**2
    inclusion = algebra.inclusion(2)
    reference = regional_purification()
    initial = np.zeros((32, 64), complex)
    initial[:16] = inclusion@reference
    records = []
    total = coins.failure_probability
    for index, label in enumerate(regional.coefficients(regional.logarithm()['matrix'])):
        p = algebra.tensor(label)
        projectors = [(np.eye(4)+sign*p)/2 for sign in (1, -1)]
        target = np.vstack([inclusion@x@reference for x in projectors])
        run, report = complete(probe_template(label), 5, initial, target)
        actual = [float(np.linalg.norm(run.success[16*j:16*(j+1)])**2) for j in range(2)]
        ideal = [float(np.linalg.norm(target[16*j:16*(j+1)])**2) for j in range(2)]
        report.update({'pauli': label, 'selector_success_probability': float(coin_probabilities[index]),
                       'pointer_success_probabilities': actual,
                       'ideal_pointer_probabilities': ideal,
                       'joint_success_probabilities': [float(coin_probabilities[index])*x for x in actual],
                       'reserved_six_level_pointer_probabilities': [0.0]*4,
                       'retained_exterior_columns': 64})
        total += float(coin_probabilities[index])*(sum(actual)+run.failure_probability)
        records.append(report)
    if abs(total-1) > 2e-11:
        raise ArithmeticError('selector and all conditional outcomes lost normalization')
    error = (regional.state()['ground_vector_error']
             + F(coin_report['complete_isometry_error_bound'])
             + max(F(r['complete_isometry_error_bound']) for r in records))
    return {'selector': coin_report, 'settings': records,
            'actual_selected_ground_complete_TV_error_upper': str(error),
            'complete_normalization_residual': abs(total-1),
            'all_32_setting_pointer_records_retained': True,
            'identity_minus_record_is_exact_zero_in_ideal_instrument': True,
            'eight_spin_preparation_reexecuted_here': False}


def modular_template(time=F(1, 64), subdivisions=8):
    if type(time) is not F or type(subdivisions) is not int or subdivisions < 1:
        raise ValueError('exact modular time and positive subdivision count required')
    log = regional.logarithm()
    coefficients = regional.coefficients(log['matrix'])
    terms = [(label, value) for label, value in coefficients.items() if value]
    half = [event for label, value in terms
            for event in cert.pauli(physical_label(label), -time*value/(2*subdivisions))]
    # Reverse the order of Pauli exponentials, preserving each exponential's circuit.
    back = [event for label, value in reversed(terms)
            for event in cert.pauli(physical_label(label), -time*value/(2*subdivisions))]
    template = (half+back)*subdivisions
    norm_sum = sum(abs(value) for _, value in terms)
    splitting = norm_sum**3*abs(time)**3/(3*subdivisions**2)
    return template, log, splitting


@lru_cache(maxsize=1)
def modular_flow():
    time, subdivisions = F(1, 64), 8
    template, log, splitting = modular_template(time, subdivisions)
    inclusion = algebra.inclusion(2)
    initial = inclusion/2  # all four code columns, with a retained Choi reference
    step = np.eye(4, dtype=complex)
    coefficients = regional.coefficients(log['matrix'])
    terms = [(label, value) for label, value in coefficients.items() if value]
    for label, value in terms+list(reversed(terms)):
        step = expm(-1j*float(time*value/(2*subdivisions))*algebra.tensor(label))@step
    ideal_word = np.linalg.matrix_power(step, subdivisions)
    run, report = complete(template, 4, initial, inclusion@ideal_word/2)
    matrix = np.asarray(log['matrix'], dtype=float)
    exact_polynomial = expm(-1j*float(time)*matrix)
    diagnostic = float(np.linalg.norm(ideal_word-exact_polynomial, ord=2))
    if diagnostic > float(splitting)+2e-13:
        raise ArithmeticError('symmetric product exceeded its analytic operator bound')
    total = (F(report['complete_isometry_error_bound'])+splitting
             + abs(time)*log['operator_error_upper'])
    report.update({'time': str(time), 'subdivisions': subdivisions,
                   'symmetric_product_operator_error_upper': str(splitting),
                   'symmetric_product_error_diagnostic': diagnostic,
                   'true_selected_regional_modular_flow_error_upper': str(total),
                   'constant_modular_phase_retained': True,
                   'retained_reference_columns': 4,
                   'arbitrary_exterior_validity_from_operator_bound': True,
                   'physical_source_clock_identified_with_modular_parameter': False})
    return report
