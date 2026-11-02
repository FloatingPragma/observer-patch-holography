"""Actual complete six-level execution of the two-site joint instrument."""
from collections import Counter
from fractions import Fraction as F
import math
import numpy as np
from model.centered_native_preparation import tree
from model.certified_reference import certificate as cert, interval as iv
from model.native_spin_observer import circuits as spin, execute
from model.native_spin_observer.precision import sqrt_interval
from . import circuits, algebra


STEPS = (32, 64)


def exact_inputs():
    a, b = sqrt_interval(5, 100)
    rlo, rhi = a - 2, b - 2
    raw = (rlo + rhi) * (1 << 47)
    rounded = (raw + F(1, 2)).numerator // (raw + F(1, 2)).denominator
    r = F(rounded, 1 << 48)
    ground_error = 2 * max(abs(r - rlo), abs(r - rhi))
    elo, ehi = F(-1, 2) - b / 4, F(-1, 2) - a / 4
    origin = (elo + ehi) / 2
    return (rounded, 1 << 48, -(1 << 48), -rounded), origin, ground_error, (ehi - elo) / 2


def plan():
    values, origin, state_error, origin_error = exact_inputs()
    word = []
    angle_error = F()
    for bit, raw, coefficients in tree.design(values, 2):
        a, b = tree.level_error(values, bit, raw, coefficients, 2)
        angle_error += a + b
        word += [dict(event, targets=[2 * j for j in event['targets']])
                 for event in tree.level_word(bit, coefficients, 2)]
    angle_error += 2 * (2 * iv.pi_angle_error(F(math.pi / 2), F(1, 2)) + iv.pi_angle_error(F(-math.pi), F(-1)))
    suffix = spin.encode(2)
    # Exact local bar-Z, including its scalar phase; this is a unitary probe.
    suffix += [('Z', (0,), cert._pi(1, 2)), ('I', (0,), cert._pi(-1, 2))]
    suffix += circuits.dimer_clock(origin, STEPS)
    actual = list(spin.emit(suffix))
    consumed, error, variables = cert.consume(actual, 0, suffix)
    if consumed != len(actual) or variables:
        raise ArithmeticError('unbound joint-instrument source word')
    word += actual
    angle_error += error
    native = tree.port_bound(Counter((e['kind'], e['angle']) for e in word))
    # ||[h,V]|| <= 2||h||||V|| = 2JB; J=1, B=1/4.
    splitting = sum((F(1, 4) * F(1 << j) ** 2 / n for j, n in enumerate(STEPS)), F())
    algebra.dimer_common_kernel()
    # After the ideal ground and Z probe, every controlled factor is the
    # identity on this common kernel. This is not a uniform operator claim.
    total = state_error + angle_error + native + 3 * origin_error
    report = {'work_bits': 7, 'complete_work_rows': 128, 'clock_size': 4,
              'steps': list(STEPS), 'pulses': len(word), 'word_sha256': cert.word_digest(word),
              'ground_vector_error': str(state_error), 'origin_error': str(origin_error),
              'literal_angle_error': str(angle_error), 'native_generator_error': str(native),
              'arbitrary_input_Hamiltonian_splitting_bound': str(splitting),
              'prepared_transverse_state_splitting_error': '0',
              'complete_isometry_error': str(total), 'unhalved_trace_error': str(2 * total),
              'success_renormalized': False, 'ground_postselected': False,
              'fixed_binary64_error_sent_to_zero': False}
    initial = np.zeros((128, 1), complex)
    initial[0, 0] = 1
    return word, initial, report


def ideal_target():
    r = math.sqrt(5) - 2
    g = np.array([r, 1, -1, -r], complex)
    g /= np.linalg.norm(g)
    psi = np.array([1, -1, 1, -1]) * g
    swap = np.array([0, 2, 1, 3])
    phase = np.array([1, -1, -1, 1])
    reflection = phase * psi[swap]
    gap = F(1, 2) + math.sqrt(5) / 4
    result = np.zeros((128, 1), complex)
    probabilities = np.zeros((2, 4))
    for k in range(2):
        q = (psi + (-1) ** k * reflection) / 2
        for j in range(4):
            amplitude = sum(np.exp(-1j * n * (gap - 2 * math.pi * j / 4)) for n in range(4)) / 4
            for mask in range(4):
                encoded = sum(1 << (2 * site + (0 if (mask >> site) & 1 else 1)) for site in range(2))
                result[encoded + (k << 4) + (j << 5), 0] = amplitude * q[mask]
            probabilities[k, j] = abs(amplitude) ** 2 * float(np.vdot(q, q).real)
    return result, probabilities


def run():
    word, initial, report = plan()
    target, probabilities = ideal_target()
    result = execute.complete(word, 7, initial)
    distance = math.sqrt(float(np.linalg.norm(result.success - target) ** 2) + result.failure_probability)
    if distance > float(F(report['complete_isometry_error'])) + execute.ADMISSION:
        raise ArithmeticError('native joint instrument exceeds its complete source bound')
    observed = np.zeros((2, 4))
    for k in range(2):
        for j in range(4):
            start = (k << 4) + (j << 5)
            observed[k, j] = float(np.linalg.norm(result.success[start:start + 16]) ** 2)
    report.update(ideal_joint_probabilities=probabilities.tolist(),
                  native_success_joint_probabilities=observed.tolist(),
                  complete_output_distance_diagnostic=distance,
                  full_first_failure_probability=result.failure_probability,
                  normalization_residual=result.normalization_residual,
                  rounded_negative_count=result.rounded_negative_count,
                  rounded_negative_mass=result.rounded_negative_mass,
                  native_failure_slots=result.native_failure_slots)
    return report, {'success': result.success, 'failure_weights': result.failure_weights}
