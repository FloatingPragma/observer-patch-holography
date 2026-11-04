"""Full-site native modular words with rational logarithm and port bounds."""
from collections import Counter
from fractions import Fraction as F
import math

import numpy as np

from ..certified_reference import certificate as cert
from ..centered_native_preparation import tree
from ..native_spin_observer import circuits, execute, precision
from .regional import rational


def _unit_log(value, terms):
    if not 1 <= value <= 2:
        raise ValueError('logarithm range reduction failed')
    z = (value-1)/(value+1)
    partial = 2*sum((z**(2*j+1)/(2*j+1) for j in range(terms+1)), F())
    remainder = 2*z**(2*terms+3)/((2*terms+3)*(1-z*z))
    return partial, partial+remainder


def log_odds_interval(q, terms=48):
    q = rational(q)
    if (not F(1, 4096) <= q <= 1-F(1, 4096)
            or type(terms) is not int or terms < 0):
        raise ValueError('certified parity interval and nonnegative integer series order required')
    ratio = (1-q)/q
    sign = 1 if ratio >= 1 else -1
    if sign < 0:
        ratio = 1/ratio
    power = 0
    while ratio > 2:
        ratio /= 2
        power += 1
    lo, hi = _unit_log(ratio, terms)
    lo2, hi2 = _unit_log(F(2), terms)
    lo, hi = lo+power*lo2, hi+power*hi2
    return (lo, hi) if sign > 0 else (-hi, -lo)


def template(angle):
    return cert.pauli('II', angle)+cert.pauli('ZZ', angle)


def word_and_bound(q=F(1, 4), time=F(1, 2), terms=48):
    q, time = rational(q), rational(time)
    if abs(time) > F(1, 2):
        raise ValueError('the certified modular-time interval is [-1/2,1/2]')
    lo, hi = log_odds_interval(q, terms)
    midpoint = (lo+hi)/2
    word = list(circuits.emit(template(-time*midpoint/2)))
    _, fixed, variables = cert.consume(word, 0, template(None))
    if len(variables) != 2 or len(word) != 8:
        raise ArithmeticError('the local modular word has unaccounted controls')
    targets = sorted((-time*lo/2, -time*hi/2))
    variable_error = sum(max(abs(v-targets[0]), abs(v-targets[1])) for v in variables)
    ports = tree.port_bound(Counter((event['kind'], event['angle']) for event in word))
    bound = fixed+variable_error+ports
    if bound >= F(1, 10**9):
        raise ArithmeticError('the modular source word exceeds its declared error budget')
    return word, {'q': str(q), 'modular_time': str(time), 'series_order': terms,
                  'log_odds_lower': str(lo), 'log_odds_upper': str(hi),
                  'native_pulses': len(word), 'word_sha256': cert.word_digest(word),
                  'fixed_pi_angle_error': str(fixed),
                  'logarithm_and_literal_angle_error': str(variable_error),
                  'native_port_error': str(ports),
                  'complete_isometry_error_upper': str(bound),
                  'complete_unhalved_trace_error_upper': str(2*bound)}


def native_execution(q=F(1, 4), time=F(1, 2)):
    word, bounds = word_and_bound(q, time)
    initial = np.eye(4, dtype=complex)/2
    target = np.diag(np.exp(-1j*float(time)*math.log(float((1-q)/q))*np.array([1, 0, 0, 1])))/2
    run = execute.complete(word, 2, initial)
    distance = float(np.sqrt(np.linalg.norm(run.success-target)**2+run.failure_probability))
    if distance > float(F(bounds['complete_isometry_error_upper']))+2e-11:
        raise ArithmeticError('complete native modular motion misses the full-site target')
    if run.failure_weights.shape != (8, 4) or run.normalization_residual > 2e-11:
        raise ArithmeticError('complete native modular outputs were lost')
    return {**bounds, 'full_site_rows': 4, 'retained_reference_columns': 4,
            'native_first_failure_slots': run.native_failure_slots,
            'complete_first_failure_probability': run.failure_probability,
            'normalization_residual': run.normalization_residual,
            'complete_output_error_diagnostic': distance,
            'all_failure_weights_retained': True,
            'success_postselected_or_renormalized': False}


def refined_commands(q=F(1, 4), time=F(1, 2), bits=80, terms=48):
    time = rational(time)
    if abs(time) > F(1, 2) or type(bits) is not int or bits < 1:
        raise ValueError('bounded modular time and positive integer port precision required')
    lo, hi = log_odds_interval(q, terms)
    angle = -time*(lo+hi)/4
    commands = [precision.dyadic_command(kind, desired, bits)
                for kind, _, desired in template(angle)]
    norm = sum(precision.generator_norms())
    error = 8*norm/F(2**bits)+abs(time)*(hi-lo)/2
    return {'q': str(q), 'modular_time': str(time), 'fraction_bits': bits,
            'series_order': terms,
            'exact_dyadic_port_commands': commands,
            'source_generator_norm_sum': str(norm),
            'complete_isometry_error_upper': str(error),
            'complete_unhalved_trace_error_upper': str(2*error),
            'legacy_binary64_API_used_for_refined_commands': False}


def certificate():
    examples = [(F(1, 4096), F(1, 2)), (F(1, 4), F(1, 2)),
                (F(1, 2), F(1, 2)), (F(3, 4), F(-1, 2)),
                (F(4095, 4096), F(1, 2))]
    return {'modular_unitary': 'exp(-i*s*log((1-q)/q)*D)',
            'clock': 'regional modular parameter s; local contact time -2*log((1-q)/q)*s/U',
            'calibrated_q_interval': ['1/4096', '4095/4096'],
            'modular_time_interval': ['-1/2', '1/2'],
            'actual_ground_q_lies_strictly_inside_calibration_domain': True,
            'actual_ground_exact_q_control_executed': False,
            'all_four_occupation_columns_and_remote_reference_extension': True,
            'native_examples': [native_execution(q, time) for q, time in examples],
            'refined_example': refined_commands(),
            'full_Hamiltonian_modular_flow_identified': False,
            'physical_precision_and_drift_control_are_inputs': True}
