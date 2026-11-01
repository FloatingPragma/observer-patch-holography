"""One identical literal source word per mode gives an exact number reference.

The exact mathematical response is exp(sum_p c_p D_p), where each c_p is
the represented binary64 control interpreted as a rational and D_p is the
retained Q(sqrt(5)) source. Numerical expm is only a diagnostic below.
Complete dephasing export yields probabilities a,b,f_2,...,f_5. Thus the
unconditioned proper-field block on M modes is diag(a,b) tensor M, without
success normalization or an assumption that the represented pulse is ideal.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
import math

import numpy as np

from ..certified_reference import interval, source
from ..observer.source_scalar import compiler


def modes(value):
    if type(value) is not int or not 1 <= value <= 1_000_000:
        raise ValueError('a positive finite mode count up to one million is required')
    return value


def source_word():
    return [compiler.pulse('Y', -math.pi/3, [0])]


def logarithm_three():
    """Exact positive atanh series, with its entire omitted tail bounded."""
    z = F(1, 2)
    lower = 2*sum((z**(2*k+1)/(2*k+1) for k in range(96)), F())
    tail = 2*z**193/(193*(1-z*z))
    return lower, lower+tail


@lru_cache(maxsize=1)
def one_mode_certificate():
    event, = source_word()
    control = source.literal_controls(event['kind'], event['angle'])
    # These are independent bounds on represented angle and full-source port
    # rounding. No numerical matrix exponential supplies their acceptance.
    angle_error = interval.pi_angle_error(F(event['angle']), F(-1, 3))
    port_error = source.generator_error('Y', F(event['angle']), control)
    error = angle_error+port_error
    delta = (2+error)*error
    if delta >= F(1, 1000):
        raise ArithmeticError('source coefficient bounds lost faithful nontracial weights')
    source.complete_export_certificate()
    low, high = logarithm_three()
    gamma_error = delta/(F(1, 4)-delta)+delta/(F(3, 4)-delta)
    return {
        'word': source_word(),
        'literal_controls': [str(x) for x in control],
        'literal_controls_sha256': hashlib.sha256(json.dumps([str(x) for x in control]).encode()).hexdigest(),
        'angle_error': angle_error, 'full_source_port_error': port_error,
        'complete_one_mode_isometry_error': error,
        'each_valid_probability_error': delta,
        'a_lower': F(1, 4)-delta, 'a_upper': F(1, 4)+delta,
        'b_lower': F(3, 4)-delta, 'b_upper': F(3, 4)+delta,
        'one_mode_abort_probability_upper': error**2,
        'gamma_minus_log3_absolute_bound': gamma_error,
        'gamma_lower': low-gamma_error, 'gamma_upper': high+gamma_error,
        'gamma_definition': 'log(b/a), a=|[exp(sum c_p D_p)]00|^2, b=|[exp(sum c_p D_p)]10|^2',
        'exact_valid_weights_need_not_sum_to_one': True,
        'all_six_dephasing_export_records_retained': True,
    }


def resource_certificate(count, region_modes=1):
    count, region_modes = modes(count), modes(region_modes)
    if region_modes > count:
        raise ValueError('the region must be contained in the mode register')
    if region_modes > 128:
        raise ValueError('the explicit local rational calculator admits at most128 modes')
    local = one_mode_certificate()
    error, delta = local['complete_one_mode_isometry_error'], local['each_valid_probability_error']
    return {
        'modes': count, 'source_pulses': count, 'complete_exports': count,
        'records_per_mode': 6, 'total_record_count_formula': f'6^{count}',
        'proper_field_dimension_formula': f'2^{count}',
        'successful_field_trace': '(a+b)^M',
        'abort_mass_upper': min(F(1), count*error**2),
        'complete_output_trace_norm_error_to_ideal_number_reference': min(F(2), 2*count*error),
        'this_error_does_not_grow_with_ideal_record_independent_history_length': True,
        'region_modes': region_modes,
        'initial_local_complete_output_trace_error': min(F(2), 2*region_modes*error),
        'initial_local_proper_block_minimum_weight': local['a_lower']**region_modes,
        'initial_local_modular_operator_error': region_modes*delta/(F(1, 4)-delta),
        'local_bound_uses_complete_output_marginal_then_local_code_projection': True,
        'global_success_projection_would_add_an_outside_region_scalar': True,
        'local_bound_uniform_in_record_independent_ideal_histories_with_abort_bypass': True,
    }


@lru_cache(maxsize=1)
def numerical_exports():
    """Execute the existing M6 pulse and every rectangular readout Kraus map."""
    event, = source_word()
    initial = np.eye(6, dtype=complex)[:, 0]
    evolved = compiler.native_pulse(event['kind'], event['angle'])@initial
    outcomes = tuple(k@evolved for k in compiler.transfer_kraus(2, coherent=False))
    for record, vector in enumerate(outcomes):
        expected = 0 if record == 0 else 1 if record == 1 else 2
        if np.any(np.delete(vector, expected)):
            raise ArithmeticError('the literal readout is not diagonal with explicit abort flags')
        vector.flags.writeable = False
    return outcomes


def numerical_branch(records):
    """Access any product record, including every local abort, without sampling."""
    if (type(records) not in (tuple, list) or not records or len(records) > 128
            or any(type(j) is not int or j not in range(6) for j in records)):
        raise ValueError('one to128 complete local record indices0..5 required')
    outputs = numerical_exports()
    address, amplitude, place = 0, 1.+0j, 1
    for record in records:
        destination = min(record, 2)
        address += destination*place
        amplitude *= outputs[record][destination]
        place *= 3
    return address, amplitude


def numerical_diagnostic():
    outcomes = numerical_exports()
    probabilities = [float(np.vdot(v, v).real) for v in outcomes]
    total = math.fsum(probabilities)
    if abs(total-1) > 1e-12:
        raise ArithmeticError('numerical source export diagnostic lost a branch')
    # A diagnostic cannot prove the exact source bounds: it is displayed and
    # independently tested, rather than used as an assumed precision premise.
    return {'probabilities_including_four_abort_records': probabilities,
            'total_probability': total,
            'valid_field_weight_sum': probabilities[0]+probabilities[1],
            'gamma_display': math.log(probabilities[1]/probabilities[0]),
            'numerical_matrix_exponential_is_not_the_certifying_arithmetic': True}


def modular_word(count, parameter=1.):
    count = modes(count)
    if count > 4096:
        raise ValueError('the explicit source tape admits at most4096 modes')
    if (type(parameter) not in (int, float) or not math.isfinite(parameter)
            or abs(parameter) > 2):
        raise ValueError('a finite real modular parameter with absolute value at most two is required')
    angle = float(parameter)*math.log(3)
    return [event for bit in range(count) for event in
            (compiler.pulse('Z', -angle/2, [bit]), compiler.pulse('I', angle/2, [bit]))]


def modular_certificate(count, parameter=1.):
    word = modular_word(count, parameter)
    local = one_mode_certificate()
    lower, upper = logarithm_three()
    t = F(float(parameter))
    first = word[:2]
    desired = sorted((t*lower, t*upper))
    actual = F(first[1]['angle'])*2
    scalar_error = max(abs(actual-x) for x in desired)+abs(t)*local['gamma_minus_log3_absolute_bound']
    port_error = sum((source.pulse_error(event) for event in first), F())
    return {'modes': count, 'modular_parameter': t, 'source_pulses': len(word),
            'centered_modular_generator': 'gamma*N, with vacuum phase fixed to zero',
            'global_reference_normalization_phase_omitted_only_from_the_automorphism': True,
            'per_mode_scalar_and_reference_error': scalar_error,
            'per_mode_full_source_control_error': port_error,
            'complete_all_input_isometry_error': count*(scalar_error+port_error),
            'all_clock_export_failures_retained': True,
            'clock_word_sha256': hashlib.sha256(json.dumps(word, sort_keys=True).encode()).hexdigest()}
