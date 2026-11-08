"""Complete coherent charge classes from the source's native phase words."""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache

import numpy as np

from model.protected_defect_catalog import charge, field
from model.native_spin_observer import algebra, circuits, execute, precision
from model.certified_reference import certificate as cert
from model.centered_native_preparation import tree
from model.centered_particle_spectrum.certificate import decimal_bounds


def _parameters(signs, modulus):
    if (not isinstance(signs, tuple) or not signs or any(type(x) is not int or x not in (-1, 1) for x in signs)
            or sum(signs) != 0 or type(modulus) is not int or modulus < 4 or modulus & (modulus-1)):
        raise ValueError('balanced source spin signs and a power-of-two modulus at least four required')
    return len(signs), modulus.bit_length()-1


def hadamard(bit):
    return [('I', (bit,), cert._pi(-1, 2)), ('H', (bit,), cert._pi(1, 2))]


def controlled_phase(first, second, multiple):
    # The native J is diag(exp(i theta), exp(i theta), exp(i theta), 1).
    # Its scalar phase is retained, including in controlled computations.
    return [('I', (first,), cert.PiMultiple(multiple)),
            ('J', (first, second), cert.PiMultiple(-multiple))]


def fourier(wires):
    """Positive discrete Fourier transform, little-endian record ordering."""
    word = []
    for j in reversed(range(len(wires))):
        word += hadamard(wires[j])
        for k in reversed(range(j)):
            word += controlled_phase(wires[j], wires[k], F(1, 1 << (j-k)))
    for j in range(len(wires)//2):
        a, b = wires[j], wires[-1-j]
        word += cert.cnot(a, b)+cert.cnot(b, a)+cert.cnot(a, b)
    return word


def template(signs, modulus):
    """Literal source program for the charge-modulo-N nondemolition pointer."""
    sites, bits = _parameters(signs, modulus)
    record = tuple(range(2*sites, 2*sites+bits))
    word = cert.inverse(circuits.encode(sites))
    for x in range(sites):
        word += hadamard(2*x)
    for wire in record:
        word += hadamard(wire)
    for j, wire in enumerate(record):
        for x, sign in enumerate(signs):
            word += controlled_phase(wire, 2*x, -sign*F(2*(1 << j), modulus))
    word += cert.inverse(fourier(record))
    for x in range(sites):
        word += hadamard(2*x)
    word += circuits.encode(sites)
    expected = 14*sites+4*bits+2*sites*bits+bits*(bits-1)+9*(bits//2)
    if len(word) != expected:
        raise ArithmeticError('native charge Fourier resource count differs from literal word')
    return word


def source_projectors(modulus):
    """Exact numerator matrices, all record classes, denominator 256."""
    op = charge.operators()
    _parameters(op['signs'], modulus)
    if modulus > 256:
        raise ValueError('the complete matrix evaluator supports at most 256 explicit records')
    w = op['walsh']
    identity = np.eye(256, dtype=np.int64)
    result = []
    for residue in range(modulus):
        columns = w[:, [i for i, q in enumerate(op['charges']) if q % modulus == residue]]
        p = columns@columns.T
        if not np.array_equal(p@p, 256*p) or not np.array_equal(p@op['H'], op['H']@p):
            raise ArithmeticError('source record is not a conserved orthogonal charge-class projector')
        result.append(p)
    if not np.array_equal(sum(result, np.zeros_like(identity)), 256*identity):
        raise ArithmeticError('charge-class records are incomplete')
    for r, p in enumerate(result):
        for q in range(r):
            if np.any(p@result[q]):
                raise ArithmeticError('distinct retained charge records overlap')
    return op, tuple(result)


def alias_bound(amplitude, modulus, original_probability_floor=None):
    """Absolute original-weight alias bound; optional conditioning is explicit."""
    if type(amplitude) not in (int, F) or amplitude == 0:
        raise ValueError('nonzero rational sine amplitude required')
    _parameters((1, -1), modulus)
    second = 6*F(amplitude)**2
    absolute = min(F(1), second/(modulus-1)**2)
    result = {'sine_amplitude': str(amplitude), 'modulus': modulus,
        'unit_annulus_configuration_variance_upper': '6',
        'unnormalized_charge_second_moment_upper': str(second),
        'absolute_original_weight_wrong_charge_probability_upper': str(absolute),
        'same_bound_after_any_commuting_contractive_energy_filter': True,
        'both_target_charge_classes': [1, modulus-1],
        'all_modulus_records_retained': True}
    if original_probability_floor is not None:
        if type(original_probability_floor) not in (int, F) or not 0 < original_probability_floor <= 1:
            raise ValueError('positive certified original preparation probability required')
        result['original_preparation_probability_lower'] = str(original_probability_floor)
        result['conditional_wrong_charge_probability_upper'] = str(min(F(1), absolute/original_probability_floor))
    return result


def prelimit_charge_bound(amplitude, scale):
    """Rational bound for the same bandlimited source test, at fixed p."""
    if any(type(x) not in (int, F) for x in (amplitude, scale)) or amplitude == 0 or not 0 < scale <= 1:
        raise ValueError('nonzero rational amplitude and scale in (0,1] required')
    a, p = abs(F(amplitude)), F(scale)
    # sqrt(6) < 5/2; this bound is deliberately not called the exact moment.
    result = a*a*(F(5, 2)+a*p/4)**2
    return {'scale': str(p), 'charge_second_moment_upper': str(result),
        'small_scale_threshold': str(min(F(1), 1/(2*a))),
        'below_seven_amplitude_squared': result < 7*a*a if p <= min(F(1), 1/(2*a)) else None}


def native_error(template_value):
    word = tuple(circuits.emit(template_value))
    _, angle, variables = cert.consume(word, 0, template_value)
    if variables:
        raise ArithmeticError('every native charge angle must be specified')
    port = tree.port_bound(Counter((e['kind'], e['angle']) for e in word))
    return word, angle+port


def precision_plan(signs, modulus, error=F(1, 10**6)):
    word = template(signs, modulus)
    if type(error) is not F or not 0 < error < 1:
        raise ValueError('complete trace error must be an exact number in (0,1)')
    bits = precision.bits_for(len(word), error)
    return {'source_pulses': len(word), 'record_bits': modulus.bit_length()-1,
        'full_occupation_and_record_bits': 2*len(signs)+modulus.bit_length()-1,
        'dyadic_port_fraction_bits': bits, 'complete_trace_error_upper': str(F(58*len(word), 1 << bits)),
        'arbitrary_untouched_reference': True, 'coherent_load_export_service_required': True,
        'fixed_binary64_port_interface_claimed_arbitrarily_precise': False}


@lru_cache(maxsize=1)
def native_two_site_execution():
    """All four logical input columns remain a coherent quantum reference."""
    signs, modulus = (1, -1), 4
    word, error = native_error(template(signs, modulus))
    inclusion = algebra.inclusion(2)
    initial = np.zeros((64, 4), complex)
    initial[:16] = inclusion/2
    run = execute.complete(word, 6, initial)
    walsh = np.array([[(-1)**((a&b).bit_count()) for b in range(4)] for a in range(4)])
    charges = [((1-2*(b&1))-(1-2*((b>>1)&1)))//2 for b in range(4)]
    target = np.zeros_like(initial)
    rows = []
    for r in range(modulus):
        columns = walsh[:, [i for i, q in enumerate(charges) if q % modulus == r]]
        p = columns@columns.T/4
        target[16*r:16*(r+1)] = inclusion@p/2
        probability = float(np.vdot(run.success[16*r:16*(r+1)], run.success[16*r:16*(r+1)]).real)
        ideal_probability = F(int(np.trace(columns@columns.T)), 16)
        if abs(probability-float(ideal_probability)) > 2*float(error)+execute.ADMISSION:
            raise ArithmeticError('a complete native charge record disagrees with its ideal value')
        rows.append({'record': r, 'ideal_Choi_probability': str(ideal_probability),
            'literal_native_probability_check_passed': True})
    distance = float(np.sqrt(np.linalg.norm(run.success-target)**2+run.failure_probability))
    if distance > float(error)+execute.ADMISSION:
        raise ArithmeticError('literal native charge instrument disagrees with complete coherent target')
    return {'sites': 2, 'modulus': 4, 'coherent_reference_dimension': 4,
        'complete_source_work_rows': 64, 'native_pulses': len(word),
        'word_sha256': cert.word_digest(word), 'complete_isometry_error_upper': str(error),
        'complete_trace_error_upper': str(2*error),
        'literal_native_complete_output_check_passed': True,
        'numerical_execution_comparison_admission': str(F(execute.ADMISSION)),
        'native_first_failure_slots': run.native_failure_slots,
        'all_first_failure_weights_evaluated_and_retained': True,
        'all_physical_pointer_outcomes': 36, 'reserved_pointer_outcomes': 32,
        'complete_class_records': rows,
        'reference_columns_dephased': False, 'successful_output_renormalized': False}


def build():
    field.configuration_charge()
    op, projectors = source_projectors(16)
    ground, error, checked = charge.ground_input()
    created = op['local_Z'].astype(object)@ground
    norm = int(ground@ground)
    records = []
    for r, p in enumerate(projectors):
        weight = F(int(created@(p.astype(object)@created)), 256*norm)
        records.append({'record': r, 'represented_charges': sorted({q for q in op['charges'] if q % 16 == r}),
            'projector_rank': int(np.trace(p))//256,
            'original_Z0_ground_probability_interval': decimal_bounds(max(F(), weight-2*error), min(F(1), weight+2*error), 12)})
    return {'same_source_charge_and_clock_replayed': True,
        'complete_8_site_modulo_16_records': records,
        'eight_site_native_program': precision_plan(op['signs'], 16),
        'literal_two_site_reference_execution': native_two_site_execution(),
        'continuum_absolute_alias_example': alias_bound(F(1), 1 << 16),
        'same_source_prelimit_charge_tail_example': prelimit_charge_bound(F(2000), F(1, 4000)),
        'uniform_odd_branch_conditional_alias_example': alias_bound(F(2000), 1 << 40, F(25, 2304*1000**2)),
        'scope': {'finite_charge_class_projectors_and_H_commutation_exact': True,
            'coherence_within_each_record_class_preserved': True,
            'all_8_site_native_pulse_tape_executed_here': False,
            'source_U1_configuration_subnet_consumed': True,
            'continuum_charge_bound_uses_actual_Stone_field_moments': True,
            'finite_p_scale_threshold_for_ground_preparation_supplied': False,
            'physical_electric_charge_or_particle_species_identified': False}}
