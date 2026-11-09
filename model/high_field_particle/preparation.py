"""A short complete native preparation using a finite residual certificate."""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt

import numpy as np

from model.certified_reference import certificate as cert
from model.centered_native_preparation import tree
from model.native_ordered_ground import source
from model.native_spin_observer import circuits, execute
from model.spin_twist_closure import field_gap


def product_bound(sites, edges, coupling, field):
    source.parameters(sites, edges, coupling, field)
    counts = Counter(edges)
    squared = sum(weight*weight for weight in counts.values())
    separation = field-coupling*len(edges)/2
    if separation <= 0:
        raise ValueError('positive product-Rayleigh to excited-spectrum separation required')
    ceil_root = isqrt(squared)
    if ceil_root*ceil_root != squared:
        ceil_root += 1
    residual_squared = coupling*coupling*squared/4
    excited_weight = min(F(1), residual_squared/separation**2)
    return {'sites': sites, 'edges': len(edges),
            'ground_Hilbert_space': 'full tensor product of two-state spins, encoded in the protected single-occupancy code',
            'squared_edge_multiplicity_sum': squared,
            'Rayleigh_value': -field*sites/2-coupling*len(edges)/2,
            'excited_spectrum_minus_Rayleigh_lower': separation,
            'residual_norm_squared': residual_squared,
            'excited_weight_upper': excited_weight,
            'ideal_local_aligned_pointer_success_lower': 1-excited_weight,
            'product_ground_half_trace_upper': min(F(1), coupling*ceil_root/(2*separation)),
            'full_finite_ground_simple': True,
            'uses_a_fitted_eigenvector': False}


def _encoded_address(mask, sites):
    # Owning inclusion: logical q=0 is down (10), q=1 is up (01).
    return sum((1 if mask >> x & 1 else 2) << (2*x) for x in range(sites))


@lru_cache(maxsize=1)
def _execution_certificate():
    sites, edges, j, b = 8, source.cube()[0], F(1), F(400)
    bound = product_bound(sites, edges, j, b)
    template = list(field_gap.preparation_word(sites, edges, j, b, F(0), 0))
    preparation_pulses = len(template)
    pointer = 2*sites
    template += circuits.probe('X', sites, 0, pointer)
    # An exact controlled logical Z on pointer outcome zero. The source J(-pi)
    # equals -CZ; I(pi) removes that scalar before the two pointer X gates.
    template += circuits.xgate(pointer)
    template += [('J', (0, pointer), cert._pi(-1)),
                 ('I', (0,), cert._pi(1))]
    template += circuits.xgate(pointer)
    word = list(circuits.emit(template))
    consumed, angle_error, variables = cert.consume(word, 0, template)
    if consumed != len(word) or variables:
        raise ArithmeticError('short source word has an unresolved angle')
    port_error = tree.port_bound(Counter((event['kind'], event['angle']) for event in word))
    native_error = angle_error+port_error
    if preparation_pulses != 136 or len(word) != 167:
        raise ArithmeticError('the literal source preparation/readout word changed')
    rows = 1 << (pointer+1)
    initial, target = np.zeros((rows, 1), complex), np.zeros((rows, 1), complex)
    initial[0] = 1
    odd = sum(1 << x for x in source.cube()[1])
    for mask in range(1 << sites):
        address = _encoded_address(mask, sites)
        # The X pointer is zero on site0's +X product; feedback Z flips its sign.
        target[address] = (-1)**((mask & odd).bit_count()+((mask >> 0) & 1))/16
    actual = execute.complete(word, pointer+1, initial)
    distance = float(np.sqrt(np.linalg.norm(actual.success-target)**2+actual.failure_probability))
    if distance > float(native_error)+execute.ADMISSION:
        raise ArithmeticError('complete source output misses the literal product/creation target')
    if max(distance, actual.failure_probability, abs(actual.normalization_residual)) >= 1e-10:
        raise ArithmeticError('short source execution exceeds its diagnostic admission')
    if native_error >= F(1, 10**13) or actual.native_failure_slots != 608:
        raise ArithmeticError('short source execution changed its admitted error or record count')
    total = bound['product_ground_half_trace_upper']+native_error
    if total >= F(1, 100):
        raise ArithmeticError('finite ground-to-created-output error exceeds one percent')
    return {'finite_source': {'coupling': str(j), 'selection_field': str(b)},
            'finite_product_ground_bound': bound,
            'ideal_exact_ground_success_probability_lower':
                bound['ideal_local_aligned_pointer_success_lower'],
            'ideal_exact_ground_success_is_particle_band_weight': False,
            'actual_native_preparation_pulses': preparation_pulses,
            'actual_native_pointer_and_feedback_pulses': len(word)-preparation_pulses,
            'actual_native_total_pulses': len(word),
            'complete_work_bits': pointer+1, 'complete_work_rows': rows,
            'native_word_sha256': actual.word_sha256,
            'native_first_failure_slots': actual.native_failure_slots,
            'complete_angle_error_bound': angle_error,
            'complete_port_error_bound': port_error,
            'complete_native_isometry_error_bound': native_error,
            'full_output_half_trace_error_from_exact_ground_instrument_upper': total,
            'complete_execution_distance_diagnostic_upper': '1/10000000000',
            'execution_diagnostic_is_a_rigorous_interval_certificate': False,
            'all_native_first_failures_and_both_pointer_outcomes_retained': True,
            'success_renormalized': False,
            'source_program_actually_executed': True,
            'literal_pointer_matches_charge_instrument_on_protected_code': True,
            'literal_pointer_equals_abstract_out_of_code_extension': False,
            'infinite_volume_Yarotsky_threshold_certified_for_this_coupling': False,
            'physical_particle_assigned': False}


def _detached(value):
    if type(value) is F:
        return F(value.numerator, value.denominator)
    if type(value) is dict:
        return {key: _detached(item) for key, item in value.items()}
    if type(value) is list:
        return [_detached(item) for item in value]
    if type(value) is tuple:
        return tuple(_detached(item) for item in value)
    return value


def execution_certificate():
    """Return detached exact evidence from the complete native execution."""
    return _detached(_execution_certificate())
