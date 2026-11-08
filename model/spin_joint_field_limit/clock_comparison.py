"""Exact selection rotation and complete clocks retaining the source field."""
from fractions import Fraction as F
from functools import lru_cache
import math

import numpy as np

from model.certified_reference import certificate as cert
from model.native_ordered_ground import source
from model.native_spin_observer import circuits, precision
from model.spin_twist_closure.field_gap import _odd_sites


def _positive(value, name):
    if type(value) not in (int, F) or value <= 0:
        raise ValueError(name + ' must be a positive rational')
    return F(value)


def selection_bounds(index, scaled_time, q_norm, p_norm, commutator):
    """Rational upper bounds in (3); supplied norms bound the actual input."""
    if type(index) is not int or index < 1:
        raise ValueError('positive integer family index required')
    values = (scaled_time, q_norm, p_norm, commutator)
    if any(type(x) not in (int, F) for x in values):
        raise ValueError('exact rational bounds required')
    if min(q_norm, p_norm, commutator) < 0:
        raise ValueError('field norms and commutator bounds must be nonnegative')
    t, q, r, c = map(F, values)
    p = F(1, index**2)
    angle = abs(t)/index**3
    cosine_error = min(F(2), angle**2/2)
    sine = min(F(1), angle)
    return {'scale': p, 'field_over_coupling': F(1, index**5),
            'selection_angle_absolute_upper': angle,
            'cosine_difference_upper': cosine_error, 'sine_absolute_upper': sine,
            'Q_gate_input_vector_error_upper': cosine_error*q+p*sine*(r+c/2),
            'P_gate_input_vector_error_upper': cosine_error*r+sine/p*(q+c/2),
            'mixed_time_history_comparison_claimed': False,
            'input_norm_bounds_are_required': True}


def clock_plan(sites, edges, coupling, field, scale, scaled_durations, error):
    """Count-only finite source clock; input preparation and instruments separate."""
    if type(sites) is not int or sites < 1:
        raise ValueError('positive integer source site count required')
    if type(edges) is not int or edges < 0 or edges > 3*sites:
        raise ValueError('cubic source edge count outside the declared bound')
    j = _positive(coupling, 'coupling')
    b = _positive(field, 'field')
    p = _positive(scale, 'scale')
    eps = _positive(error, 'error')
    if eps >= 1:
        raise ValueError('error must be less than one')
    durations = tuple(scaled_durations)
    if any(type(t) not in (int, F) for t in durations):
        raise ValueError('exact rational scaled durations required')
    physical = tuple(abs(F(t))/(j*p) for t in durations if t)
    norm_sum = j*edges+b*sites/2
    squared_time_sum = sum((t*t for t in physical), F())
    slices = math.ceil(norm_sum**2*squared_time_sum/eps) if physical else 0
    pulses = len(physical)*slices*(62*edges+11*sites)
    split = norm_sum**2*squared_time_sum/(2*slices) if slices else F()
    bits = 1
    while 29*pulses*F(1, 1 << bits) > eps/2:
        bits += 1
    native = 29*pulses*F(1, 1 << bits)
    if split+native > eps:
        raise ArithmeticError('complete retained-field clock budget exceeded')
    return {'sites': sites, 'edges': edges, 'coupling': j, 'selection_field': b,
            'scale': p, 'scaled_durations': list(map(F, durations)),
            'absolute_physical_durations': list(physical),
            'term_norm_sum_upper': norm_sum, 'nonzero_clock_slots': len(physical),
            'product_slices_per_nonzero_slot': slices,
            'native_pulse_count_upper': pulses, 'dyadic_port_bits': bits,
            'product_operator_error_upper': split,
            'complete_native_isometry_error_upper': native,
            'complete_clock_half_trace_error_upper': split+native,
            'requested_error': eps, 'all_failure_records_retained': True,
            'source_clock': 'the same H_B as the prepared finite source ground',
            'input_preparation_and_instrument_errors_included': False,
            'coherent_time_control_compilation_included': False,
            'source_scalar_exchange_phase_retained': True,
            'admitted_coherent_services_required': True,
            'clock_sequence_executed': False,
            'H_B_and_H_0_histories_identified': False}


def cube_clock_plan(radius, index, scaled_durations, error, coupling=F(1)):
    if type(radius) is not int or radius < 0:
        raise ValueError('nonnegative integer box radius required')
    if type(index) is not int or index < 1:
        raise ValueError('positive integer family index required')
    j = _positive(coupling, 'coupling')
    side = 2*radius+1
    row = clock_plan(side**3, 3*(side-1)*side**2, j, j/index**5,
                     F(1, index**2), scaled_durations, error)
    return {'source_radius': radius, 'family_index': index, **row}


def clock_word(sites, edges, coupling, field, scale, scaled_duration, slices):
    """Lazy physical-frame H_B clock on an encoded native spin register."""
    source.parameters(sites, edges, coupling, field)
    p = _positive(scale, 'scale')
    if type(scaled_duration) not in (int, F):
        raise ValueError('exact rational scaled duration required')
    if type(slices) is not int or slices < 0:
        raise ValueError('nonnegative integer slice count required')
    if bool(scaled_duration) != bool(slices):
        raise ValueError('nonzero duration and slices, or both zero, required')
    odd = set(_odd_sites(sites, edges))
    if not slices:
        return
    delta = F(scaled_duration)/(coupling*p*slices)
    for _ in range(slices):
        for edge in edges:
            yield from circuits.exchange(delta, coupling, sites, edge)
        for site in range(sites):
            sign = -1 if site in odd else 1
            yield from circuits.spin_rotation('X', delta*field*sign/2, site, sites)


def dyadic_clock_word(sites, edges, coupling, field, scale, scaled_duration, slices, bits):
    for kind, targets, angle in clock_word(sites, edges, coupling, field,
                                           scale, scaled_duration, slices):
        yield {'targets': list(targets),
               'port_command': precision.dyadic_command(kind, angle, bits)}


def _spin_operators(sites):
    pauli = (np.array([[0, 1], [1, 0]], complex),
             np.array([[0, -1j], [1j, 0]], complex),
             np.diag([1., -1.]).astype(complex))
    output = []
    for site in range(sites):
        axes = []
        for axis in pauli:
            value = np.ones((1, 1), complex)
            for index in range(sites):
                value = np.kron(axis if index == site else np.eye(2), value)
            axes.append(value/2)
        output.append(axes)
    return output


@lru_cache(maxsize=1)
def source_rotation_certificate():
    """Binary-rational exact finite identities against the actual source matrix."""
    rows = []
    j, b, p = F(1), F(1, 4), F(1, 4)
    for sites, edges in ((1, ()), (2, ((0, 1),)), (3, ((0, 1), (1, 2)))):
        ops = _spin_operators(sites)
        dim = 1 << sites
        odd = set(_odd_sites(sites, edges))
        signs = [1-2*(i in odd) for i in range(sites)]
        h0 = sum((sum((ops[x][a]@ops[y][a] for a in range(3)),
                      np.zeros((dim, dim), complex))-np.eye(dim)/4
                  for x, y in edges), np.zeros((dim, dim), complex))
        sel = -float(b)*sum((signs[x]*ops[x][0] for x in range(sites)),
                            np.zeros((dim, dim), complex))
        marshall = np.diag([(-1)**sum((mask >> site) & 1 for site in odd)
                           for mask in range(dim)])
        integer, scale = source.matrix(sites, edges, j, b)
        if not np.array_equal(scale*marshall@(h0+sel)@marshall, integer):
            raise ArithmeticError('physical frame differs from actual source')
        profile = [F(x+1, 2) for x in range(sites)]
        q = float(p*p)*sum((float(profile[x])*signs[x]*ops[x][2]
                            for x in range(sites)), np.zeros((dim, dim), complex))
        momentum = -float(p)*sum((float(profile[x])*ops[x][1]
                                  for x in range(sites)), np.zeros((dim, dim), complex))
        generator = sel/float(j*p)
        qdot = 1j*(generator@q-q@generator)
        pdot = 1j*(generator@momentum-momentum@generator)
        if not np.array_equal(qdot, float(b/j)*momentum):
            raise ArithmeticError('selection Q drift sign or normalization differs')
        if not np.array_equal(pdot, -float(b/(j*p*p))*q):
            raise ArithmeticError('selection P drift sign or normalization differs')
        rows.append({'sites': sites, 'all_rows': dim, 'actual_source_frame_exact': True,
                     'selection_Q_drift_exact': True, 'selection_P_drift_exact': True})
    return {'finite_source_checks': rows,
            'exact_selection_angle': 'B T/(J p)',
            'joint_family_angle': 'T/n^3',
            'joint_family_P_to_Q_absolute_coefficient_upper': '|T|/n',
            'joint_family_Q_to_P_absolute_coefficient_upper': '|T|/n^5',
            'mixed_history_Duhamel_estimate_proved': False,
            'equal_time_drift_used_to_identify_two_clocks': False}


@lru_cache(maxsize=1)
def small_native_clock():
    """A complete reference-entangled clock execution, including every failure."""
    from collections import Counter
    from scipy.linalg import expm
    from model.centered_native_preparation import tree
    from model.native_spin_observer import algebra, execute

    sites, edges = 2, ((0, 1),)
    j, b, p, time, slices = F(1), F(1, 4), F(1, 2), F(1, 16), 2
    template = list(clock_word(sites, edges, j, b, p, time, slices))
    word = list(circuits.emit(template))
    consumed, angle_error, variables = cert.consume(word, 0, template)
    if consumed != len(word) or variables:
        raise ArithmeticError('source clock contains unresolved commands')
    port_error = tree.port_bound(Counter((event['kind'], event['angle']) for event in word))
    native_error = angle_error+port_error
    ops = _spin_operators(sites)
    h0 = sum((ops[0][a]@ops[1][a] for a in range(3)), np.zeros((4, 4), complex))-np.eye(4)/4
    sel = -float(b)*(ops[0][0]-ops[1][0])
    delta = float(time/(j*p*slices))
    unitary = np.eye(4, dtype=complex)
    for _ in range(slices):
        unitary = expm(-1j*delta*float(j)*h0)@unitary
        unitary = expm(-1j*delta*sel)@unitary
    inclusion = algebra.inclusion(sites)
    initial = inclusion/2  # Reference dimension four; no work row is omitted.
    result = execute.complete(word, 4, initial)
    target = inclusion@unitary/2
    distance = float(np.sqrt(np.linalg.norm(result.success-target)**2+
                             result.failure_probability))
    if distance > float(native_error)+execute.ADMISSION:
        raise ArithmeticError('complete native clock differs from the literal source product')
    if max(distance, result.failure_probability, abs(result.normalization_residual)) >= 1e-10:
        raise ArithmeticError('native clock witness exceeds fixed numerical enclosure')
    exact = inclusion@expm(-1j*float(time/(j*p))*(float(j)*h0+sel))/2
    exact_distance = float(np.sqrt(np.linalg.norm(result.success-exact)**2+
                                   result.failure_probability))
    split = (j+b*sites/2)**2*(time/(j*p))**2/(2*slices)
    if exact_distance > float(split+native_error)+execute.ADMISSION:
        raise ArithmeticError('actual retained-field clock exceeds product/native bound')
    return {'sites': sites, 'all_work_rows': 16, 'reference_dimension': 4,
            'slices': slices, 'scaled_duration': str(time), 'scale': str(p),
            'physical_duration': str(time/(j*p)), 'pulses': len(word),
            'word_sha256': cert.word_digest(word), 'failure_slots': result.native_failure_slots,
            'native_complete_error_upper': str(native_error),
            'product_operator_error_upper': str(split),
            'literal_product_distance_interval': ['0', '1/10000000000'],
            'failure_probability_interval': ['0', '1/10000000000'],
            'normalization_residual_absolute_upper': '1/10000000000',
            'exact_H_B_clock_distance_within_sum_bound': True,
            'all_failure_records_retained': True,
            'large_joint_family_schedule_executed': False}


def build():
    return {'selection_rotation': source_rotation_certificate(),
            'equal_time_unit_norm_examples': [
                selection_bounds(n, F(1), F(1), F(1), F(1)) for n in (2, 4, 8)],
            'retained_field_clock_examples': [
                cube_clock_plan(1, n, (F(1), F(-1, 2), F()), F(1, 100))
                for n in (2, 4)],
            'complete_native_execution': small_native_clock(),
            'prior_ordered_limit_identified': False,
            'H_B_and_H_0_multitime_equivalence_proved': False}
