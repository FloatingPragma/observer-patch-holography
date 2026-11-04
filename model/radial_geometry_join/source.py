"""Complete native parity/energy instrument on the actual field occupations."""
from collections import Counter
from fractions import Fraction as F
import math

import numpy as np

from model.certified_reference import certificate as cert, interval as iv
from model.centered_native_preparation import tree
from model.common_source import scalar_coupling
from model.matter_geometry_response import source as parent
from model.native_spin_observer import circuits, execute
from model.observer.source_scalar import compiler


def parity_copy(up=0, down=1, probe=2):
    if any(type(i) is not int or i < 0 for i in (up, down, probe)) or len({up, down, probe}) != 3:
        raise ValueError('three distinct nonnegative occupation/probe wires required')
    return circuits.xgate(probe)+cert.cnot(up, probe)+cert.cnot(down, probe)


def interaction_template():
    """Independent exact routing of the retained dimension-two subword."""
    pi = cert.PiMultiple
    cnot = cert.cnot(0, 2)
    return [('X', (2,), pi(F(-1, 4))), ('Z', (2,), pi(F(-1, 4))),
            ('X', (2,), pi(F(1, 4))), ('Z', (2,), F(0)),
            ('I', (2,), F(0)), ('I', (0,), F(0)), ('Z', (0,), F(0)),
            ('Z', (2,), F(0))]+cnot+[
                ('Z', (2,), None)]+cert.inverse(cnot)+[
            ('I', (2,), F(0)), ('Z', (2,), F(0)),
            ('X', (2,), pi(F(-1, 4))), ('Z', (2,), pi(F(1, 4))),
            ('X', (2,), pi(F(1, 4)))]


def word_and_bound():
    """Bound literal controls and the nonrational angle independently."""
    interaction = scalar_coupling.scalar_subwords(
        2, 4*parent.ANGLE, coupling=F(1, 4))['interaction']
    template = interaction_template()
    offset, fixed_error, variables = cert.consume(interaction, 0, template)
    if offset != len(interaction) or len(variables) != 1:
        raise ArithmeticError('the interaction has unaccounted angles or routing')
    sine, cosine = iv.sine_cosine(variables[0])
    rotation_error = (sine-F(3, 5)).absolute_upper()+(cosine-F(4, 5)).absolute_upper()
    copy = parity_copy()
    prepare = [('X', (4,), cert.PiMultiple(F(1, 4)))]
    pointer = cert.cnot(4, 6)
    before, after = copy+prepare, pointer+cert.inverse(copy)
    before_word, after_word = list(circuits.emit(before)), list(circuits.emit(after))
    _, before_error, before_variables = cert.consume(before_word, 0, before)
    _, after_error, after_variables = cert.consume(after_word, 0, after)
    if before_variables or after_variables:
        raise ArithmeticError('unaccounted copy, preparation or pointer angle')
    word = before_word+compiler.remap(interaction, [2, 3, 4, 5])+after_word
    port_error = tree.port_bound(Counter((event['kind'], event['angle']) for event in word))
    error = fixed_error+rotation_error+before_error+after_error+port_error
    coin_template = [('X', (0,), cert.PiMultiple(F(1, 4)))]+cert.cnot(0, 1)
    coin_word = list(circuits.emit(coin_template))
    _, coin_angle, coin_variables = cert.consume(coin_word, 0, coin_template)
    if coin_variables:
        raise ArithmeticError('unaccounted fair-coin angle')
    coin_error = coin_angle+tree.port_bound(Counter(
        (event['kind'], event['angle']) for event in coin_word))
    if error+9*coin_error >= F(1, 10**9):
        raise ArithmeticError('the complete fixed-control error exceeds its declared budget')
    return word, {
        'native_pulses': len(word), 'word_sha256': cert.word_digest(word),
        'interaction_trig_error': str(rotation_error),
        'fixed_angle_error': str(fixed_error+before_error+after_error),
        'native_port_error': str(port_error),
        'complete_quantum_isometry_error_upper': str(error),
        'one_complete_fair_coin_isometry_error_upper': str(coin_error),
        'nine_coin_attempt_trace_distance_upper': str(error+9*coin_error),
        'nine_coin_attempt_trace_distance_strict_bound': '1/1000000000',
        'trace_distance_convention': 'one_half_trace_norm',
        'all_native_failure_outputs_included_in_bound': True,
        'nonrational_angle_checked_by_integer_trigonometric_enclosure': True,
    }


def complete_target():
    initial = np.zeros((128, 4), complex)
    initial[:4] = np.eye(4)/2
    target = np.zeros_like(initial)
    for occupation in range(4):
        even = int(occupation.bit_count() % 2 == 0)
        target[occupation, occupation] = (7-6*even)/(2*math.sqrt(50))
        target[occupation+16+64, occupation] = 1j*(1+6*even)/(2*math.sqrt(50))
    return initial, target


def quantum_execution():
    word, error = word_and_bound()
    initial, target = complete_target()
    run = execute.complete(word, 7, initial)
    distance = float(np.sqrt(np.linalg.norm(run.success-target)**2+run.failure_probability))
    if distance > float(F(error['complete_quantum_isometry_error_upper']))+2e-11:
        raise ArithmeticError('complete native parity/energy execution fails its exact target')
    if run.normalization_residual > 2e-11:
        raise ArithmeticError('complete native parity/energy execution loses probability')
    return {**error,
        'logical_work_qubits': 7, 'full_work_rows': 128,
        'retained_quantum_reference_columns': 4,
        'central_field_occupation_columns': 4,
        'complete_logical_pointer_outcomes': 2, 'physical_pointer_outcomes': 6,
        'ideal_zero_reserved_pointer_outcomes': 4,
        'native_first_failure_slots': run.native_failure_slots,
        'complete_output_error_diagnostic': distance,
        'complete_first_failure_probability': run.failure_probability,
        'normalization_residual': run.normalization_residual,
        'success_postselected_or_renormalized': False,
        'full_27_site_Fock_evolution_executed': False,
        'local_identity_extends_with_arbitrary_remote_field_and_reference': True,
        'probe_charge_bit_uncomputed_to_its_blank': True,
        'scalar_and_pointer_bit_retained': True,
        'readout_field_Kraus': ['(7 I-6 D0)/sqrt(50)', 'i(I+6 D0)/sqrt(50)'],
        'D0': '(n_center-1)^2',
        'record1_probability': '1/50+(24/25)<D0>',
        'unconditional_parity_offdiagonal_multiplier': '7/25',
        'field_energy_injection': '-(18/25)<H_origin_incident_hopping>',
    }


def kraus(even_projector):
    if np.ma.isMaskedArray(even_projector):
        raise ValueError('masked field operators are not complete projectors')
    value = np.asarray(even_projector)
    if (value.ndim != 2 or value.shape[0] != value.shape[1]
            or not np.isfinite(value).all()
            or np.linalg.norm(value-value.conj().T) > 1e-12
            or np.linalg.norm(value@value-value) > 1e-12):
        raise ValueError('a complete Hermitian parity projector is required')
    identity = np.eye(len(value))
    return (7*identity-6*value)/np.sqrt(50), 1j*(identity+6*value)/np.sqrt(50)


def mean_matrix(parity_probability):
    if type(parity_probability) is not F or not 0 <= parity_probability <= 1:
        raise ValueError('an exact parity probability in [0,1] is required')
    from model.matter_geometry_response.geometry import laplacian
    theta = F(3, 25)*(1-2*parity_probability)
    matrix = laplacian(theta)
    return [[F(int(i == j))-matrix[i][j]/64 for j in range(12)] for i in range(12)]
