"""Preparation and readout by contracted twelve-port observer source words.

The classical thermal draw, blank register and coherent proper-code services
are supplied apparatus. Every branch's Fock preparation, displacement,
quadrature rotation and two-outcome detector is nevertheless executed from
native source pulses. No fitted density matrix or analytic characteristic
function is inserted in the prediction path.

The d=64 basis has 1,162,054 native pulses per direction. A complete cold
run can take several minutes; repeated settings share already contracted
operation matrices. The original pulse counts and full-phase guards remain
unchanged. Single-threaded BLAS avoids excessive overhead on these small
matrices (for example OPENBLAS_NUM_THREADS=1).
"""
from __future__ import annotations

from copy import deepcopy
from functools import lru_cache
import hashlib
import json
import math
from numbers import Integral, Real

import numpy as np
from scipy.linalg import expm

from ...observer.source_scalar import compiler as c
from ..source import _basis
from . import canonical_source


MAXIMUM_DIMENSION = 64
# Fixed host-contraction grouping, chosen before execution or target checks.
# Every pulse remains in the original word and every block uses the unchanged
# compiler guard. This does not alter the requested physical source operations.
MAXIMUM_CONTRACTION_BLOCK_PULSES = 4096
COUNT_KEYS = ('native_pulses', 'single_code_pulses', 'pair_code_pulses',
              'identity_phase_pulses', 'proper_code_load_events',
              'proper_code_export_events', 'proper_code_service_events')


def _dimension(value):
    if (isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral)
            or value < 2 or value > MAXIMUM_DIMENSION or value & (value - 1)):
        raise ValueError('a power-of-two occupation dimension from 2 through 64 is required')
    return int(value)


def _real(value, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real) or not math.isfinite(value):
        raise ValueError(f'{name} must be a finite real number')
    return float(value)


def _word_hash(word):
    payload = json.dumps(word, sort_keys=True, separators=(',', ':'), allow_nan=False)
    return hashlib.sha256(payload.encode()).hexdigest()


def _contract(word, qubits, target, *, name):
    """Execute fixed pulse blocks with separately retained floating drift.

    R_k is an executed raw block, Q_k its guarded polar representation.
    For exact matrix products, E_k <= ||Q_k|| E_(k-1)
      + ||Q_k-R_k|| product_(i<k)||R_i||
    bounds the complete difference. The evaluated expression below, like
    all the norms and products, is a floating diagnostic, not an interval
    certificate. Raw and corrected full-phase target gates remain mandatory.
    """
    raw = corrected = None
    guards, counts = [], []
    norm_product, telescope = 1., 0.
    for start in range(0, max(1, len(word)), MAXIMUM_CONTRACTION_BLOCK_PULSES):
        block = word[start:start+MAXIMUM_CONTRACTION_BLOCK_PULSES]
        executed = c.execute(block, qubits)
        represented, guard = c.stabilize_executed_unitary(executed)
        telescope = (float(np.linalg.norm(represented, 2))*telescope
                     + guard['projection_change']*norm_product)
        norm_product *= float(np.linalg.norm(executed, 2))
        raw = executed if raw is None else executed @ raw
        corrected = represented if corrected is None else represented @ corrected
        guards.append(guard)
        counts.append(len(block))
    error = float(np.linalg.norm(raw - target, 2))
    if not math.isfinite(error) or error > 2e-10:
        raise ArithmeticError(f'{name}: native contraction differs from its independent target')
    if len(guards) == 1:
        # This one block already passed the final compiler guard. Preserve
        # the existing short-word contraction, including its numerical data.
        matrix, stabilization, final_change = corrected, guards[0], 0.
    else:
        matrix, stabilization = c.stabilize_executed_unitary(corrected)
        final_change = stabilization['projection_change']
    corrected_error = float(np.linalg.norm(matrix - target, 2))
    if not math.isfinite(corrected_error) or corrected_error > 2e-10:
        raise ArithmeticError(f'{name}: corrected native contraction differs from its independent target')
    raw_defect = float(np.linalg.norm(raw.conj().T @ raw-np.eye(len(raw)), 2))
    complete_change = float(np.linalg.norm(matrix-raw, 2))
    telescope += final_change
    if not all(math.isfinite(value) for value in (raw_defect, complete_change, telescope)):
        raise ArithmeticError(f'{name}: nonfinite complete contraction diagnostic')
    matrix.flags.writeable = False
    return matrix, {
        'name': name, 'word_sha256': _word_hash(word),
        'resources': c.resources(word), 'register_qubits': qubits,
        'full_phase_target_residual': error,
        'corrected_full_phase_target_residual': corrected_error,
        'numerical_stabilization': stabilization,
        'blockwise_contraction': {
            'method': 'fixed chronological pulse blocks; raw and guarded products retained',
            'maximum_native_pulses_per_block': MAXIMUM_CONTRACTION_BLOCK_PULSES,
            'block_pulse_counts': counts,
            'block_count': len(guards),
            'maximum_block_raw_unitarity_defect': max(g['raw_unitarity_defect'] for g in guards),
            'maximum_block_projection_change': max(g['projection_change'] for g in guards),
            'maximum_block_restored_unitarity_defect': max(g['restored_unitarity_defect'] for g in guards),
            'sum_block_projection_changes': math.fsum(g['projection_change'] for g in guards),
            'final_projection_change': final_change,
            'raw_block_product_unitarity_defect': raw_defect,
            'complete_corrected_minus_raw_norm': complete_change,
            'algebraic_telescoping_expression_evaluated': telescope,
            'raw_product_grouping': 'product of executed blocks, not monolithic pulse-by-pulse contraction',
            'floating_matrix_product_roundoff_in_telescoping_expression': False,
            'outward_roundoff_or_physical_error_certificate': False,
            'source_word_or_resource_counts_changed': False,
        },
        'matrix_origin': 'contracted_native_source_pulses',
        'physical_service_error_certified': False,
    }


def _ladder(dimension):
    return np.diag(np.sqrt(np.arange(1, dimension)), 1)


@lru_cache(maxsize=32)
def _displacement(dimension, alpha_real):
    ladder = _ladder(dimension)
    target = expm(alpha_real * (ladder.T - ladder))
    word = canonical_source.displacement_word(dimension, alpha_real)
    return _contract(word, dimension.bit_length() - 1, target, name='real_displacement')


@lru_cache(maxsize=8)
def _coordinate_basis(dimension):
    ladder = _ladder(dimension)
    coordinate = (ladder + ladder.T) / 2
    values, basis = _basis(coordinate)
    word = canonical_source.coordinate_word(dimension)
    bits = dimension.bit_length() - 1
    forward, metadata = _contract(word, bits, basis, name='coordinate_basis')
    backward, inverse_metadata = _contract(c.inverse(word), bits, basis.T, name='coordinate_inverse_basis')
    return values, forward, backward, metadata, inverse_metadata


@lru_cache(maxsize=8)
def _fock_preparations(dimension):
    bits = dimension.bit_length() - 1
    columns = []
    records = []
    blank = np.zeros(dimension, complex)
    blank[0] = 1
    for occupation in range(dimension):
        # exp(-i*pi/2*I) exp(i*pi/2*X)=X, with the scalar phase retained.
        word = []
        for bit in range(bits):
            if occupation >> bit & 1:
                word.extend([c.pulse('X', math.pi / 2, [bit]),
                             c.pulse('I', -math.pi / 2, [bit])])
        state = c.execute(word, bits, blank)
        target = np.zeros(dimension, complex)
        target[occupation] = 1
        residual = float(np.linalg.norm(state - target))
        if residual > 2e-12:
            raise ArithmeticError('native Fock preparation failed')
        columns.append(state)
        records.append({'occupation': occupation, 'word_sha256': _word_hash(word),
                        'resources': c.resources(word), 'full_phase_target_residual': residual})
    matrix = np.column_stack(columns)
    matrix.flags.writeable = False
    return matrix, records


def thermal_weights(dimension, thermal):
    """A supplied capped geometric draw, with no rejected preparation branch.

    Occupations >= dimension-1 are routed to the last code state. The trace
    distance to the unbounded thermal mixture before displacement is r**d.
    The finite classical weights are defined directly, not by conditioning
    on an accepted draw or renormalizing a quantum state.
    """
    dimension = _dimension(dimension)
    thermal = _real(thermal, 'thermal occupation')
    if thermal < 0:
        raise ValueError('thermal occupation must be nonnegative')
    ratio = thermal / (1 + thermal)
    weights = [(1 - ratio) * ratio ** n for n in range(dimension - 1)]
    weights.append(ratio ** (dimension - 1))
    return np.array(weights), {
        'rule': 'capped_geometric_classical_preparation',
        'mean_thermal_occupation_parameter': thermal,
        'ratio': ratio,
        'classical_branch_probabilities': weights,
        'initial_thermal_trace_distance': ratio ** dimension,
        'thermal_draw_supplied': True,
        'no_preparation_branch_discarded': True,
        'state_or_outcome_renormalization': False,
    }


def preparation(dimension, alpha_real, thermal):
    """Return actually prepared state columns and supplied classical weights."""
    dimension = _dimension(dimension)
    alpha_real = _real(alpha_real, 'real displacement')
    weights, thermal_metadata = thermal_weights(dimension, thermal)
    fock, fock_metadata = _fock_preparations(dimension)
    displacement, displacement_metadata = _displacement(dimension, alpha_real)
    states = displacement @ fock
    metadata = {
        'dimension': dimension, 'alpha_real': alpha_real,
        'thermal': thermal_metadata,
        'fock_branches': deepcopy(fock_metadata),
        'displacement': deepcopy(displacement_metadata),
        'weighted_prepared_mass': float(np.dot(weights, np.sum(abs(states) ** 2, axis=0))),
        'blank_register_and_classical_draw_supplied': True,
        'initialization_storage_readout_physical_costs_included': False,
        'postselection': False,
    }
    return states, weights, metadata


def _setting(setting):
    if not isinstance(setting, dict) or set(setting) != {'u', 'v', 'theta'}:
        raise ValueError('each setting must contain exactly u, v and theta')
    return {key: _real(setting[key], key) for key in ('u', 'v', 'theta')}


def source_run(dimension, alpha_real, thermal, settings):
    """Emit binary probabilities from source-executed preparation and detector.

    ``u+i*v`` denotes D(u+i*v), and the bit difference is
    Re(exp(i*theta)*chi(u+i*v)); theta=0 and -pi/2 yield real and imaginary
    components. The quadrature rotation angle is a setting, not model time.
    This function accepts no observed means or withheld response values.
    """
    dimension = _dimension(dimension)
    if not isinstance(settings, (tuple, list)) or not settings:
        raise ValueError('a nonempty list of detector settings is required')
    settings = [_setting(setting) for setting in settings]
    states, weights, preparation_metadata = preparation(dimension, alpha_real, thermal)
    values, basis, inverse_basis, basis_metadata, inverse_metadata = _coordinate_basis(dimension)
    bits = dimension.bit_length() - 1
    pointer = bits
    inverse_axis_word = [c.pulse('X', -math.pi / 4, [pointer])]
    axis_word = [c.pulse('X', math.pi / 4, [pointer])]
    inverse_axis, inverse_axis_metadata = _contract(
        inverse_axis_word, bits + 1,
        np.kron(expm(-1j * math.pi / 4 * np.array([[0, 1], [1, 0]])), np.eye(dimension)),
        name='pointer_inverse_axis')
    axis, axis_metadata = _contract(
        axis_word, bits + 1,
        np.kron(expm(1j * math.pi / 4 * np.array([[0, 1], [1, 0]])), np.eye(dimension)),
        name='pointer_axis')
    records = []
    for setting in settings:
        u, v, theta = (setting[key] for key in ('u', 'v', 'theta'))
        phi, radius, rotation_word, conditional_word = canonical_source.quadrature_words(
            dimension, u, v, theta)
        # Independent binary64 targets remain checks, never the word producer.
        target_radius = math.hypot(u, v)
        target_phi = 0. if target_radius == 0 else math.atan2(-u, v)
        phases = -target_phi * np.arange(dimension)
        rotation, rotation_metadata = _contract(rotation_word, bits, np.diag(np.exp(1j * phases)),
                                                name='quadrature_rotation')
        rotated = inverse_basis @ rotation @ states
        joint = np.vstack([rotated, np.zeros_like(rotated)])
        angles = target_radius * values + theta / 2
        pointer_phases = np.r_[-angles, angles]
        conditional, conditional_metadata = _contract(
            conditional_word, bits + 1, np.diag(np.exp(1j * pointer_phases)),
            name='conditional_detector_phase')
        joint = axis @ conditional @ inverse_axis @ joint
        branch_zero, branch_one = basis @ joint[:dimension], basis @ joint[dimension:]
        branch_probabilities = np.array([np.sum(abs(branch_zero) ** 2, axis=0),
                                        np.sum(abs(branch_one) ** 2, axis=0)])
        probabilities = branch_probabilities @ weights
        stages = [preparation_metadata['displacement'], rotation_metadata,
                  inverse_metadata, inverse_axis_metadata, conditional_metadata,
                  axis_metadata, basis_metadata]
        common_counts = {key: sum(stage['resources'][key] for stage in stages)
                         for key in COUNT_KEYS}
        branch_counts = [{key: common_counts[key] + branch['resources'][key]
                          for key in COUNT_KEYS}
                         for branch in preparation_metadata['fock_branches']]
        resources = {key: {
            'maximum_per_shot': max(branch[key] for branch in branch_counts),
            'classical_weighted_mean_per_shot': float(sum(weight * branch[key]
                       for weight, branch in zip(weights, branch_counts, strict=True))),
        } for key in COUNT_KEYS}
        records.append({
            **setting, 'quadrature_angle': phi, 'detector_strength_eta': 2 * radius,
            'probabilities': probabilities.tolist(),
            'mean': float(probabilities[0] - probabilities[1]),
            'probability_mass': float(sum(probabilities)),
            'probabilities_by_initial_occupation': branch_probabilities.T.tolist(),
            'rotation': rotation_metadata, 'conditional_detector_phase': conditional_metadata,
            'chronological_word_hashes_after_fock_preparation': [stage['word_sha256'] for stage in stages],
            'resources': resources,
            'native_pulses_by_initial_occupation': [branch['native_pulses'] for branch in branch_counts],
            'both_detector_outcomes_retained': True,
            'postselection': False,
        })
    return {
        'schema': 'source_executed_oscillator_readout_v1',
        'dimension': dimension,
        'preparation': preparation_metadata,
        'shared_detector_operations': {
            'coordinate_basis': deepcopy(basis_metadata),
            'coordinate_inverse_basis': deepcopy(inverse_metadata),
            'pointer_inverse_axis': inverse_axis_metadata, 'pointer_axis': axis_metadata,
        },
        'records': records,
        'prediction_origin': 'norms_of_source_executed_pointer_branches',
        'canonical_source_compilation': canonical_source.policy(),
        'analytic_characteristic_used_for_prediction': False,
        'source_summary': c.source_summary(),
        'supplied_apparatus': ['blank register', 'classical thermal draw',
                              'coherent proper-code load/export', 'pointer readout'],
        'service_failure_outcomes_retained': True,
        'valid_code_failure_probability_in_exact_arithmetic': 0,
        'physical_service_error_certified': False,
        'occupation_error_certified_here': False,
        'time_parameter': 'none; quadrature rotation is an apparatus setting',
    }
