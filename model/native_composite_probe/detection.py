"""Full native instrument with the unobserved field kept as a quantum reference."""
from functools import lru_cache
from fractions import Fraction as F
from pathlib import Path

import numpy as np

from ..centered_native_preparation import packet as prepared, target, tree
from ..charged_composite_search import packet as spectral
from ..native_spectral_detector import execution
from . import circuits

ROOT = Path(__file__).resolve().parents[2]


@lru_cache(maxsize=1)
def field_indices():
    """Occupation-qubit permutation, with the CAR signs inside the operator."""
    selected = circuits.FIELD_MODES
    complement = tuple(j for j in range(16) if j not in selected)
    masks = np.arange(65536, dtype=np.int64)
    local = sum(((masks >> bit) & 1) << j for j, bit in enumerate(selected))
    reference = sum(((masks >> bit) & 1) << j for j, bit in enumerate(complement))
    if len(np.unique(local+256*reference)) != 65536:
        raise ArithmeticError('field/reference permutation is not complete')
    return local, reference


def as_reference_columns(field):
    value = np.asarray(field)
    if np.ma.isMaskedArray(field) or value.shape not in ((65536,), (65536, 1)) or not np.isfinite(value).all():
        raise ValueError('all 65536 field amplitudes are required')
    local, reference = field_indices()
    result = np.empty((256, 256), complex)
    result[local, reference] = value.reshape(-1)
    return result


def selector_permutation(selector):
    """Independent signed-Pauli truth table; no compiled pulse is consumed."""
    if type(selector) is not int or not 0 <= selector < 32:
        raise ValueError('one of all 32 pointer labels required')
    masks = np.arange(256, dtype=np.int64)
    output = masks.copy()
    sign = 1-2*(masks & 1)
    # X+iY is represented by X then Z; X-iY by Z then X.
    output ^= 2
    if selector & 1:
        sign *= 1-2*((output >> 1) & 1)
    neighbor = selector >> 3
    if neighbor < 3:
        up = (2, 4, 6)[neighbor]
        if selector & 2:
            sign *= 1-2*((output >> up) & 1)
        output ^= 1 << up
        output ^= 1 << (up+1)
        if selector & 4:
            sign *= 1-2*((output >> (up+1)) & 1)
    return output, sign


def selector_preparation():
    h = np.array([[1., 1.], [1., -1.]])/np.sqrt(2.)
    r = np.array([[1., -1.], [1., 1.]])/np.sqrt(2.)
    high = np.array([[np.sqrt(2/3), -np.sqrt(1/3)],
                     [np.sqrt(1/3), np.sqrt(2/3)]])
    low = np.eye(4)
    low[:2, :2] = r
    return np.kron(low @ np.kron(high, np.eye(2)), np.kron(np.kron(h, h), h))


def exact_selector_identity():
    """Check every occupation matrix entry against independent signed CAR."""
    car = np.zeros((256, 256), np.int64)
    for mask in range(256):
        for up in (2, 4, 6):
            current, sign = mask, 1
            for mode, creation in ((1, False), (up+1, False), (up, True)):
                if bool(current & (1 << mode)) == creation:
                    sign = 0
                    break
                sign *= (-1)**((current & ((1 << mode)-1)).bit_count())
                current ^= 1 << mode
            car[current, mask] += sign
    selected = np.zeros_like(car)
    for label in range(24):
        output, sign = selector_permutation(label)
        selected[output, np.arange(256)] += sign
    if not np.array_equal(selected, 8*car):
        raise ArithmeticError('coherent signed selection does not equal the cubic CAR probe')
    return car


def native_pointer_probabilities(logical):
    if np.shape(logical) != (32,) or not np.isfinite(logical).all() or np.min(logical) < 0:
        raise ValueError('all nonnegative logical pointer probabilities required')
    return np.array([logical[sum(v << j for j, v in enumerate(raw))]
                     if all(v < 2 for v in raw) else 0.
                     for raw in np.ndindex(*(6,)*5)], np.float64)


def ideal_output(field):
    columns = as_reference_columns(field)
    prepare = selector_preparation()
    intermediate = prepare[:, 0, None, None]*columns[None, :, :]
    routed = np.empty_like(intermediate)
    for selector in range(32):
        output, sign = selector_permutation(selector)
        routed[selector, output, :] = sign[:, None]*intermediate[selector]
    return (prepare.conj().T @ routed.reshape(32, -1)).reshape(32, 256, 256)


def complete_native(field, *, chunk_size=32):
    """Partition only an untouched reference basis; preserve its coherences."""
    if type(chunk_size) is not int or not 1 <= chunk_size <= 32:
        raise ValueError('one through 32 coherent reference columns per worker')
    columns = as_reference_columns(field)
    result = np.empty((8192, 256), complex)
    failures = None
    maximum_accounting_residual = 0.
    negative_slots = 0
    negative_mass = 0.
    blocks = circuits.blocks()
    for start in range(0, 256, chunk_size):
        stop = min(256, start+chunk_size)
        initial = np.zeros((8192, stop-start), complex)
        initial[:256] = columns[:, start:stop]
        run = execution.complete(blocks, circuits.BITS, initial)
        result[:, start:stop] = run.success
        failures = (np.array(run.failure_weights) if failures is None
                    else failures+run.failure_weights)
        maximum_accounting_residual = max(maximum_accounting_residual, run.normalization_error)
        negative_slots += run.rounded_negative_weight_count
        negative_mass += run.rounded_negative_weight_total
    output = result.reshape(32, 256, 256)
    pointers = np.sum(abs(output)**2, axis=(1, 2))
    original = float(np.vdot(columns, columns).real)
    total = float(pointers.sum()+failures.sum())
    if abs(total-original) > execution.ADMISSION:
        raise ArithmeticError('complete composite native output loses probability')
    native_slots = sum(sum(len(branches) for branches in instrument.failures)
                       for instrument, _ in blocks)
    arrays = {'complete_pointer_field_reference_amplitudes': output,
              'first_export_failure_weights': failures, 'pointer_probabilities': pointers,
              'native_pointer_probabilities': native_pointer_probabilities(pointers)}
    result = {'complete_joint_amplitudes': int(output.size), 'logical_pointer_labels': 32,
              'native_pointer_outcomes': 7776, 'native_pointer_exact_zero_slots': 7744,
              'local_field_rows': 256, 'coherent_reference_columns': 256,
              'coherent_reference_chunk_size': chunk_size,
              'native_pulses': len(failures), 'native_failure_slots': native_slots,
              'padding_slots': int(failures.size-native_slots),
              'first_export_failure_probability': float(failures.sum()),
              'all_pointer_probability': float(pointers.sum()),
              'complete_probability_residual': abs(total-original),
              'largest_chunk_probability_residual': maximum_accounting_residual,
              'rounded_negative_slot_count': negative_slots, 'rounded_negative_slot_mass': negative_mass,
              'success_pointer_probability': float(pointers[0]),
              'every_reference_column_is_coherent_not_a_classical_sample': True,
              'success_postselection_or_renormalization': False}
    return result, arrays


def inherited_input():
    """Authenticate the retained literal tree; do not rerun its 262137 pulses."""
    receipt, arrays = prepared.read_pair(prepared.OUTPUT)
    values = arrays['hidden_success']
    if values.shape != (65536, 1) or not np.isfinite(values).all():
        raise ValueError('literal full hidden-frame preparation is missing')
    tree_proof = tree.certificate()
    if receipt['complete_error']['preparation']['Gray_tree'] != tree_proof:
        raise ValueError('inherited tree certificate differs from the literal source')
    target_proof = target.summary()
    error = F(tree_proof['complete_tree_error'])+F(target_proof['normalized_target_to_exact_ground_error'])
    return receipt, values, error, float(arrays['tree_failure_weights'].sum())


def failure_reference_amplitudes(field, block_index, event, slot):
    """Replay one actual abort, preserving all untouched reference coherences."""
    blocks = circuits.blocks()
    if type(block_index) is not int or not 0 <= block_index < len(blocks):
        raise ValueError('actual source block required')
    instrument, _ = blocks[block_index]
    if (type(event) is not int or not 0 <= event < len(instrument.tape)
            or type(slot) is not int or not 0 <= slot < len(instrument.failures[event])):
        raise ValueError('actual native failure required; padding is excluded')
    columns = as_reference_columns(field)
    parts = []
    for start in range(0, 256, 32):
        initial = np.zeros((8192, 32), complex)
        initial[:256] = columns[:, start:start+32]
        run = execution.complete(blocks[:block_index+1], circuits.BITS, initial)
        parts.append(run.failure(block_index, event, slot))
    return np.concatenate(parts, axis=1)


def report():
    receipt, field, input_error, tree_abort = inherited_input()
    exact_selector_identity()
    controls = circuits.certificate()
    result, arrays = complete_native(field)
    comparison = ideal_output(field)
    diagnostic = float(np.sqrt(np.linalg.norm(arrays['complete_pointer_field_reference_amplitudes']-comparison)**2
                               +result['first_export_failure_probability']))
    if diagnostic > execution.ADMISSION:
        raise ArithmeticError('native composite execution exceeds its numerical admission')
    ideal_ground = ideal_output(target.ideal_hidden())
    source = spectral.strict_json.load(spectral.OUTPUT)['spectrum']
    lo, hi = map(F, source['scaled_probe_success_interval'])
    exact_probability_diagnostic = float(np.vdot(ideal_ground[0], ideal_ground[0]).real)
    if not float(lo)-1e-12 <= exact_probability_diagnostic <= float(hi)+1e-12:
        raise ArithmeticError('independent cubic spectral measure and LCU filter disagree')
    error = input_error+F(controls['complete_operator_and_untouched_reference_error'])
    result.update({'native_to_ideal_complete_output_error_diagnostic': diagnostic,
                   'ideal_filter_probability_diagnostic': exact_probability_diagnostic,
                   'inherited_tree_abort_probability': tree_abort,
                   'preparation_and_detector_probability_residual':
                       abs(result['all_pointer_probability']+result['first_export_failure_probability']+tree_abort-1),
                   'source_success_probability_interval': source['scaled_probe_success_interval'],
                   'complete_state_and_probability_error_upper': str(error),
                   'inherited_preparation_error_upper': str(input_error),
                   'inherited_preparation_pulses': 262137,
                   'inherited_preparation_reexecuted': False,
                   'inherited_source_custody_verified': receipt['dependency_custody_pass'],
                   'physical_species_or_mass_identification': False,
                   'native_energy_clock_or_resolved_138_lines_executed': False})
    return {'control_certificate': controls, 'execution': result}, arrays
