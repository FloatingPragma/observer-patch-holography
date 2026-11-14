"""Source-compiled native quantum repairs around the canonical field history.

Factorized Kraus maps retain the full nine-qubit native output, the disturbed
field and every classical outcome. Factoring tensor products is a simulation
optimization, not a measurement of the native input. The terminal interference
read measures all nine code bits after one source rotation; it is deliberately
not the earlier two-outcome nondemolition instrument.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import json
import math
from pathlib import Path

import numpy as np

from tools import strict_json
from ..observer import native
from ..observer.source_discrimination import circuit
from ..observer.source_scalar import compiler
from ..observer.scalar_response.experiment import canonical, compare_receipt
from . import composition, continuation_bounds as bounds, protocol
from . import native_routing as routing


STATES = native.states()
INDEX = {state: index for index, state in enumerate(STATES)}
CODES = circuit.SOURCE_CODES
TOLERANCE = 2e-8
STEPS = 512
NATIVE_MODES = ('full_carrier', 'protected_qubit', 'passive_qubit')


def _mode(value):
    if type(value) is not str or value not in NATIVE_MODES:
        raise ValueError('an explicit native input and recovery mode required')
    return value


def input_dimension(mode):
    return 78 if _mode(mode) == 'full_carrier' else 2


def final_seam(mode, scalar_bit):
    _mode(mode)
    if type(scalar_bit) is not int or scalar_bit not in (0, 1):
        raise ValueError('a committed binary smooth record required')
    return ((0, 5) if mode == 'full_carrier' else (2, 3)) if scalar_bit == 0 else (0, 1)


def _law(value):
    if type(value) is not str or value not in ('C', 'E'):
        raise ValueError('native instrument C or E required')
    return value


def _array(value, shape, name):
    routing._original_scalars(value)
    if np.ma.isMaskedArray(value):
        raise ValueError(f'{name} cannot be masked')
    raw = np.asarray(value)
    if (raw.shape != shape or raw.dtype.kind not in 'iufc'
            or not np.isfinite(raw).all()):
        raise ValueError(f'{name} must be a finite numeric array of shape {shape}')
    result = np.array(raw, dtype=complex, copy=True)
    if not np.isfinite(result).all():
        raise ValueError(f'{name} must remain finite in the execution representation')
    return result


def _immutable(value):
    raw = np.ascontiguousarray(value, dtype=complex)
    return np.frombuffer(raw.tobytes(), dtype=complex).reshape(raw.shape)


def _programs(dimension, steps, backend):
    if type(backend) is not str or backend not in protocol.BACKENDS:
        raise ValueError('an admitted field execution backend required')
    return {current: composition._executed(protocol.continuation_program(
                dimension, steps, native_current=current), backend)
            for current in (-1, 0, 1)}


@dataclass(frozen=True)
class JointBranch:
    # (first coin, first raw code, quantum record, final coin, final raw code)
    # An invalid first raw code halts: quantum record=(), final fields=None.
    record: tuple
    abort: bool
    native_map: np.ndarray
    work_map: np.ndarray
    work_state: np.ndarray


@dataclass(frozen=True)
class JointHistory:
    law: str
    dimension: int
    steps: int
    field_backend: str
    branches: tuple[JointBranch, ...]
    native_mode: str = 'full_carrier'


def _injection():
    result = np.zeros((512, 78), complex)
    result[list(CODES), np.arange(78)] = 1
    return result


@lru_cache(maxsize=6)
def _native_tree(law, mode='full_carrier'):
    """Execute every native outcome on the declared input columns."""
    _law(law)
    _mode(mode)
    if mode == 'full_carrier':
        repair, initial = routing.routed_repair_branches, _injection()
    else:
        from . import encoded_memory as memory
        initial = memory.code_matrix()

        def repair(candidate, seam, state):
            return memory.repair_branches(candidate, seam, state,
                                         recover=mode == 'protected_qubit')
    first = repair(law, (0, 1), initial)
    rows = []
    for start in first:
        a = _immutable(start['amplitude'])
        endings = {}
        if not start['abort']:
            for scalar_bit in (0, 1):
                seam = final_seam(mode, scalar_bit)
                endings[scalar_bit] = tuple(
                    (end['coin'], end['record_code'], end['abort'],
                     _immutable(end['amplitude']))
                    for end in repair(law, seam, a))
        rows.append((start['coin'], start['record_code'], start['current'],
                     start['abort'], a, endings))
    return tuple(rows)


def execute(law='C', dimension=2, steps=STEPS, *, field_backend='source',
            native_mode='full_carrier'):
    """Complete finite channel on native input and independently prepared work.

    Every ``native_map tensor work_map`` also defines the finite channel on
    arbitrary native/work/reference input; that override omits vacuum
    preparation and does not inherit its occupation error bound.
    """
    law, native_mode = _law(law), _mode(native_mode)
    histories = _programs(dimension, steps, field_backend)
    prepared = histories[0].initial_state
    identity = _immutable(np.eye(4*dimension**2))
    branches = []
    for c1, raw1, current, abort, a, endings in _native_tree(law, native_mode):
        if abort:
            branches.append(JointBranch((c1, raw1, (), None, None), True,
                                        a, identity, prepared))
            continue
        for middle in histories[current].branches:
            for c2, raw2, stopped, combined in endings[middle.record[2]]:
                branches.append(JointBranch((c1, raw1, middle.record, c2, raw2),
                                            stopped, combined, middle.kraus, middle.state))
    result = JointHistory(law, dimension, steps, field_backend, tuple(branches), native_mode)
    _validate_menu(result)
    return result


def _record_current(law, coin, raw):
    """Independent raw-record interpretation from the native endpoint law.

    This checker must not reuse the execution decoder: a changed raw/current
    dictionary must not redefine which physical preparation the oracle expects.
    E retains parity p and contrast 2*j+(2*c-1)*p, bounded by total load two.
    """
    if type(coin) is not int or coin not in (0, 1) or type(raw) is not int:
        return None
    current = {0: 0, 1: 1, 2: -1}.get(raw % 4)
    if current is None or not 0 <= raw < (4 if law == 'C' else 8):
        return None
    if law == 'E' and abs(2*current+(2*coin-1)*(raw//4)) > 2:
        return None
    return current


def _menu(law):
    records = 1 << circuit.RECORD_BITS[_law(law)]
    result = {}
    for c1, raw1 in product((0, 1), range(records)):
        if _record_current(law, c1, raw1) is None:
            result[(c1, raw1, (), None, None)] = True
        else:
            for bits in product((0, 1), repeat=4):
                for c2, raw2 in product((0, 1), range(records)):
                    result[(c1, raw1, bits, c2, raw2)] = _record_current(law, c2, raw2) is None
    return result


def _validate_menu(history):
    if type(history) is not JointHistory:
        raise ValueError('an executed JointHistory required')
    protocol.continuation_program(history.dimension, history.steps)
    if history.field_backend not in protocol.BACKENDS:
        raise ValueError('unknown field backend')
    expected = _menu(history.law)
    native_size = input_dimension(history.native_mode)
    seen = set()
    size = 4*history.dimension**2
    if type(history.branches) is not tuple:
        raise ValueError('immutable complete branch sequence required')
    for branch in history.branches:
        if type(branch) is not JointBranch:
            raise ValueError('typed joint branches required')
        key = branch.record
        if (type(key) is not tuple or len(key) != 5 or type(key[2]) is not tuple
                or any(type(x) is not int for x in (*key[:2], *key[2]))
                or any(x is not None and type(x) is not int for x in key[3:])
                or key not in expected or key in seen or type(branch.abort) is not bool
                or branch.abort != expected[key]):
            raise ValueError('invalid, duplicate, or relabelled joint outcome')
        seen.add(key)
        _array(branch.native_map, (512, native_size), 'native branch')
        _array(branch.work_map, (size, size), 'work branch')
        _array(branch.work_state, (size, 1), 'prepared work branch')
    if seen != set(expected):
        raise ValueError('every raw outcome, including aborts, must be retained')


def _input(value, size):
    if np.ma.isMaskedArray(value):
        raise ValueError('native input cannot be masked')
    raw = np.asarray(value)
    if raw.ndim not in (1, 2) or raw.shape[0] != size:
        raise ValueError('native vector or native-by-reference columns required')
    result = _array(value, raw.shape, 'native input').reshape(size, -1)
    norm_squared = float(np.vdot(result, result).real)
    if not result.shape[1] or not math.isfinite(norm_squared) or abs(norm_squared-1) > 2e-10:
        raise ValueError('normalized nonempty native/reference input required')
    return result


def output(history, initial):
    """Return factorized full native/reference and prepared work amplitudes."""
    _validate_menu(history)
    initial = _input(initial, input_dimension(history.native_mode))
    return {branch.record: (branch.native_map@initial, branch.work_state)
            for branch in history.branches}


def output_on_joint_input(history, initial):
    """Iterate the full native/work/reference output without a vacuum claim.

    Streaming bounds host storage to one branch and retains every outcome.
    """
    _validate_menu(history)
    size = 4*history.dimension**2
    if np.ma.isMaskedArray(initial):
        raise ValueError('joint input cannot be masked')
    raw = np.asarray(initial)
    if raw.ndim != 3 or raw.shape[:2] != (input_dimension(history.native_mode), size) or not raw.shape[2]:
        raise ValueError('native-by-work-by-reference amplitudes required')
    value = _array(initial, raw.shape, 'joint input')
    norm_squared = float(np.vdot(value, value).real)
    if not math.isfinite(norm_squared) or abs(norm_squared-1) > 2e-10:
        raise ValueError('normalized joint input required')
    return ((b.record, np.einsum('ai,bj,ijr->abr', b.native_map, b.work_map, value,
                                 optimize=True)) for b in history.branches)


def _endpoint(state, seam, coin):
    # Independent floor/ceiling rule; does not call native.native or compiler.
    a, b = seam
    total = state[a]+state[b]
    result = list(state)
    result[a] = total//2+coin*(total % 2)
    result[b] = total-result[a]
    return tuple(result), state[a]-result[a], state[a]-state[b]


def _raw(law, current, contrast):
    code = {0: 0, 1: 1, -1: 2}[current]
    return code if law == 'C' else code+4*(abs(contrast) % 2)


@lru_cache(maxsize=2)
def _ideal_native_tree(law):
    """Independent exact half-amplitude partial permutations on every input."""
    _law(law)
    maps = {}
    for c1, c2, selected in product((0, 1), repeat=3):
        for column, initial in enumerate(STATES):
            middle, j1, d1 = _endpoint(initial, (0, 1), c1)
            end, j2, d2 = _endpoint(middle, (0, 5) if selected == 0 else (0, 1), c2)
            key = (c1, _raw(law, j1, d1), selected, c2, _raw(law, j2, d2))
            matrix = maps.setdefault(key, np.zeros((512, 78), complex))
            matrix[circuit.encode(end), column] = .5
    return {key: _immutable(value) for key, value in maps.items()}


def _logical_supports():
    # Independent statement of the existing four-configuration code. The
    # scientific oracle must not obtain its columns from the source encoder.
    return tuple(tuple(tuple(int(port in pair) for port in range(12))
                       for pair in pairs)
                 for pairs in (((0, 2), (1, 3)), ((0, 3), (1, 2))))


def _logical_code():
    value = np.zeros((512, 2), complex)
    for column, supports in enumerate(_logical_supports()):
        for state in supports:
            value[circuit.encode(state), column] = 1/math.sqrt(2)
    return value


def _code_weights(law, seam):
    counts = {}
    for column, supports in enumerate(_logical_supports()):
        for state in supports:
            for coin in (0, 1):
                _, current, contrast = _endpoint(state, seam, coin)
                key = (coin, _raw(law, current, contrast))
                counts.setdefault(key, [0, 0])[column] += 1
    if any(left != right for left, right in counts.values()):
        raise ArithmeticError('this native seam does not satisfy the declared logical code condition')
    return {key: F(values[0], 4) for key, values in counts.items()}


@lru_cache(maxsize=4)
def _ideal_encoded_tree(law, mode):
    _law(law)
    if _mode(mode) == 'full_carrier':
        raise ValueError('an encoded native mode required')
    maps = {}
    if mode == 'protected_qubit':
        code = _logical_code()
        for (c1, raw1), first in _code_weights(law, (0, 1)).items():
            for selected in (0, 1):
                # Pin the declared selector independently of the executor.
                seam = (2, 3) if selected == 0 else (0, 1)
                for (c2, raw2), last in _code_weights(law, seam).items():
                    maps[c1, raw1, selected, c2, raw2] = math.sqrt(float(first*last))*code
    else:
        for c1, c2, selected in product((0, 1), repeat=3):
            for column, supports in enumerate(_logical_supports()):
                for state in supports:
                    middle, j1, d1 = _endpoint(state, (0, 1), c1)
                    end, j2, d2 = _endpoint(middle, (2, 3) if selected == 0 else (0, 1), c2)
                    key = c1, _raw(law, j1, d1), selected, c2, _raw(law, j2, d2)
                    matrix = maps.setdefault(key, np.zeros((512, 2), complex))
                    matrix[circuit.encode(end), column] += 1/(2*math.sqrt(2))
    return {key: _immutable(value) for key, value in maps.items()}


def _oracle_branch(branch, law, histories, mode='full_carrier'):
    c1, raw1, bits, c2, raw2 = branch.record
    zero = np.zeros((512, input_dimension(mode)), complex)
    if not bits:
        size = histories[0].program.work_dimension
        return zero, np.eye(size), histories[0].initial_state
    current = _record_current(law, c1, raw1)
    if current is None:
        raise ValueError('an admitted raw first record is required for continued work')
    middle = next(b for b in histories[current].branches if b.record == bits)
    tree = _ideal_native_tree(law) if mode == 'full_carrier' else _ideal_encoded_tree(law, mode)
    a = tree.get((c1, raw1, bits[2], c2, raw2), zero)
    return a, middle.kraus, middle.state


def audit(history):
    """Validate all factored branch maps, not only diagonal probabilities.

    Frobenius bounds on the stacked Kraus differences upper-bound its operator
    norm in exact arithmetic. These evaluated residuals are floating diagnostics,
    not outward numerical certificates. Independent exact targets fix each phase.
    """
    _validate_menu(history)
    if history.field_backend != 'source':
        raise ValueError('the source implementation is required for this audit')
    targets = _programs(history.dimension, history.steps, 'product')
    native_size = input_dimension(history.native_mode)
    prepared_effect = np.zeros((native_size, native_size), complex)
    native_max = work_max = 0.
    errors, prepared_errors = [], []
    abort_mass = padding_mass = 0.
    padding = [i for i in range(512) if i not in set(CODES)]
    for branch in history.branches:
        a, b, v = _oracle_branch(branch, history.law, targets, history.native_mode)
        dn = float(np.linalg.norm(branch.native_map-a, 'fro'))
        dw = float(np.linalg.norm(branch.work_map-b, 'fro'))
        dv = float(np.linalg.norm(branch.work_state-v))
        native_max, work_max = max(native_max, dn), max(work_max, dw)
        errors.append((dn*np.linalg.norm(branch.work_map, 'fro')+np.linalg.norm(a, 'fro')*dw)**2)
        prepared_errors.append((dn*np.linalg.norm(branch.work_state)+np.linalg.norm(a, 'fro')*dv)**2)
        weight = float(np.vdot(branch.work_state, branch.work_state).real)
        prepared_effect += weight*(branch.native_map.conj().T@branch.native_map)
        mass = weight*float(np.linalg.norm(branch.native_map, 'fro')**2)/native_size
        abort_mass += mass*branch.abort
        padding_mass += weight*float(np.linalg.norm(branch.native_map[padding], 'fro')**2)/native_size
    result = {
        'full_outcome_count_including_abort': len(history.branches),
        'first_abort_outcome_count': sum(not b.record[2] for b in history.branches),
        'final_abort_outcome_count': sum(b.abort and bool(b.record[2]) for b in history.branches),
        'maximum_native_map_frobenius_error': native_max,
        'maximum_work_map_frobenius_error': work_max,
        'complete_stacked_map_frobenius_error_upper': math.sqrt(math.fsum(errors)),
        'prepared_complete_stacked_map_error_upper': math.sqrt(math.fsum(prepared_errors)),
        'prepared_native_input_completeness_operator_error': float(np.linalg.norm(prepared_effect-np.eye(native_size), 2)),
        'maximally_mixed_native_abort_probability': abort_mass,
        'maximally_mixed_native_padding_probability': padding_mass,
    }
    if any(not math.isfinite(x) or x > TOLERANCE for key, x in result.items()
           if key.endswith(('error', 'upper', 'probability'))):
        raise ArithmeticError('complete source-native channel differs from independent targets')
    return result


@lru_cache(maxsize=1)
def _read_rotation():
    # Full measured basis read, not a nondemolition two-outcome pointer.
    return _immutable(compiler.execute([compiler.pulse('Y', math.pi/4, [0])], 9))


def prepared_native(setting, *, native_mode='full_carrier'):
    _mode(native_mode)
    if type(setting) is not str or setting not in circuit.SETTINGS[:3]:
        raise ValueError('one of the three declared coherent preparations required')
    blank = np.zeros(512)
    blank[0] = 1
    result = compiler.execute(circuit.preparation_word(setting), 9, blank)
    retained = list(CODES) if native_mode == 'full_carrier' else [0, 1]
    result = _array(result, (512,), 'source native preparation')
    if np.count_nonzero(np.delete(result, retained)):
        raise ArithmeticError('native preparation escaped its input code')
    return result[retained]


def readout_outputs(history, initial):
    """Stream every final quantum/reference output, including stopped branches.

    Each entry is (old record, final code or None, native/reference factor,
    work factor). For a measured code the native factor is its reference row;
    the explicit code labels the native basis ket. An abort retains all 512
    native rows. No complementary outcome is removed or renormalized.
    """
    blocks = output(history, initial)
    rotation = _read_rotation()
    for branch in history.branches:
        a, v = blocks[branch.record]
        if branch.abort:
            yield branch.record, None, a, v
            continue
        if history.native_mode != 'full_carrier':
            from . import encoded_memory as memory
            a = memory.decoding(a)
        measured = rotation@a
        if not np.isfinite(measured).all():
            raise ArithmeticError('nonfinite complete native readout')
        for code in range(512):
            yield branch.record, code, measured[code:code+1], v


def readout(history, initial):
    """Summarize the complete measured output; no outcome-dependent selection."""
    masses, plus, abort = [], [], []
    for _, code, a, v in readout_outputs(history, initial):
        native_mass, work_mass = float(np.vdot(a, a).real), float(np.vdot(v, v).real)
        if not all(math.isfinite(mass) and mass >= 0 for mass in (native_mass, work_mass)):
            raise ArithmeticError('nonfinite native or work probability')
        weight = native_mass*work_mass
        if not math.isfinite(weight):
            raise ArithmeticError('nonfinite joint probability')
        masses.append(weight)
        if code is None:
            abort.append(weight)
        elif code == 0:
            plus.append(weight)
    total = math.fsum(masses)
    if not math.isfinite(total) or abs(total-1) > TOLERANCE:
        raise ArithmeticError('native readout discarded part of the attempted history')
    return {'probability_mass': total, 'plus_code_zero_probability': math.fsum(plus),
            'abort_probability': math.fsum(abort), 'complete_outcome_count': len(masses)}


def coherence_witness(histories):
    """A terminal observer read distinguishes quantum laws with the same loads."""
    if type(histories) is not dict or set(histories) != {'C', 'E'}:
        raise ValueError('both native law histories required')
    signature = {(h.dimension, h.steps, h.field_backend) for h in histories.values()}
    if (len(signature) != 1
            or any(h.law != law or h.native_mode != 'full_carrier' for law, h in histories.items())):
        raise ValueError('C/E must use the same field program and execution backend')
    d, steps, backend = next(iter(signature))
    middle = _programs(d, steps, backend)[0]
    q = math.fsum(b.probability for b in middle.branches if b.record[2] == 0)
    results = {}
    for law, history in histories.items():
        results[law] = {}
        for setting in circuit.SETTINGS[:3]:
            row = readout(history, prepared_native(setting))
            results[law][setting] = {key: value for key, value in row.items() if key != 'rows'}
    contrast = {law: rows['phase_0']['plus_code_zero_probability']
                     -rows['phase_pi']['plus_code_zero_probability'] for law, rows in results.items()}
    if (abs(contrast['C']-(1+q)/4) > TOLERANCE or abs(contrast['E']) > TOLERANCE
            or contrast['C'] < .24):
        raise ArithmeticError('quantum native coherence or its scalar-controlled routing was lost')
    for law, rows in results.items():
        midpoint = (rows['phase_0']['plus_code_zero_probability']+rows['phase_pi']['plus_code_zero_probability'])/2
        if abs(rows['phase_pi_over_2']['plus_code_zero_probability']-midpoint) > TOLERANCE:
            raise ArithmeticError('complex native phase control differs from independent prediction')
    return {'second_smooth_plus_probability_at_current_zero': q,
            'independent_C_phase_contrast_formula': (1+q)/4,
            'executed_phase_contrasts': contrast, 'executed_settings': results,
            'same_classical_load_law_selects_unique_quantum_instrument': False}


def resources(law, dimension, steps, *, native_mode='full_carrier'):
    law = _law(law)
    if _mode(native_mode) != 'full_carrier':
        from . import encoded_observer
        return encoded_observer.resources(law, dimension, steps, native_mode=native_mode)
    histories = _programs(dimension, steps, 'source')
    first = routing.routed_repair_resources(law, (0, 1))
    last = routing.routed_repair_resources(law, (0, 5))
    work = {key: max(dict(branch.resources)[key] for h in histories.values() for branch in h.branches)
            for key in protocol.RESOURCE_KEYS}
    # Three native preparation pulses and one terminal read rotation.
    total_pulses = work['native_pulses']+first['native_pulses']+last['native_pulses']+4
    return {'maximum_native_pulses_per_attempt': total_pulses,
            'maximum_proper_code_load_events_per_attempt': total_pulses,
            'maximum_proper_code_export_events_per_attempt': total_pulses,
            'maximum_proper_code_service_events_per_attempt': 2*total_pulses,
            'field_path_maxima': work, 'first_native_repair': first,
            'routed_last_native_repair': last,
            'native_input_qubits': 9, 'field_qubits': 2+2*(dimension.bit_length()-1),
            'peak_code_qubits_excluding_six_level_processor': 9+2+2*(dimension.bit_length()-1)+circuit.RECORD_BITS[law],
            'fixed_width_terminal_record_bits_including_two_status_bits': 2*(1+circuit.RECORD_BITS[law])+4+9+2,
            'native_register_initializations': 1, 'native_coin_initializations': 2,
            'native_repair_record_initializations': 2, 'native_classical_readout_registers': 5,
            'field_preparations': 1, 'field_pointer_initializations_and_reads': 4,
            'host_tree_storage_or_physical_durations_counted': False}


def run(steps=STEPS):
    rows = []
    for dimension in (2, 4):
        histories = {law: execute(law, dimension, steps) for law in ('C', 'E')}
        program = protocol.continuation_program(dimension, steps)
        ideal = bounds.program_bound(program)
        keys = (*bounds.required_service_keys(program.operations), 'native_preparation',
                'first_native_coin', 'first_native_repair', 'last_native_coin',
                'last_native_routing_and_repair', 'native_terminal_rotation_and_readout')
        rows.append({'dimension': dimension,
                     'law_audits': {law: audit(history) for law, history in histories.items()},
                     'coherence_witness': coherence_witness(histories),
                     'resources': {law: resources(law, dimension, steps) for law in histories},
                     'ideal_prepared_work_budget': ideal,
                     'physical_service_coverage': bounds.service_error_budget(
                         F(ideal['complete_ideal_product_to_infinite_instrument_error_upper']), required_stage_keys=keys)})
    return {'schema': 'oph.common_source.quantum_native.v1', 'steps': steps,
            'routing': routing.routing_word_summary(),
            'routing_execution': routing.routing_execution_diagnostics(), 'cutoffs': rows,
            'scope': {'full_native_input_dimension': 78, 'retained_native_output_dimension': 512,
                      'complete_native_reference_outputs_retained': True,
                      'cross_record_coherence_preserved': False,
                      'declared_native_integer_primitive_used_in_execution': False,
                      'native_terminal_read_is_full_computational_read_after_rotation': True,
                      'physical_services_or_clock_calibration_established': False,
                      'new_physical_law_or_empirical_data': False}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    files = parser.add_mutually_exclusive_group()
    files.add_argument('--write', type=Path)
    files.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.write:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(canonical(result))
    if args.check:
        compare_receipt(strict_json.loads(args.check.read_text(encoding='utf-8')), result)
    print(json.dumps({'complete_native_quantum_observer': True, 'cutoffs': [2, 4]}))


if __name__ == '__main__':
    main()
