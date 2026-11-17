"""Source-pulse Ramsey control of the finite native density walk.

The field is the existing 720-arc representation with five clean scratch
bits, plus one Ramsey bit. Raw source routines are reused on the full
register; no logical target matrix substitutes for an executed routine.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from model.observer import coherent_continuation as wave
from model.observer.source_scalar import compiler as c
from . import clock_services, continuation, joint_clock
from .joins_contract import source_bindings, check_receipt
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
CAMPAIGN = HERE / 'ramsey_campaign.json'
OUTPUT = HERE / 'ramsey_service_results.json'
QUBITS = 16
DIMENSION = 1 << QUBITS
ARM = 15
GUARDS = {'prepare': F(1, 10**8), 'controlled_reflection': F(1, 10**8),
          'ccx': F(1, 10**10), 'flip': F(1, 10**10),
          'ramsey': F(1, 10**10), 'site_prepare': F(1, 10**9)}


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                     allow_nan=False).encode()).hexdigest()


def _quantize(tape):
    return [c.pulse(p['kind'], round(p['angle']*2**44)/2**44, p['targets']) for p in tape]


def port_probe():
    """The normalized centered port-zero vector, embedded into 16 addresses."""
    result = np.zeros(16)
    result[:12] = -1/math.sqrt(132)
    result[0] = 11/math.sqrt(132)
    return result


def _site_preparation():
    """The inherited Gray elimination, with exact integer input ratios.

    Fixed high precision before dyadic rounding keeps the retained word
    independent of platform libm hypot/atan2 last-bit differences.
    """
    import mpmath as mp
    gray = [j ^ (j >> 1) for j in range(16)]
    integers = [11]+[-1]*11+[0]*4
    rotations = []
    with mp.workdps(80):
        residual = [mp.mpf(integers[j]) for j in gray]
        for row in range(15, 0, -1):
            x, y = residual[row-1], residual[row]
            theta = mp.mpf(0) if x == y == 0 else mp.atan2(y, x)
            rotations.append((gray[row-1], gray[row], float(theta)))
            residual[row-1], residual[row] = mp.sqrt(x*x+y*y), mp.mpf(0)
    return _quantize([pulse for a, b, angle in reversed(rotations)
                      for pulse in c._gray_rotation(a, b, -angle, 4)])


@lru_cache(maxsize=1)
def _program():
    base = wave.density_walk_program()
    prepare = base['routines']['prepare']
    phase = c.compile_diagonal([0.]*64+[0.]+[math.pi]*63)
    reflection = c.inverse(prepare)+phase+prepare
    hadamard = np.array([[1., 1.], [1., -1.]])/math.sqrt(2)
    routines = {name: tape for name, tape in base['routines'].items() if name != 'reflection'}
    routines.update(controlled_reflection=_quantize(reflection),
                    ramsey=_quantize(c.compile_real_orthogonal(hadamard)),
                    site_prepare=_site_preparation())
    shift = []
    # At each endpoint-target gate, scratch13 is blank. Computing its two
    # controls, acting only on arm1, then uncomputing adds one control without
    # adding scratch. Event/zero-bit computations cancel on the reference arm.
    for op in base['shift']:
        bits = op['targets']
        if op['routine'] == 'ccx' and bits[2] < 4:
            shift.extend({'routine': 'ccx', 'targets': list(t)} for t in
                         ((bits[0], bits[1], 13), (13, ARM, bits[2]),
                          (bits[0], bits[1], 13)))
        else:
            shift.append(op)
    return {'routines': routines, 'shift': shift,
            'walk': [{'routine': 'controlled_reflection', 'targets': [4, 5, 6, 7, 8, 9, ARM]}]+shift,
            'address_qubits': 10, 'clean_scratch_qubits': 5, 'Ramsey_qubits': 1}


def program():
    # Public mutable schedules cannot poison a later verification.
    return json.loads(json.dumps(_program()))


def expanded_word(schedule):
    p = _program()
    return [pulse for op in schedule for pulse in c.remap(p['routines'][op['routine']], op['targets'])]


def shift_permutation():
    """Exact execution on every full-register basis state, including padding."""
    value = np.arange(DIMENSION, dtype=np.int64)
    for op in _program()['shift']:
        bits = op['targets']
        if op['routine'] == 'flip':
            value ^= 1 << bits[0]
        else:
            value ^= (((value >> bits[0]) & 1) & ((value >> bits[1]) & 1)) << bits[2]
    expected = np.arange(DIMENSION, dtype=np.int64)
    for address in range(1024):
        site, event = address & 15, address >> 4
        if site < 12 and event < 60 and event & 1:
            a, b = wave.quantum_mixing.seams()[event//2]
            if site in (a, b):
                expected[(1 << ARM)+address] ^= a ^ b
    clean = np.r_[np.arange(1024), (1 << ARM)+np.arange(1024)]
    if (len(np.unique(value)) != DIMENSION or not np.array_equal(value[clean], expected[clean])
            or not np.array_equal(value[:1 << ARM], np.arange(1 << ARM))):
        raise ArithmeticError('controlled shift loses its reference arm, clean scratch, or padding')
    return value


def _targets():
    uniform = np.r_[np.full(60, 1/math.sqrt(60)), np.zeros(4)]
    reflection = np.eye(128)
    reflection[64:, 64:] = 2*np.outer(uniform, uniform)-np.eye(64)
    ccx = np.zeros((8, 8))
    for i in range(8):
        ccx[i ^ (4 if i & 3 == 3 else 0), i] = 1
    return {'prepare': uniform, 'controlled_reflection': reflection, 'ccx': ccx,
            'flip': np.array([[0., 1.], [1., 0.]]),
            'ramsey': np.array([[1., 1.], [1., -1.]])/math.sqrt(2),
            'site_prepare': port_probe()}


@lru_cache(maxsize=1)
def _raw():
    dims = {'prepare': 6, 'controlled_reflection': 7, 'ccx': 3, 'flip': 1,
            'ramsey': 1, 'site_prepare': 4}
    result = {}
    for name, tape in _program()['routines'].items():
        matrix = c.execute(tape, dims[name])
        result[name] = matrix
    certify_routines(result)
    return {name: (matrix.dtype.str, matrix.tobytes(), matrix.shape)
            for name, matrix in result.items()}


def certify_routines(arrays):
    """All-column guards for operations, first-column guards for blank inputs."""
    targets = _targets()
    if set(arrays) != set(targets):
        raise ValueError('every source routine required')
    for name, target in targets.items():
        matrix = arrays[name]
        size = len(target)
        if (not isinstance(matrix, np.ndarray) or matrix.shape != (size, size)
                or not np.isfinite(matrix).all()):
            raise ValueError('original finite square routine arrays required')
        residual = np.linalg.norm(matrix[:, 0]-target) if target.ndim == 1 else np.linalg.norm(matrix-target, 'fro')
        if (residual > float(GUARDS[name]) or
                np.linalg.norm(matrix.conj().T@matrix-np.eye(len(matrix)), 'fro') > 2*float(GUARDS[name])):
            raise ArithmeticError(f'raw source {name} exceeds its numerical guard')
    return True


def raw_routines():
    return {name: np.frombuffer(payload, dtype=np.dtype(dtype)).reshape(shape).copy()
            for name, (dtype, payload, shape) in _raw().items()}


def execute(schedule, columns, *, routines=None):
    """Full-register raw-array action with arbitrary spectator columns."""
    value = np.asarray(columns, dtype=complex)
    if value.ndim != 2 or value.shape[0] != DIMENSION or not np.isfinite(value).all():
        raise ValueError('finite full 16-qubit reference columns required')
    value = value.copy()
    arrays = raw_routines() if routines is None else routines
    for op in schedule:
        value = c._apply(value, arrays[op['routine']], list(reversed(op['targets'])), QUBITS)
    return value


def budget(ticks, *, fixed_probe=True, read=True):
    if type(ticks) is not int or not 0 <= ticks <= 4:
        raise ValueError('declared integer tick count 0 through 4 required')
    if type(fixed_probe) is not bool or type(read) is not bool:
        raise ValueError('explicit preparation/read booleans required')
    uses = Counter(op['routine'] for op in _program()['walk'])
    uses = Counter({name: count*ticks for name, count in uses.items()})
    uses.update(prepare=1, ramsey=1+int(read), site_prepare=int(fixed_probe))
    total = sum((GUARDS[name]*count for name, count in uses.items()), F())
    if total >= 1:
        raise ArithmeticError('source product guard is not perturbative')
    delta = total/(1-total)
    pulses = sum(len(_program()['routines'][name])*count for name, count in uses.items())
    return {'ticks': ticks, 'routine_uses': dict(sorted(uses.items())), 'source_pulses': pulses,
            'proper_code_load_events': pulses, 'proper_code_export_events': pulses,
            'source_operator_error_upper': str(delta),
            'complete_halftrace_error_upper': str(delta*(2+delta)/2),
            'source_channel_norm_upper': str((1+delta)**2),
            'guard_sum': str(total), 'guards_are_interval_or_physical_certificates': False}


def _field_columns(ticks, *, read=True, sites=None):
    """Complete field output; scratch and address padding remain present.

    sites=None executes the fixed source-prepared port probe. Otherwise the
    supplied 16-row site columns retain arbitrary reference amplitudes.
    Event, scratch and arm inputs are blank. A read splits the arm into two
    complete outputs; an abort retains the unmeasured full Ramsey register.
    """
    budget(ticks, fixed_probe=sites is None, read=read)
    schedule = []
    if sites is None:
        initial = np.zeros((DIMENSION, 1), complex)
        initial[0, 0] = 1
        schedule.append({'routine': 'site_prepare', 'targets': [0, 1, 2, 3]})
    else:
        sites = np.asarray(sites, dtype=complex)
        if sites.ndim != 2 or sites.shape[0] != 16 or not np.isfinite(sites).all():
            raise ValueError('finite 16-row site/reference columns required')
        initial = np.zeros((DIMENSION, sites.shape[1]), complex)
        initial[:16] = sites
    schedule.extend(({'routine': 'prepare', 'targets': [4, 5, 6, 7, 8, 9]},
                     {'routine': 'ramsey', 'targets': [ARM]}))
    schedule.extend(_program()['walk']*ticks)
    if read:
        schedule.append({'routine': 'ramsey', 'targets': [ARM]})
    output = execute(schedule, initial)
    return (output[:1 << ARM], output[1 << ARM:]) if read else output


@lru_cache(maxsize=8)
def _fixed_field(ticks, read):
    return _field_columns(ticks, read=read)


def field_columns(ticks, *, read=True, sites=None):
    """Fresh complete outputs, with cached fixed-probe source execution."""
    budget(ticks, fixed_probe=sites is None, read=read)
    if sites is not None:
        return _field_columns(ticks, read=read, sites=sites)
    value = _fixed_field(ticks, read)
    return tuple(v.copy() for v in value) if read else value.copy()


def ideal_field_columns(ticks, *, read=True, sites=None):
    """Independent 720-arc reflection/permutation oracle on clean scratch."""
    site = port_probe()[:, None] if sites is None else np.asarray(sites, dtype=complex)
    initial = np.zeros((1 << 15, site.shape[1]), complex)
    for event in range(60):
        initial[16*event:16*(event+1)] = site/math.sqrt(60)
    live = initial.copy()
    reversal = np.arange(1 << 15)
    for event in range(1, 60, 2):
        a, b = wave.quantum_mixing.seams()[event//2]
        reversal[16*event+a], reversal[16*event+b] = 16*event+b, 16*event+a
    uniform = np.r_[np.full(60, 1/math.sqrt(60)), np.zeros(4)]
    reflection = 2*np.outer(uniform, uniform)-np.eye(64)
    for _ in range(ticks):
        reflected = c._apply(live, reflection, list(range(9, 3, -1)), 15)
        live = np.empty_like(reflected)
        live[reversal] = reflected
    return ((initial+live)/2, (initial-live)/2) if read else np.vstack((initial, live))/math.sqrt(2)


def _gram(columns):
    return columns.conj().T@columns


def joint_reference(frame, source_columns, field_columns):
    """Reduced untouched-reference output from explicit source/field factors.

    Input sum_k C[a,r,k] tensor f_k; a source abort traces padding rather
    than deleting the field or its reference. This formula also retains
    arbitrary off-diagonal reference entries and unnormalized masses.
    """
    value = np.asarray(source_columns, dtype=complex)
    fields = np.asarray(field_columns, dtype=complex)
    if (value.ndim != 3 or value.shape[0] != 66 or fields.ndim != 2
            or value.shape[2] != fields.shape[1] or not np.isfinite(value).all()
            or not np.isfinite(fields).all()):
        raise ValueError('compatible finite source/reference/field Schmidt columns required')
    source = np.einsum('ab,bri->ari', frame.matrix, value)
    return np.einsum('ari,asj,ji->rs', source, source.conj(), _gram(fields))


def processor_export(code_dimension, columns):
    """Every proper-code success/complement output, including references."""
    value = np.asarray(columns, dtype=complex)
    if value.ndim != 2 or value.shape[0] != 6 or not np.isfinite(value).all():
        raise ValueError('six-level processor/reference columns required')
    return tuple(k@value for k in c.transfer_kraus(code_dimension, coherent=True))


@lru_cache(maxsize=1)
def _clock_budget():
    # Reuse the source-bound predecessor certificate; the new packet freshly
    # executes the actual clock frames but does not recertify every unchanged
    # native/probe matrix twice. Its strict source closure includes the proof,
    # independent tests and retained source-word certificate dependencies.
    from tools import strict_json
    receipt = strict_json.load(HERE/'joint_clock_results.json')
    for path, digest in bound_pins(receipt['source_sha256']).items():
        if hashlib.sha256((ROOT/path).read_bytes()).hexdigest() != digest:
            raise ValueError('predecessor joint-clock source binding changed')
    return receipt['source_and_service_budget']


def complete_joint(law):
    """Replace every admitted wave/Ramsey factor in the complete joint law."""
    fields = {n: field_columns(n) for n in (1, 2)}
    waits = {n: field_columns(n, read=False) for n in (0, 1, 2)}
    reference = np.diag(list(map(float, joint_clock.clock.transport.reference()))).astype(complex)
    frames = joint_clock.clock_frames(law, raw=True)
    actual, ideal = {}, joint_clock.factor_probabilities(frames, reference)
    for frame in frames:
        mass = float(np.einsum('ai,ai->', frame.matrix@reference, frame.matrix.conj()).real)
        if frame.abort:
            actual[frame.key] = mass*float(np.linalg.norm(waits[frame.ticks])**2)
        else:
            for sign, field in zip((1, -1), fields[frame.ticks], strict=True):
                actual[(*frame.key, sign)] = mass*float(np.linalg.norm(field)**2)
    error = sum(abs(actual[key]-value) for key, value in ideal.items())/2
    bound = max(F(budget(n)['complete_halftrace_error_upper']) for n in (1, 2))
    source_gain = F(_clock_budget()['source_channel_norm_upper'])
    if error > float(source_gain*bound):
        raise ArithmeticError('compiled field joint law exceeds its complete product bound')
    if error > 1e-6:
        raise ArithmeticError('compiled joint-law numerical diagnostic failed')
    preparation = joint_clock.clock.preparation.population_maps(law=law)
    success = np.diag(preparation['success'][:, 0]).astype(complex)
    failures = {flag: float(np.sum(matrix[:, 0])) for flag, matrix in preparation.items() if flag != 'success'}
    finite_mass = sum(failures.values())
    for frame in frames:
        mass = float(np.einsum('ai,ai->', frame.matrix@success, frame.matrix.conj()).real)
        field_mass = (np.linalg.norm(waits[frame.ticks])**2 if frame.abort
                      else sum(np.linalg.norm(v)**2 for v in fields[frame.ticks]))
        finite_mass += mass*field_mass
    if abs(finite_mass-1) > float(source_gain*bound)+2e-9:
        raise ArithmeticError('field compilation dropped a preparation or export flag')
    return {'law': law, 'complete_read_and_abort_records': len(actual),
            'source_export_aborts_retained': sum(frame.abort for frame in frames),
            'wave_halftrace_bound': str(source_gain*bound),
            'compiled_vs_admitted_wave_TV_below': '1/1000000',
            'complete_clock_record_digest': _digest(sorted(map(repr, actual))),
            'all_predecessor_record_keys_preserved': set(actual) == set(ideal),
            'finite_native_preparation_flags': sorted(preparation),
            'failed_preparations_bypass_field_services': True,
            'finite_preparation_complete_mass_within_bound': True}


def continuing_history(law):
    """Four actual native repairs, interleaved clocks and live wave ticks.

    This retained path is unnormalized. The enclosing complete instrument
    includes every branch from continuation.branches and one explicit source
    export failure per requested clock. Field evolution is never restarted.
    """
    records = clock_services.exported_history(law)['complete_native_tape']
    params = clock_services.FamilyParameters(c_weight=F(int(law == 'C')))
    codes = np.array(clock_services.full_flow.support()[0])
    outside = np.array([i for i in range(512) if i not in set(codes)])
    source = clock_services._probe()
    ideal = source.copy()
    memory, core, exports, weight, cores, aborts = (), (0, 1), [], F(1), [], []
    for position, fields in enumerate(records):
        record = continuation.observer.Record(**fields)
        branch = next(b for b in continuation.branches(memory, core, params) if b.record == record)
        source = clock_services._native_vector(source, branch)
        ideal = clock_services._native_vector(ideal, branch)
        word = clock_services.matrix(branch.core)
        escaped = word[np.ix_(outside, codes)]@source
        source = word[np.ix_(codes, codes)]@source
        phase = np.exp(-4j*np.log(np.array(list(map(float, continuation.reference(branch.core))))*60))
        ideal = phase[:, None]*ideal
        # The clock abort occurs before this attempt's live wave increment.
        wave_wait = field_columns(position, read=False)
        aborts.append({'position': position, 'source_padding_columns': escaped,
                       'unmeasured_wave_columns': wave_wait})
        weight *= branch.weight
        memory, core = branch.memory, branch.core
        if branch.exported is not None:
            exports.append(branch.exported)
        cores.append(core)
    fields = field_columns(4)
    targets = ideal_field_columns(4)
    scale = math.sqrt(float(weight))
    source_error = np.linalg.norm(source-ideal)/scale
    wave_error = math.sqrt(sum(np.linalg.norm(a-b)**2 for a, b in zip(fields, targets, strict=True)))
    source_bound = (1+clock_services.OPERATOR_GATE)**4-1
    if (source_error > float(source_bound) or wave_error > float(F(budget(4)['source_operator_error_upper']))
            or core in continuation.compiled_cores() or len(exports) != 3):
        raise ArithmeticError('continuing clock/field path failed after FIFO export')
    return {'law': law, 'native_tape': records, 'transported_cores': [list(x) for x in cores],
            'native_path_weight': str(weight), 'retained_records': len(memory), 'exported_records': len(exports),
            'clock_requests': 4, 'wave_ticks': 4, 'preparations_during_continuation': 0,
            'clock_export_abort_slots': len(aborts),
            'abort_wave_ticks': [row['position'] for row in aborts],
            'final_clock_core_outside_original_four': True,
            'clock_native_pulses': 4*2042, 'field_budget': budget(4),
            'complete_selected_history_outputs_within_product_bound': True,
            'selected_path_is_complete_history_distribution': False}


def run():
    permutation = shift_permutation()
    arrays = raw_routines()
    cases = []
    for n in (0, 1, 2, 4):
        actual, target = field_columns(n), ideal_field_columns(n)
        difference = math.sqrt(sum(np.linalg.norm(a-b)**2 for a, b in zip(actual, target, strict=True)))
        bound = F(budget(n)['source_operator_error_upper'])
        if difference > float(bound):
            raise ArithmeticError('complete raw Ramsey output exceeds its operator product bound')
        if n in (1, 2):
            expected = joint_clock.field_probabilities(n)
            if max(abs(np.linalg.norm(a)**2-float(p)) for a, p in zip(actual, expected, strict=True)) > 1e-6:
                raise ArithmeticError('density representation lost the original native field probability')
        cases.append({'ticks': n, 'budget': budget(n), 'complete_raw_output_error_within_bound': True})
    word = expanded_word(_program()['walk'])
    paths = set(wave.SOURCE_PATHS) | set(joint_clock.clock.transport.source_bindings())
    paths.update(bound_pins(json.loads((HERE/'clock_services_results.json').read_text())['source_sha256']))
    paths.update({'model/clock_defect_energy/ramsey_campaign.json',
                  'model/clock_defect_energy/ramsey_service.py',
                  'model/clock_defect_energy/test_ramsey_service.py',
                  'model/clock_defect_energy/RAMSEY_SERVICE.md',
                  'model/clock_defect_energy/joint_clock.py',
                  'model/clock_defect_energy/joint_clock_results.json',
                  'model/clock_defect_energy/clock_services_results.json',
                  'Lean/Dynamics/ChoiCPTP.lean', 'Lean/Time/ClockReadout.lean'})
    campaign = json.loads(CAMPAIGN.read_text())
    return {'schema': 'oph.clock_defect_energy.ramsey_service.v1',
            'campaign_id': campaign['campaign_id'],
            'contract_sha256': hashlib.sha256(CAMPAIGN.read_bytes()).hexdigest(),
            'predecessor_git_commit': campaign['predecessor_git_commit'],
            'controlled_walk': {'valid_field_arcs': 720, 'full_register_basis_states': len(permutation),
                'clean_addresses_both_arms_checked': 2048, 'Ramsey_reference_arm_identity_all_dirty_inputs': True,
                'native_pulses': len(word), 'expanded_word_sha256': _digest(word),
                'routine_uses': dict(Counter(op['routine'] for op in _program()['walk'])),
                'source_routine_dimensions': {name: len(a) for name, a in arrays.items()},
                'all_raw_routine_columns_executed': True, 'target_arrays_used_as_execution': False,
                'full_register_execution': 'factorized reuse of raw source arrays, with complete scratch and padding',
                'controlled_walk_dirty_scratch_extension_claimed_equal_to_predecessor': False},
            'resource_accounting': {
                'field_qubits_including_Ramsey_and_scratch': 16,
                'clock_address_qubits': 9, 'clock_plus_field_register_qubits': 25,
                'shared_source_processor_levels': 6,
                'compact_program_utf8_bytes': len(json.dumps(_program(), sort_keys=True, separators=(',', ':')).encode()),
                'expanded_walk_utf8_bytes': len(json.dumps(word, sort_keys=True, separators=(',', ':')).encode()),
                'expanded_walk_position_counter_bits': len(word).bit_length(),
                'core_plus_repair_policy_states': 132,
                'source_word_interpreter_memory_counted_separately': True,
                'host_workspace_or_physical_storage_derived': False},
            'proper_code_failure_contract': {
                'per_pulse_success_and_every_processor_complement_retained': True,
                'exact_valid_code_complement_probability': '0',
                'reason': 'invariant proper-code source generators, checked by the predecessor compiler',
                'nonzero_invalid_processor_export_tested': True,
                'address_padding_is_not_processor_failure': True,
                'physical_pulse_error_formula': 'P*epsilon/(1-P*epsilon) when P*epsilon<1; transfer errors additional',
                'coherent_receiver_overwrite_and_sender_reset_are_supplied': True},
            'fixed_field_cases': cases, 'joint_laws': [complete_joint(law) for law in ('C', 'E')],
            'continuing_histories': [continuing_history(law) for law in ('C', 'E')],
            'predecessor_clock_bound_reuse': 'hash-verified joint_clock receipt; actual clock arrays execute, inherited exact native/probe certificates are not re-proved',
            'density_sector_intertwiner': wave.density_sector_audit(),
            'physical_seconds_calibrated': False, 'full_3960_arc_unitary_compiled': False,
            'integrated_cp02_accepted': False, 'source_sha256': bound_pins(source_bindings(paths))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', type=Path, nargs='?', const=OUTPUT)
    parser.add_argument('--check', type=Path, nargs='?', const=OUTPUT)
    args = parser.parse_args()
    result = run()
    if args.write:
        args.write.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    if args.check:
        check_receipt(args.check, result)
    if not args.write and not args.check:
        print(json.dumps(result, indent=2, sort_keys=True))
    else:
        print('Source Ramsey controller, complete joint outputs and finite density map verified.')


if __name__ == '__main__':
    main()
