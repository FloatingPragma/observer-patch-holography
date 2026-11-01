"""A complete four-mode observer of the two spectator-compatible lifts.

The finite native endpoint rule is compiled reversibly, retaining the coin,
every raw record and every first failed proper-code export. The phase feedback
is an explicit choice of coherent instrument, not a population-law inference.
"""
from __future__ import annotations

from fractions import Fraction
import math
from types import MappingProxyType

import numpy as np

from ..observer.source_scalar import compiler
from ..reference_recovery import recovery
from . import many_body


SOURCE_CEILING = Fraction(1, 10**8)
SEAM = (0, 2)
NATIVE_PORTS = (0, 1, 5, 6)
NATIVE_SEAM = (0, 5)


def _domain(law, character):
    if type(law) is not str or law not in ('C', 'E'):
        raise ValueError('a literal C/E law required')
    if type(character) is not int or character not in (-1, 1):
        raise ValueError('an exact positive or negative character required')
    return 7 if law == 'C' else 8


def _flip(controls, values, target):
    """Controlled X from a phase projector, including its scalar phase."""
    bits = len(controls)+1
    phases = np.zeros(1 << bits)
    address = sum(value << i for i, value in enumerate(values))+(1 << (bits-1))
    phases[address] = -math.pi
    diagonal = compiler.remap(compiler.compile_diagonal(phases), [*controls, target])
    return ([compiler.pulse('Y', math.pi/4, [target])]+diagonal
            +[compiler.pulse('Y', -math.pi/4, [target])])


def repair_word(law, character):
    _domain(law, character)
    word = [compiler.pulse('Y', -math.pi/4, [4])]
    word += _flip((0, 2, 4), (1, 0, 0), 5)
    word += _flip((0, 2, 4), (0, 1, 1), 6)
    if law == 'E':
        word += compiler.cnot(0, 7)+compiler.cnot(2, 7)
    if character == -1:
        # A moving fermion crosses exactly mode 1; mode 3 is outside the seam.
        phases = np.array([0., 0., 0., math.pi])
        for flag in (5, 6):
            word += compiler.remap(compiler.compile_diagonal(phases), [1, flag])
    for flag in (5, 6):
        word += compiler.cnot(flag, 0)+compiler.cnot(flag, 2)
    return word


def target_columns(law, character):
    """Independent endpoint arithmetic, not circuit-generated target data."""
    bits = _domain(law, character)
    target = np.zeros((1 << bits, 16), complex)
    for mask in range(16):
        left, right = mask & 1, (mask >> 2) & 1
        for coin in (0, 1):
            low, odd = divmod(left+right, 2)
            new_left = low+coin*odd
            current = left-new_left
            output = mask if current == 0 else mask ^ 5
            raw = {0: 0, 1: 1, -1: 2}[current]
            if law == 'E':
                raw += 4*odd
            sign = character if current and mask & 2 else 1
            target[output+(coin << 4)+(raw << 5), mask] = sign/math.sqrt(2)
    return target


def _check_embedding(law, character, target):
    """Bind the four-mode field chart to the actual twelve-port instrument."""
    bits = _domain(law, character)
    if (type(NATIVE_PORTS) is not tuple or len(NATIVE_PORTS) != 4
            or any(type(p) is not int or not 0 <= p < 12 for p in NATIVE_PORTS)
            or tuple(sorted(set(NATIVE_PORTS))) != NATIVE_PORTS
            or type(SEAM) is not tuple or SEAM != (0, 2)
            or any(type(p) is not int for p in SEAM)
            or type(NATIVE_SEAM) is not tuple
            or any(type(p) is not int for p in NATIVE_SEAM)
            or tuple(NATIVE_PORTS[p] for p in SEAM) != NATIVE_SEAM):
        raise ValueError('the ordered active-port chart and actual native seam must agree')
    expected = np.zeros((1 << bits, 16), complex)
    for mask in range(16):
        embedded = sum(((mask >> i) & 1) << p for i, p in enumerate(NATIVE_PORTS))
        for coin in (0, 1):
            raw, moved, phase, _ = many_body.column(embedded, law, NATIVE_SEAM, coin, character)
            packed = sum(((moved >> p) & 1) << i for i, p in enumerate(NATIVE_PORTS))
            if moved != sum(((packed >> i) & 1) << p for i, p in enumerate(NATIVE_PORTS)):
                raise ArithmeticError('native movement escaped the declared four-mode chart')
            expected[packed+(coin << 4)+(raw << 5), mask] = phase/math.sqrt(2)
    if not np.array_equal(_plain_complex(target, expected.shape), expected):
        raise ArithmeticError('four-mode target differs from the actual embedded native instrument')
    return 32


def _plain_complex(value, shape):
    if (type(value) is not np.ndarray or value.shape != shape
            or value.dtype.kind not in 'fc' or not np.isfinite(value).all()):
        raise ValueError('plain finite floating arrays on the entire output domain required')
    return np.array(value, dtype=np.complex128, copy=True)


def _archive(run, *, word, bits, initial):
    """Bind and replay every prefix once; aggregate masses cannot replace it."""
    if type(bits) is not int or not 1 <= bits <= 9:
        raise ValueError('one through nine exact work bits required')
    if (type(initial) is not np.ndarray or initial.ndim != 2
            or initial.shape[0] != 1 << bits or not 1 <= initial.shape[1] <= 16):
        raise ValueError('complete initial work register/reference columns required')
    initial = _plain_complex(initial, initial.shape)
    tape = compiler._validate(word, bits)
    expected = tuple((event['kind'], event['angle'], tuple(event['targets'])) for event in tape)
    if (type(run) is not recovery.CompleteRun or type(run.bits) is not int
            or run.bits != bits or type(run.word) is not tuple
            or any(type(event) is not tuple or len(event) != 3
                   or type(event[0]) is not str or type(event[1]) is not float
                   or not math.isfinite(event[1]) or type(event[2]) is not tuple
                   or any(type(p) is not int for p in event[2]) for event in run.word)
            or run.word != expected):
        raise ValueError('complete archive must retain the exact requested source word and work bits')
    archived_initial = _plain_complex(run.initial, initial.shape)
    if not np.array_equal(archived_initial, initial):
        raise ArithmeticError('complete archive changed its requested input')
    keys = {(kind, angle) for kind, angle, _ in expected}
    if (type(run.primitives) is not MappingProxyType or set(run.primitives) != keys
            or any(type(key) is not tuple or len(key) != 2
                   or type(key[0]) is not str or type(key[1]) is not float
                   for key in run.primitives)):
        raise ValueError('complete immutable six-level primitive catalogue required')
    primitives = {}
    for key in keys:
        actual = _plain_complex(run.primitives[key], (6, 6))
        fresh = _plain_complex(compiler.native_pulse(*key), (6, 6))
        if not np.array_equal(actual, fresh):
            raise ArithmeticError('archive primitive differs from the requested native source response')
        if np.linalg.norm(fresh.conj().T@fresh-np.eye(6), 2) > 2e-11:
            raise ArithmeticError('complete source primitive is not unitary')
        primitives[key] = fresh
    if (type(run.failure_weights) is not np.ndarray
            or run.failure_weights.shape != (len(expected), 4)
            or run.failure_weights.dtype.kind != 'f'
            or not np.isfinite(run.failure_weights).all()
            or np.any(run.failure_weights < 0)):
        raise ValueError('one complete four-slot failure row per source pulse required')
    weights = np.asarray(run.failure_weights, dtype=np.float64)
    gram = _plain_complex(run.failure_gram, (initial.shape[1], initial.shape[1]))
    success = _plain_complex(run.success, initial.shape)
    state = initial.copy()
    rebuilt_weights = np.zeros((len(expected), 4))
    rebuilt_gram = np.zeros_like(gram)
    for event_index, (kind, angle, targets) in enumerate(expected):
        full = primitives[kind, angle]
        failed = recovery._failure_vectors(state, full, targets, bits)
        for slot, vector in enumerate(failed):
            rebuilt_weights[event_index, slot] = np.vdot(vector, vector).real
            rebuilt_gram += vector.conj().T@vector
        dimension = 1 << len(targets)
        state = compiler._apply(state, full[:dimension, :dimension], targets, bits)
    # The same frozen input and native matrices define a deterministic linear
    # replay. Exact equality here introduces no probability floor: even tiny
    # retained complements and the zero padding in pair-code rows are bound.
    if (not np.array_equal(weights, rebuilt_weights)
            or not np.array_equal(gram, rebuilt_gram)
            or not np.array_equal(success, state)):
        raise ArithmeticError('complete archive outputs differ from its literal source prefixes')
    return success, gram, weights


def admit(run, target, *, word, bits, initial):
    """Complete Stinespring error upper bound, including all failure outputs."""
    actual, gram, weights = _archive(run, word=word, bits=bits, initial=initial)
    target = _plain_complex(target, actual.shape)
    failure = float(np.sum(weights, dtype=np.float64))
    # execute_complete owns the immutable prefix archive. Its complete input
    # Gram is checked as well as the amplitude distance to the independent map.
    roundoff = 64*np.finfo(float).eps*failure+64*float(np.spacing(failure))
    if (not math.isfinite(failure) or np.linalg.norm(gram-gram.conj().T) > roundoff
            or np.linalg.norm(gram) > failure+roundoff
            or abs(float(np.trace(gram).real)-failure) > roundoff):
        raise ArithmeticError('inconsistent complete failure archive')
    error = float(np.linalg.norm(actual-target))+math.sqrt(failure)
    if not math.isfinite(error) or error > float(SOURCE_CEILING):
        raise ArithmeticError('complete source output exceeds its admission ceiling')
    return error


def _resources(word):
    values = compiler.resources(word)
    # Floating processor diagnostics are freshly checked, not retained as
    # historical interval bounds. Counts and the actual word remain exact.
    return {key: values[key] for key in (
        'native_pulses', 'single_code_pulses', 'pair_code_pulses',
        'identity_phase_pulses', 'proper_code_load_events',
        'proper_code_export_events', 'proper_code_service_events')} | {
        'first_failure_Kraus_slots': sum(2 if e['kind'] == 'J' else 4 for e in word)}


def _beam():
    # The two commuting Pauli products rotate only the one-occupation sector.
    word = compiler.pauli_rotation(['X', 'Y'], -math.pi/8)
    word += compiler.pauli_rotation(['Y', 'X'], math.pi/8)
    return compiler.remap(word, [1, 3])


def experiment(law, character, *, which_path=False):
    if type(which_path) is not bool:
        raise ValueError('an explicit which-path switch required')
    work_bits = _domain(law, character)
    bits = work_bits+int(which_path)
    # Fixed total number two: one endpoint particle and one coherent spectator.
    # Apparatus control pulses and blank preparation are explicit resources.
    seed = _flip((), (), 0)+_flip((), (), 1)
    prepare = seed+_beam()
    record = compiler.cnot(1, work_bits) if which_path else []
    word = prepare+record+repair_word(law, character)+compiler.inverse(_beam())
    initial = np.zeros((1 << bits, 1), complex)
    initial[0, 0] = 1
    run = recovery.execute_complete(word, bits, initial)
    # A source-independent ideal encoder has columns |3> and |9>. The
    # emitted beam's sign convention is fixed analytically by its generator.
    state = np.zeros(16, complex)
    state[3] = state[9] = 1/math.sqrt(2)
    target = target_columns(law, character)
    _check_embedding(law, character, target)
    repaired = target@state
    ideal = np.zeros(1 << bits, complex)
    for address, amplitude in enumerate(repaired):
        output = address | (((address >> 1) & 1) << work_bits) if which_path else address
        ideal[output] = amplitude
    # The inverse beam maps the symmetric spectator to mode1 and the
    # antisymmetric spectator to mode3, with a harmless overall minus sign.
    expected = np.zeros_like(ideal)
    for base in range(1 << bits):
        if base & 10:
            continue
        a, b = ideal[base | 2], ideal[base | 8]
        expected[base | 2] = (a+b)/math.sqrt(2)
        expected[base | 8] = (b-a)/math.sqrt(2)
    admit(run, expected[:, None], word=word, bits=bits, initial=initial)
    probabilities = np.sum(abs(run.success)**2, axis=1)
    # Every pointer address, including forbidden addresses, participates in
    # this gate. No coarse marginal or normalization can conceal leakage.
    amplitude_lower = np.linalg.norm(np.sqrt(probabilities)-abs(expected))
    if amplitude_lower > float(SOURCE_CEILING):
        raise ArithmeticError('complete terminal law exceeds source amplitude budget')
    nonzero = {str(i): '1/8' if which_path else '1/2'
               for i in np.flatnonzero(expected != 0)}
    if len(nonzero) != (8 if which_path else 2):
        raise ArithmeticError('unexpected ideal terminal support')
    return {'law': law, 'character': character, 'total_input_occupation': 2,
            'active_native_ports_in_mode_order': list(NATIVE_PORTS),
            'native_seam': list(NATIVE_SEAM),
            'field_modes': 4, 'work_qubits': bits, 'blank_work_bits': bits,
            'retained_which_path_record': which_path,
            'complete_terminal_addresses': 1 << bits,
            'ideal_nonzero_terminal_probabilities': nonzero,
            'every_other_ideal_terminal_probability': '0',
            'ideal_proper_code_failure_probability': '0',
            'complete_output_amplitude_error_ceiling': str(SOURCE_CEILING),
            'physical_initialization_storage_and_readout_cost_calibrated': False,
            'resources': _resources(word)}


def certificate():
    audits, experiments = [], []
    for law in ('C', 'E'):
        bits = _domain(law, 1)
        for character in (1, -1):
            word = repair_word(law, character)
            initial = np.eye(1 << bits, dtype=complex)[:, :16]
            target = target_columns(law, character)
            embedded = _check_embedding(law, character, target)
            run = recovery.execute_complete(word, bits, initial)
            admit(run, target, word=word, bits=bits, initial=initial)
            audits.append({'law': law, 'character': character,
                           'all_field_input_columns': 16,
                           'actual_native_embedding_columns_checked': embedded,
                           'complete_terminal_addresses': 1 << bits,
                           'raw_slots_per_coin': 4 if law == 'C' else 8,
                           'complete_output_amplitude_error_ceiling': str(SOURCE_CEILING),
                           'resources': _resources(word)})
            experiments.append(experiment(law, character))
            experiments.append(experiment(law, character, which_path=True))
    return {'all_column_audits': audits, 'fixed_number_phase_readouts': experiments,
            'numerical_ceiling_is_interval_or_hardware_certificate': False,
            'native_12_port_full_register_source_execution_claimed': False,
            'source_primitive': 'existing six-level proper-code compiler',
            'phase_feedback': 'explicit coherent positive or CAR-negative lift',
            'all_failure_outputs': 'immutable first-failure prefix archive from execute_complete'}
