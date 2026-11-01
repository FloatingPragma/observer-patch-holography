"""Complete source probes of the full-occupation noise-algebra boundary.

Twelve native occupation bits and every coin/raw record are explicit.  The
extreme-sector cat requires a supplied cross-number phase apparatus and is
outside the fixed-global-number V2 observer class.  Its fixed-number control
uses six supplied particles and number-preserving preparation macros.

The bounded executor retains the same lazy, literal first-export-failure
archive as reference_recovery, without enlarging that predecessor's domain.
"""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
import math

import numpy as np

from ..exchange_sector import many_body
from ..observer.source_scalar import compiler
from ..reference_observation.observation import frozen
from ..reference_recovery import recovery


MODES = 12
FULL = (1 << MODES)-1
SEAM = (0, 5)
COIN, PLUS, MINUS, PARITY = 12, 13, 14, 15
PULSE_CEILING = Fraction(1, 10**11)
IDEAL_WORD_CEILING = Fraction(1, 10**11)
MAX_BITS, MAX_PULSES = 17, 1024
PAIRS = ((1, 7), (2, 8), (3, 9), (4, 10), (6, 11))
SIX_LEFT = sum(1 << j for j in (0, 1, 2, 3, 4, 6))
SIX_RIGHT = SIX_LEFT ^ (1 << 0) ^ (1 << 5)


@dataclass(frozen=True, slots=True)
class _FrozenPrimitives(Mapping):
    """Tuple-owned catalogue: no caller-held dictionary can change a prefix."""
    entries: tuple

    def __getitem__(self, key):
        for candidate, value in self.entries:
            if candidate == key:
                return value
        raise KeyError(key)

    def __iter__(self):
        return (key for key, _ in self.entries)

    def __len__(self):
        return len(self.entries)


def _immutable_array(value):
    """Accept the bytes-backed snapshots emitted by frozen(), not flags alone."""
    current = value
    if type(current) is not np.ndarray or current.dtype.kind not in 'fc':
        raise ValueError('immutable bytes-backed source archive arrays required')
    while type(current) is np.ndarray:
        if current.flags.writeable:
            raise ValueError('immutable bytes-backed source archive arrays required')
        current = current.base
    if type(current) is not bytes:
        raise ValueError('immutable bytes-backed source archive arrays required')


def _domain(law, character):
    if type(law) is not str or law not in ('C', 'E'):
        raise ValueError('literal C or E law required')
    if type(character) is not int or character not in (-1, 1):
        raise ValueError('exact positive or negative character required')
    return 15 if law == 'C' else 16


def _flip(controls, values, target):
    bits = len(controls)+1
    phases = np.zeros(1 << bits)
    phases[sum(v << j for j, v in enumerate(values))+(1 << (bits-1))] = -math.pi
    return ([compiler.pulse('Y', math.pi/4, [target])]
            +compiler.remap(compiler.compile_diagonal(phases), [*controls, target])
            +[compiler.pulse('Y', -math.pi/4, [target])])


def repair_word(law, character):
    """Literal complete twelve-mode endpoint rule, including spectator signs."""
    _domain(law, character)
    word = [compiler.pulse('Y', -math.pi/4, [COIN])]
    word += _flip((0, 5, COIN), (1, 0, 0), PLUS)
    word += _flip((0, 5, COIN), (0, 1, 1), MINUS)
    if law == 'E':
        word += compiler.cnot(0, PARITY)+compiler.cnot(5, PARITY)
    if character == -1:
        for spectator in range(1, 5):
            for flag in (PLUS, MINUS):
                word += compiler.remap(compiler.compile_diagonal(np.array([0., 0., 0., math.pi])),
                                       [spectator, flag])
    for flag in (PLUS, MINUS):
        word += compiler.cnot(flag, 0)+compiler.cnot(flag, 5)
    return word


def _ideal_pulse(event):
    kind, angle = event['kind'], event['angle']
    if kind == 'J':
        block = np.diag([np.exp(1j*angle)]*3+[1.])
    else:
        x = np.array([[0., 1.], [1., 0.]], complex)
        y = np.array([[0., -1j], [1j, 0.]])
        z = np.diag([1., -1.])
        p = {'X': x, 'Y': y, 'Z': z, 'I': np.eye(2), 'H': (x+z)/math.sqrt(2)}[kind]
        block = math.cos(angle)*np.eye(2)+1j*math.sin(angle)*p
    out = np.zeros((6, len(block)), complex)
    out[:len(block)] = block
    return out


def pulse_budget(word):
    """Complete primitive admission and arbitrary-input telescoping bound."""
    tape = compiler._validate(word, MAX_BITS)
    if len(tape) > MAX_PULSES:
        raise ValueError('source word exceeds bounded pulse domain')
    maximum = 0.
    for kind, angle in {(e['kind'], e['angle']) for e in tape}:
        actual = compiler.native_pulse(kind, angle)
        target = _ideal_pulse(dict(kind=kind, angle=angle))
        if (type(actual) is not np.ndarray or actual.shape != (6, 6)
                or actual.dtype.kind not in 'fc' or not np.isfinite(actual).all()):
            raise ArithmeticError('finite complete six-level primitive required')
        if np.linalg.norm(actual.conj().T@actual-np.eye(6), 2) > 2e-11:
            raise ArithmeticError('complete six-level primitive is not unitary')
        error = float(np.linalg.norm(actual[:, :target.shape[1]]-target, 2))
        if not math.isfinite(error) or error > float(PULSE_CEILING):
            raise ArithmeticError('primitive complete-output admission failed')
        maximum = max(maximum, error)
    total = len(tape)*PULSE_CEILING
    # Full success and complement operator norms are each bounded by B.
    # Their orthogonal direct sum has norm at most 2B.
    return 2*total/(1-total), maximum


def _input(word, bits, initial):
    if type(bits) is not int or not 1 <= bits <= MAX_BITS:
        raise ValueError('one through seventeen bounded work bits required')
    if (type(initial) is not np.ndarray or initial.ndim != 2
            or initial.shape[0] != 1 << bits or not 1 <= initial.shape[1] <= 2
            or initial.dtype.kind not in 'fc' or not np.isfinite(initial).all()):
        raise ValueError('complete finite work register with one or two reference columns required')
    tape = compiler._validate(word, bits)
    if len(tape) > MAX_PULSES:
        raise ValueError('source word exceeds bounded pulse domain')
    return tape, np.array(initial, dtype=complex, copy=True)


def execute_complete(word, bits, initial):
    """Stream all first failures, retaining immutable reconstructible outputs."""
    tape, initial = _input(word, bits, initial)
    pulse_budget(tape)
    word = tuple((e['kind'], e['angle'], tuple(e['targets'])) for e in tape)
    primitives = {key: frozen(compiler.native_pulse(*key)) for key in
                  {(kind, angle) for kind, angle, _ in word}}
    state = initial.copy()
    weights = np.zeros((len(word), 4))
    gram = np.zeros((initial.shape[1],)*2, complex)
    for i, (kind, angle, targets) in enumerate(word):
        full = primitives[kind, angle]
        for slot, branch in enumerate(recovery._failure_vectors(state, full, targets, bits)):
            gram += branch.conj().T@branch
            weights[i, slot] = np.vdot(branch, branch).real
        d = 1 << len(targets)
        state = compiler._apply(state, full[:d, :d], targets, bits)
    if (not np.isfinite(state).all() or not np.isfinite(gram).all()
            or np.linalg.norm(state.conj().T@state+gram-initial.conj().T@initial, 2) > 1e-8):
        raise ArithmeticError('complete source/reference probability accounting failed')
    return recovery.CompleteRun(word, bits, frozen(initial), _FrozenPrimitives(tuple(sorted(primitives.items()))),
                                frozen(state), frozen(weights), frozen(gram))


def admit(run, target, *, word, bits, initial):
    """Bind every output and failure prefix to the requested complete history."""
    tape, initial = _input(word, bits, initial)
    expected = tuple((e['kind'], e['angle'], tuple(e['targets'])) for e in tape)
    if (type(run) is not recovery.CompleteRun or type(run.bits) is not int or run.bits != bits
            or type(run.word) is not tuple
            or any(type(e) is not tuple or len(e) != 3 or type(e[0]) is not str
                   or type(e[1]) is not float or not math.isfinite(e[1])
                   or type(e[2]) is not tuple or any(type(p) is not int for p in e[2])
                   for e in run.word)
            or run.word != expected
            or type(run.primitives) is not _FrozenPrimitives
            or type(run.primitives.entries) is not tuple
            or any(type(entry) is not tuple or len(entry) != 2 or type(entry[0]) is not tuple
                   or len(entry[0]) != 2 or type(entry[0][0]) is not str
                   or type(entry[0][1]) is not float or not math.isfinite(entry[0][1])
                   for entry in run.primitives.entries)
            or set(run.primitives) != {(kind, angle) for kind, angle, _ in expected}):
        raise ValueError('literal complete source archive required')
    if len(run.primitives) != len(set(run.primitives)):
        raise ValueError('distinct immutable source primitive entries required')
    for name in ('initial', 'success', 'failure_weights', 'failure_gram'):
        _immutable_array(getattr(run, name))
    for actual in run.primitives.values():
        _immutable_array(actual)
    if not np.array_equal(run.initial, initial):
        raise ArithmeticError('source archive changed its requested initial state')
    for key, actual in run.primitives.items():
        if not np.array_equal(actual, compiler.native_pulse(*key)):
            raise ArithmeticError('source archive changed its complete primitive')
    replay = execute_complete(tape, bits, initial)
    for name in ('success', 'failure_weights', 'failure_gram'):
        if not np.array_equal(getattr(run, name), getattr(replay, name)):
            raise ArithmeticError('source archive differs from its literal failure-prefix replay')
    if (type(target) is not np.ndarray or target.shape != initial.shape
            or target.dtype.kind not in 'fc' or not np.isfinite(target).all()):
        raise ValueError('finite target on every output address required')
    bound, _ = pulse_budget(tape)
    error = float(np.linalg.norm(run.success-target, 2))+math.sqrt(max(0., float(np.linalg.norm(run.failure_gram, 2))))
    if not math.isfinite(error) or error > float(bound):
        raise ArithmeticError('complete source output exceeds the primitive-composition budget')
    return error, bound


def _mask_tuple(masks, maximum=2):
    if (type(masks) is not tuple or not 1 <= len(masks) <= maximum
            or any(type(mask) is not int or not 0 <= mask <= FULL for mask in masks)
            or len(set(masks)) != len(masks)):
        raise ValueError('distinct bounded native occupation masks required')
    return masks


def target_columns(law, character, masks=(0, FULL)):
    bits = _domain(law, character)
    masks = _mask_tuple(masks)
    out = np.zeros((1 << bits, len(masks)), complex)
    for column, mask in enumerate(masks):
        for coin in (0, 1):
            raw, moved, sign, _ = many_body.column(mask, law, SEAM, coin, character)
            out[moved+(coin << COIN)+(raw << PLUS), column] = sign/math.sqrt(2)
    return out


def certify_repair(law, character):
    """All active columns bind the word; other six native ports are untouched."""
    bits = _domain(law, character)
    word = repair_word(law, character)
    active = tuple(range(6))+tuple(range(12, bits))
    chart = {p: j for j, p in enumerate(active)}
    if any(p not in chart for e in word for p in e['targets']):
        raise ArithmeticError('repair word touches an unaccounted spectator')
    compact = [compiler.pulse(e['kind'], e['angle'], [chart[p] for p in e['targets']]) for e in word]
    out = np.zeros((1 << len(active), 64), complex)
    out[:64] = np.eye(64)
    for event in compact:
        columns = _ideal_pulse(event)
        out = compiler._apply(out, columns[:columns.shape[1]], event['targets'], len(active))
    expected = np.zeros_like(out)
    for mask in range(64):
        for coin in (0, 1):
            raw, moved, sign, _ = many_body.column(mask, law, SEAM, coin, character)
            expected[moved+(coin << 6)+(raw << 7), mask] = sign/math.sqrt(2)
    error = float(np.linalg.norm(out-expected, 2))
    if not math.isfinite(error) or error > float(IDEAL_WORD_CEILING):
        raise ArithmeticError('literal repair word fails general-input native continuation')
    # Independently bind every untouched-spectator placement, not just vacuum.
    for mask in range(1 << MODES):
        for coin in (0, 1):
            raw, moved, sign, current = many_body.column(mask, law, SEAM, coin, character)
            base = many_body.column(mask & 63, law, SEAM, coin, character)
            if (raw, moved & 63, sign, current) != base or moved & ~63 != mask & ~63:
                raise ArithmeticError('native rule depends on an omitted spectator')
    bound, primitive = pulse_budget(word)
    return dict(law=law, character=character, all_native_input_columns=4096,
                compact_active_input_columns=64, complete_raw_slots_per_coin=4 if law == 'C' else 8,
                full_input_complete_amplitude_bound=str(bound+IDEAL_WORD_CEILING),
                ideal_word_amplitude_ceiling=str(IDEAL_WORD_CEILING), raw_primitive_residual=primitive,
                raw_ideal_word_residual=error, resources=_resources(word))


def _controlled_swap(control, a, b):
    return compiler.cnot(b, a)+_flip((control, a), (1, 1), b)+compiler.cnot(b, a)


def _beam():
    return compiler.remap(compiler.pauli_rotation(['X', 'Y'], -math.pi/8)
                          +compiler.pauli_rotation(['Y', 'X'], math.pi/8), [0, 5])


def _preparation(family):
    if family == 'extreme':
        word = [compiler.pulse('Y', -math.pi/4, [0])]
        for j in range(1, MODES):
            word += compiler.cnot(0, j)
        return word, (0, 1), (0, FULL)
    if family == 'fixed_number':
        word = _beam()
        for a, b in PAIRS:
            word += _controlled_swap(5, a, b)
        return word, (SIX_LEFT, SIX_RIGHT), (SIX_LEFT, FULL ^ SIX_LEFT)
    raise ValueError('extreme or fixed_number preparation family required')


def _decode(mask, family):
    """Independent inverse encoder on an occupation address, including phases."""
    if family == 'extreme':
        if mask & 1:
            mask ^= FULL ^ 1
        base, occupied = mask & ~1, mask & 1
        return ((base, 1/math.sqrt(2)), (base | 1, (1 if occupied else -1)/math.sqrt(2)))
    for a, b in reversed(PAIRS):
        if mask & (1 << 5) and ((mask >> a) ^ (mask >> b)) & 1:
            mask ^= (1 << a) | (1 << b)
    left, right = mask & 1, (mask >> 5) & 1
    if left == right:
        return ((mask, 1.),)
    base = mask & ~33
    return ((base | 1, 1/math.sqrt(2)),
            (base | 32, (1 if right else -1)/math.sqrt(2)))


def _resources(word):
    values = compiler.resources(word)
    return {key: values[key] for key in ('native_pulses', 'single_code_pulses', 'pair_code_pulses',
            'identity_phase_pulses', 'proper_code_load_events', 'proper_code_export_events',
            'proper_code_service_events')} | {
                'first_failure_Kraus_slots': sum(2 if e['kind'] == 'J' else 4 for e in word)}


def preparation_certificate():
    """Independent whole-macro number conservation, before any code promise."""
    beam = compiler.remap(_beam(), list(range(6)))
    # Compress the two active endpoints; other four coordinates are untouched.
    beam = [compiler.pulse(e['kind'], e['angle'], [0 if p == 0 else 1 for p in e['targets']])
            for e in beam]
    swap = _controlled_swap(0, 1, 2)
    for bits, word, name in ((2, beam, 'beam'), (3, swap, 'controlled_swap')):
        out = np.eye(1 << bits, dtype=complex)
        for event in word:
            columns = _ideal_pulse(event)
            out = compiler._apply(out, columns[:columns.shape[1]], event['targets'], bits)
        expected = np.eye(1 << bits, dtype=complex)
        if name == 'beam':
            expected[np.ix_([1, 2], [1, 2])] = np.array([[1., -1.], [1., 1.]])/math.sqrt(2)
        else:
            expected = np.zeros_like(expected)
            for mask in range(8):
                moved = mask ^ 6 if mask & 1 and ((mask >> 1) ^ (mask >> 2)) & 1 else mask
                expected[moved, mask] = 1
        if np.linalg.norm(out-expected, 2) > float(IDEAL_WORD_CEILING):
            raise ArithmeticError('preparation macro differs from its number-preserving target')
        numbers = np.array([mask.bit_count() for mask in range(1 << bits)])
        if np.any(expected[numbers[:, None] != numbers[None, :]]):
            raise ArithmeticError('fixed-number preparation macro changed number')
    return dict(fixed_number_preparation_macro_columns=12,
                ideal_macros_preserve_total_number_on_all_inputs=True,
                numerical_source_leakage_retained_within_complete_amplitude_budget=True,
                extreme_cross_number_preparation_outside_fixed_number_V2=True)


def experiment(law, character, *, control='none', family='extreme', return_run=False):
    """Two logical columns, all records, and a phase-sensitive inverse encoder."""
    base_bits = _domain(law, character)
    if type(control) is not str or control not in ('none', 'local_z', 'which_path'):
        raise ValueError('literal none, local_z or which_path control required')
    if type(family) is not str or family not in ('extreme', 'fixed_number'):
        raise ValueError('extreme or fixed_number family required')
    if type(return_run) is not bool:
        raise ValueError('explicit return_run Boolean required')
    if family == 'fixed_number' and control != 'none':
        raise ValueError('fixed-number comparison already retains the native distinguishing record')
    bits = base_bits+int(control == 'which_path')
    prepare, seeds, codes = _preparation(family)
    intervention = []
    if control == 'local_z':
        intervention = [compiler.pulse('Z', math.pi/2, [6]), compiler.pulse('I', -math.pi/2, [6])]
    elif control == 'which_path':
        intervention = compiler.cnot(6, base_bits)
    word = prepare+repair_word(law, character)+intervention+compiler.inverse(prepare)
    if family == 'fixed_number':
        preparation_certificate()
    initial = np.zeros((1 << bits, 2), complex)
    initial[seeds[0], 0] = initial[seeds[1], 1] = 1
    target = np.zeros_like(initial)
    for column in range(2):
        for k, mask in enumerate(codes):
            prepared = (-1 if column == 1 and k == 0 else 1)/math.sqrt(2)
            for coin in (0, 1):
                raw, moved, sign, _ = many_body.column(mask, law, SEAM, coin, character)
                amplitude = prepared*sign/math.sqrt(2)
                if control == 'local_z':
                    amplitude *= (-1)**((moved >> 6) & 1)
                record = (coin << COIN)+(raw << PLUS)
                if control == 'which_path':
                    record += ((moved >> 6) & 1) << base_bits
                for decoded, coefficient in _decode(moved, family):
                    target[decoded+record, column] += amplitude*coefficient
    run = execute_complete(word, bits, initial)
    error, bound = admit(run, target, word=word, bits=bits, initial=initial)
    probabilities = np.abs(run.success)**2
    ideal_probabilities = np.abs(target)**2
    if np.linalg.norm(np.sqrt(probabilities)-np.sqrt(ideal_probabilities), 2) > float(bound):
        raise ArithmeticError('complete terminal law exceeds amplitude budget')
    # Retain the ideal law sparsely with an explicit dimension/default. Every
    # actual amplitude remains in the returned/replayed run; this report stores
    # its full mass and forbidden mass without thresholding small entries.
    laws = []
    for column in range(2):
        nonzero = {}
        for i in np.flatnonzero(ideal_probabilities[:, column]):
            probability = Fraction(float(ideal_probabilities[i, column])).limit_denominator(8)
            if probability <= 0 or abs(float(probability)-ideal_probabilities[i, column]) > 1e-14:
                raise ArithmeticError('ideal law is not the declared exact dyadic distribution')
            nonzero[str(i)] = str(probability)
        laws.append(dict(ideal_nonzero_terminal_probabilities=nonzero,
                         every_other_ideal_probability='0',
                         actual_probability_mass=float(probabilities[:, column].sum()),
                         actual_forbidden_probability_mass=float(probabilities[ideal_probabilities[:, column] == 0, column].sum())))
    report = dict(law=law, character=character, family=family, control=control,
                  native_occupation_bits=12, work_bits=bits, logical_reference_columns=2,
                  complete_terminal_addresses=1 << bits, supplied_initial_field_occupations=[s.bit_count() for s in seeds],
                  blank_apparatus_bits=bits-12, encoded_occupation_numbers=[s.bit_count() for s in codes],
                  complete_output_amplitude_bound=str(bound), raw_complete_output_residual=error,
                  total_first_failure_weight=float(run.failure_weights.sum()),
                  terminal_laws=laws, resources=_resources(word),
                  cat_cross_number_phase_apparatus_supplied=family == 'extreme',
                  encoded_preparation_preserves_declared_total=family == 'fixed_number',
                  encoded_preparation_parity_preservation_claimed=family == 'fixed_number',
                  primitive_encoded_number_conservation_claimed=False,
                  program_class_v2_execution_claimed=False,
                  physical_particle_or_clock_identified=False,
                  all_failure_outputs='immutable literal first-failure prefix archive',
                  source_word_sha256=hashlib.sha256(json.dumps(word, sort_keys=True, separators=(',', ':')).encode()).hexdigest(),
                  postselection=False)
    return (report, run, target) if return_run else report


def certificate():
    """Deterministic receipt; fresh execution admits every unrounded output.

    Floating observations remain available from experiment(return_run=True).
    This receipt binds fixed numerical ceilings, exact ideal laws and source
    words, rather than comparing platform-dependent residuals with a floor.
    """
    repairs, experiments = [], []
    for law in ('C', 'E'):
        for sign in (1, -1):
            repair = certify_repair(law, sign)
            repairs.append({key: value for key, value in repair.items() if not key.startswith('raw_')})
            for family, control in (('extreme', 'none'), ('extreme', 'local_z'),
                                    ('extreme', 'which_path'), ('fixed_number', 'none')):
                report = experiment(law, sign, family=family, control=control)
                entry = {key: report[key] for key in ('law', 'character', 'family', 'control',
                         'work_bits', 'complete_terminal_addresses', 'complete_output_amplitude_bound',
                         'source_word_sha256', 'resources')}
                entry['ideal_terminal_laws'] = [row['ideal_nonzero_terminal_probabilities']
                                                for row in report['terminal_laws']]
                experiments.append(entry)
    return dict(native_occupation_bits=12, maximum_work_bits=MAX_BITS,
                logical_reference_columns=2, every_unlisted_ideal_terminal_probability='0',
                all_actual_complex_columns_and_failure_prefixes_freshly_admitted=True,
                all_failure_outputs='immutable literal first-failure prefix archive',
                numerical_primitive_amplitude_ceiling=str(PULSE_CEILING),
                preparations=preparation_certificate(), general_input_repairs=repairs,
                experiments=experiments,
                supplied_initial_field_occupations={'extreme': [0, 1], 'fixed_number': [6, 6]},
                encoded_occupation_numbers={'extreme': [0, 12], 'fixed_number': [6, 6]},
                encoded_preparation_preserves_declared_total={'extreme': False, 'fixed_number': True},
                encoded_preparation_parity_preservation_claimed={'extreme': False, 'fixed_number': True},
                extreme_cross_number_phase_apparatus_supplied=True,
                primitive_encoded_number_conservation_claimed=False,
                numerical_admissions_are_interval_or_hardware_certificates=False,
                physical_particle_or_clock_identified=False, postselection=False,
                program_class_v2_execution_claimed=False,
                extreme_sector_cat_is_fixed_number_V2_observer=False)
