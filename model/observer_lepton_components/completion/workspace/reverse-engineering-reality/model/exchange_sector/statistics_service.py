"""Source-executed exchange controls on the full padded native register.

The signs are supplied coherent controls, not inferred from C/E populations.
Every executed gate comes from the existing proper-code source compiler.
Integer Walsh coefficients certify the diagonal words on all addresses;
Gray permutations use the executed source macro, never target matrices.
The optional path bit is an explicitly prepared and read observer register.
"""
from __future__ import annotations

from functools import lru_cache
from itertools import combinations_with_replacement
import math

import numpy as np

from ..common_source import native_routing as routing, record_code
from ..observer.source_discrimination import circuit, experiment as recorded
from ..observer.source_scalar import compiler
from . import statistics as exact


SIZE = 512
LOOP = ((0, 5), (0, 11), (5, 11))
LOOP_PAIR = (0, 11)
COUNT_KEYS = routing.COUNT_KEYS
NUMERICAL_LIMIT = 2e-7


def _address(pair):
    return circuit.encode(tuple(int(p in pair) for p in range(12)))


def _state(value, bits):
    routing._original_scalars(value)
    raw = np.asarray(value)
    if (raw.ndim not in (1, 2) or raw.shape[0] != 1 << bits
            or (raw.ndim == 2 and not raw.shape[1])
            or raw.dtype.kind not in 'iufc' or not np.isfinite(raw).all()):
        raise ValueError('finite full-register amplitudes or nonempty reference columns required')
    out = np.array(raw, dtype=complex, copy=True)
    if not np.isfinite(out).all():
        raise ValueError('amplitudes must remain finite after conversion')
    return out


def _bool(value):
    if type(value) is not bool:
        raise ValueError('an explicit Boolean control switch required')
    return value


def _quarters(values):
    if (type(values) is not tuple or len(values) not in (512, 1024)
            or any(type(x) is not int for x in values)):
        raise ValueError('512 or 1024 integer quarter-turn phases required')
    return values


def walsh_coefficients(values):
    """Unnormalized integer Walsh transform; no floating threshold trimming."""
    out = list(_quarters(values))
    stride = 1
    while stride < len(out):
        for start in range(0, len(out), 2*stride):
            for i in range(start, start+stride):
                a, b = out[i], out[i+stride]
                out[i], out[i+stride] = a+b, a-b
        stride *= 2
    return tuple(out)


def diagonal_word(values):
    """Emit exp(i*pi*values/2), sharing parity along cyclic Gray walks.

    Each group restores its parity wires.  Coefficients, including the
    scalar coefficient, are exact integers divided by the register size.
    This is the same native parity compilation as the routing macro, with
    arbitrary coefficients and a possible explicit tenth path bit.
    """
    values = _quarters(values)
    coefficients = walsh_coefficients(values)
    bits, size = len(values).bit_length()-1, len(values)
    word = []
    if coefficients[0]:
        word.append(compiler.pulse('I', math.pi*coefficients[0]/(2*size), [0]))
    for pivot in range(bits):
        if not any(coefficients[1 << pivot:1 << (pivot+1)]):
            continue
        count = 1 << pivot
        for index in range(count):
            gray = index ^ (index >> 1)
            coefficient = coefficients[(1 << pivot) | gray]
            if coefficient:
                word.append(compiler.pulse('Z', math.pi*coefficient/(2*size), [pivot]))
            following = (index+1) % count
            change = gray ^ (following ^ (following >> 1))
            if change:
                word.extend(compiler.cnot(change.bit_length()-1, pivot))
    return word


def _audit_diagonal_word(values, word):
    """Read the emitted pulses independently as an exact parity polynomial.

    Only literal phase-exact native CNOT triples and rational-quarter I/Z
    pulses are admitted.  A transform agreeing with itself is insufficient:
    this verifies the actual emitted order, angles, wires and scalar phase.
    """
    size, bits = len(values), len(values).bit_length()-1
    parity, accumulated = [1 << bit for bit in range(bits)], [0]*size
    position = 0
    while position < len(word):
        event = word[position]
        if event['kind'] in ('I', 'Z'):
            if len(event['targets']) != 1 or not 0 <= event['targets'][0] < bits:
                raise ArithmeticError('diagonal word has an invalid phase wire')
            coefficient = round(event['angle']*2*size/math.pi)
            if event['angle'] != math.pi*coefficient/(2*size):
                raise ArithmeticError('diagonal phase is outside its exact rational-angle menu')
            mask = 0 if event['kind'] == 'I' else parity[event['targets'][0]]
            accumulated[mask] += coefficient
            position += 1
        else:
            if position+3 > len(word) or word[position+1]['kind'] != 'J':
                raise ArithmeticError('literal native CNOT triple required')
            targets = word[position+1]['targets']
            if len(targets) != 2 or any(not 0 <= t < bits for t in targets):
                raise ArithmeticError('invalid native CNOT wires')
            control, target = targets
            if control == target or word[position:position+3] != compiler.cnot(control, target):
                raise ArithmeticError('phase-inexact or incomplete native CNOT triple')
            parity[target] ^= parity[control]
            position += 3
    if parity != [1 << bit for bit in range(bits)]:
        raise ArithmeticError('diagonal word did not restore its parity wires')
    if tuple(accumulated) != walsh_coefficients(values):
        raise ArithmeticError('emitted diagonal phases disagree on the full address register')


def phase_quarters(seam, current=1, character=-1):
    seam, character = exact._seam(seam), exact._character(character)
    if type(current) is not int or current not in (-1, 0, 1):
        raise ValueError('an integer signed native current required')
    phases = [0]*SIZE
    for pair in exact.PAIRS:
        if exact.phase_on_output(seam, current, pair, character) == -1:
            phases[_address(pair)] = 2
    return tuple(phases)


def phase_word(seam, current=1, character=-1):
    return diagonal_word(phase_quarters(seam, current, character))


def apply_phase(seam, current, state, character=-1):
    word = phase_word(seam, current, character)
    return routing._immutable(compiler.execute(word, 9, _state(state, 9)))


def _word_key(word):
    return tuple((p['kind'], p['angle'], tuple(p['targets'])) for p in word)


@lru_cache(maxsize=2)
def _executed_gray(key):
    # Cache is keyed by the emitted word itself.  No polar fit or ideal
    # permutation replaces the raw source-executed matrix.
    word = [compiler.pulse(kind, angle, list(targets)) for kind, angle, targets in key]
    actual = compiler.execute(word, 9)
    if actual.shape != (512, 512) or not np.isfinite(actual).all():
        raise ArithmeticError('finite complete source Gray matrix required')
    target = np.eye(512, dtype=complex)
    target[[0, 256]] = target[[256, 0]]
    error = float(np.linalg.norm(actual-target))
    _guard(error, 'full512 executed Gray macro operator upper bound')
    return routing._immutable(actual), error


def _permutation(images, state, inverse=False):
    images, inverse = routing._permutation(images), _bool(inverse)
    vector = state.ndim == 1
    columns = state.reshape(SIZE, -1)
    macro = _executed_gray(_word_key(routing._canonical_word(inverse)))[0]
    charts = {}
    for target, pattern in routing.permutation_edges(images, inverse=inverse):
        if target not in charts:
            charts[target] = routing._chart(macro, target)
        before = routing._conjugation_word(pattern)
        columns = compiler.execute(before, 9, columns)
        columns = charts[target]@columns
        columns = compiler.execute(compiler.inverse(before), 9, columns)
    return columns[:, 0] if vector else columns


def _source_permutation(images, state, inverse=False):
    """Nine source wires, with an optional untouched path/reference factor."""
    if state.shape[0] == SIZE:
        return _permutation(images, state, inverse)
    vector = state.ndim == 1
    count = 1 if vector else state.shape[1]
    columns = state.reshape(2, SIZE, count).transpose(1, 0, 2).reshape(SIZE, 2*count)
    out = _permutation(images, columns, inverse)
    out = out.reshape(SIZE, 2, count).transpose(1, 0, 2).reshape(2*SIZE, count)
    return out[:, 0] if vector else out


def _bare_repair(law, seam, state):
    """Fresh input execution, complete native raw slots and unchanged references."""
    law, seam = routing._candidate(law), exact._seam(seam)
    route = record_code._route_permutation(seam)
    routed = _permutation(route, state)
    rows = recorded.repair_branches(law, routed)
    _validate_grid(law, rows, shape=state.shape)
    # Apply inverse routing to all raw branches, including abort amplitudes.
    width = 1 if state.ndim == 1 else state.shape[1]
    joined = np.concatenate([r['amplitude'].reshape(SIZE, width) for r in rows], axis=1)
    joined = _permutation(route, joined, True)
    for i, row in enumerate(rows):
        value = joined[:, i*width:(i+1)*width]
        row['amplitude'] = value[:, 0] if state.ndim == 1 else value
    return rows


def repair_branches(law, seam, state, character=-1):
    """Native repair followed by explicitly supplied record-selected feedback.

    Bare C/E populations and records are unchanged.  The negative-character
    quantum instrument differs coherently on moving blocks.  No slot,
    padding row, small amplitude or untouched reference is discarded.
    """
    law, seam, character = routing._candidate(law), exact._seam(seam), exact._character(character)
    state = _state(state, 9)
    rows = _bare_repair(law, seam, state)
    _validate_grid(law, rows, shape=state.shape)
    for row in rows:
        if row['current'] in (-1, 1):
            row['amplitude'] = compiler.execute(
                phase_word(seam, row['current'], character), 9, row['amplitude'])
        row['amplitude'] = routing._immutable(row['amplitude'])
    _validate_grid(law, rows, shape=state.shape)
    return rows


def _witness_permutation(pairs):
    return routing.complete_permutation(tuple((j, _address(pair)) for j, pair in enumerate(pairs)))


def witness_columns(law, character=-1):
    """Execute both unknown logical input columns through every raw branch.

    Final amplitudes are full512 x 2 Kraus columns after inverse branch
    transport and X read rotation.  Computational readout retains every
    address; applying these columns to any reference state is linear.
    """
    law, character = routing._candidate(law), exact._character(character)
    seam = exact.WITNESS_SEAM
    prepared = _permutation(_witness_permutation(exact.WITNESS_INPUT_PAIRS),
                            np.eye(SIZE, 2, dtype=complex))
    rows = repair_branches(law, seam, prepared, character)
    for row in rows:
        pairs = ((1, 5), (5, 6)) if row['current'] == 1 else exact.WITNESS_INPUT_PAIRS
        decoded = _permutation(_witness_permutation(pairs), row['amplitude'], True)
        row['amplitude'] = routing._immutable(compiler.execute(
            [compiler.pulse('Y', math.pi/4, [0])], 9, decoded))
    return rows


def witness_events(law, character=-1):
    """Actual source preparation of |+>, followed by the complete observer."""
    law, character = routing._candidate(law), exact._character(character)
    # Native preparation precedes the encoder.  Matrix columns are replayed
    # from fresh unknown inputs, not generated from declared ideal branches.
    logical = compiler.execute([compiler.pulse('Y', -math.pi/4, [0])], 1, [1., 0.])
    return [{**row, 'amplitude': routing._immutable(row['amplitude']@logical)}
            for row in witness_columns(law, character)]


def hop_layout(seam, character=-1, *, controlled=False):
    """A full-register diagonalization of hard-core/CAR seam hopping.

    The ten disjoint blocks are {a,s}<->{b,s}.  Collisions, empty seams
    and padding have generator zero.  For a path control only its arm1
    carries the generator.  No phase is attached to a blocked double.
    """
    a, b = exact._seam(seam)
    character, controlled = exact._character(character), _bool(controlled)
    assignments, diagonal = [], [0]*(1024 if controlled else 512)
    offset = 512 if controlled else 0
    for k, spectator in enumerate(p for p in range(12) if p not in (a, b)):
        first, second = tuple(sorted((a, spectator))), tuple(sorted((b, spectator)))
        sign = exact.sign(first, second, character)
        assignments.extend(((2*k, _address(first)), (2*k+1, _address(second))))
        diagonal[offset+2*k], diagonal[offset+2*k+1] = -sign, sign
    return routing.complete_permutation(tuple(assignments)), tuple(diagonal)


def apply_hop(seam, state, character=-1, *, controlled=False):
    """Source-execute exp(-i*pi*H_seam/2), with all padding retained."""
    controlled = _bool(controlled)
    permutation, quarters = hop_layout(seam, character, controlled=controlled)
    bits = 10 if controlled else 9
    out = _source_permutation(permutation, _state(state, bits), True)
    out = compiler.execute([compiler.pulse('Y', math.pi/4, [0])], bits, out)
    out = compiler.execute(diagonal_word(quarters), bits, out)
    out = compiler.execute([compiler.pulse('Y', -math.pi/4, [0])], bits, out)
    return routing._immutable(_source_permutation(permutation, out))


def _path_read_word():
    # Z(pi/4) turns (|0>+i|1>)/sqrt2 into |+> up to a common phase.
    return [compiler.pulse('Z', math.pi/4, [9]), compiler.pulse('Y', math.pi/4, [9])]


def loop_events(law, character=-1, *, dephase=False):
    """Controlled triangular exchange, then bare C/E and a Ramsey Y read.

    Dephasing means a retained computational measurement of the supplied
    path bit after the first hop.  Both archive outcomes survive; this is
    an explicit which-arm measurement, not an inferred source projection.
    Terminal amplitudes retain both path slots and all512 source addresses.
    """
    law, character, dephase = routing._candidate(law), exact._character(character), _bool(dephase)
    blank = np.zeros(1024, complex)
    blank[0] = 1
    preparation = _witness_permutation((LOOP_PAIR,))
    state = _source_permutation(preparation, blank)
    state = compiler.execute([compiler.pulse('Y', -math.pi/4, [9])], 10, state)
    state = apply_hop(LOOP[0], state, character, controlled=True)
    components = [(None, state)]
    if dephase:
        components = []
        for arm in (0, 1):
            projected = np.zeros_like(state)
            projected[arm*SIZE:(arm+1)*SIZE] = state[arm*SIZE:(arm+1)*SIZE]
            components.append((arm, projected))
    result = []
    for archive, state in components:
        for seam in LOOP[1:]:
            state = apply_hop(seam, state, character, controlled=True)
        # Bare continuation acts identically on both path arms; the path is
        # an untouched reference, and every native raw slot is retained.
        columns = state.reshape(2, SIZE).T
        for row in _bare_repair(law, exact.WITNESS_SEAM, columns):
            value = row['amplitude'].T.reshape(1024)
            value = compiler.execute(_path_read_word(), 10, value)
            result.append({**row, 'dephasing_record': archive,
                           'amplitude': routing._immutable(value)})
    return result


def _counts(word):
    resources = compiler.resources(word)
    return {key: resources[key] for key in COUNT_KEYS}


def _add(*rows):
    return {key: sum(row[key] for row in rows) for key in COUNT_KEYS}


def _failure_account(counts):
    # The source export has one coherent success Kraus map and one retained
    # complement map per unused M6 level: four for d2 and two for d4.
    return {**counts, 'proper_code_complement_outcome_slots':
            4*counts['single_code_pulses']+2*counts['pair_code_pulses']}


def phase_resources(seam, current=1, character=-1):
    return _counts(phase_word(seam, current, character))


def hop_resources(seam, character=-1, *, controlled=False):
    permutation, quarters = hop_layout(seam, character, controlled=controlled)
    route = record_code._permutation_resources(permutation)
    word = ([compiler.pulse('Y', math.pi/4, [0])]+diagonal_word(quarters)
            +[compiler.pulse('Y', -math.pi/4, [0])])
    return {**_add(route, route, _counts(word)),
            'work_qubits': 10 if controlled else 9,
            'gray_macro_invocations': 2*route['controlled_flip_invocations']}


def witness_resources(law, character=-1):
    law, character = routing._candidate(law), exact._character(character)
    seam = exact.WITNESS_SEAM
    encoding = record_code._permutation_resources(_witness_permutation(exact.WITNESS_INPUT_PAIRS))
    prepare = _counts([compiler.pulse('Y', -math.pi/4, [0])])
    read = _counts([compiler.pulse('Y', math.pi/4, [0])])
    repair = record_code.repair_resources(law, seam)
    branches = []
    for coin in (0, 1):
        labels = circuit.record_values(law, coin)
        for raw in range(1 << circuit.RECORD_BITS[law]):
            value = labels.get(raw)
            current = value if law == 'C' or value is None else (
                value-(2*coin-1)*(abs(value) % 2))//2
            pairs = ((1, 5), (5, 6)) if current == 1 else exact.WITNESS_INPUT_PAIRS
            decoder = record_code._permutation_resources(_witness_permutation(pairs))
            phase = phase_resources(seam, current, character) if current in (-1, 1) else dict.fromkeys(COUNT_KEYS, 0)
            branches.append({'coin': coin, 'raw': raw, 'abort': raw not in labels,
                             **_failure_account(_add(prepare, encoding, repair, phase, decoder, read))})
    return {'blank_source_qubits': 9, 'fresh_coin_qubits': 1,
            'unknown_logical_column_execution_supplied_qubits': 1,
            'unknown_logical_column_execution_blank_qubits': 8,
            'fresh_record_qubits': circuit.RECORD_BITS[law],
            'native_coin_readout_bits': 1,
            'native_record_readout_bits': circuit.RECORD_BITS[law],
            'terminal_computational_source_readout_bits': 9,
            'retained_classical_bits_per_run': 10+circuit.RECORD_BITS[law],
            'terminal_source_outcome_slots_per_raw_record': 512,
            'raw_record_slots': len(branches), 'branches': branches,
            'physical_initialization_storage_readout_costs_calibrated': False,
            'proper_code_failure_outcomes_retained': True,
            'proper_code_failure_probability_exact': '0'}


def loop_resources(law, character=-1, *, dephase=False):
    law, character, dephase = routing._candidate(law), exact._character(character), _bool(dephase)
    prepare = record_code._permutation_resources(_witness_permutation((LOOP_PAIR,)))
    path = _counts([compiler.pulse('Y', -math.pi/4, [9])]+_path_read_word())
    hops = [hop_resources(seam, character, controlled=True) for seam in LOOP]
    repair = record_code.repair_resources(law, exact.WITNESS_SEAM)
    return {**_failure_account(_add(prepare, path, *hops, repair)),
            'blank_source_qubits': 9, 'blank_coherent_path_qubits': 1,
            'fresh_coin_qubits': 1, 'fresh_record_qubits': circuit.RECORD_BITS[law],
            'native_coin_readout_bits': 1,
            'native_record_readout_bits': circuit.RECORD_BITS[law],
            'terminal_path_readout_bits': 1,
            'terminal_source_computational_measurements': 0,
            'retained_classical_bits_per_run': 2+circuit.RECORD_BITS[law]+int(dephase),
            'retained_intermediate_path_measurements': int(dephase),
            'dephasing_archive_slots': 2 if dephase else 1,
            'raw_native_record_slots': 2 << circuit.RECORD_BITS[law],
            'terminal_path_outcomes_per_raw_record': 2,
            'terminal_source_addresses_per_path_outcome': 512,
            'physical_initialization_storage_readout_costs_calibrated': False,
            'proper_code_failure_outcomes_retained': True,
            'proper_code_failure_probability_exact': '0'}


def _word_error(word):
    # Numerical algebraic product bound on arbitrary inputs/reference
    # extensions.  It is not an interval or hardware error certificate.
    value = math.expm1(record_code._word_log_error(word))
    if not math.isfinite(value) or value < 0:
        raise ArithmeticError('finite nonnegative primitive-product error required')
    return value


def _permutation_error(images, inverse=False):
    log_error = 0.
    macro = math.log1p(_executed_gray(_word_key(routing._canonical_word(inverse)))[1])
    for _, pattern in routing.permutation_edges(images, inverse=inverse):
        before = routing._conjugation_word(pattern)
        log_error += macro+record_code._word_log_error(before)+record_code._word_log_error(compiler.inverse(before))
    return math.expm1(log_error)


def hop_error_bound(seam, character=-1, *, controlled=False):
    """Complete source product bound on any full-register/reference input.

    The macro's full512 Frobenius residual bounds its operator residual;
    wire charts and untouched path/reference extensions preserve that bound.
    Both permutation directions, both Y rotations and every diagonal pulse
    are included.  Matrix arithmetic remains ordinary floating arithmetic.
    """
    permutation, quarters = hop_layout(seam, character, controlled=controlled)
    parts = (_permutation_error(permutation, True),
             _word_error([compiler.pulse('Y', math.pi/4, [0])]),
             _word_error(diagonal_word(quarters)),
             _word_error([compiler.pulse('Y', -math.pi/4, [0])]),
             _permutation_error(permutation))
    value = float(math.expm1(math.fsum(math.log1p(x) for x in parts)))
    _guard(value, 'complete hopping arbitrary-input product bound')
    return value


def _validate_grid(law, rows, *, loop=False, dephase=False, columns=False, shape=None):
    expected = {(archive, coin, raw) for archive in ((0, 1) if dephase else (None,))
                for coin in (0, 1) for raw in range(1 << circuit.RECORD_BITS[law])}
    if type(rows) is not list or len(rows) != len(expected):
        raise ArithmeticError('complete unique native outcome grid required')
    fields = {'coin', 'record_code', 'current', 'contrast', 'abort', 'amplitude'}
    if loop:
        fields.add('dephasing_record')
    seen = set()
    for row in rows:
        if type(row) is not dict or set(row) != fields:
            raise ArithmeticError('exact source event fields required')
        coin, raw = row['coin'], row['record_code']
        archive = row.get('dephasing_record')
        if (type(coin) is not int or type(raw) is not int or type(row['abort']) is not bool
                or (dephase and type(archive) is not int)):
            raise ArithmeticError('native event labels require exact types')
        key = archive, coin, raw
        if key not in expected or key in seen:
            raise ArithmeticError('invalid, duplicate or missing native event')
        seen.add(key)
        value = circuit.record_values(law, coin).get(raw)
        current = value if law == 'C' or value is None else (
            value-(2*coin-1)*(abs(value) % 2))//2
        contrast = value if law == 'E' else None
        if (type(row['current']) is not type(current) or row['current'] != current
                or type(row['contrast']) is not type(contrast) or row['contrast'] != contrast
                or row['abort'] != (value is None)):
            raise ArithmeticError('raw flags and native decoded labels disagree')
        amplitude = row['amplitude']
        expected_shape = shape if shape is not None else (1024,) if loop else (512, 2) if columns else (512,)
        if (np.ma.isMaskedArray(amplitude) or not isinstance(amplitude, np.ndarray)
                or amplitude.shape != expected_shape or amplitude.dtype.kind not in 'iufc'
                or not np.isfinite(amplitude).all()):
            raise ArithmeticError('finite complete source amplitudes required')


def _guard(value, name):
    if type(value) is not float or not math.isfinite(value) or not 0 <= value <= NUMERICAL_LIMIT:
        raise ArithmeticError(f'{name} exceeds the fixed finite numerical admission')


def _source_address_certificate():
    """Independent literal78-state chart, not a self-consistent encode round trip."""
    local = ((0, 0), (0, 1), (1, 0), (0, 2), (1, 1), (2, 0))
    shells = []
    for total in range(3):
        entries = []
        for particles in combinations_with_replacement(range(10), total):
            entries.append(tuple(particles.count(p) for p in range(10)))
        entries.sort()
        if total:
            preferred = (1,)+(0,)*9 if total == 1 else (1, 1)+(0,)*8
            entries.remove(preferred)
            entries.insert(0, preferred)
        shells.append(tuple(entries))
    checked = set()
    for pair in combinations_with_replacement(range(12), 2):
        state = tuple(pair.count(p) for p in range(12))
        expected = local.index(state[:2])+8*shells[sum(state[2:])].index(state[2:])
        if circuit.encode(state) != expected or circuit.decode(expected) != state:
            raise ArithmeticError('source address chart differs from the literal carrier binding')
        if pair[0] != pair[1] and _address(pair) != expected:
            raise ArithmeticError('exchange source address changes the bound amplitude chart')
        checked.add(expected)
    if len(checked) != 78:
        raise ArithmeticError('source carrier address binding is not injective')
    return len(checked)


def _proper_code_certificate():
    # Loading/exact invariant-code controls are supplied services already
    # checked against the bound full-M6 generators by compiler._source().
    # Keep the complete export instrument, including every zero failure map
    # on admitted proper-code inputs, as an explicit terminal-halt family.
    source = compiler.source_summary()
    _guard(float(source['native_generator_antihermitian_residual']), 'native source generator invariant')
    coefficient = source['maximum_base_coefficient']
    if (type(coefficient) is not float or not math.isfinite(coefficient) or coefficient < 0
            or source['source_generators'] != 12 or source['exact_control_rows_checked'] != 10):
        raise ArithmeticError('finite bound native source control data required')
    rows = {}
    for dimension in (2, 4):
        kraus = compiler.transfer_kraus(dimension)
        expected = np.zeros((7-dimension, dimension+1, 6), complex)
        expected[0, :dimension, :dimension] = np.eye(dimension)
        for level in range(dimension, 6):
            expected[1+level-dimension, dimension, level] = 1
        if (len(kraus) != 7-dimension
                or any(k.shape != (dimension+1, 6) or not np.isfinite(k).all() for k in kraus)
                or any(not np.array_equal(k, target) for k, target in zip(kraus, expected, strict=True))
                or not np.array_equal(sum((k.conj().T@k for k in kraus), np.zeros((6, 6))), np.eye(6))
                or any(np.any(k[:, :dimension]) for k in kraus[1:])):
            raise ArithmeticError('complete proper-code export/failure instrument changed')
        rows[str(dimension)] = {'success_outcomes': 1, 'retained_complement_outcomes': 6-dimension,
                                'literal_success_and_complement_columns_checked': True,
                                'admitted_input_complement_probability_exact': '0',
                                'complement_disposition': 'retained terminal halt'}
    return rows


def _audit_record_phase(seam, current, character, quarters):
    # Independent endpoint oracle, including all collision and padding rows.
    # Walsh inversion alone only certifies the supplied phases, even if wrong.
    expected = [0]*SIZE
    a, b = seam
    if character == -1 and current:
        destination = b if current == 1 else a
        for spectator in range(a+1, b):
            expected[_address(tuple(sorted((destination, spectator))))] = 2
    if quarters != tuple(expected):
        raise ArithmeticError('record phase disagrees with native endpoint action')


def _layout_certificate():
    addresses = _source_address_certificate()
    phases, hops, largest_word = 0, 0, 0
    for seam in exact.SEAMS:
        for current in (-1, 0, 1):
            for character in (-1, 1):
                quarters = phase_quarters(seam, current, character)
                _audit_record_phase(seam, current, character, quarters)
                recovered = walsh_coefficients(walsh_coefficients(quarters))
                if recovered != tuple(SIZE*x for x in quarters):
                    raise ArithmeticError('full512 integer Walsh inversion failed')
                word = diagonal_word(quarters)
                _audit_diagonal_word(quarters, word)
                _guard(float(_word_error(word)), 'phase product bound')
                largest_word = max(largest_word, len(word))
                phases += 1
        for character in (-1, 1):
            for controlled in (False, True):
                permutation, quarters = hop_layout(seam, character, controlled=controlled)
                schedule = record_code._integer_schedule_action(routing.permutation_edges(permutation))
                inverse = record_code._integer_schedule_action(routing.permutation_edges(permutation, inverse=True))
                if schedule != permutation or any(inverse[schedule[i]] != i for i in range(SIZE)):
                    raise ArithmeticError('complete hopping diagonalizer schedule changed')
                # Independent endpoint oracle verifies actual block columns.
                a, b = seam
                for k, s in enumerate(p for p in range(12) if p not in seam):
                    sign = character if a < s < b else 1
                    offset = SIZE if controlled else 0
                    if (permutation[2*k] != _address(tuple(sorted((a, s))))
                            or permutation[2*k+1] != _address(tuple(sorted((b, s))))
                            or quarters[offset+2*k:offset+2*k+2] != (-sign, sign)):
                        raise ArithmeticError('source hopping columns disagree with native endpoints')
                allowed = range(SIZE if controlled else 0, (SIZE if controlled else 0)+20)
                if any(value for i, value in enumerate(quarters) if i not in allowed):
                    raise ArithmeticError('hop changes its zero-generator complement')
                if walsh_coefficients(walsh_coefficients(quarters)) != tuple(len(quarters)*x for x in quarters):
                    raise ArithmeticError('controlled hopping Walsh inversion failed')
                word = diagonal_word(quarters)
                _audit_diagonal_word(quarters, word)
                _guard(float(_word_error(word)), 'hop diagonal product bound')
                hop_error_bound(seam, character, controlled=controlled)
                hops += 1
    return {'literal_native_source_addresses_checked': addresses,
            'record_selected_diagonals_checked': phases,
            'record_phase_full_register_endpoint_oracle': True,
            'full_register_hopping_layouts_checked': hops,
            'largest_record_selected_diagonal_pulse_count': largest_word,
            'complete_hop_arbitrary_input_product_error_upper': '2e-7',
            'complete_hop_unknown_input_register_dimensions': [512, 1024],
            'identity_collision_and_padding_extension': True,
            'full_address_integer_Walsh_inversion': True}


def certificate():
    """Fresh executable admissions, then a deterministic exact receipt body.

    Floating residuals are checked, never serialized as platform-dependent
    receipt bytes.  Cached source matrices are regenerated from their words.
    """
    _executed_gray.cache_clear()
    compiler._native.cache_clear()
    compiler._source.cache_clear()
    recorded._repair.cache_clear()
    recorded._coin.cache_clear()
    record_code._primitive_error.cache_clear()
    record_code._macro_resource_counts.cache_clear()
    record_code._native_repair_counts.cache_clear()
    exports = _proper_code_certificate()
    layout = _layout_certificate()
    witnesses, loops = {}, {}
    for law in ('C', 'E'):
        witnesses[law], loops[law] = {}, {}
        for character in (-1, 1):
            columns = witness_columns(law, character)
            _validate_grid(law, columns, columns=True)
            errors = []
            for row in columns:
                target = np.zeros((512, 2), complex)
                if row['record_code'] == (1 if law == 'C' else 5) and row['coin'] == 0:
                    target[:2] = np.array([[character, 1], [-character, 1]])/2
                elif row['record_code'] == (0 if law == 'C' else 4) and row['coin'] == 1:
                    target[:2] = np.array([[1, 1], [-1, 1]])/2
                errors.append(row['amplitude']-target)
            _guard(float(np.linalg.norm(np.concatenate(errors, axis=0))),
                   'complete unknown-input witness Stinespring columns')
            # Source preparation itself is replayed, preserving its phase.
            logical = compiler.execute([compiler.pulse('Y', -math.pi/4, [0])], 1, [1., 0.])
            _guard(float(np.linalg.norm(logical-np.ones(2)/math.sqrt(2))), 'source plus preparation')
            witnesses[law][str(character)] = {
                'positive_read_probability_exact': '1' if character == 1 else '1/2',
                'negative_read_probability_exact': '0' if character == 1 else '1/2',
                'logical_input_domain_dimension': 2,
                'untouched_reference_dimension': 'arbitrary',
                'combined_raw_record_Stinespring_error_checked': True,
                'all_raw_full512_Kraus_columns_checked': True,
                'resources': witness_resources(law, character)}
            loops[law][str(character)] = {}
            for dephase in (False, True):
                rows = loop_events(law, character, dephase=dephase)
                _validate_grid(law, rows, loop=True, dephase=dephase)
                # Independent endpoint continuation: final pair is (0,11),
                # coin0 moves to (5,11), coin1 holds, each probability1/2.
                expected_support = {0: _address((5, 11)), 1: _address((0, 11))}
                errors = []
                for row in rows:
                    target = np.zeros(1024, complex)
                    wanted = (1 if law == 'C' else 5) if row['coin'] == 0 else (0 if law == 'C' else 4)
                    if row['record_code'] == wanted:
                        address = expected_support[row['coin']]
                        if dephase:
                            scale = np.exp(1j*math.pi/4)/(2*math.sqrt(2))
                            if row['dephasing_record'] == 0:
                                target[address], target[512+address] = scale, -scale
                            else:
                                target[address] = target[512+address] = character*scale
                        else:
                            target[address+(512 if character == -1 else 0)] = character*np.exp(1j*math.pi/4)/math.sqrt(2)
                    errors.append(row['amplitude']-target)
                _guard(float(np.linalg.norm(np.concatenate(errors))),
                       'complete phase-sensitive triangular Ramsey amplitudes')
                loops[law][str(character)]['dephased' if dephase else 'coherent'] = {
                    'positive_path_read_probability_exact': '1/2' if dephase else '1' if character == 1 else '0',
                    'negative_path_read_probability_exact': '1/2' if dephase else '0' if character == 1 else '1',
                    'all_raw_full1024_terminal_amplitudes_retained': True,
                    'combined_complex_record_and_archive_amplitudes_checked': True,
                    'resources': loop_resources(law, character, dephase=dephase)}
    return {'schema': 'oph.exchange_statistics.source_service.v1',
            'layout': layout, 'proper_code_export_instruments': exports,
            'one_attempt_witnesses': witnesses,
            'triangular_exchange': {'seams': [list(s) for s in LOOP],
                                    'initial_and_final_pair': list(LOOP_PAIR),
                                    'common_dynamical_phase': 'i', 'observations': loops},
            'source_admission': {'fresh_source_execution_passed': True,
                                 'arbitrary_reference_columns_retained': True,
                                 'all_zero_and_abort_slots_retained': True,
                                 'fixed_numerical_limit': '2e-7',
                                 'interval_or_hardware_certificate': False},
            'scope': 'Supplied coherent controls distinguish the two characters; bare C/E records do not select either.'}
