"""Source native instruments marked by a local occupation-changing event.

No-tick steps include remote moves. The absorbing resolvent therefore acts
on the retained source state; instantaneous escape normalization is invalid.
All matrices use output rows and input columns.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
import math
from types import MappingProxyType

import numpy as np
import sympy as sp

from ..observer import native
from ..native_quantum_history import construction as quantum
from ..primitive_observer import runtime


CLOCK_PORT = 2  # first source port outside the existing phase pair (0,1)
WAIT = (-1, 0)
ZERO = (0,) * 12
SOURCE_PATHS = tuple(sorted(set((*quantum.SOURCE_PATHS,
    'model/relational_event_clock/kernel.py',
    'Lean/Time/TimeOrderLedger.lean',
    'Lean/Time/ClockReadout.lean',
    'Lean/Screen/RegionalContinuity.lean',
))))


@lru_cache(maxsize=1)
def marks():
    neighbours = sorted(b if a == CLOCK_PORT else a
                        for a, b in quantum.mixing.seams() if CLOCK_PORT in (a, b))
    if len(neighbours) != 5:
        raise ArithmeticError('local source clock does not have five incident seams')
    return (WAIT,) + tuple((v, sign) for v in neighbours for sign in (-1, 1))


def intervention_seam():
    """First remote source seam touching a local clock neighbour."""
    neighbours = {v for v, _ in marks()[1:]}
    return next(seam for seam in quantum.mixing.seams()
                if CLOCK_PORT not in seam and neighbours.intersection(seam))


@dataclass(frozen=True)
class NativeClock:
    law: str
    # Groups retain the native seam, coin and committed C/E label internally.
    # Each item is (visible mark index, seam index, coin, source label, mapping).
    groups: tuple
    instrument: runtime.Instrument


@lru_cache(maxsize=2)
def _source(law):
    law = quantum.candidate(law)
    grouped = defaultdict(dict)
    assigned = {}
    for edge, seam in enumerate(quantum.mixing.seams()):
        a, b = seam
        for coin in (0, 1):
            for i, state in enumerate(quantum.STATES):
                final, record = native.native(state, seam, coin)
                if record != (coin, state[a]-final[a]):
                    raise ArithmeticError('native current grammar changed')
                current = record[1]
                delta = current * (int(b == CLOCK_PORT)-int(a == CLOCK_PORT))
                mark = WAIT if not delta else (b if a == CLOCK_PORT else a, delta)
                label = current if law == 'C' else state[a]-state[b]
                key = (edge, coin, label)
                index = marks().index(mark)
                if key in assigned and assigned[key] != index:
                    raise ArithmeticError('native committed label does not determine clock record')
                assigned[key] = index
                grouped[key][i] = quantum.INDEX[final]
    groups = []
    matrices = [[] for _ in marks()]
    for (edge, coin, label), mapping in sorted(grouped.items()):
        if len(set(mapping.values())) != len(mapping):
            raise ArithmeticError('source native branch is not a partial injection')
        mark = assigned[edge, coin, label]
        groups.append((mark, edge, coin, label, tuple(sorted(mapping.items()))))
        matrix = np.zeros((66, 66), complex)
        for i, j in mapping.items():
            matrix[j, i] = 1/math.sqrt(60)
        matrices[mark].append(matrix)
    counts = [0]*66
    for _, _, _, _, mapping in groups:
        for i, _ in mapping:
            counts[i] += 1
    if set(counts) != {60}:
        raise ArithmeticError('source local clock instrument is incomplete')
    instrument = runtime.Instrument(kraus=tuple(tuple(group) for group in matrices), visible=True,
        resources=(('native_attempts', 1), ('local_clock_reads', 1)),
        binding={'candidate': law, 'clock_port': CLOCK_PORT,
                 'outcomes': [list(mark) for mark in marks()],
                 'tick': 'actual nonzero occupation change at clock port',
                 'zero_tick_contains_remote_repairs': True,
                 'scheduler': 'uniform thirty source seams and fair native coin',
                 'hidden_source_labels_probability_summed': True})
    return NativeClock(law, tuple(groups), instrument)


def source(law='C'):
    cached = _source(quantum.candidate(law))
    instrument = runtime.Instrument(kraus=cached.instrument.kraus,
        visible=cached.instrument.visible, resources=cached.instrument.resources,
        binding=cached.instrument.binding)
    return NativeClock(cached.law, cached.groups, instrument)


def sparse(state):
    if type(state) is not dict or any(
            type(key) is not tuple or len(key) != 2
            or any(type(i) is not int or not 0 <= i < 66 for i in key)
            or type(value) is not F for key, value in state.items()):
        raise ValueError('native matrix units with exact real Fraction entries required')
    return {key: value for key, value in state.items() if value}


def trace(state):
    return sum((value for (i, j), value in sparse(state).items() if i == j), F())


def dense(state):
    result = np.zeros((66, 66), complex)
    for key, value in sparse(state).items():
        result[key] = float(value)
    return result


def exact_step(law, state, *, activity=F(1)):
    """All eleven complete outcomes, without normalizing successful states."""
    law = quantum.candidate(law)
    state = sparse(state)
    if type(activity) is not F or not 0 < activity <= 1:
        raise ValueError('a positive Fraction activity at most one is required')
    outputs = [defaultdict(F) for _ in marks()]
    for key, value in state.items():
        outputs[0][key] = (1-activity)*value
    for mark, _, _, _, pairs in _source(law).groups:
        mapping = dict(pairs)
        for (i, j), value in state.items():
            if i in mapping and j in mapping:
                outputs[mark][mapping[i], mapping[j]] += activity*value/60
    return tuple(sparse(dict(out)) for out in outputs)


@lru_cache(maxsize=1)
def sectors():
    groups = defaultdict(list)
    for i, a in enumerate(quantum.STATES):
        for j, b in enumerate(quantum.STATES):
            groups[tuple(x-y for x, y in zip(a, b, strict=True))].append((i, j))
    return MappingProxyType({delta: tuple(pairs) for delta, pairs in groups.items()})


@lru_cache(maxsize=16, typed=True)
def block_counts(law, difference=ZERO, bias=False):
    """Exact native block with all local marks; denominator60 (or62 biased)."""
    law = quantum.candidate(law)
    if (type(difference) is not tuple or len(difference) != 12
            or any(type(v) is not int for v in difference)
            or difference not in sectors() or type(bias) is not bool):
        raise ValueError('an attained signed-difference sector and Boolean bias required')
    pairs = sectors()[difference]
    index = {pair: k for k, pair in enumerate(pairs)}
    arrays = [[[0]*len(pairs) for _ in pairs] for _ in marks()]
    for mark, edge, _, _, entries in _source(law).groups:
        mapping = dict(entries)
        weight = 2 if bias and edge == 0 else 1
        for column, (i, j) in enumerate(pairs):
            if i in mapping and j in mapping:
                image = mapping[i], mapping[j]
                if image not in index:
                    raise ArithmeticError('local clock source map changed signed-difference sector')
                arrays[mark][index[image]][column] += weight
    return tuple(tuple(tuple(row) for row in matrix) for matrix in arrays)


@lru_cache(maxsize=16, typed=True)
def resolvent(law, difference=ZERO, bias=False):
    """Immutable exact (dI-Q_count)^(-1); not a fractional channel power."""
    counts = block_counts(law, difference, bias)
    denominator = 62 if bias else 60
    matrix = denominator*sp.eye(len(counts[0]))-sp.Matrix(counts[0])
    inverse = matrix.inv(method='DM')
    if matrix*inverse != sp.eye(matrix.rows):
        raise ArithmeticError('exact local-clock resolvent inverse failed')
    return tuple(tuple(F(value) for value in row) for row in inverse.tolist())


@lru_cache(maxsize=128, typed=True)
def event_block(law, difference=ZERO, mark=None, bias=False):
    """R_mark (dI-Q_count)^(-1), preserving the source operation order."""
    if mark is not None and (type(mark) is not int or not 1 <= mark < len(marks())):
        raise ValueError('a positive local clock outcome index required')
    matrices = block_counts(law, difference, bias)
    moving = sp.zeros(len(matrices[0]))
    for chosen in range(1, len(marks())) if mark is None else (mark,):
        moving += sp.Matrix(matrices[chosen])
    result = moving*sp.Matrix(resolvent(law, difference, bias))
    return tuple(tuple(F(value) for value in row) for row in result.tolist())


def event_apply(law, state, *, mark=None, bias=False):
    """Exact eventual event CP series on only the populated difference blocks."""
    law = quantum.candidate(law)
    if (type(bias) is not bool or mark is not None
            and (type(mark) is not int or not 1 <= mark < len(marks()))):
        raise ValueError('a local event mark and Boolean bias required')
    state = sparse(state)
    groups = defaultdict(dict)
    for key, value in state.items():
        i, j = key
        delta = tuple(a-b for a, b in zip(quantum.STATES[i], quantum.STATES[j], strict=True))
        groups[delta][key] = value
    result = {}
    for delta, values in groups.items():
        pairs = sectors()[delta]
        matrix = event_block(law, delta, mark, bias)
        vector = tuple(values.get(key, F()) for key in pairs)
        for key, row in zip(pairs, matrix, strict=True):
            value = sum((entry*coefficient for entry, coefficient in zip(row, vector, strict=True)
                         if entry and coefficient), F())
            if value:
                result[key] = value
    return result


def diagonal(vector):
    return {(i, i): F(value) for i, value in enumerate(vector) if value}


def source_palm():
    """Actual one-attempt local-movement herald from the supplied uniform state."""
    uniform = diagonal((F(1, 66),)*66)
    outcomes = exact_step('C', uniform)
    success = defaultdict(F)
    for state in outcomes[1:]:
        for key, value in state.items():
            success[key] += value
    rate = trace(dict(success))
    palm = {key: value/rate for key, value in success.items() if value}
    if rate != F(5, 198) or event_apply('C', palm) != palm:
        raise ArithmeticError('source local-event Palm preparation failed')
    return outcomes[0], dict(success), palm


def tail_certificate(horizons=(1, 2, 4, 16, 64, 128)):
    """Exact all-input/reference tail effects of the no-local-tick CP map."""
    if not horizons or any(type(n) is not int or not 1 <= n <= 256 for n in horizons):
        raise ValueError('bounded positive tail horizons required')
    no = block_counts('C')[0]
    # Every native Kraus branch is a partial injection. A history's effect
    # is diagonal, so these basis survival probabilities certify the effect
    # operator and every external reference, not only diagonal inputs.
    survival = (1,)*66
    rows = []
    for n in range(1, max(horizons)+1):
        survival = tuple(sum(no[j][i]*survival[j] for j in range(66)) for i in range(66))
        if n in horizons:
            rows.append({'attempt_cap': n, 'maximum_timeout_probability': str(F(max(survival), 60**n)),
                         'minimum_timeout_probability': str(F(min(survival), 60**n))})
    if F(rows[-1]['maximum_timeout_probability']) >= 1:
        raise ArithmeticError('local source clock has no uniform eventual-tick certificate')
    return rows


def fixed_intervention(law='C'):
    """One extra remote native attempt, with signed current retained."""
    groups = [[], [], []]
    chosen = quantum.mixing.seams().index(intervention_seam())
    a, _ = intervention_seam()
    for _, edge, _, label, entries in _source(quantum.candidate(law)).groups:
        if edge != chosen:
            continue
        currents = {quantum.STATES[i][a]-quantum.STATES[j][a] for i, j in entries}
        if len(currents) != 1:
            raise ArithmeticError('intervention label does not determine its current')
        current = currents.pop()
        matrix = np.zeros((66, 66), complex)
        for i, j in entries:
            matrix[j, i] = 1/math.sqrt(2)
        groups[current+1].append(matrix)
    return runtime.Instrument(kraus=tuple(tuple(group) for group in groups), visible=True,
        resources=(('native_attempts', 1), ('intervention_current_reads', 1)),
        binding={'designated_source_seam': list(intervention_seam()), 'current_outcomes': [-1, 0, 1],
                 'event_trigger': 'immediately after first local clock tick',
                 'source_seam_is_disjoint_from_clock_port': True})


def exact_intervention(law, state):
    outputs = [defaultdict(F) for _ in range(3)]
    chosen = quantum.mixing.seams().index(intervention_seam())
    a, _ = intervention_seam()
    for _, edge, _, _, entries in _source(quantum.candidate(law)).groups:
        if edge != chosen:
            continue
        mapping = dict(entries)
        current = next(quantum.STATES[i][a]-quantum.STATES[j][a] for i, j in entries)
        for (i, j), value in sparse(state).items():
            if i in mapping and j in mapping:
                outputs[current+1][mapping[i], mapping[j]] += value/2
    return tuple(sparse(dict(out)) for out in outputs)
