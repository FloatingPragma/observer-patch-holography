"""Exact Gaussian-rational execution of the original binary source arrays.

The work register is the complete 512-address chart. Native repair is extended
by a retained code-complement abort. No amplitude, failed record or roundoff
component is projected away or normalized. The old ideal experiment is an
independent oracle, never an execution input.
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from sympy.polys.domains import QQ_I

from ..observer import native, quantum_mixing
from ..quantum_observer_correspondence import construction as ideal
from . import source
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

ZERO = QQ_I.zero
CODE = frozenset(source.codes())
MAX_CAP = 3
_GAUSSIAN = type(QQ_I(1))
_RATIONAL = type(QQ_I(1).x)
_INTEGER = type(QQ_I(1).x.numerator)


def rational(value):
    """Exact value of one validated binary64 source component."""
    return QQ_I(F(float(value.real)), F(float(value.imag)))


def conjugate(value):
    _gaussian(value)
    return QQ_I(value.x, -value.y)


def fraction(value):
    if type(value) is not _RATIONAL:
        raise ValueError('a canonical exact rational component is required')
    try:
        numerator, denominator = value.numerator, value.denominator
    except AttributeError as exc:
        raise ValueError('a canonical exact rational component is required') from exc
    if (type(numerator) is not _INTEGER
            or type(denominator) is not _INTEGER or denominator <= 0
            or math.gcd(int(numerator), int(denominator)) != 1):
        raise ValueError('canonical exact integer numerator and positive denominator required')
    return F(int(value.numerator), int(value.denominator))


def _gaussian(value):
    if type(value) is not _GAUSSIAN:
        raise ValueError('an original exact Gaussian-rational value is required')
    try:
        fraction(value.x)
        fraction(value.y)
    except AttributeError as exc:
        raise ValueError('complete exact Gaussian-rational components required') from exc


def _clean(state):
    return {key: value for key, value in state.items() if value}


def add(left, right):
    result = dict(left)
    for key, value in right.items():
        result[key] = result.get(key, ZERO) + value
    return _clean(result)


def trace(state):
    for value in state.values():
        _gaussian(value)
    result = sum((value for (i, j, r, s), value in state.items()
                  if i == j and r == s), ZERO)
    if result.y:
        raise ArithmeticError('an exact density has a nonreal trace')
    return fraction(result.x)


def prepare(reference=False):
    matrix = source.preparation(reference)
    ket = {(i, r): rational(matrix[i, r]) for i in range(512)
           for r in range(matrix.shape[1]) if matrix[i, r] != 0}
    return {(i, j, r, s): x*conjugate(y)
            for (i, r), x in ket.items() for (j, s), y in ket.items()}


def local(matrix, state):
    """Apply all original local components with spectator identity exactly."""
    if (type(matrix) is not np.ndarray or matrix.shape != (8, 8)
            or matrix.dtype.kind not in 'fc'
            or (matrix.dtype.kind == 'f' and matrix.dtype.itemsize > 8)
            or (matrix.dtype.kind == 'c' and matrix.dtype.itemsize > 16)
            or not np.isfinite(matrix).all()):
        raise ValueError('a finite original 8 by 8 binary source array required')
    columns = {i: tuple((j, rational(matrix[j, i])) for j in range(8)
                        if matrix[j, i] != 0) for i in range(8)}
    right_columns = {i: tuple((j, conjugate(value)) for j, value in entries)
                     for i, entries in columns.items()}
    result = defaultdict(lambda: ZERO)
    for (i, j, r, s), value in state.items():
        _gaussian(value)
        for a, x in columns[i % 8]:
            for b, y in right_columns[j % 8]:
                result[8*(i//8)+a, 8*(j//8)+b, r, s] += x*value*y
    return _clean(result)


def complement(state):
    return {key: value for key, value in state.items()
            if key[0] not in CODE and key[1] not in CODE}


@lru_cache(maxsize=1)
def _native_maps():
    """Independent native endpoint enumeration; no ideal transfer is used."""
    states = source.STATES
    index = {state: code for state, code in zip(states, source.codes(), strict=True)}
    result = []
    for seam in quantum_mixing.seams():
        groups = {}
        for coin in (0, 1):
            for before in states:
                after, (_, current) = native.native(before, seam, coin)
                groups.setdefault((coin, current), {})[index[before]] = index[after]
        result.append((seam, tuple((current, tuple(sorted(mapping.items())))
                                  for (_, current), mapping in sorted(groups.items()))))
    return tuple(result)


def _mark(seam, current):
    return (tuple(int(bool(set(seam).intersection(region))) for region in ideal.REGIONS),
            tuple(current*(int(v == seam[1])-int(v == seam[0])) for v in source.UNION))


@lru_cache(maxsize=3)
def _counts(seam=None):
    remote, marked = defaultdict(lambda: defaultdict(int)), {}
    for edge, branches in _native_maps():
        if seam is not None and edge != seam:
            continue
        incident = bool(set(edge).intersection(source.UNION))
        for current, mapping in branches:
            target = marked.setdefault(_mark(edge, current), defaultdict(lambda: defaultdict(int))) if incident else remote
            for i, a in mapping:
                for j, b in mapping:
                    target[i, j][a, b] += 1
    freeze = lambda table: {key: tuple(sorted(values.items())) for key, values in table.items()}
    return freeze(remote), {key: freeze(table) for key, table in sorted(marked.items())}


def _apply(counts, state, denominator):
    result = defaultdict(lambda: ZERO)
    for (i, j, r, s), value in state.items():
        for (a, b), count in counts.get((i, j), ()):
            result[a, b, r, s] += value*QQ_I(F(count, denominator))
    return _clean(result)


def interval(state, cap):
    remote, maps = _counts()
    completed = {mark: {} for mark in maps}
    aborts = []
    waiting = state
    for _ in range(cap):
        aborts.append(complement(waiting))
        for mark, table in maps.items():
            completed[mark] = add(completed[mark], _apply(table, waiting, 60))
        waiting = _apply(remote, waiting, 60)
    return completed, waiting, tuple(aborts)


def _arguments(gamma, cap, reference):
    if type(gamma) is not int or gamma not in (0, 1):
        raise ValueError('binary integer interaction switch required')
    if type(cap) is not int or not 1 <= cap <= MAX_CAP:
        raise ValueError('integer cap between one and three required')
    if type(reference) is not bool:
        raise ValueError('Boolean reference choice required')


def first_stage(*, gamma=1, cap=2, reference=False):
    _arguments(gamma, cap, reference)
    completed, waiting, aborts = interval(prepare(reference), cap)
    return {'gamma': gamma, 'cap': cap, 'reference': reference,
            'completed': completed, 'waiting': waiting, 'aborts': aborts}


def _finish(stage):
    gamma, cap = stage['gamma'], stage['cap']
    gate, readout = source.local_apparatus(gamma), source.local_readout()
    result = {}

    def finish(word, status, state):
        for bit, pointer in enumerate(readout):
            result[word, status, bit] = local(pointer, state)

    finish((), 'timeout_first', stage['waiting'])
    for step, state in enumerate(stage['aborts']):
        finish((), f'code_abort_first_{step}', state)
    for mark, state in stage['completed'].items():
        prefix = (('scheduled', mark),)
        finish(prefix, 'code_abort_feedback', complement(state))
        seam = (0, 1) if mark[1][0] > 0 else (0, 5)
        for feedback, table in _counts(seam)[1].items():
            word = prefix+(('feedback', feedback),)
            controlled = local(gate, _apply(table, state, 2))
            completed, waiting, aborts = interval(controlled, cap)
            finish(word, 'timeout_second', waiting)
            for step, aborted in enumerate(aborts):
                finish(word, f'code_abort_second_{step}', aborted)
            for last, terminal in completed.items():
                finish(word+(('scheduled', last),), 'complete', terminal)
    return result


def execute(*, gamma=1, cap=2, reference=False, checkpoint=None):
    _arguments(gamma, cap, reference)
    stage = first_stage(gamma=gamma, cap=cap, reference=reference)
    if checkpoint is not None:
        expected = checkpoint_json(stage)
        if type(checkpoint) is not str or checkpoint != expected:
            raise ValueError('checkpoint differs from exact source-prefix replay')
        # The canonical serialized payload contains every exact state entry.
        # Equality to replay binds source, parameters, record alphabet and costs.
        # Continue from the serialized arrays after verification, not a reset
        # or the freshly recomputed prefix.
        payload = json.loads(checkpoint)
        def state(rows):
            return {tuple(key): QQ_I(F(real), F(imag)) for key, real, imag in rows}
        def mark(value):
            return tuple(tuple(part) for part in value)
        stage['completed'] = {mark(key): state(rows) for key, rows in payload['completed']}
        stage['waiting'] = state(payload['waiting'])
        stage['aborts'] = tuple(state(rows) for rows in payload['aborts'])
    result = _finish(stage)
    return validate_history(result, cap=cap, reference=reference)


def record_menu(cap=2):
    """Independent endpoint menu, including zero-mass code-complement aborts."""
    from ..quantum_observer_correspondence import clock_response as oracle
    _arguments(0, cap, False)
    menu = set()
    def finish(word, status):
        menu.update((word, status, bit) for bit in (0, 1))
    finish((), 'timeout_first')
    for step in range(cap):
        finish((), f'code_abort_first_{step}')
    for mark in oracle.mark_menu():
        if mark is None:
            continue
        prefix = (('scheduled', mark),)
        finish(prefix, 'code_abort_feedback')
        seam = (0, 1) if mark[1][0] > 0 else (0, 5)
        for feedback in oracle.mark_menu(seam):
            word = prefix+(('feedback', feedback),)
            finish(word, 'timeout_second')
            for step in range(cap):
                finish(word, f'code_abort_second_{step}')
            for last in oracle.mark_menu():
                if last is not None:
                    finish(word+(('scheduled', last),), 'complete')
    return menu


def _state_snapshot(state, reference_dimension=None):
    if type(state) is not dict:
        raise ValueError('original sparse exact state dictionaries required')
    copied = {}
    for key, value in state.items():
        if (type(key) is not tuple or len(key) != 4
                or any(type(index) is not int for index in key)
                or not 0 <= key[0] < 512 or not 0 <= key[1] < 512
                or key[2] < 0 or key[3] < 0
                or (reference_dimension is not None
                    and (key[2] >= reference_dimension or key[3] >= reference_dimension))):
            raise ValueError('original exact source/reference entries required')
        _gaussian(value)
        copied[key] = QQ_I(fraction(value.x), fraction(value.y))
    return copied


def _label(value):
    # Distance comparisons also admit simple string labels for matrix-unit
    # controls, but never dispatch equality/hash on caller-defined objects.
    if type(value) in (str, int):
        return
    if type(value) is tuple:
        for item in value:
            _label(item)
        return
    raise ValueError('original tuple/string/integer record fields required')


def _snapshot(result):
    if type(result) is not dict:
        raise ValueError('an original exact output dictionary is required')
    snapshot = {}
    for record, state in result.items():
        _label(record)
        snapshot[record] = _state_snapshot(state)
    return snapshot


def validate_history(result, *, cap=2, reference=False, native_target=False):
    """Validate original fields, then check and return a detached exact copy.

    Positivity follows from the executed CP construction, not from this
    structural validator. This is not authentication of a caller's channel.
    """
    _arguments(0, cap, reference)
    if type(native_target) is not bool:
        raise ValueError('Boolean native-target menu choice required')
    if type(result) is not dict:
        raise ValueError('the complete clock/timeout/abort record menu is required')
    for key in result:
        if (type(key) is not tuple or len(key) != 3 or type(key[0]) is not tuple
                or type(key[1]) is not str or type(key[2]) is not int):
            raise ValueError('original tuple/string/integer record fields required')
        for event in key[0]:
            if (type(event) is not tuple or len(event) != 2 or type(event[0]) is not str
                    or type(event[1]) is not tuple or len(event[1]) != 2
                    or any(type(part) is not tuple for part in event[1])
                    or any(type(item) is not int for part in event[1] for item in part)):
                raise ValueError('original integer current/clock marks required')
    # Validate every scalar before any Hermitian comparison can call equality
    # on its transpose. QQ_I and its rational slots are mutable: retain copies.
    snapshot = {record: _state_snapshot(state, 1+int(reference))
                for record, state in result.items()}
    menu = record_menu(cap)
    if native_target:
        menu = {record for record in menu if not record[1].startswith('code_abort')}
    if set(snapshot) != menu:
        raise ValueError('the complete clock/timeout/abort record menu is required')
    for state in snapshot.values():
        for key, value in state.items():
            i, j, r, s = key
            if state.get((j, i, s, r), ZERO) != conjugate(value):
                raise ValueError('exact output must be Hermitian without repair')
        if trace(state) < 0:
            raise ValueError('output record has negative mass')
    return snapshot


def _entries(state):
    return [[list(key), str(fraction(value.x)), str(fraction(value.y))]
            for key, value in sorted(state.items())]


def checkpoint_json(stage):
    gamma, cap, reference = (stage[name] for name in ('gamma', 'cap', 'reference'))
    _arguments(gamma, cap, reference)
    root = Path(__file__).resolve().parents[2]
    bindings = source.bindings()
    bindings.update(bound_pins({path: hashlib.sha256((root/path).read_bytes()).hexdigest()
                     for path in (*ideal.SOURCE_PATHS, 'model/compiled_clock_history/history.py')}))
    payload = {'schema': 'oph.compiled-clock-prefix.v1',
               'parameters': [gamma, cap, reference],
               'source_binding': bindings,
               'preparation_resources': source.resources('preparation', reference=reference),
               'scheduled_attempt_cap': cap,
               'completed': [[mark, _entries(state)] for mark, state in sorted(stage['completed'].items())],
               'waiting': _entries(stage['waiting']),
               'aborts': [_entries(state) for state in stage['aborts']]}
    return json.dumps(payload, sort_keys=True, separators=(',', ':'), allow_nan=False)


def ideal_history(*, gamma=1, cap=2, reference=False):
    """Lift the independently executed old exact history only for comparison."""
    _arguments(gamma, cap, reference)
    codes = source.codes()
    return {word: {(codes[i], codes[j], r, s): QQ_I(value)
                   for (i, j, r, s), value in state.items()}
            for word, state in ideal.execute(gamma=gamma, cap=cap, reference=reference).items()}


def distance_upper(actual, target):
    """Exact entrywise nuclear-norm bound, including all missing/abort blocks."""
    actual, target = _snapshot(actual), _snapshot(target)
    error = F(0)
    for word in set(actual) | set(target):
        a, b = actual.get(word, {}), target.get(word, {})
        for key in set(a) | set(b):
            left, right = a.get(key, ZERO), b.get(key, ZERO)
            _gaussian(left)
            _gaussian(right)
            value = left-right
            error += abs(fraction(value.x))+abs(fraction(value.y))
    return error/2


def digest(history):
    stream = hashlib.sha256()
    for word, state in sorted(history.items()):
        stream.update(json.dumps([word, _entries(state)], separators=(',', ':')).encode())
    return stream.hexdigest()
