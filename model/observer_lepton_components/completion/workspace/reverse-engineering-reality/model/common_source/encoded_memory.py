"""Source-encoded logical memory corrected after complete native observations.

The code is the four-support logical qubit of RECORD_RECOVERY equation (8).
Each committed native record selects a full 512-dimensional source unitary.
It restores the code's branch image without changing or coherently erasing
the classical archive. Raw aborts and numerical code leakage remain present.
This protects the declared code against this instrument, not arbitrary noise.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
import math

import numpy as np

from ..observer.source_discrimination import circuit
from ..observer.source_scalar import compiler
from . import native_routing as routing


LOGICAL_PAIRS = (((0, 2), (1, 3)), ((0, 3), (1, 2)))
SUPPORTS = tuple(tuple(tuple(int(port in pair) for port in range(12))
                       for pair in logical) for logical in LOGICAL_PAIRS)
ENCODER_ASSIGNMENTS = ((0, 2), (2, 73), (1, 74), (3, 1))


def _label(law, seam, coin, raw):
    law, seam = routing._candidate(law), routing._seam(seam)
    if type(coin) is not int or coin not in (0, 1):
        raise ValueError('an integer native coin zero or one required')
    if type(raw) is not int or not 0 <= raw < 1 << circuit.RECORD_BITS[law]:
        raise ValueError('a raw integer record in the complete native menu required')
    return law, seam, coin, raw


def _admitted(law, coin, raw):
    current = {0: 0, 1: 1, 2: -1}.get(raw % 4)
    return current is not None and (law == 'C' or abs(2*current+(2*coin-1)*(raw//4)) <= 2)


def branch_images(law, seam, coin, raw):
    """Exact basis supports of this code's native branch, before correction.

    This table is derived from the supplied endpoint law, not source execution
    or observed probabilities. Empty supports retain their raw output labels.
    """
    law, seam, coin, raw = _label(law, seam, coin, raw)
    images = [[], []]
    a, b = seam
    for logical, states in enumerate(SUPPORTS):
        for state in states:
            total = state[a]+state[b]
            left = total//2+coin*(total % 2)
            current = state[a]-left
            code = {0: 0, 1: 1, -1: 2}[current]
            if law == 'E':
                code += 4*(total % 2)
            if code == raw:
                end = list(state)
                end[a], end[b] = left, total-left
                images[logical].append(circuit.encode(tuple(end)))
    if (len(images[0]) != len(images[1]) or set(images[0]).intersection(images[1])
            or any(len(set(row)) != len(row) for row in images)):
        raise ArithmeticError('recorded native code branch violates the logical Gram identity')
    return tuple(tuple(row) for row in images)


def branch_weight(law, seam, coin, raw):
    """Exact code-input outcome probability, including coin and code amplitudes."""
    return F(len(branch_images(law, seam, coin, raw)[0]), 4)


def _encoder_permutation():
    return routing.complete_permutation(ENCODER_ASSIGNMENTS)


def _branch_permutation(law, seam, coin, raw):
    images = branch_images(law, seam, coin, raw)
    assignments = tuple((logical+2*term, image)
                        for logical, row in enumerate(images) for term, image in enumerate(row))
    return routing.complete_permutation(assignments)


def _mixing_word(inverse=False):
    return [compiler.pulse('Y', math.pi/4 if inverse else -math.pi/4, [1])]


@lru_cache(maxsize=2)
def _encoder(inverse):
    permutation = routing.permutation_unitary(_encoder_permutation(), inverse=inverse)
    if inverse:
        result = compiler.execute(_mixing_word(True), 9, permutation)
    else:
        mixed = compiler.execute(_mixing_word(), 9)
        result = permutation@mixed
    if not np.isfinite(result).all():
        raise ArithmeticError('nonfinite source encoder execution')
    return routing._immutable(result)


def code_matrix():
    """Actual source encoder on both logical input columns; no target insertion."""
    return _encoder(False)[:, :2]


def ideal_code_matrix():
    """Independent declared logical code for diagnostics, not state preparation."""
    result = np.zeros((512, 2), complex)
    for logical, states in enumerate(SUPPORTS):
        for state in states:
            result[circuit.encode(state), logical] = 1/math.sqrt(2)
    return routing._immutable(result)


def _logical_state(value):
    routing._original_scalars(value)
    raw = np.asarray(value)
    if (raw.ndim not in (1, 2) or raw.shape[0] != 2
            or (raw.ndim == 2 and raw.shape[1] < 1)
            or raw.dtype.kind not in 'iufc' or not np.isfinite(raw).all()):
        raise ValueError('finite two-component logical amplitudes or reference columns required')
    result = np.array(raw, dtype=complex, copy=True)
    if not np.isfinite(result).all():
        raise ValueError('logical amplitudes must remain finite in the execution representation')
    return result


def _finite_output(value):
    if not np.isfinite(value).all():
        raise ArithmeticError('native memory operation produced nonfinite amplitudes')
    return routing._immutable(value)


def encoding(state):
    """Attach eight blank code qubits and execute the full source encoder.

    Logical/reference amplitudes may be unnormalized. The input logical
    qubit is supplied; this operation does not prepare its unknown state.
    """
    value = _logical_state(state)
    return _finite_output(code_matrix()@value)


def decoding(state):
    """Execute the source inverse, retaining every syndrome and padding row."""
    value = routing._state(state)
    return _finite_output(_encoder(True)@value)


@lru_cache(maxsize=96)
def _recovery(law, seam, coin, raw):
    count = len(branch_images(law, seam, coin, raw)[0])
    if not count:
        # A zero-code-probability but admitted record still has a complete
        # specified action on any numerical leakage or off-code input.
        return routing._immutable(compiler.execute([], 9))
    result = routing.permutation_unitary(_branch_permutation(law, seam, coin, raw), inverse=True)
    if count == 2:
        result = compiler.execute(_mixing_word(True), 9, result)
    return _finite_output(_encoder(False)@result)


def recovery_unitary(law, seam, coin, raw):
    """Record-selected unitary on all 512 dimensions; invalid raw labels abort."""
    law, seam, coin, raw = _label(law, seam, coin, raw)
    if not _admitted(law, coin, raw):
        raise ValueError('invalid raw native records halt instead of selecting recovery')
    return _recovery(law, seam, coin, raw)


def _counts(*items):
    return {key: sum(item[key] for item in items) for key in routing.COUNT_KEYS}


def encoding_resources():
    counts = _counts(routing.permutation_resources(_encoder_permutation()),
                     compiler.resources(_mixing_word()))
    return {**counts, 'supplied_logical_qubits': 1, 'blank_code_qubits': 8,
            'new_classical_records': 0, 'physical_initialization_error_included': False}


def decoding_resources():
    counts = _counts(routing.permutation_resources(_encoder_permutation()),
                     compiler.resources(_mixing_word(True)))
    return {**counts, 'blank_code_qubits': 0, 'new_classical_records': 0,
            'output_rows_retained': 512, 'syndrome_postselection_applied': False}


def recovery_resources(law, seam, coin, raw):
    law, seam, coin, raw = _label(law, seam, coin, raw)
    count = len(branch_images(law, seam, coin, raw)[0])
    admitted = _admitted(law, coin, raw)
    counts = dict.fromkeys(routing.COUNT_KEYS, 0)
    if admitted and count:
        counts = _counts(encoding_resources(),
                         routing.permutation_resources(_branch_permutation(law, seam, coin, raw)))
        if count == 2:
            counts = _counts(counts, compiler.resources(_mixing_word(True)))
    return {**counts, 'code_branch_weight_exact': str(F(count, 4)),
            'raw_abort': not admitted, 'identity_on_valid_zero_weight_branch': admitted and not count,
            'blank_code_qubits': 0, 'new_classical_records': 0,
            'archive_record_rewrites': 0, 'input_dependent_normalization': False}


def repair_branches(law, seam, state, *, recover=True):
    """Observe and optionally actively correct, retaining every raw branch.

    Input may occupy all 512 rows with arbitrary reference columns. Only
    encoded inputs inherit the logical-memory theorem. No code projection,
    syndrome selection, state normalization or quantum reset is performed.
    """
    law, seam = routing._candidate(law), routing._seam(seam)
    if type(recover) is not bool:
        raise ValueError('active recovery must be an explicit Boolean')
    value = routing._state(state)
    branches = routing.routed_repair_branches(law, seam, value)
    for branch in branches:
        coin, raw = branch['coin'], branch['record_code']
        if branch['abort'] != (not _admitted(law, coin, raw)):
            raise ArithmeticError('source abort flag differs from independent raw-record law')
        applied = recover and not branch['abort']
        if applied:
            branch['amplitude'] = _finite_output(recovery_unitary(law, seam, coin, raw)@branch['amplitude'])
        branch['code_branch_weight'] = branch_weight(law, seam, coin, raw)
        branch['recovery_applied'] = applied
        branch['recovery_resources'] = (recovery_resources(law, seam, coin, raw) if applied
                                        else {**dict.fromkeys(routing.COUNT_KEYS, 0),
                                              'blank_code_qubits': 0, 'new_classical_records': 0,
                                              'archive_record_rewrites': 0})
    return branches
