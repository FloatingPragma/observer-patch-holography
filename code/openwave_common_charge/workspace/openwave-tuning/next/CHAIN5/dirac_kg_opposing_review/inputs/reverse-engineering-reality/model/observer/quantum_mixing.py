"""Exact signed-difference blocks of the uniform native C/E channels.

The seam and placement coin are independently uniform. Their records are
forgotten for this channel calculation; no claim is made about a growing
retained record tape, the E/F controller, physical time or source execution.
All public matrices use output rows and input columns of matrix units.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json

from . import native


STATES = native.states()
INDEX = {state: i for i, state in enumerate(STATES)}
BINARY = tuple(max(state) == 1 for state in STATES)
ZERO = (0,) * 12


@dataclass(frozen=True)
class Sector:
    difference: tuple[int, ...]
    pairs: tuple[tuple[int, int], ...]


def _candidate(value):
    if type(value) is not str or value not in ("C", "E"):
        raise ValueError("candidate C or E required")
    return value


def _binary(value):
    if type(value) is not bool:
        raise ValueError("binary-shell switch must be Boolean")
    return value


@lru_cache(maxsize=1)
def seams():
    """Reuse the existing hash-bound carrier, without selecting a new chart."""
    from ..common_source.slow_geometry import source
    result = tuple(tuple(edge) for edge in source().seams())
    if len(result) != 30 or len(set(result)) != 30:
        raise ArithmeticError("the declared carrier must have thirty distinct seams")
    return result


@lru_cache(maxsize=2)
def _sectors(binary):
    groups = defaultdict(list)
    for i, left in enumerate(STATES):
        for j, right in enumerate(STATES):
            if not binary or (BINARY[i] and BINARY[j]):
                groups[tuple(a-b for a, b in zip(left, right, strict=True))].append((i, j))
    return tuple(Sector(delta, tuple(pairs)) for delta, pairs in sorted(groups.items()))


def sectors(*, binary=False):
    """Immutable partition of all 6,084 (or 4,356 binary) matrix units."""
    return _sectors(_binary(binary))


@lru_cache(maxsize=2)
def _lookup(binary):
    return {sector.difference: sector for sector in _sectors(binary)}


def _sector(difference, binary):
    _binary(binary)
    if (type(difference) is not tuple or len(difference) != 12
            or any(type(x) is not int for x in difference)
            or difference not in _lookup(binary)):
        raise ValueError("an attained twelve-integer signed load difference required")
    return _lookup(binary)[difference]


@lru_cache(maxsize=1)
def _images():
    # Native primitive calls, not an assumed swap or target superoperator.
    images = []
    for a, b in seams():
        for coin in (0, 1):
            event = []
            for state in STATES:
                output, record = native.native(state, (a, b), coin)
                if (type(record) is not tuple or len(record) != 2
                        or any(type(value) is not int for value in record)
                        or record != (coin, state[a]-output[a])):
                    raise ArithmeticError("native repair record violates the coin/current grammar")
                event.append((INDEX[output], record[1], state[a]-state[b]))
            images.append(tuple(event))
    return tuple(images)


def native_block(candidate, difference, *, binary=False):
    """Assemble one exact block from equality of actual Kraus records."""
    candidate = _candidate(candidate)
    sector = _sector(difference, binary)
    positions = {pair: i for i, pair in enumerate(sector.pairs)}
    counts = [[0] * len(positions) for _ in positions]
    record_index = 1 if candidate == "C" else 2
    for column, (i, j) in enumerate(sector.pairs):
        for event in _images():
            left, right = event[i], event[j]
            if left[record_index] == right[record_index]:
                pair = left[0], right[0]
                if pair not in positions:
                    raise ArithmeticError("a Kraus branch changed signed difference or shell")
                counts[positions[pair]][column] += 1
    return tuple(tuple(F(value, 60) for value in row) for row in counts)


def killing_rate(candidate, difference):
    """Binary-sector loss from seams touching unequal signed differences."""
    candidate = _candidate(candidate)
    _sector(difference, True)
    if candidate == "C":
        return F(sum(abs(difference[a]-difference[b]) for a, b in seams()), 60)
    return F(sum(difference[a] != difference[b] for a, b in seams()), 30)


def binary_reduction(candidate, difference):
    """(1-killing) I minus the induced exclusion Laplacian divided by 60.

    The occupied positive/negative difference sites are fixed. Common
    particles swap only along seams inside the zero-difference vertices.
    This construction never calls the native primitive or its Kraus images.
    """
    sector = _sector(difference, True)
    loss = killing_rate(candidate, difference)
    positions = {pair: i for i, pair in enumerate(sector.pairs)}
    out = [[F(0)] * len(positions) for _ in positions]
    for column, (i, j) in enumerate(sector.pairs):
        out[column][column] = 1-loss
        for a, b in seams():
            if difference[a] or difference[b] or STATES[i][a] == STATES[i][b]:
                continue
            left, right = list(STATES[i]), list(STATES[j])
            left[a], left[b] = left[b], left[a]
            right[a], right[b] = right[b], right[a]
            row = positions[INDEX[tuple(left)], INDEX[tuple(right)]]
            out[column][column] -= F(1, 60)
            out[row][column] += F(1, 60)
    return tuple(map(tuple, out))


def audit():
    """Check all exact blocks, their triangular extension and decay bounds.

    The spectral-gap proof additionally uses the connected-graph exclusion
    gap theorem on the diagonal block; this audit does not replace that
    theorem with floating diagonalization or assert physical calibration.
    """
    shell = _lookup(True)
    reports = {}
    for candidate in ("C", "E"):
        digest = hashlib.sha256()
        binary_max = transient_max = F(0)
        transient_offdiagonal_max = F(0)
        gershgorin_min = F(1)
        transient_units = transfer_entries = 0
        for sector in sectors():
            delta, pairs = sector.difference, sector.pairs
            actual = native_block(candidate, delta)
            digest.update(json.dumps([delta, [[str(x) for x in row] for row in actual]],
                                     separators=(",", ":")).encode())
            inside = [i for i, (a, b) in enumerate(pairs) if BINARY[a] and BINARY[b]]
            outside = [i for i in range(len(pairs)) if i not in inside]
            transient_units += len(outside)
            if inside:
                expected = binary_reduction(candidate, delta)
                restricted = tuple(tuple(actual[i][j] for j in inside) for i in inside)
                if tuple(pairs[i] for i in inside) != shell[delta].pairs or restricted != expected:
                    raise ArithmeticError("native Kraus block differs from induced exclusion reduction")
                if expected != tuple(zip(*expected, strict=True)):
                    raise ArithmeticError("a binary block is not self-adjoint")
                for i, row in enumerate(expected):
                    lower = row[i]-sum(abs(x) for j, x in enumerate(row) if j != i)
                    gershgorin_min = min(gershgorin_min, lower)
                    if sum(row, F(0)) != 1-killing_rate(candidate, delta):
                        raise ArithmeticError("the constant sector vector has the wrong decay")
                if delta != ZERO:
                    binary_max = max(binary_max, 1-killing_rate(candidate, delta))
            for i in outside:
                if any(actual[i][j] for j in inside):
                    raise ArithmeticError("the binary algebra leaks into transient matrix units")
                if any(actual[i][j] for j in outside if j != i):
                    raise ArithmeticError("the transient block is not diagonal")
                transient_max = max(transient_max, actual[i][i])
                if delta != ZERO:
                    transient_offdiagonal_max = max(transient_offdiagonal_max, actual[i][i])
            transfer_entries += sum(bool(actual[i][j]) for i in inside for j in outside)
        expected_max = F(5, 6) if candidate == "C" else F(7, 10)
        if binary_max != expected_max or transient_max != F(5, 6) or gershgorin_min < 0:
            raise ArithmeticError("the fixed carrier's exact decay/positivity controls failed")
        reports[candidate] = {
            "block_sha256": digest.hexdigest(),
            "binary_offdiagonal_spectral_radius": str(binary_max),
            "binary_gershgorin_lower_bound": str(gershgorin_min),
            "transient_spectral_radius": str(transient_max),
            "offdiagonal_transient_spectral_radius": str(transient_offdiagonal_max),
            "transient_matrix_units": transient_units,
            "nonzero_transient_to_binary_entries": transfer_entries,
            "all_blocks_exact": True,
        }
    return {
        "schema": "oph.observer.quantum_mixing_blocks.v1",
        "full_matrix_units": sum(len(s.pairs) for s in sectors()),
        "binary_matrix_units": sum(len(s.pairs) for s in sectors(binary=True)),
        "binary_diagonal_dimension": len(shell[ZERO].pairs),
        "binary_centered_diagonal_dimension": len(shell[ZERO].pairs)-1,
        "binary_offdiagonal_dimension": sum(len(s.pairs) for s in sectors(binary=True)
                                            if s.difference != ZERO),
        "full_block_sizes": {str(k): v for k, v in sorted(Counter(len(s.pairs) for s in sectors()).items())},
        "binary_block_sizes": {str(k): v for k, v in sorted(Counter(len(s.pairs) for s in sectors(binary=True)).items())},
        "candidates": reports,
        "assumptions": {"seams": "independent uniform choices on the existing thirty-seam carrier",
                        "coins": "independent fair placement coins",
                        "channel": "classical seam, coin and repair records forgotten after each event",
                        "sector": "existing nonnegative integer total-two carrier"},
        "new_source_execution_or_physical_clock": False,
    }
