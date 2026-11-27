"""Recover what the scalar graph's encoded even sector identifies."""

from dataclasses import dataclass
from math import isfinite
from numbers import Real

import numpy as np

from .. import compiler, instruments


@dataclass
class OneBodyReadout:
    matrix: np.ndarray | None
    trace: float
    pauli_residual: float


def _decode_on_cycles(graph, terms, cycles):
    def representative(word):
        # A fundamental cycle contains its own chord and no other chord.
        # Right multiplication fixes the action on every cycle-+1 input,
        # including terms that send that input outside the prepared sector.
        for chord, _, stabilizer in cycles:
            if word[1] & (1 << chord):
                word = compiler.mul(word, stabilizer)
        return instruments.canonical(word)

    reduced = {}
    for word, value in terms.items():
        key, phase = representative(word)
        reduced[key] = reduced.get(key, 0) + phase * value

    # Distinct compiler words stay distinct on the code for n >= 3. Two
    # different graph edges have X masks differing on two edges, which cannot
    # form a cycle of a simple graph. With equal X masks, a nontrivial cycle
    # always changes at least one chord, so different Z masks stay distinct.
    words = {(0, 0, 0), *graph["B"]}
    for edge in range(graph["M"]):
        for coefficient in (1, 1j):
            words.update(word for _, word in compiler.hopping_words(graph, edge, coefficient))
    lifted = {}
    for word in words:
        key, phase = representative(word)
        if key in reduced:
            lifted[word] = (reduced.pop(key) * phase.conjugate()).real
    matrix, residual = instruments.decode_paulis(graph, lifted)
    return matrix, max(residual, max(map(abs, reduced.values()), default=0.))


def read_one_body(graph, terms):
    """Read an encoded number-conserving one-body operator on a connected graph.

    The residual is the largest Pauli coefficient discrepancy on cycle-+1
    inputs. On two modes only the trace is identifiable.
    """
    if graph["n"] != 2:
        matrix, residual = instruments.decode_paulis(graph, terms)
        if residual and graph["M"] >= graph["n"]:
            cycles, _ = compiler.fundamental_cycles(graph)
            matrix, residual = _decode_on_cycles(graph, terms, cycles)
        trace = float(np.trace(matrix).real)
        if not np.isfinite(matrix).all() or not isfinite(trace) or not isfinite(residual):
            raise ValueError("The one-body readout exceeds floating-point range")
        return OneBodyReadout(matrix, trace, residual)

    # The single encoded qubit has vacuum at Z=+1 and double occupancy
    # at Z=-1. A valid one-body operator there is trace(h) * (I-Z)/2.
    identity, z = (0, 0, 0), (0, 0, 1)
    x, y = (0, 1, 0), (1, 1, 1)
    if not isinstance(terms, dict) or any(
        not isinstance(word, tuple) or not all(type(v) is int for v in word)
        or word not in (identity, x, y, z)
        or not isinstance(value, Real) or not isfinite(value)
        for word, value in terms.items()
    ):
        raise ValueError("Expected finite real coefficients on canonical single-qubit Paulis")
    scalar = float(terms.get(identity, 0))
    diagonal = float(terms.get(z, 0))
    trace = scalar - diagonal
    if not isfinite(trace):
        raise ValueError("The identifiable two-mode trace exceeds floating-point range")
    residual = max(abs(scalar - trace/2), abs(diagonal + trace/2),
                   abs(terms.get(x, 0)), abs(terms.get(y, 0)))
    return OneBodyReadout(None, trace, float(residual))
