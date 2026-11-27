"""Recover what the scalar graph's encoded even sector identifies."""

from dataclasses import dataclass
from math import isfinite
from numbers import Real

import numpy as np

from .. import instruments


@dataclass
class OneBodyReadout:
    matrix: np.ndarray | None
    trace: float
    pauli_residual: float


def read_one_body(graph, terms):
    """Read an encoded number-conserving one-body operator on a connected graph.

    The residual is the largest Pauli coefficient discrepancy. On two modes
    only the trace is identifiable and the returned matrix is None.
    """
    if graph["n"] != 2:
        matrix, residual = instruments.decode_paulis(graph, terms)
        return OneBodyReadout(matrix, float(np.trace(matrix).real), residual)

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
