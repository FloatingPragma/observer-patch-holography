#!/usr/bin/env python3
"""Exact finite controls for the fixed-unitary asymptotic-limit no-go.

This is not a numerical proof of the general Lean theorem.  It is a small,
dependency-free adversarial replay on integer orthogonal matrices.  It checks
the algebraic mechanism behind the theorem and two scope controls: the
identity evolution has zero displacement, and the relative comparison of
identical nontrivial evolutions is constant.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from typing import Sequence

Matrix = tuple[tuple[int, ...], ...]


def _matrix(rows: Sequence[Sequence[int]]) -> Matrix:
    if any(not isinstance(value, int) for row in rows for value in row):
        raise ValueError("matrix entries must be exact integers")
    result = tuple(tuple(value for value in row) for row in rows)
    if not result or any(len(row) != len(result) for row in result):
        raise ValueError("matrix must be nonempty and square")
    return result


def identity(dimension: int) -> Matrix:
    if dimension <= 0:
        raise ValueError("dimension must be positive")
    return tuple(
        tuple(int(row == column) for column in range(dimension))
        for row in range(dimension)
    )


def transpose(matrix: Matrix) -> Matrix:
    return tuple(tuple(matrix[column][row] for column in range(len(matrix)))
                 for row in range(len(matrix)))


def multiply(left: Matrix, right: Matrix) -> Matrix:
    if len(left) != len(right):
        raise ValueError("matrix dimensions differ")
    dimension = len(left)
    return tuple(
        tuple(
            sum(left[row][inner] * right[inner][column]
                for inner in range(dimension))
            for column in range(dimension)
        )
        for row in range(dimension)
    )


def subtract(left: Matrix, right: Matrix) -> Matrix:
    if len(left) != len(right):
        raise ValueError("matrix dimensions differ")
    return tuple(
        tuple(a - b for a, b in zip(left_row, right_row, strict=True))
        for left_row, right_row in zip(left, right, strict=True)
    )


def power(matrix: Matrix, exponent: int) -> Matrix:
    if exponent < 0:
        raise ValueError("only natural powers are used in this replay")
    result = identity(len(matrix))
    factor = matrix
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = multiply(result, factor)
        factor = multiply(factor, factor)
        remaining >>= 1
    return result


def frobenius_squared(matrix: Matrix) -> int:
    return sum(value * value for row in matrix for value in row)


def is_orthogonal(matrix: Matrix) -> bool:
    return multiply(transpose(matrix), matrix) == identity(len(matrix))


def relative_comparison(left: Matrix, right: Matrix, exponent: int) -> Matrix:
    """Return `(left^n)^T right^n` for real orthogonal controls."""

    return multiply(transpose(power(left, exponent)), power(right, exponent))


@dataclass(frozen=True)
class ExactReplay:
    name: str
    dimension: int
    nonidentity: bool
    one_step_displacement_squared: int
    adjacent_displacements_squared: tuple[int, ...]
    factorization_holds: bool
    identical_relative_comparison_is_constant: bool


def audit_step(name: str, step_rows: Sequence[Sequence[int]], horizon: int = 12) -> ExactReplay:
    """Audit a finite real-unitary (orthogonal) one-step evolution exactly."""

    if horizon <= 0:
        raise ValueError("horizon must be positive")
    step = _matrix(step_rows)
    if not is_orthogonal(step):
        raise ValueError("step must be an exact orthogonal matrix")

    unit = identity(len(step))
    generator_gap = subtract(step, unit)
    gaps: list[int] = []
    factorization_holds = True
    relative_constant = True

    for exponent in range(horizon):
        current = power(step, exponent)
        following = power(step, exponent + 1)
        adjacent_gap = subtract(following, current)
        factored_gap = multiply(current, generator_gap)
        factorization_holds &= adjacent_gap == factored_gap
        gaps.append(frobenius_squared(adjacent_gap))
        relative_constant &= relative_comparison(step, step, exponent) == unit

    return ExactReplay(
        name=name,
        dimension=len(step),
        nonidentity=step != unit,
        one_step_displacement_squared=frobenius_squared(generator_gap),
        adjacent_displacements_squared=tuple(gaps),
        factorization_holds=factorization_holds,
        identical_relative_comparison_is_constant=relative_constant,
    )


def default_replays(horizon: int = 12) -> tuple[ExactReplay, ...]:
    """Identity boundary plus two nontrivial exact-unitary attacks."""

    return (
        audit_step("identity", ((1, 0), (0, 1)), horizon),
        audit_step("swap", ((0, 1), (1, 0)), horizon),
        audit_step("quarter_turn", ((0, -1), (1, 0)), horizon),
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--horizon", type=int, default=12)
    args = parser.parse_args()
    payload = {
        "claim_boundary": (
            "Exact finite controls only; the general topological-group theorem "
            "is proved in Lean. Relative/Moller limits are not excluded."
        ),
        "replays": [asdict(replay) for replay in default_replays(args.horizon)],
    }
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
