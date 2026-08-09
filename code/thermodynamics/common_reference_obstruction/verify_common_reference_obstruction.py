#!/usr/bin/env python3
"""Independent exact-arithmetic check of the B12 common-object no-go.

This verifier deliberately reconstructs the recurrent kernel from integer
transition counts and uses only :class:`fractions.Fraction`.  It does not
parse Lean output and it does not use floating-point tolerances.
"""

from __future__ import annotations

import json
from fractions import Fraction
from math import gcd
from typing import Sequence


CANONICAL_COUNTS = ((1324, 107), (5, 503))
CANONICAL_STATIONARY = (Fraction(7155, 61511), Fraction(54356, 61511))
CANONICAL_MODE = (Fraction(54356), Fraction(-7155))
CANONICAL_EIGENVALUE = Fraction(665437, 726948)
CANONICAL_SAMPLE_TOTAL = 16384


def normalize_rows(
    counts: Sequence[Sequence[int]],
) -> tuple[tuple[Fraction, ...], ...]:
    """Return the row-normalized exact kernel."""

    if len(counts) != 2 or any(len(row) != 2 for row in counts):
        raise ValueError("the audited recurrent kernel must be 2 by 2")
    if any(entry < 0 for row in counts for entry in row):
        raise ValueError("transition counts must be nonnegative")
    row_sums = tuple(sum(row) for row in counts)
    if any(total == 0 for total in row_sums):
        raise ValueError("every transition row must have positive mass")
    return tuple(
        tuple(Fraction(entry, total) for entry in row)
        for row, total in zip(counts, row_sums, strict=True)
    )


def matrix_vector(
    matrix: Sequence[Sequence[Fraction]],
    vector: Sequence[Fraction],
) -> tuple[Fraction, ...]:
    return tuple(
        sum((entry * value for entry, value in zip(row, vector, strict=True)),
            Fraction(0))
        for row in matrix
    )


def row_vector_matrix(
    vector: Sequence[Fraction],
    matrix: Sequence[Sequence[Fraction]],
) -> tuple[Fraction, ...]:
    return tuple(
        sum((vector[i] * matrix[i][j] for i in range(2)), Fraction(0))
        for j in range(2)
    )


def verify_packet(
    *,
    counts: Sequence[Sequence[int]] = CANONICAL_COUNTS,
    stationary: Sequence[Fraction] = CANONICAL_STATIONARY,
    mode: Sequence[Fraction] = CANONICAL_MODE,
    eigenvalue: Fraction = CANONICAL_EIGENVALUE,
    sample_total: int = CANONICAL_SAMPLE_TOTAL,
) -> dict[str, object]:
    """Verify every arithmetic premise used by the Lean obstruction.

    The deterministic-pushforward check is equivalent to asking whether the
    reduced denominator of the stationary mass divides ``sample_total``.
    The adjacent-grid bracket is reported as a second, independent witness.
    """

    kernel = normalize_rows(counts)
    stationary = tuple(Fraction(value) for value in stationary)
    mode = tuple(Fraction(value) for value in mode)
    eigenvalue = Fraction(eigenvalue)
    if len(stationary) != 2 or len(mode) != 2:
        raise ValueError("stationary law and mode must each have length two")
    if sample_total <= 0:
        raise ValueError("sample total must be positive")

    row_stochastic = all(sum(row, Fraction(0)) == 1 for row in kernel)
    stationary_normalized = sum(stationary, Fraction(0)) == 1
    stationary_exact = row_vector_matrix(stationary, kernel) == stationary
    eigenpair_exact = matrix_vector(kernel, mode) == tuple(
        eigenvalue * value for value in mode
    )
    nonzero_mode = any(value != 0 for value in mode)
    genuinely_mixing_mode = 0 < eigenvalue < 1

    # If H^2 = H and HT = TP, an eigenvector Pv = lambda v obeys
    # lambda(lambda - 1) T(v) = 0.  A nonzero coefficient forces T(v)=0.
    intertwiner_coefficient = eigenvalue * (eigenvalue - 1)
    intertwiner_obstruction = intertwiner_coefficient != 0

    first_mass = stationary[0]
    lower_numerator = (first_mass.numerator * sample_total) // first_mass.denominator
    lower_mass = Fraction(lower_numerator, sample_total)
    upper_mass = Fraction(lower_numerator + 1, sample_total)
    adjacent_grid_bracket = lower_mass < first_mass < upper_mass
    denominator_divides_sample_total = sample_total % first_mass.denominator == 0
    coprime_denominators = gcd(sample_total, first_mass.denominator) == 1
    deterministic_pushforward_obstruction = not denominator_divides_sample_total

    checks = {
        "row_stochastic": row_stochastic,
        "stationary_normalized": stationary_normalized,
        "stationary_exact": stationary_exact,
        "eigenpair_exact": eigenpair_exact,
        "nonzero_mode": nonzero_mode,
        "genuinely_mixing_mode": genuinely_mixing_mode,
        "intertwiner_obstruction": intertwiner_obstruction,
        "adjacent_grid_bracket": adjacent_grid_bracket,
        "coprime_denominators": coprime_denominators,
        "deterministic_pushforward_obstruction": (
            deterministic_pushforward_obstruction
        ),
    }
    return {
        "ok": all(checks.values()),
        "checks": checks,
        "kernel": [[str(value) for value in row] for row in kernel],
        "stationary": [str(value) for value in stationary],
        "mode": [str(value) for value in mode],
        "eigenvalue": str(eigenvalue),
        "intertwiner_coefficient": str(intertwiner_coefficient),
        "empirical_bracket": [str(lower_mass), str(upper_mass)],
        "sample_total": sample_total,
    }


def main() -> int:
    result = verify_packet()
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
