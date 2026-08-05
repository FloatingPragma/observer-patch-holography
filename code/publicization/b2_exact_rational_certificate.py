"""Exact-rational finite certificate for the B2 publicization packet.

The certificate uses the nontrivial coordinate partition ``2 + 1`` of a
three-dimensional matrix algebra.  Kraus operators are stored as an unscaled
rational matrix together with their rational squared scale.  This avoids
floating-point square roots while still checking every quadratic Kraus
identity exactly.

The nine rational matrix units form a basis of the full matrix algebra, so
checking the map identities on all nine rows is an exhaustive finite linear
certificate, not random sampling.
"""

from __future__ import annotations

from fractions import Fraction
import json
from typing import Iterable

Q = Fraction
N = 3
Matrix = tuple[tuple[Q, ...], ...]


def matrix(rows: Iterable[Iterable[int | Q]]) -> Matrix:
    return tuple(tuple(Q(x) for x in row) for row in rows)


ZERO = matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
ONE = matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
P0 = matrix([[1, 0, 0], [0, 1, 0], [0, 0, 0]])
P1 = matrix([[0, 0, 0], [0, 0, 0], [0, 0, 1]])
PARTITION = (P0, P1)


def add(*terms: Matrix) -> Matrix:
    return tuple(
        tuple(sum((term[i][j] for term in terms), Q(0)) for j in range(N))
        for i in range(N)
    )


def neg(a: Matrix) -> Matrix:
    return scale(Q(-1), a)


def sub(a: Matrix, b: Matrix) -> Matrix:
    return add(a, neg(b))


def scale(c: Q, a: Matrix) -> Matrix:
    return tuple(tuple(c * a[i][j] for j in range(N)) for i in range(N))


def mul(a: Matrix, b: Matrix) -> Matrix:
    return tuple(
        tuple(sum((a[i][r] * b[r][j] for r in range(N)), Q(0)) for j in range(N))
        for i in range(N)
    )


def transpose(a: Matrix) -> Matrix:
    return tuple(tuple(a[j][i] for j in range(N)) for i in range(N))


def trace(a: Matrix) -> Q:
    return sum((a[i][i] for i in range(N)), Q(0))


def unit(a: int, b: int) -> Matrix:
    return tuple(
        tuple(Q(1) if (i, j) == (a, b) else Q(0) for j in range(N))
        for i in range(N)
    )


BASIS = tuple(unit(a, b) for a in range(N) for b in range(N))


def kraus_rows() -> tuple[tuple[int, int, int, Q, Matrix], ...]:
    """Rows ``(block,a,b,scale_squared,unscaled_operator)``."""

    rows: list[tuple[int, int, int, Q, Matrix]] = []
    for block, p in enumerate(PARTITION):
        rank = trace(p)
        assert rank > 0
        for a in range(N):
            for b in range(N):
                rows.append((block, a, b, Q(1, rank), mul(mul(p, unit(a, b)), p)))
    return tuple(rows)


def kraus_complete() -> Matrix:
    return add(
        *(scale(c2, mul(transpose(a), a)) for _, _, _, c2, a in kraus_rows())
    )


def kraus_average(x: Matrix) -> Matrix:
    return add(
        *(scale(c2, mul(mul(a, x), transpose(a))) for _, _, _, c2, a in kraus_rows())
    )


def partition_average(x: Matrix) -> Matrix:
    return add(
        *(scale(trace(mul(x, p)) / trace(p), p) for p in PARTITION)
    )


def partition_pinching(x: Matrix) -> Matrix:
    return add(*(mul(mul(p, x), p) for p in PARTITION))


def public_relax(amplitude: Q, x: Matrix) -> Matrix:
    ex = partition_average(x)
    return add(ex, scale(amplitude, sub(x, ex)))


def projector_gksl(rate: Q, x: Matrix) -> Matrix:
    half = Q(1, 2)
    return add(
        *(
            scale(
                rate,
                sub(
                    mul(mul(p, x), p),
                    scale(half, add(mul(p, x), mul(x, p))),
                ),
            )
            for p in PARTITION
        )
    )


def _row(a: Matrix) -> list[str]:
    return [str(a[i][j]) for i in range(N) for j in range(N)]


def build_certificate() -> dict[str, object]:
    assert add(*PARTITION) == ONE
    assert all(mul(p, p) == p for p in PARTITION)
    assert mul(P0, P1) == ZERO and mul(P1, P0) == ZERO
    assert kraus_complete() == ONE

    map_rows = []
    for index, x in enumerate(BASIS):
        kx = kraus_average(x)
        ex = partition_average(x)
        assert kx == ex
        assert trace(kx) == trace(x)
        assert public_relax(Q(1, 2), public_relax(Q(1, 3), x)) == public_relax(Q(1, 6), x)
        assert projector_gksl(Q(2), x) == scale(Q(2), sub(partition_pinching(x), x))
        map_rows.append(
            {
                "basis_index": index,
                "partition_average": _row(ex),
                "pinching_generator_rate_2": _row(projector_gksl(Q(2), x)),
            }
        )

    nonzero_kraus_rows = [
        {
            "block": block,
            "a": a,
            "b": b,
            "scale_squared": str(c2),
            "unscaled_operator": _row(op),
        }
        for block, a, b, c2, op in kraus_rows()
        if op != ZERO
    ]

    return {
        "schema": "oph.b2-publicization.exact-rational.v1",
        "dimension": N,
        "partition_ranks": [str(trace(p)) for p in PARTITION],
        "basis_rows_checked": len(BASIS),
        "trace_rows_checked": len(BASIS),
        "kraus_rows_total": len(kraus_rows()),
        "kraus_rows_nonzero": nonzero_kraus_rows,
        "kraus_completeness": _row(kraus_complete()),
        "linear_map_rows": map_rows,
        "relaxation_composition": {"a": "1/2", "b": "1/3", "composite": "1/6"},
        "gksl_rate": "2",
        "claim_boundary": (
            "Exact rational Kraus-normalization, Kraus-map, relaxation-composition, "
            "and projector-GKSL identities on a complete matrix-unit basis; no formal "
            "CP/CPTP predicate and no physical clock or rate attachment."
        ),
    }


if __name__ == "__main__":
    print(json.dumps(build_certificate(), indent=2, sort_keys=True))
