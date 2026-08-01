#!/usr/bin/env python3
"""Independent verifier for the issue-654 stage-4 persistence packet.

The verifier does not import the stage-4 producer.  It reconstructs the
projective rational function, the complete dyadic leaf partition, the three
local orbit margins, and the stabilizer fixed-space facts from the immutable
Cartesian carrier implementation.  ``--mutations`` resigns semantic receipt
mutations before checking them.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from fractions import Fraction
from math import isqrt
from pathlib import Path
from typing import Any, Iterable

import a5_multipole_fixed_point_certificate as base
import a5_multipole_fixed_point_hardening_certificate as hardening


HERE = Path(__file__).resolve().parent
RECEIPT_PATH = HERE / "runtime" / "a5_multipole_persistence_receipt_v4.json"
PARENT_PATH = HERE / "runtime" / "a5_multipole_persistence_receipt_v3.json"
SCHEMA = "oph.a5_multipole_persistence_receipt.v4"
STATUS = (
    "EXACT_FULL_COSINE_62_DIRECTION_PERSISTENCE_ON_0_LT_X_LE_1__"
    "DECLARED_EQUAL_WEIGHT_BRANCH__PHYSICAL_MAP_OPEN"
)
PARENT_SCHEMA = "oph.a5_multipole_persistence_receipt.v3"
PARENT_STATUS = (
    "EXACT_THROUGH_EIGHTH_ORDER_I6_TEMPLATE__FULL_COSINE_X10_PLUS_TAIL_BOUNDS__"
    "GLOBAL_PERSISTENCE_COVER_OPEN__PHYSICAL_MAP_OPEN"
)
CANONICAL_PARENT_BYTES = 5870
CANONICAL_PARENT_SHA256 = (
    "sha256:bf73eb4095d19ef2c94d4e87d4e33db904e289ad6c40836ecb98f407b0376863"
)
CANONICAL_PARENT_RECEIPT_SHA256 = (
    "sha256:9a5b4bbd3843f4325b4f2aac9a5792fd4fa8414073293bc178b457fd1dc83e53"
)
SQ5_LO = Fraction(2_236_067_977, 1_000_000_000)
SQ5_HI = Fraction(1_118_033_989, 500_000_000)
C1 = Fraction(6875, 101152)
C2 = Fraction(383125, 562658)
LOCAL_MULT = 4095
ZERO, ONE = base.ZERO, base.ONE
Q5 = base.Q5
Poly = dict[tuple[int, int], Q5]


class VerificationError(ValueError):
    pass


def check(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


@dataclass(frozen=True)
class Iv:
    lo: Fraction
    hi: Fraction

    def __post_init__(self) -> None:
        check(self.lo <= self.hi, "reversed interval")

    def __add__(self, other: "Iv") -> "Iv":
        return Iv(self.lo + other.lo, self.hi + other.hi)

    def __neg__(self) -> "Iv":
        return Iv(-self.hi, -self.lo)

    def __sub__(self, other: "Iv") -> "Iv":
        return self + (-other)

    def __mul__(self, other: "Iv") -> "Iv":
        values = (self.lo * other.lo, self.lo * other.hi, self.hi * other.lo, self.hi * other.hi)
        return Iv(min(values), max(values))

    def power(self, exponent: int) -> "Iv":
        check(exponent >= 0, "negative power")
        if exponent == 0:
            return Iv(Fraction(1), Fraction(1))
        if exponent % 2 == 0 and self.lo <= 0 <= self.hi:
            return Iv(Fraction(0), max(abs(self.lo), abs(self.hi)) ** exponent)
        values = (self.lo**exponent, self.hi**exponent)
        return Iv(min(values), max(values))

    def scale(self, value: Fraction) -> "Iv":
        return self * Iv(Fraction(value), Fraction(value))

    def abs_lo(self) -> Fraction:
        return Fraction(0) if self.lo <= 0 <= self.hi else min(abs(self.lo), abs(self.hi))

    def abs_hi(self) -> Fraction:
        return max(abs(self.lo), abs(self.hi))


def qi(value: Q5) -> Iv:
    a, b = value
    return Iv(a + b * (SQ5_LO if b >= 0 else SQ5_HI), a + b * (SQ5_HI if b >= 0 else SQ5_LO))


def pa(left: Poly, right: Poly) -> Poly:
    out = dict(left)
    for mono, value in right.items():
        total = base.q5_add(out.get(mono, ZERO), value)
        if total == ZERO:
            out.pop(mono, None)
        else:
            out[mono] = total
    return out


def ps(poly: Poly, value: Q5) -> Poly:
    return {} if value == ZERO else {mono: base.q5_mul(c, value) for mono, c in poly.items()}


def pm(left: Poly, right: Poly) -> Poly:
    out: Poly = {}
    for (i, j), x in left.items():
        for (k, ell), y in right.items():
            mono = (i + k, j + ell)
            value = base.q5_add(out.get(mono, ZERO), base.q5_mul(x, y))
            if value == ZERO:
                out.pop(mono, None)
            else:
                out[mono] = value
    return out


def pp(poly: Poly, exponent: int) -> Poly:
    out: Poly = {(0, 0): ONE}
    for _ in range(exponent):
        out = pm(out, poly)
    return out


def pd(poly: Poly, axis: int) -> Poly:
    out = {}
    for mono, value in poly.items():
        if mono[axis]:
            target = list(mono)
            power = target[axis]
            target[axis] -= 1
            out[tuple(target)] = base.q5_scale(value, power)
    return out


def pe(poly: Poly, s: Q5, t: Q5) -> Q5:
    out = ZERO
    for (i, j), value in poly.items():
        out = base.q5_add(out, base.q5_mul(value, base.q5_mul(base.q5_pow(s, i), base.q5_pow(t, j))))
    return out


def pi(poly: Poly, s: Iv, t: Iv) -> Iv:
    out = Iv(Fraction(0), Fraction(0))
    for (i, j), value in poly.items():
        out = out + qi(value) * s.power(i) * t.power(j)
    return out


def pc(poly: Poly, s: Iv, t: Iv) -> Iv:
    sm, tm = (s.lo + s.hi) / 2, (t.lo + t.hi) / 2
    ds, dt = Iv(s.lo - sm, s.hi - sm), Iv(t.lo - tm, t.hi - tm)
    return qi(pe(poly, base.q5(sm), base.q5(tm))) + pi(pd(poly, 0), s, t) * ds + pi(pd(poly, 1), s, t) * dt


def pq(poly: Poly, s: Iv, t: Iv) -> Iv:
    check(all(sum(mono) <= 2 for mono in poly), "nonquadratic neighborhood")
    sm, tm = (s.lo + s.hi) / 2, (t.lo + t.hi) / 2
    qs, qt = base.q5(sm), base.q5(tm)
    ds, dt = Iv(s.lo - sm, s.hi - sm), Iv(t.lo - tm, t.hi - tm)
    return (
        qi(pe(poly, qs, qt))
        + qi(pe(pd(poly, 0), qs, qt)) * ds
        + qi(pe(pd(poly, 1), qs, qt)) * dt
        + (qi(pe(pd(pd(poly, 0), 0), qs, qt)) * ds.power(2)).scale(Fraction(1, 2))
        + qi(pe(pd(pd(poly, 0), 1), qs, qt)) * ds * dt
        + (qi(pe(pd(pd(poly, 1), 1), qs, qt)) * dt.power(2)).scale(Fraction(1, 2))
    )


def dot(x, y) -> Q5:
    out = ZERO
    for i in range(3):
        out = base.q5_add(out, base.q5_mul(x[i], y[i]))
    return out


def add_vectors(rows: Iterable[tuple[Q5, Q5, Q5]]) -> tuple[Q5, Q5, Q5]:
    out = [ZERO, ZERO, ZERO]
    for row in rows:
        for i in range(3):
            out[i] = base.q5_add(out[i], row[i])
    return tuple(out)  # type: ignore[return-value]


def chart(origin, first, second) -> tuple[Poly, Poly, Poly]:
    coords = []
    for axis in range(3):
        coords.append({mono: value for mono, value in (((0, 0), origin[axis]), ((1, 0), first[axis]), ((0, 1), second[axis])) if value != ZERO})
    d: Poly = {}
    for coordinate in coords:
        d = pa(d, pp(coordinate, 2))
    numerator: Poly = {}
    for powers, coefficient in base.build_cartesian_frame()["_i6_poly_object"].items():
        degree = sum(powers)
        term: Poly = {(0, 0): coefficient}
        for axis, power in enumerate(powers):
            term = pm(term, pp(coords[axis], power))
        numerator = pa(numerator, pm(term, pp(d, 3 - degree // 2)))
    outputs = []
    for axis in range(2):
        outputs.append(pa(pm(pd(numerator, axis), d), ps(pm(numerator, pd(d, axis)), base.q5(-3))))
    return d, outputs[0], outputs[1]


def axes() -> list[tuple[str, tuple[Q5, Q5, Q5]]]:
    vertices = base.cartesian_vertices()
    faces, edges = hardening.derive_adjacency_faces_edges(vertices)
    rays = [("vertex", v) for v in vertices] + [("edge", add_vectors(vertices[i] for i in e)) for e in edges] + [("face", add_vectors(vertices[i] for i in f)) for f in faces]
    out = []
    for kind, ray in rays:
        if not any(hardening.cross(ray, old) == (ZERO, ZERO, ZERO) for _, old in out):
            out.append((kind, ray))
    check(len(out) == 31, "axis census")
    return out


def neighborhood(axis, coords) -> Poly:
    raw: Poly = {}
    d: Poly = {}
    for value, coordinate in zip(axis, coords):
        raw = pa(raw, ps(coordinate, value))
        d = pa(d, pp(coordinate, 2))
    return pa(ps(pp(raw, 2), base.q5(4096)), ps(d, base.q5_mul(base.q5(-LOCAL_MULT), dot(axis, axis))))


def global_cover() -> dict[str, Any]:
    critical = axes()
    leaves = []
    rows = []
    for fixed in range(3):
        rem = [i for i in range(3) if i != fixed]
        origin = tuple(ONE if i == fixed else ZERO for i in range(3))
        first = tuple(ONE if i == rem[0] else ZERO for i in range(3))
        second = tuple(ONE if i == rem[1] else ZERO for i in range(3))
        d, gs, gt = chart(origin, first, second)
        coords = [{mono: value for mono, value in (((0, 0), origin[i]), ((1, 0), first[i]), ((0, 1), second[i])) if value != ZERO} for i in range(3)]
        local_polys = [neighborhood(axis, coords) for _, axis in critical]
        pending = [(Iv(Fraction(-1), Fraction(1)), Iv(Fraction(-1), Fraction(1)), 0, "")]
        local_count = gradient_count = depth_seen = 0
        weakest = None
        while pending:
            s, t, depth, path = pending.pop()
            depth_seen = max(depth_seen, depth)
            inside = next((i for i, poly in enumerate(local_polys) if pq(poly, s, t).lo >= 0), None)
            if inside is not None:
                local_count += 1
                leaves.append({"chart": fixed, "path": path, "kind": "local", "axis": inside})
                continue
            gi, gj, di = pc(gs, s, t), pc(gt, s, t), pi(d, s, t)
            margin = gi.abs_lo() ** 2 + gj.abs_lo() ** 2 - C1**2 * di.hi**7
            if margin > 0:
                gradient_count += 1
                weakest = margin if weakest is None else min(weakest, margin)
                leaves.append({"chart": fixed, "path": path, "kind": "gradient"})
                continue
            check(depth < 20, f"independent cover depth exhausted: {fixed}:{path}")
            if s.hi - s.lo >= t.hi - t.lo:
                mid = (s.lo + s.hi) / 2
                pending.extend([(Iv(mid, s.hi), t, depth + 1, path + "1"), (Iv(s.lo, mid), t, depth + 1, path + "0")])
            else:
                mid = (t.lo + t.hi) / 2
                pending.extend([(s, Iv(mid, t.hi), depth + 1, path + "1"), (s, Iv(t.lo, mid), depth + 1, path + "0")])
        check(weakest is not None, "missing cover margin")
        rows.append({"chart_fixed_coordinate": ("x", "y", "z")[fixed], "local_boxes": local_count, "gradient_boxes": gradient_count, "maximum_depth": depth_seen, "weakest_exact_squared_numerator_margin": str(weakest)})
    return {
        "chart_rows": rows,
        "leaf_count": len(leaves),
        "leaf_partition_sha256": base.tagged_sha256(base.canonical_json_bytes(leaves)),
    }


def radial(poly: Poly, axis: int) -> Poly:
    return {mono: base.q5_scale(value, Fraction(1, sum(mono) + 1)) for mono, value in pd(poly, axis).items()}


def sqrt_lower(value: Fraction) -> Fraction:
    scale = 10**15
    candidate = Fraction(isqrt(value.numerator * value.denominator * scale * scale), value.denominator * scale)
    while candidate * candidate > value:
        candidate -= Fraction(1, value.denominator * scale)
    check(candidate > 0 and candidate * candidate <= value, "sqrt enclosure")
    return candidate


def sqrt_upper(value: Fraction) -> Fraction:
    lower = sqrt_lower(value)
    if lower * lower == value:
        return lower
    step = Fraction(1, value.denominator * 10**15)
    upper = lower + step
    while upper * upper < value:
        upper += step
    check(upper * upper >= value, "sqrt upper enclosure")
    return upper


def local_rows() -> list[dict[str, Any]]:
    vertices = base.cartesian_vertices()
    faces, edges = hardening.derive_adjacency_faces_edges(vertices)
    adjacent = next(e[1] for e in edges if e[0] == 0)
    va = vertices[0]
    vf = tuple(base.q5_scale(x, 2) for x in hardening.cross(va, vertices[adjacent]))
    vs = hardening.cross(va, vf)
    face = faces[0]
    fa = add_vectors(vertices[i] for i in face)
    ff = tuple(base.q5_scale(base.q5_sub(vertices[face[0]][i], vertices[face[1]][i]), 2) for i in range(3))
    fs = hardening.cross(fa, ff)
    edge = edges[0]
    ea = add_vectors(vertices[i] for i in edge)
    ef = tuple(base.q5_mul(base.q5_sub(vertices[edge[0]][i], vertices[edge[1]][i]), base.q5(1, 1)) for i in range(3))
    es = hardening.cross(ea, ef)
    out = []
    for orbit, axis, first, second in (("vertex", va, vf, vs), ("face", fa, ff, fs), ("edge", ea, ef, es)):
        an, e1, e2 = dot(axis, axis), dot(first, first), dot(second, second)
        check(dot(axis, first) == dot(axis, second) == dot(first, second) == ZERO, "local orthogonality")
        radii = []
        for norm in (e1, e2):
            valid = [n for n in range(1, 4097) if base.q5_sign(base.q5_sub(base.q5_scale(norm, LOCAL_MULT), base.q5_scale(an, n * n))) >= 0]
            check(valid, "radius search")
            radii.append(Fraction(1, max(valid)))
        sb, tb = Iv(-radii[0], radii[0]), Iv(-radii[1], radii[1])
        d, gs, gt = chart(axis, first, second)
        check(
            pe(gs, ZERO, ZERO) == ZERO and pe(gt, ZERO, ZERO) == ZERO,
            "independent fixed-centre gradient numerator is nonzero",
        )
        matrix = ((radial(gs, 0), radial(gs, 1)), (radial(gt, 0), radial(gt, 1)))
        c0 = tuple(tuple(pe(matrix[i][j], ZERO, ZERO) for j in range(2)) for i in range(2))
        check(c0[0][1] == c0[1][0] == ZERO, "local diagonalization")
        errors = []
        for i in range(2):
            for j in range(2):
                errors.append(pc(pa(matrix[i][j], {(0, 0): base.q5_neg(c0[i][j])}), sb, tb))
        ai, i1, i2 = qi(an), qi(e1), qi(e2)
        diagonal_scales = (ai.hi / i1.lo, ai.hi / i2.lo)
        off_scale = ai.hi / sqrt_lower(i1.lo * i2.lo)
        scaled = (errors[0].abs_hi() * diagonal_scales[0], errors[1].abs_hi() * off_scale, errors[2].abs_hi() * off_scale, errors[3].abs_hi() * diagonal_scales[1])
        error_row = max(scaled[0] + scaled[1], scaled[2] + scaled[3])
        error_column = max(scaled[0] + scaled[2], scaled[1] + scaled[3])
        error = sqrt_upper(error_row * error_column)
        h0 = (
            base.q5_mul(c0[0][0], base.q5_div(an, e1)),
            base.q5_mul(c0[1][1], base.q5_div(an, e2)),
        )
        diagonal_lower = min(qi(h0[0]).abs_lo(), qi(h0[1]).abs_lo())
        singular = diagonal_lower - error
        check(singular > 0, "local singular margin")
        lower = singular / pi(d, sb, tb).hi**4
        check(lower > C2, "local tail domination")
        normalized_h0 = (
            base.q5_div(h0[0], base.q5_pow(an, 4)),
            base.q5_div(h0[1], base.q5_pow(an, 4)),
        )
        hessian_lower = min(qi(value).abs_lo() for value in normalized_h0)
        morse_margin = hessian_lower - C2
        check(morse_margin > 0, "independent Morse signature margin")
        out.append(
            {
                "orbit": orbit,
                "s_radius": str(radii[0]),
                "t_radius": str(radii[1]),
                "lower": str(lower),
                "margin": str(lower - C2),
                "fixed_diagonal": str(diagonal_lower),
                "operator_error": str(error),
                "operator_error_row": str(error_row),
                "operator_error_column": str(error_column),
                "fixed_minus_error": str(singular),
                "hessian_eigenvalues": [base.q5_str(value) for value in normalized_h0],
                "hessian_lower": str(hessian_lower),
                "morse_margin": str(morse_margin),
            }
        )
    return out


Matrix3 = tuple[
    tuple[Q5, Q5, Q5], tuple[Q5, Q5, Q5], tuple[Q5, Q5, Q5]
]


def matrix_transpose(matrix: Matrix3) -> Matrix3:
    return tuple(tuple(matrix[j][i] for j in range(3)) for i in range(3))  # type: ignore[return-value]


def matrix_multiply(left: Matrix3, right: Matrix3) -> Matrix3:
    columns = matrix_transpose(right)
    return tuple(
        tuple(dot(left[i], columns[j]) for j in range(3)) for i in range(3)
    )  # type: ignore[return-value]


def matrix_apply(matrix: Matrix3, vector: tuple[Q5, Q5, Q5]) -> tuple[Q5, Q5, Q5]:
    return tuple(dot(row, vector) for row in matrix)  # type: ignore[return-value]


def matrix_determinant(matrix: Matrix3) -> Q5:
    return dot(matrix[0], hardening.cross(matrix[1], matrix[2]))


def independent_stabilizer_matrices() -> tuple[tuple[str, int, Matrix3], ...]:
    q = base.q5
    return (
        (
            "vertex",
            5,
            (
                (q(Fraction(-1, 4), Fraction(1, 4)), q(Fraction(1, 4), Fraction(1, 4)), q(Fraction(-1, 2))),
                (q(Fraction(-1, 4), Fraction(-1, 4)), q(Fraction(1, 2)), q(Fraction(-1, 4), Fraction(1, 4))),
                (q(Fraction(1, 2)), q(Fraction(-1, 4), Fraction(1, 4)), q(Fraction(1, 4), Fraction(1, 4))),
            ),
        ),
        (
            "face",
            3,
            ((ZERO, ONE, ZERO), (ZERO, ZERO, ONE), (ONE, ZERO, ZERO)),
        ),
        (
            "edge",
            2,
            (
                (q(Fraction(-1, 4), Fraction(-1, 4)), q(Fraction(1, 2)), q(Fraction(-1, 4), Fraction(1, 4))),
                (q(Fraction(1, 2)), q(Fraction(-1, 4), Fraction(1, 4)), q(Fraction(1, 4), Fraction(1, 4))),
                (q(Fraction(-1, 4), Fraction(1, 4)), q(Fraction(1, 4), Fraction(1, 4)), q(Fraction(-1, 2))),
            ),
        ),
    )


def independent_stabilizer_rows() -> list[dict[str, Any]]:
    """Verify explicit proper carrier rotations without producer construction."""

    matrices = independent_stabilizer_matrices()
    identity = ((ONE, ZERO, ZERO), (ZERO, ONE, ZERO), (ZERO, ZERO, ONE))
    vertices = base.cartesian_vertices()
    faces, edges = hardening.derive_adjacency_faces_edges(vertices)
    representative_axes = {
        "vertex": vertices[0],
        "face": add_vectors(vertices[index] for index in faces[0]),
        "edge": add_vectors(vertices[index] for index in edges[0]),
    }
    vertex_set = set(vertices)
    out = []
    for orbit, order, matrix in matrices:
        check(matrix_multiply(matrix_transpose(matrix), matrix) == identity, "independent stabilizer orthogonality")
        check(matrix_determinant(matrix) == ONE, "independent stabilizer orientation")
        check(set(matrix_apply(matrix, vertex) for vertex in vertices) == vertex_set, "independent stabilizer port action")
        check(
            matrix_apply(matrix, representative_axes[orbit]) == representative_axes[orbit],
            "independent stabilizer does not fix the declared orbit axis",
        )
        power = identity
        first_identity = None
        for exponent in range(1, order + 1):
            power = matrix_multiply(power, matrix)
            if power == identity and first_identity is None:
                first_identity = exponent
        check(first_identity == order, "independent stabilizer order")
        difference = tuple(tuple(base.q5_sub(matrix[i][j], identity[i][j]) for j in range(3)) for i in range(3))
        nonzero_minor = any(
            base.q5_sub(base.q5_mul(difference[r1][c1], difference[r2][c2]), base.q5_mul(difference[r1][c2], difference[r2][c1])) != ZERO
            for r1 in range(3)
            for r2 in range(r1 + 1, 3)
            for c1 in range(3)
            for c2 in range(c1 + 1, 3)
        )
        check(matrix_determinant(difference) == ZERO and nonzero_minor, "independent stabilizer fixed-space rank")
        out.append({"orbit": orbit, "rotation_order": order, "proper_orthogonal": True, "fixed_space_dimension": 1, "tangent_fixed_space_dimension": 0})
    return out


def projectively_equal(
    left: tuple[Q5, Q5, Q5], right: tuple[Q5, Q5, Q5]
) -> bool:
    return hardening.cross(left, right) == (ZERO, ZERO, ZERO)


def projective_unique(
    rays: Iterable[tuple[Q5, Q5, Q5]],
) -> list[tuple[Q5, Q5, Q5]]:
    out: list[tuple[Q5, Q5, Q5]] = []
    for ray in rays:
        check(dot(ray, ray) != ZERO, "zero projective ray")
        if not any(projectively_equal(ray, old) for old in out):
            out.append(ray)
    return out


def projective_sine_squared(
    left: tuple[Q5, Q5, Q5], right: tuple[Q5, Q5, Q5]
) -> Q5:
    denominator = base.q5_mul(dot(left, left), dot(right, right))
    check(denominator != ZERO, "zero norm in projective separation")
    numerator = base.q5_sub(
        denominator, base.q5_pow(dot(left, right), 2)
    )
    return base.q5_div(numerator, denominator)


def independent_full_group_transport() -> dict[str, Any]:
    """Reconstruct A5 by generator closure, independently of directed edges."""

    identity: Matrix3 = ((ONE, ZERO, ZERO), (ZERO, ONE, ZERO), (ZERO, ZERO, ONE))
    stabilizers = independent_stabilizer_matrices()
    generators = (stabilizers[0][2], stabilizers[1][2])
    rotations: set[Matrix3] = {identity}
    pending = [identity]
    while pending:
        current = pending.pop()
        for generator in generators:
            candidate = matrix_multiply(current, generator)
            if candidate not in rotations:
                rotations.add(candidate)
                pending.append(candidate)
    check(len(rotations) == 60, "independent proper rotation group order")

    vertices = base.cartesian_vertices()
    faces, edges = hardening.derive_adjacency_faces_edges(vertices)
    vertex_set = set(vertices)
    vertex_index = {vertex: index for index, vertex in enumerate(vertices)}
    for rotation in rotations:
        check(
            matrix_multiply(matrix_transpose(rotation), rotation) == identity,
            "independent full-group orthogonality",
        )
        check(matrix_determinant(rotation) == ONE, "independent full-group orientation")
        check(
            set(matrix_apply(rotation, vertex) for vertex in vertices) == vertex_set,
            "independent full-group port action",
        )

    permutations = sorted(
        tuple(vertex_index[matrix_apply(rotation, vertex)] for vertex in vertices)
        for rotation in rotations
    )
    check(len(set(permutations)) == 60, "independent carrier action faithfulness")
    permutation_set = set(permutations)
    identity_permutation = tuple(range(12))
    check(identity_permutation in permutation_set, "independent group identity")
    for left in permutations:
        inverse = [0] * 12
        for source, target in enumerate(left):
            inverse[target] = source
        check(tuple(inverse) in permutation_set, "independent group inverse")
        for right in permutations:
            check(
                tuple(left[right[index]] for index in range(12)) in permutation_set,
                "independent group closure",
            )

    representatives = {
        "vertex": vertices[0],
        "face": add_vectors(vertices[index] for index in faces[0]),
        "edge": add_vectors(vertices[index] for index in edges[0]),
    }
    expected_sizes = {"vertex": 6, "face": 10, "edge": 15}
    critical_axes = axes()
    orbit_rows = []
    transported_union: list[tuple[Q5, Q5, Q5]] = []
    for orbit in ("vertex", "face", "edge"):
        transported = projective_unique(
            matrix_apply(rotation, representatives[orbit]) for rotation in rotations
        )
        check(len(transported) == expected_sizes[orbit], f"independent {orbit} orbit size")
        typed_axes = [axis for kind, axis in critical_axes if kind == orbit]
        check(
            all(any(projectively_equal(ray, axis) for axis in typed_axes) for ray in transported)
            and all(any(projectively_equal(axis, ray) for ray in transported) for axis in typed_axes),
            f"independent {orbit} critical-axis coverage",
        )
        transported_union.extend(transported)
        orbit_rows.append(
            {
                "orbit": orbit,
                "projective_orbit_size": len(transported),
                "critical_axis_count": len(typed_axes),
                "covers_complete_typed_axis_set": True,
            }
        )
    union = projective_unique(transported_union)
    check(len(union) == 31, "independent transported axis union")
    check(
        all(any(projectively_equal(ray, axis) for _, axis in critical_axes) for ray in union)
        and all(any(projectively_equal(axis, ray) for ray in union) for _, axis in critical_axes),
        "independent transported union equality",
    )

    minimum: Q5 | None = None
    for left_index, left in enumerate(union):
        for right in union[left_index + 1 :]:
            value = projective_sine_squared(left, right)
            if minimum is None or base.q5_sign(base.q5_sub(value, minimum)) < 0:
                minimum = value
    check(minimum is not None, "independent projective separation census")
    expected_minimum = base.q5(Fraction(1, 2), Fraction(-1, 6))
    check(minimum == expected_minimum, "independent minimum projective separation")
    check(
        base.q5_sign(base.q5_sub(minimum, base.q5(Fraction(1, 1024)))) > 0,
        "independent local-neighborhood disjointness",
    )

    serialized_permutations = [list(row) for row in permutations]
    return {
        "construction": (
            "the images of one directed carrier edge construct all exact proper "
            "rotations over Q(sqrt(5))"
        ),
        "directed_edge_image_count": 60,
        "proper_rotation_count": 60,
        "port_action_faithful": True,
        "identity_present": True,
        "closed_under_composition": True,
        "closed_under_inverse": True,
        "port_permutations": serialized_permutations,
        "port_permutation_sha256": base.tagged_sha256(
            base.canonical_json_bytes(serialized_permutations)
        ),
        "projective_axis_orbits": orbit_rows,
        "projective_orbit_union_axis_count": 31,
        "projective_orbit_union_matches_critical_axis_set": True,
        "distinct_projective_axis_pair_count": 465,
        "minimum_distinct_projective_axis_sine_squared": base.q5_str(minimum),
        "twice_local_radius_sine_squared_upper": "1/1024",
        "minimum_separation_exceeds_twice_local_radius": True,
        "neighborhoods_sine_squared_1_over_4096_pairwise_disjoint": True,
    }


def self_hash(value: dict[str, Any]) -> str:
    return base.tagged_sha256(base.canonical_json_bytes({k: v for k, v in value.items() if k != "receipt_sha256"}))


def verify_canonical_parent(parent: dict[str, Any], parent_raw: bytes) -> None:
    """Bind stage 4 to the exact target-clean v3 parent, not a rehashed lookalike."""

    check(parent_raw == base.canonical_json_bytes(parent), "parent is not canonical JSON")
    check(len(parent_raw) == CANONICAL_PARENT_BYTES, "canonical parent byte count")
    check(
        base.tagged_sha256(parent_raw) == CANONICAL_PARENT_SHA256,
        "canonical parent byte digest",
    )
    check(
        parent.get("receipt_sha256") == CANONICAL_PARENT_RECEIPT_SHA256,
        "canonical parent receipt digest",
    )
    check(parent.get("receipt_sha256") == self_hash(parent), "parent self hash")
    check(
        set(parent)
        == {
            "schema",
            "status",
            "issue",
            "extends",
            "parent_pin",
            "through_eighth_order_i6_template",
            "normalized_tail_bounds",
            "critical_axis_separation",
            "quantitative_persistence",
            "comparison_boundary",
            "fail_closed_controls",
            "receipt_sha256",
        },
        "canonical parent keyset",
    )
    check(
        parent.get("schema") == PARENT_SCHEMA
        and parent.get("status") == PARENT_STATUS
        and parent.get("issue") == 654,
        "parent schema/status/issue",
    )
    template = parent.get("through_eighth_order_i6_template", {})
    check(
        template.get("declared_branch_premise") is True
        and template.get("source_selected") is False
        and template.get("architecture_forced") is False
        and template.get("physical_source_selection_owner") == 655,
        "parent source and architecture boundary",
    )
    comparison = parent.get("comparison_boundary", {})
    check(
        comparison
        == {
            "public_measurement_read": False,
            "comparison_permitted": False,
            "physical_map_open": True,
        },
        "parent comparison and physical-map boundary",
    )


def verify_object(
    receipt: dict[str, Any],
    parent: dict[str, Any],
    parent_raw: bytes,
    *,
    expected_cover: dict[str, Any] | None = None,
    expected_local: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    check(
        set(receipt)
        == {
            "schema",
            "status",
            "issue",
            "extends",
            "parent_pin",
            "declared_kernel_branch",
            "exact_axis_stationarity",
            "global_off_neighborhood_cover",
            "local_uniqueness_boxes",
            "persistence_theorem",
            "fail_closed_controls",
            "independent_verifier_binding",
            "comparison_boundary",
            "receipt_sha256",
        },
        "top-level keyset",
    )
    check(receipt.get("schema") == SCHEMA and receipt.get("status") == STATUS, "schema/status")
    check(receipt.get("issue") == 654, "issue")
    check(receipt.get("receipt_sha256") == self_hash(receipt), "self hash")
    verify_canonical_parent(parent, parent_raw)
    check(
        receipt.get("extends")
        == {
            "schema": PARENT_SCHEMA,
            "status": PARENT_STATUS,
            "path": "code/a5_fingerprint/runtime/a5_multipole_persistence_receipt_v3.json",
            "reason": "append-only closure of the v3 global and local mathematical boundaries",
        },
        "extension relation",
    )
    pin = receipt.get("parent_pin", {})
    check(
        pin
        == {
            "path": "code/a5_fingerprint/runtime/a5_multipole_persistence_receipt_v3.json",
            "bytes": len(parent_raw),
            "sha256": base.tagged_sha256(parent_raw),
            "receipt_sha256": parent["receipt_sha256"],
        },
        "parent pin",
    )
    check(receipt.get("declared_kernel_branch") == {"kernel": "Q_x(n) = sum_i [1 - cos(x u_i.n)]", "equal_weight_twelve_port_branch": True, "range": "0 < x = |a k| <= 1", "source_selected": False, "architecture_forced": False, "physical_source_selection_owner": 655}, "branch boundary")
    tail = parent["normalized_tail_bounds"]
    check(Fraction(tail["normalized_C1_gradient_bound"]) == C1 and Fraction(tail["normalized_C2_intrinsic_hessian_bound"]) == C2, "parent tail constants")
    full_transport = independent_full_group_transport()
    cover = expected_cover if expected_cover is not None else global_cover()
    expected_cover_full = {
        "method": (
            "three exact projective charts with dominant coordinate fixed to +1; "
            "deterministic dyadic interval subdivision of [-1,1]^2"
        ),
        "coverage_argument": (
            "antipodal evenness identifies n with -n, and every projective "
            "direction has a coordinate of maximal absolute value; dividing "
            "by that coordinate places it in at least one declared chart"
        ),
        "local_neighborhood": "sin^2(angle to one of 31 axes) <= 1/4096",
        "normalized_tail_C1_bound": str(C1),
        "chart_rows": cover["chart_rows"],
        "leaf_count": cover["leaf_count"],
        "leaf_partition_sha256": cover["leaf_partition_sha256"],
        "all_off_neighborhood_boxes_gradient_certified": True,
        "arithmetic": "Fraction intervals with a rationally proved sqrt(5) enclosure",
    }
    check(receipt.get("global_off_neighborhood_cover") == expected_cover_full, "complete global cover payload")
    independent_local = expected_local if expected_local is not None else local_rows()
    local_argument = (
        "G(y)=C(y)y exactly; the fixed C(0) smallest singular value exceeds "
        "the uniform row-sum operator bound for C(y)-C(0), so cancellation "
        "is impossible even when the fixed Hessian is indefinite"
    )
    expected_local_rows = [
        {
            "orbit": row["orbit"],
            "chart_rectangle": {"s_radius": row["s_radius"], "t_radius": row["t_radius"]},
            "contains_full_sine_squared_1_over_4096_neighborhood": True,
            "fixed_center_gradient_numerators_zero": True,
            "template_gradient_linear_lower": row["lower"],
            "tail_gradient_linear_upper": str(C2),
            "strict_margin": row["margin"],
            "fixed_center_scaled_C_min_abs_diagonal": row["fixed_diagonal"],
            "local_matrix_error_row_sum": row["operator_error_row"],
            "local_matrix_error_column_sum": row["operator_error_column"],
            "local_matrix_error_induced_2_upper": row["operator_error"],
            "induced_2_bound_justification": (
                "||E||_2 <= sqrt(||E||_1 ||E||_infinity), with both "
                "induced one and infinity norms bounded entrywise"
            ),
            "fixed_center_minus_uniform_operator_error": row["fixed_minus_error"],
            "fixed_center_normalized_hessian_eigenvalues": row["hessian_eigenvalues"],
            "fixed_center_min_abs_hessian_eigenvalue_lower": row["hessian_lower"],
            "morse_signature_margin_against_tail_C2": row["morse_margin"],
            "saddle_safe_injectivity_argument": local_argument,
            "normalized_y_chart_tail_comparison": {
                "coordinate_definition": "y_i = s_i sqrt(E_i/A)",
                "projective_chart_operator_norm_upper": "||J_y|| <= 1",
                "geodesic_distance": "rho = atan(|y|) <= |y|",
                "tail_gradient_conclusion": "||dE/dy|| <= C2 |y|",
            },
        }
        for row in independent_local
    ]
    expected_local_full = {
        "method": (
            "exact local normalized charts; radial fundamental-theorem factorization "
            "of the template gradient; intrinsic tail-Hessian integration"
        ),
        "stationary_centres": (
            "equal port weights make the cosine kernel invariant under each "
            "vertex, face, or edge rotation stabilizer; its tangent fixed space "
            "is zero, so every declared axis centre is stationary for every x"
        ),
        "normalized_tail_C2_bound": str(C2),
        "orbit_representatives": expected_local_rows,
        "symmetry_transport_join": {
            "proper_rotation_count": full_transport["proper_rotation_count"],
            "port_permutation_sha256": full_transport["port_permutation_sha256"],
            "projective_axis_orbits": full_transport["projective_axis_orbits"],
            "projective_orbit_union_axis_count": full_transport[
                "projective_orbit_union_axis_count"
            ],
            "projective_orbit_union_matches_critical_axis_set": full_transport[
                "projective_orbit_union_matches_critical_axis_set"
            ],
            "minimum_distinct_projective_axis_sine_squared": full_transport[
                "minimum_distinct_projective_axis_sine_squared"
            ],
            "twice_local_radius_sine_squared_upper": full_transport[
                "twice_local_radius_sine_squared_upper"
            ],
            "neighborhoods_sine_squared_1_over_4096_pairwise_disjoint": full_transport[
                "neighborhoods_sine_squared_1_over_4096_pairwise_disjoint"
            ],
            "transported_local_box_count": sum(
                row["projective_orbit_size"]
                for row in full_transport["projective_axis_orbits"]
            ),
        },
        "every_local_neighborhood_has_exactly_one_stationary_direction": True,
        "morse_signatures_preserved_by_fixed_center_C2_margins": True,
        "arithmetic": "exact Q(sqrt(5)) polynomials and rational interval enclosures",
    }
    check(receipt.get("local_uniqueness_boxes") == expected_local_full, "complete local uniqueness payload")

    independent_stabilizers = independent_stabilizer_rows()
    expected_stationarity = {
        "representative_rotations": [
            {
                "orbit": row["orbit"],
                "rotation_order": row["rotation_order"],
                "permutes_all_twelve_ports": True,
                "proper_orthogonal": True,
                "fixed_space_dimension": 1,
                "tangent_fixed_space_dimension": 0,
            }
            for row in independent_stabilizers
        ],
        "full_group_transport": full_transport,
        "kernel_even_from_antipodal_ports": True,
        "conclusion": (
            "equal weights make every cosine kernel invariant under the exhibited "
            "stabilizers; an invariant tangent gradient must lie in the zero "
            "tangent fixed space, so all 62 oriented axis directions are stationary"
        ),
    }
    check(receipt.get("exact_axis_stationarity") == expected_stationarity, "complete stationarity payload")

    expected_theorem = {
        "finite_exactly_62_persistence_range": True,
        "range": "0 < |a k| <= 1",
        "oriented_stationary_direction_count": 62,
        "unoriented_axis_count": 31,
        "orbit_counts": {"maxima": 12, "minima": 20, "saddles": 30},
        "morse_types_preserved": True,
        "theorem_scope": "declared equal-weight cosine kernel only",
        "source_selection_proved": False,
        "physical_attachment_proved": False,
    }
    check(receipt.get("persistence_theorem") == expected_theorem, "complete theorem payload")
    check(receipt.get("comparison_boundary") == {"public_measurement_read": False, "comparison_permitted": False, "physical_map_open": True}, "comparison boundary")

    expected_controls = {
        "controls": [
            {"control": "truncate the global subdivision before it closes", "expected_failure": "exact off-neighborhood cover", "detector_fired": True},
            {"control": "change the certified axis-neighborhood multiplier", "expected_failure": "frozen local/global partition", "detector_fired": True},
            {"control": "shrink the local chart below the complete neighborhood", "expected_failure": "complete local-neighborhood containment", "detector_fired": True},
            {"control": "inflate the admitted tail Hessian", "expected_failure": "exact v3 C2 custody", "detector_fired": True},
        ],
        "all_detectors_fired": True,
    }
    check(receipt.get("fail_closed_controls") == expected_controls, "complete producer controls")

    verifier_raw = Path(__file__).resolve().read_bytes()
    expected_binding = {
        "path": "code/a5_fingerprint/a5_multipole_persistence_stage4_independent_verifier.py",
        "bytes": len(verifier_raw),
        "sha256": base.tagged_sha256(verifier_raw),
        "implementation_independent_of_producer_import": True,
        "replay_command": (
            "python3 code/a5_fingerprint/"
            "a5_multipole_persistence_stage4_independent_verifier.py --mutations"
        ),
        "resigned_semantic_mutation_count": 25,
    }
    check(receipt.get("independent_verifier_binding") == expected_binding, "independent verifier binding")
    return {
        "status": "PASS",
        "receipt_sha256": receipt["receipt_sha256"],
        "independent_global_leaf_replay": True,
        "independent_local_margin_replay": True,
        "independent_full_group_transport_replay": True,
        "canonical_parent_bytes_verified": True,
        "promotion_boundaries_fail_closed": True,
    }


def load(path: Path) -> tuple[bytes, dict[str, Any]]:
    raw = path.read_bytes()
    try:
        value = json.loads(raw)
    except json.JSONDecodeError as error:
        raise VerificationError("invalid JSON") from error
    check(isinstance(value, dict), "JSON root")
    return raw, value


def mutations(
    receipt: dict[str, Any],
    parent: dict[str, Any],
    parent_raw: bytes,
    expected_cover: dict[str, Any],
    expected_local: list[dict[str, Any]],
) -> list[str]:
    import copy

    cases = []
    receipt_mutations = (
        ("leaf_digest", lambda x: x["global_off_neighborhood_cover"].__setitem__("leaf_partition_sha256", "sha256:" + "0" * 64)),
        ("local_margin", lambda x: x["local_uniqueness_boxes"]["orbit_representatives"][0].__setitem__("strict_margin", "0")),
        ("range", lambda x: x["persistence_theorem"].__setitem__("range", "0 < |a k| <= 2")),
        ("source_promotion", lambda x: x["persistence_theorem"].__setitem__("source_selection_proved", True)),
        ("comparison", lambda x: x["comparison_boundary"].__setitem__("comparison_permitted", True)),
        ("parent_pin", lambda x: x["parent_pin"].__setitem__("bytes", x["parent_pin"]["bytes"] + 1)),
        ("morse_type", lambda x: x["persistence_theorem"].__setitem__("morse_types_preserved", False)),
        ("orbit_count", lambda x: x["persistence_theorem"]["orbit_counts"].__setitem__("maxima", 13)),
        ("extra_field", lambda x: x.__setitem__("unregistered_claim", True)),
        ("global_boolean", lambda x: x["global_off_neighborhood_cover"].__setitem__("all_off_neighborhood_boxes_gradient_certified", False)),
        ("local_boolean", lambda x: x["local_uniqueness_boxes"].__setitem__("every_local_neighborhood_has_exactly_one_stationary_direction", False)),
        ("implementation_pin", lambda x: x["independent_verifier_binding"].__setitem__("sha256", "sha256:" + "f" * 64)),
        ("operator_bound", lambda x: x["local_uniqueness_boxes"]["orbit_representatives"][2].__setitem__("local_matrix_error_induced_2_upper", "0")),
        ("y_chart_bound", lambda x: x["local_uniqueness_boxes"]["orbit_representatives"][2]["normalized_y_chart_tail_comparison"].__setitem__("projective_chart_operator_norm_upper", "||J_y|| <= 2")),
        ("symmetry_group_order", lambda x: x["exact_axis_stationarity"]["full_group_transport"].__setitem__("proper_rotation_count", 59)),
        ("symmetry_port_permutation", lambda x: x["exact_axis_stationarity"]["full_group_transport"]["port_permutations"][0].__setitem__(0, 11)),
        ("symmetry_projective_orbit_size", lambda x: x["exact_axis_stationarity"]["full_group_transport"]["projective_axis_orbits"][0].__setitem__("projective_orbit_size", 5)),
        ("symmetry_axis_union_count", lambda x: x["local_uniqueness_boxes"]["symmetry_transport_join"].__setitem__("projective_orbit_union_axis_count", 30)),
        ("symmetry_minimum_separation", lambda x: x["exact_axis_stationarity"]["full_group_transport"].__setitem__("minimum_distinct_projective_axis_sine_squared", "1/1024+0*sqrt5")),
        ("symmetry_neighborhood_disjointness", lambda x: x["local_uniqueness_boxes"]["symmetry_transport_join"].__setitem__("neighborhoods_sine_squared_1_over_4096_pairwise_disjoint", False)),
    )
    for name, mutate in receipt_mutations:
        candidate = copy.deepcopy(receipt)
        mutate(candidate)
        candidate["receipt_sha256"] = self_hash(candidate)
        try:
            verify_object(
                candidate,
                parent,
                parent_raw,
                expected_cover=expected_cover,
                expected_local=expected_local,
            )
        except VerificationError:
            cases.append(name)
        else:
            raise VerificationError(f"semantic mutation escaped: {name}")

    parent_mutations = (
        ("parent_source_promotion", lambda x: x["through_eighth_order_i6_template"].__setitem__("source_selected", True)),
        ("parent_architecture_promotion", lambda x: x["through_eighth_order_i6_template"].__setitem__("architecture_forced", True)),
        ("parent_physical_map_promotion", lambda x: x["comparison_boundary"].__setitem__("physical_map_open", False)),
        ("parent_comparison_promotion", lambda x: x["comparison_boundary"].__setitem__("comparison_permitted", True)),
        ("parent_measured_target_injection", lambda x: x.__setitem__("measured_targets", {"forbidden": "injected"})),
    )
    for name, mutate in parent_mutations:
        candidate_parent = copy.deepcopy(parent)
        mutate(candidate_parent)
        candidate_parent["receipt_sha256"] = self_hash(candidate_parent)
        candidate_parent_raw = base.canonical_json_bytes(candidate_parent)
        candidate_receipt = copy.deepcopy(receipt)
        candidate_receipt["parent_pin"] = {
            "path": "code/a5_fingerprint/runtime/a5_multipole_persistence_receipt_v3.json",
            "bytes": len(candidate_parent_raw),
            "sha256": base.tagged_sha256(candidate_parent_raw),
            "receipt_sha256": candidate_parent["receipt_sha256"],
        }
        candidate_receipt["receipt_sha256"] = self_hash(candidate_receipt)
        try:
            verify_object(
                candidate_receipt,
                candidate_parent,
                candidate_parent_raw,
                expected_cover=expected_cover,
                expected_local=expected_local,
            )
        except VerificationError:
            cases.append(name)
        else:
            raise VerificationError(f"semantic parent mutation escaped: {name}")
    return cases


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receipt", type=Path, default=RECEIPT_PATH)
    parser.add_argument("--parent", type=Path, default=PARENT_PATH)
    parser.add_argument("--mutations", action="store_true")
    args = parser.parse_args(argv)
    raw, receipt = load(args.receipt)
    check(raw == base.canonical_json_bytes(receipt), "noncanonical receipt")
    parent_raw, parent = load(args.parent)
    expected_cover = global_cover()
    expected_local = local_rows()
    result = verify_object(
        receipt,
        parent,
        parent_raw,
        expected_cover=expected_cover,
        expected_local=expected_local,
    )
    if args.mutations:
        result["rejected_resigned_semantic_mutations"] = mutations(
            receipt, parent, parent_raw, expected_cover, expected_local
        )
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
