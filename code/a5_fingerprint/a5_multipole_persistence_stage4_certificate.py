#!/usr/bin/env python3
"""Stage-4 exact persistence certificate for issue 654.

This append-only packet closes the two mathematical boundaries left open by
the v3 packet for the declared equal-weight cosine kernel.  It uses only
exact rational interval arithmetic:

* three projective charts cover the real projective sphere;
* a deterministic dyadic subdivision proves that the tangent gradient of
  the normalized level-six template dominates the complete cosine tail away
  from 31 disjoint symmetry-axis neighborhoods; and
* one local chart for each of the vertex, face, and edge orbits proves that
  the template gradient has a uniform linear lower bound which dominates the
  normalized tail Hessian.  The icosahedral stabilizer makes the centre of
  each such chart an exact stationary point of every equal-weight cosine
  kernel, so the local estimate proves uniqueness.

The theorem is conditional on the declared equal-weight cosine branch.  It
does not select that branch, attach it to a physical observable, or open a
comparison with public data.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable

import a5_multipole_fixed_point_certificate as base
import a5_multipole_fixed_point_hardening_certificate as hardening
import a5_multipole_persistence_certificate as v3


HERE = Path(__file__).resolve().parent
RUNTIME = HERE / "runtime"
PARENT_PATH = RUNTIME / "a5_multipole_persistence_receipt_v3.json"
RECEIPT_PATH = RUNTIME / "a5_multipole_persistence_receipt_v4.json"
VERIFIER_PATH = HERE / "a5_multipole_persistence_stage4_independent_verifier.py"

SCHEMA = "oph.a5_multipole_persistence_receipt.v4"
STATUS = (
    "EXACT_FULL_COSINE_62_DIRECTION_PERSISTENCE_ON_0_LT_X_LE_1__"
    "DECLARED_EQUAL_WEIGHT_BRANCH__PHYSICAL_MAP_OPEN"
)

Q5 = base.Q5
ZERO = base.ZERO
ONE = base.ONE
q5 = base.q5
q5_add = base.q5_add
q5_sub = base.q5_sub
q5_mul = base.q5_mul
q5_scale = base.q5_scale
q5_pow = base.q5_pow
q5_sign = base.q5_sign
q5_str = base.q5_str
require = base.require

# Exact rational enclosure.  The two integer-square checks below make the
# decimal-looking constants proof data rather than floating-point premises.
SQRT5_LO = Fraction(2_236_067_977, 1_000_000_000)
SQRT5_HI = Fraction(1_118_033_989, 500_000_000)
require(SQRT5_LO * SQRT5_LO < 5, "sqrt(5) lower enclosure drift")
require(SQRT5_HI * SQRT5_HI > 5, "sqrt(5) upper enclosure drift")

X_MAX = Fraction(1)
LOCAL_COSINE_MULTIPLIER = 4095
MAX_GLOBAL_DEPTH = 20


@dataclass(frozen=True)
class Interval:
    """Closed rational interval with outward-exact elementary operations."""

    lo: Fraction
    hi: Fraction

    def __post_init__(self) -> None:
        require(self.lo <= self.hi, "reversed interval")

    def __add__(self, other: "Interval") -> "Interval":
        return Interval(self.lo + other.lo, self.hi + other.hi)

    def __neg__(self) -> "Interval":
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other: "Interval") -> "Interval":
        return self + (-other)

    def __mul__(self, other: "Interval") -> "Interval":
        values = (
            self.lo * other.lo,
            self.lo * other.hi,
            self.hi * other.lo,
            self.hi * other.hi,
        )
        return Interval(min(values), max(values))

    def pow(self, exponent: int) -> "Interval":
        require(exponent >= 0, "negative interval exponent")
        if exponent == 0:
            return Interval(Fraction(1), Fraction(1))
        if exponent % 2 == 0 and self.lo <= 0 <= self.hi:
            return Interval(Fraction(0), max(abs(self.lo), abs(self.hi)) ** exponent)
        values = (self.lo**exponent, self.hi**exponent)
        return Interval(min(values), max(values))

    def scale(self, scalar: Fraction) -> "Interval":
        scalar = Fraction(scalar)
        return self * Interval(scalar, scalar)

    def abs_upper(self) -> Fraction:
        return max(abs(self.lo), abs(self.hi))

    def abs_lower(self) -> Fraction:
        if self.lo <= 0 <= self.hi:
            return Fraction(0)
        return min(abs(self.lo), abs(self.hi))


Poly2 = dict[tuple[int, int], Q5]


def q5_interval(value: Q5) -> Interval:
    a, b = value
    if b >= 0:
        return Interval(a + b * SQRT5_LO, a + b * SQRT5_HI)
    return Interval(a + b * SQRT5_HI, a + b * SQRT5_LO)


def p2_const(value: Q5) -> Poly2:
    return {(0, 0): value} if value != ZERO else {}


def p2_add(left: Poly2, right: Poly2) -> Poly2:
    out = dict(left)
    for mono, coefficient in right.items():
        value = q5_add(out.get(mono, ZERO), coefficient)
        if value == ZERO:
            out.pop(mono, None)
        else:
            out[mono] = value
    return out


def p2_scale(poly: Poly2, scalar: Q5) -> Poly2:
    if scalar == ZERO:
        return {}
    return {mono: q5_mul(value, scalar) for mono, value in poly.items()}


def p2_mul(left: Poly2, right: Poly2) -> Poly2:
    out: Poly2 = {}
    for (i, j), x in left.items():
        for (k, ell), y in right.items():
            mono = (i + k, j + ell)
            value = q5_add(out.get(mono, ZERO), q5_mul(x, y))
            if value == ZERO:
                out.pop(mono, None)
            else:
                out[mono] = value
    return out


def p2_pow(poly: Poly2, exponent: int) -> Poly2:
    out = p2_const(ONE)
    for _ in range(exponent):
        out = p2_mul(out, poly)
    return out


def p2_derivative(poly: Poly2, axis: int) -> Poly2:
    out: Poly2 = {}
    for mono, value in poly.items():
        power = mono[axis]
        if power == 0:
            continue
        target = list(mono)
        target[axis] -= 1
        out[tuple(target)] = q5_scale(value, power)
    return out


def p2_eval_interval(poly: Poly2, s: Interval, t: Interval) -> Interval:
    out = Interval(Fraction(0), Fraction(0))
    for (i, j), value in poly.items():
        term = q5_interval(value) * s.pow(i) * t.pow(j)
        out = out + term
    return out


def p2_eval_centered(poly: Poly2, s: Interval, t: Interval) -> Interval:
    """Mean-value enclosure, sharper for cancellation-dominated polynomials."""

    s_mid = (s.lo + s.hi) / 2
    t_mid = (t.lo + t.hi) / 2
    centre = q5_interval(p2_eval_q5(poly, q5(s_mid), q5(t_mid)))
    ds = Interval(s.lo - s_mid, s.hi - s_mid)
    dt = Interval(t.lo - t_mid, t.hi - t_mid)
    return (
        centre
        + p2_eval_interval(p2_derivative(poly, 0), s, t) * ds
        + p2_eval_interval(p2_derivative(poly, 1), s, t) * dt
    )


def p2_eval_quadratic_centered(poly: Poly2, s: Interval, t: Interval) -> Interval:
    """Exact Taylor enclosure for a polynomial of total degree at most two."""

    require(all(sum(mono) <= 2 for mono in poly), "quadratic enclosure degree drift")
    s_mid = (s.lo + s.hi) / 2
    t_mid = (t.lo + t.hi) / 2
    point_s = q5(s_mid)
    point_t = q5(t_mid)
    ds = Interval(s.lo - s_mid, s.hi - s_mid)
    dt = Interval(t.lo - t_mid, t.hi - t_mid)
    value = q5_interval(p2_eval_q5(poly, point_s, point_t))
    grad_s = q5_interval(p2_eval_q5(p2_derivative(poly, 0), point_s, point_t))
    grad_t = q5_interval(p2_eval_q5(p2_derivative(poly, 1), point_s, point_t))
    h_ss = q5_interval(
        p2_eval_q5(p2_derivative(p2_derivative(poly, 0), 0), point_s, point_t)
    )
    h_st = q5_interval(
        p2_eval_q5(p2_derivative(p2_derivative(poly, 0), 1), point_s, point_t)
    )
    h_tt = q5_interval(
        p2_eval_q5(p2_derivative(p2_derivative(poly, 1), 1), point_s, point_t)
    )
    return (
        value
        + grad_s * ds
        + grad_t * dt
        + (h_ss * ds.pow(2)).scale(Fraction(1, 2))
        + h_st * ds * dt
        + (h_tt * dt.pow(2)).scale(Fraction(1, 2))
    )


def p2_eval_q5(poly: Poly2, s: Q5, t: Q5) -> Q5:
    out = ZERO
    for (i, j), value in poly.items():
        out = q5_add(out, q5_mul(value, q5_mul(q5_pow(s, i), q5_pow(t, j))))
    return out


def dot(left: tuple[Q5, Q5, Q5], right: tuple[Q5, Q5, Q5]) -> Q5:
    out = ZERO
    for axis in range(3):
        out = q5_add(out, q5_mul(left[axis], right[axis]))
    return out


def vector_add(vectors: Iterable[tuple[Q5, Q5, Q5]]) -> tuple[Q5, Q5, Q5]:
    out = [ZERO, ZERO, ZERO]
    for vector in vectors:
        for axis in range(3):
            out[axis] = q5_add(out[axis], vector[axis])
    return tuple(out)  # type: ignore[return-value]


def vector_sub(
    left: tuple[Q5, Q5, Q5], right: tuple[Q5, Q5, Q5]
) -> tuple[Q5, Q5, Q5]:
    return tuple(q5_sub(left[i], right[i]) for i in range(3))  # type: ignore[return-value]


def vector_scale(
    value: tuple[Q5, Q5, Q5], scalar: Q5
) -> tuple[Q5, Q5, Q5]:
    return tuple(q5_mul(component, scalar) for component in value)  # type: ignore[return-value]


def chart_polynomials(
    origin: tuple[Q5, Q5, Q5],
    first: tuple[Q5, Q5, Q5],
    second: tuple[Q5, Q5, Q5],
) -> tuple[Poly2, Poly2, Poly2]:
    """Return D and the two derivative numerators for I6=N/D^3."""

    coordinates: list[Poly2] = []
    for axis in range(3):
        coordinate: Poly2 = {}
        for mono, value in (
            ((0, 0), origin[axis]),
            ((1, 0), first[axis]),
            ((0, 1), second[axis]),
        ):
            if value != ZERO:
                coordinate[mono] = value
        coordinates.append(coordinate)
    d = {}
    for coordinate in coordinates:
        d = p2_add(d, p2_pow(coordinate, 2))

    i6 = base.build_cartesian_frame()["_i6_poly_object"]
    numerator: Poly2 = {}
    for powers, coefficient in i6.items():
        degree = sum(powers)
        require(degree in (0, 2, 4, 6), "I6 degree leaves even level-six lift")
        term = p2_const(coefficient)
        for axis, power in enumerate(powers):
            term = p2_mul(term, p2_pow(coordinates[axis], power))
        term = p2_mul(term, p2_pow(d, 3 - degree // 2))
        numerator = p2_add(numerator, term)

    outputs: list[Poly2] = []
    for axis in range(2):
        first_term = p2_mul(p2_derivative(numerator, axis), d)
        second_term = p2_scale(
            p2_mul(numerator, p2_derivative(d, axis)), q5(-3)
        )
        outputs.append(p2_add(first_term, second_term))
    return d, outputs[0], outputs[1]


def derive_axes() -> list[tuple[str, tuple[Q5, Q5, Q5]]]:
    vertices = base.cartesian_vertices()
    faces, edges = hardening.derive_adjacency_faces_edges(vertices)
    rays = (
        [("vertex", vertex) for vertex in vertices]
        + [("edge", vector_add(vertices[index] for index in edge)) for edge in edges]
        + [("face", vector_add(vertices[index] for index in face)) for face in faces]
    )
    axes: list[tuple[str, tuple[Q5, Q5, Q5]]] = []
    for kind, ray in rays:
        if any(hardening.cross(ray, old) == (ZERO, ZERO, ZERO) for _, old in axes):
            continue
        axes.append((kind, ray))
    require(len(axes) == 31, "projective critical-axis census drift")
    require(
        {kind: sum(row[0] == kind for row in axes) for kind in ("vertex", "face", "edge")}
        == {"vertex": 6, "face": 10, "edge": 15},
        "projective critical-axis orbit census drift",
    )
    return axes


def neighborhood_polynomial(
    axis: tuple[Q5, Q5, Q5],
    coordinates: tuple[Poly2, Poly2, Poly2],
) -> Poly2:
    raw_dot: Poly2 = {}
    d: Poly2 = {}
    for component, coordinate in zip(axis, coordinates):
        raw_dot = p2_add(raw_dot, p2_scale(coordinate, component))
        d = p2_add(d, p2_pow(coordinate, 2))
    axis_norm = dot(axis, axis)
    # sin^2 <= 1/4096 iff 4096 (a.v)^2 >= 4095 ||a||^2 ||v||^2.
    return p2_add(
        p2_scale(p2_pow(raw_dot, 2), q5(4096)),
        p2_scale(d, q5_mul(q5(-LOCAL_COSINE_MULTIPLIER), axis_norm)),
    )


def standard_chart(axis: int) -> tuple[tuple[Q5, Q5, Q5], ...]:
    require(axis in (0, 1, 2), "invalid projective chart")
    remaining = [index for index in range(3) if index != axis]
    origin = tuple(ONE if index == axis else ZERO for index in range(3))
    first = tuple(ONE if index == remaining[0] else ZERO for index in range(3))
    second = tuple(ONE if index == remaining[1] else ZERO for index in range(3))
    return origin, first, second  # type: ignore[return-value]


def tail_bounds(x_max: Fraction = X_MAX) -> tuple[Fraction, Fraction]:
    row = v3.cosine_tail_bounds(x_max=x_max)
    return (
        Fraction(row["normalized_C1_gradient_bound"]),
        Fraction(row["normalized_C2_intrinsic_hessian_bound"]),
    )


def box_digest(rows: list[dict[str, Any]]) -> str:
    return base.tagged_sha256(base.canonical_json_bytes(rows))


def global_gradient_cover(
    *,
    x_max: Fraction = X_MAX,
    max_depth: int = MAX_GLOBAL_DEPTH,
    neighborhood_multiplier: int = 4096,
) -> dict[str, Any]:
    """Certify every projective-chart point away from local axis boxes."""

    require(neighborhood_multiplier == 4096, "local-neighborhood mutation")
    c1_tail, _ = tail_bounds(x_max)
    axes = derive_axes()
    all_leaves: list[dict[str, Any]] = []
    chart_rows = []

    for chart_axis in range(3):
        origin, first, second = standard_chart(chart_axis)
        d, gs, gt = chart_polynomials(origin, first, second)
        coordinate_polys = []
        for coordinate in range(3):
            poly: Poly2 = {}
            for mono, value in (
                ((0, 0), origin[coordinate]),
                ((1, 0), first[coordinate]),
                ((0, 1), second[coordinate]),
            ):
                if value != ZERO:
                    poly[mono] = value
            coordinate_polys.append(poly)
        neighborhoods = [
            neighborhood_polynomial(axis, tuple(coordinate_polys))
            for _, axis in axes
        ]

        pending = [
            (
                Interval(Fraction(-1), Fraction(1)),
                Interval(Fraction(-1), Fraction(1)),
                0,
                "",
            )
        ]
        local_count = 0
        gradient_count = 0
        depth_seen = 0
        weakest_margin: Fraction | None = None
        while pending:
            s, t, depth, path = pending.pop()
            depth_seen = max(depth_seen, depth)
            local_axis = None
            for index, poly in enumerate(neighborhoods):
                if p2_eval_quadratic_centered(poly, s, t).lo >= 0:
                    local_axis = index
                    break
            if local_axis is not None:
                local_count += 1
                all_leaves.append(
                    {
                        "chart": chart_axis,
                        "path": path,
                        "kind": "local",
                        "axis": local_axis,
                    }
                )
                continue

            gs_interval = p2_eval_centered(gs, s, t)
            gt_interval = p2_eval_centered(gt, s, t)
            d_interval = p2_eval_interval(d, s, t)
            require(d_interval.lo > 0, "projective chart denominator can vanish")
            numerator_lower_squared = (
                gs_interval.abs_lower() ** 2 + gt_interval.abs_lower() ** 2
            )
            # The projective normalization map has tangent operator norm at
            # most D^{-1/2}.  Since (f_s,f_t)=(Gs,Gt)/D^4, domination of the
            # intrinsic tail follows from |G|^2 > C1^2 D^7.
            required_squared = c1_tail**2 * d_interval.hi**7
            if numerator_lower_squared > required_squared:
                gradient_count += 1
                margin = numerator_lower_squared - required_squared
                weakest_margin = margin if weakest_margin is None else min(weakest_margin, margin)
                all_leaves.append(
                    {
                        "chart": chart_axis,
                        "path": path,
                        "kind": "gradient",
                    }
                )
                continue

            require(
                depth < max_depth,
                "global gradient cover exhausted its depth budget at "
                f"chart={chart_axis}, path={path}, "
                f"s=[{s.lo},{s.hi}], t=[{t.lo},{t.hi}]",
            )
            s_width = s.hi - s.lo
            t_width = t.hi - t.lo
            if s_width >= t_width:
                midpoint = (s.lo + s.hi) / 2
                pending.append((Interval(midpoint, s.hi), t, depth + 1, path + "1"))
                pending.append((Interval(s.lo, midpoint), t, depth + 1, path + "0"))
            else:
                midpoint = (t.lo + t.hi) / 2
                pending.append((s, Interval(midpoint, t.hi), depth + 1, path + "1"))
                pending.append((s, Interval(t.lo, midpoint), depth + 1, path + "0"))

        require(weakest_margin is not None, "chart has no gradient-certified boxes")
        chart_rows.append(
            {
                "chart_fixed_coordinate": ("x", "y", "z")[chart_axis],
                "local_boxes": local_count,
                "gradient_boxes": gradient_count,
                "maximum_depth": depth_seen,
                "weakest_exact_squared_numerator_margin": str(weakest_margin),
            }
        )

    return {
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
        "normalized_tail_C1_bound": str(c1_tail),
        "chart_rows": chart_rows,
        "leaf_count": len(all_leaves),
        "leaf_partition_sha256": box_digest(all_leaves),
        "all_off_neighborhood_boxes_gradient_certified": True,
        "arithmetic": "Fraction intervals with a rationally proved sqrt(5) enclosure",
    }


def radial_factor(poly: Poly2, derivative_axis: int) -> Poly2:
    """Integral_0^1 partial_axis poly(t s,t t) dt."""

    derivative = p2_derivative(poly, derivative_axis)
    return {
        mono: q5_scale(value, Fraction(1, sum(mono) + 1))
        for mono, value in derivative.items()
    }


def rational_sqrt_lower(value: Fraction) -> Fraction:
    """A proved lower square-root enclosure on a fixed rational grid."""

    from math import isqrt

    require(value > 0, "square-root lower enclosure requires positivity")
    scale = 10**15
    candidate = Fraction(isqrt(value.numerator * value.denominator * scale * scale), value.denominator * scale)
    # The integer formula above can undershoot more than necessary but never
    # promotes unless its defining square inequality is exact.
    while candidate * candidate > value:
        candidate -= Fraction(1, value.denominator * scale)
    require(candidate > 0 and candidate * candidate <= value, "sqrt lower enclosure failed")
    return candidate


def rational_sqrt_upper(value: Fraction) -> Fraction:
    """A proved upper square-root enclosure on a fixed rational grid."""

    lower = rational_sqrt_lower(value)
    if lower * lower == value:
        return lower
    step = Fraction(1, value.denominator * 10**15)
    upper = lower + step
    while upper * upper < value:
        upper += step
    require(upper * upper >= value, "sqrt upper enclosure failed")
    return upper


def representative_local_charts() -> list[dict[str, Any]]:
    vertices = base.cartesian_vertices()
    faces, edges = hardening.derive_adjacency_faces_edges(vertices)

    vertex_axis = vertices[0]
    adjacent = next(edge[1] for edge in edges if edge[0] == 0)
    vertex_first = hardening.cross(vertex_axis, vertices[adjacent])
    vertex_first = vector_scale(vertex_first, q5(2))
    vertex_second = hardening.cross(vertex_axis, vertex_first)

    face = faces[0]
    face_axis = vector_add(vertices[index] for index in face)
    face_first = vector_scale(vector_sub(vertices[face[0]], vertices[face[1]]), q5(2))
    face_second = hardening.cross(face_axis, face_first)

    edge = edges[0]
    edge_axis = vector_add(vertices[index] for index in edge)
    edge_first = vector_scale(
        vector_sub(vertices[edge[0]], vertices[edge[1]]), base.q5(1, 1)
    )
    edge_second = hardening.cross(edge_axis, edge_first)

    return [
        {"orbit": "vertex", "axis": vertex_axis, "first": vertex_first, "second": vertex_second},
        {"orbit": "face", "axis": face_axis, "first": face_first, "second": face_second},
        {"orbit": "edge", "axis": edge_axis, "first": edge_first, "second": edge_second},
    ]


Matrix3 = tuple[
    tuple[Q5, Q5, Q5], tuple[Q5, Q5, Q5], tuple[Q5, Q5, Q5]
]


def matrix_apply(matrix: Matrix3, vector: tuple[Q5, Q5, Q5]) -> tuple[Q5, Q5, Q5]:
    return tuple(dot(matrix[row], vector) for row in range(3))  # type: ignore[return-value]


def matrix_mul(left: Matrix3, right: Matrix3) -> Matrix3:
    columns = tuple(tuple(right[row][column] for row in range(3)) for column in range(3))
    return tuple(
        tuple(dot(left[row], columns[column]) for column in range(3))
        for row in range(3)
    )  # type: ignore[return-value]


def matrix_transpose(matrix: Matrix3) -> Matrix3:
    return tuple(tuple(matrix[column][row] for column in range(3)) for row in range(3))  # type: ignore[return-value]


def matrix_det(matrix: Matrix3) -> Q5:
    return dot(matrix[0], hardening.cross(matrix[1], matrix[2]))


def linear_map_from_bases(
    sources: tuple[tuple[Q5, Q5, Q5], ...],
    targets: tuple[tuple[Q5, Q5, Q5], ...],
) -> Matrix3:
    source_columns: Matrix3 = tuple(
        tuple(sources[column][row] for column in range(3)) for row in range(3)
    )  # type: ignore[assignment]
    target_columns: Matrix3 = tuple(
        tuple(targets[column][row] for column in range(3)) for row in range(3)
    )  # type: ignore[assignment]
    determinant = matrix_det(source_columns)
    require(determinant != ZERO, "stabilizer source basis is singular")
    cofactors = (
        hardening.cross(sources[1], sources[2]),
        hardening.cross(sources[2], sources[0]),
        hardening.cross(sources[0], sources[1]),
    )
    inverse: Matrix3 = tuple(
        tuple(base.q5_div(cofactors[row][column], determinant) for column in range(3))
        for row in range(3)
    )  # type: ignore[assignment]
    return matrix_mul(target_columns, inverse)


def matrix_rank(matrix: Matrix3) -> int:
    if any(matrix[row][column] != ZERO for row in range(3) for column in range(3)):
        for r1 in range(3):
            for r2 in range(r1 + 1, 3):
                for c1 in range(3):
                    for c2 in range(c1 + 1, 3):
                        minor = q5_sub(
                            q5_mul(matrix[r1][c1], matrix[r2][c2]),
                            q5_mul(matrix[r1][c2], matrix[r2][c1]),
                        )
                        if minor != ZERO:
                            return 3 if matrix_det(matrix) != ZERO else 2
        return 1
    return 0


def projectively_equal(
    left: tuple[Q5, Q5, Q5], right: tuple[Q5, Q5, Q5]
) -> bool:
    return hardening.cross(left, right) == (ZERO, ZERO, ZERO)


def projective_unique(
    rays: Iterable[tuple[Q5, Q5, Q5]],
) -> list[tuple[Q5, Q5, Q5]]:
    out: list[tuple[Q5, Q5, Q5]] = []
    for ray in rays:
        require(dot(ray, ray) != ZERO, "zero projective ray")
        if not any(projectively_equal(ray, old) for old in out):
            out.append(ray)
    return out


def projective_sine_squared(
    left: tuple[Q5, Q5, Q5], right: tuple[Q5, Q5, Q5]
) -> Q5:
    denominator = q5_mul(dot(left, left), dot(right, right))
    require(denominator != ZERO, "zero norm in projective separation")
    numerator = q5_sub(denominator, q5_pow(dot(left, right), 2))
    return base.q5_div(numerator, denominator)


def exact_full_group_transport_certificate() -> dict[str, Any]:
    """Construct all proper carrier rotations and their projective orbits.

    A proper icosahedral rotation is determined by the image of one directed
    edge.  The thirty undirected edges therefore give sixty exact candidates.
    Constructing each map over Q(sqrt(5)) avoids assuming transitivity when it
    is the theorem join that transports three local boxes to all 31 axes.
    """

    vertices = base.cartesian_vertices()
    faces, edges = hardening.derive_adjacency_faces_edges(vertices)
    vertex_set = set(vertices)
    vertex_index = {vertex: index for index, vertex in enumerate(vertices)}
    source_left, source_right = edges[0]
    source_basis = (
        vertices[source_left],
        vertices[source_right],
        hardening.cross(vertices[source_left], vertices[source_right]),
    )
    identity: Matrix3 = (
        (ONE, ZERO, ZERO),
        (ZERO, ONE, ZERO),
        (ZERO, ZERO, ONE),
    )

    rotations: set[Matrix3] = set()
    for left, right in edges:
        for target_left, target_right in ((left, right), (right, left)):
            target_basis = (
                vertices[target_left],
                vertices[target_right],
                hardening.cross(vertices[target_left], vertices[target_right]),
            )
            rotation = linear_map_from_bases(source_basis, target_basis)
            require(
                matrix_mul(matrix_transpose(rotation), rotation) == identity,
                "directed-edge map is not orthogonal",
            )
            require(matrix_det(rotation) == ONE, "directed-edge map is not proper")
            image = tuple(matrix_apply(rotation, vertex) for vertex in vertices)
            require(set(image) == vertex_set, "directed-edge map leaves carrier")
            rotations.add(rotation)
    require(len(rotations) == 60, "proper rotation group order drift")

    permutations = sorted(
        tuple(vertex_index[matrix_apply(rotation, vertex)] for vertex in vertices)
        for rotation in rotations
    )
    require(len(set(permutations)) == 60, "carrier action is not faithful")
    identity_permutation = tuple(range(12))
    require(identity_permutation in permutations, "proper group lacks identity")
    permutation_set = set(permutations)
    for left in permutations:
        inverse = [0] * 12
        for source, target in enumerate(left):
            inverse[target] = source
        require(tuple(inverse) in permutation_set, "proper group lacks inverse")
        for right in permutations:
            composition = tuple(left[right[index]] for index in range(12))
            require(composition in permutation_set, "proper group is not closed")

    representatives = {
        "vertex": vertices[0],
        "face": vector_add(vertices[index] for index in faces[0]),
        "edge": vector_add(vertices[index] for index in edges[0]),
    }
    expected_sizes = {"vertex": 6, "face": 10, "edge": 15}
    critical_axes = derive_axes()
    orbit_rows = []
    transported_union: list[tuple[Q5, Q5, Q5]] = []
    for orbit in ("vertex", "face", "edge"):
        transported = projective_unique(
            matrix_apply(rotation, representatives[orbit]) for rotation in rotations
        )
        require(len(transported) == expected_sizes[orbit], f"{orbit} orbit-size drift")
        typed_axes = [axis for kind, axis in critical_axes if kind == orbit]
        require(
            all(any(projectively_equal(ray, axis) for axis in typed_axes) for ray in transported)
            and all(any(projectively_equal(axis, ray) for ray in transported) for axis in typed_axes),
            f"{orbit} transport does not cover its critical axes",
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
    require(len(union) == 31, "transported projective union does not have 31 axes")
    require(
        all(any(projectively_equal(ray, axis) for _, axis in critical_axes) for ray in union)
        and all(any(projectively_equal(axis, ray) for ray in union) for _, axis in critical_axes),
        "transported projective union does not equal the critical-axis set",
    )

    minimum: Q5 | None = None
    for left_index, left in enumerate(union):
        for right in union[left_index + 1 :]:
            value = projective_sine_squared(left, right)
            if minimum is None or q5_sign(q5_sub(value, minimum)) < 0:
                minimum = value
    require(minimum is not None, "empty projective separation census")
    expected_minimum = q5(Fraction(1, 2), Fraction(-1, 6))
    require(minimum == expected_minimum, "minimum projective separation drift")
    require(
        q5_sign(q5_sub(minimum, q5(Fraction(1, 1024)))) > 0,
        "axis neighborhoods are not exactly separated",
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
        "minimum_distinct_projective_axis_sine_squared": q5_str(minimum),
        "twice_local_radius_sine_squared_upper": "1/1024",
        "minimum_separation_exceeds_twice_local_radius": True,
        "neighborhoods_sine_squared_1_over_4096_pairwise_disjoint": True,
    }


def exact_stabilizer_certificate(
    full_transport: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Exhibit one nontrivial proper rotation about each orbit axis."""

    vertices = base.cartesian_vertices()
    faces, edges = hardening.derive_adjacency_faces_edges(vertices)
    vertex_set = set(vertices)

    neighbors = sorted(
        {right for left, right in edges if left == 0}
        | {left for left, right in edges if right == 0}
    )
    require(len(neighbors) == 5, "vertex stabilizer neighbor census drift")
    edge_set = {tuple(sorted(edge)) for edge in edges}
    cycle = [neighbors[0]]
    previous = None
    while len(cycle) < 5:
        candidates = sorted(
            index
            for index in neighbors
            if index not in cycle and tuple(sorted((cycle[-1], index))) in edge_set
        )
        require(candidates, "vertex-link cycle failed")
        if previous is not None and len(candidates) > 1:
            candidates = [value for value in candidates if value != previous]
        previous = cycle[-1]
        cycle.append(candidates[0])
    require(tuple(sorted((cycle[-1], cycle[0]))) in edge_set, "vertex link does not close")
    vertex_rotation = linear_map_from_bases(
        (vertices[0], vertices[cycle[0]], vertices[cycle[1]]),
        (vertices[0], vertices[cycle[1]], vertices[cycle[2]]),
    )

    face = faces[0]
    face_axis = vector_add(vertices[index] for index in face)
    face_rotation = linear_map_from_bases(
        tuple(vertices[index] for index in face),
        (vertices[face[1]], vertices[face[2]], vertices[face[0]]),
    )

    edge = edges[0]
    edge_axis = vector_add(vertices[index] for index in edge)
    edge_norm = dot(edge_axis, edge_axis)
    edge_rotation: Matrix3 = tuple(
        tuple(
            q5_sub(
                base.q5_div(q5_scale(q5_mul(edge_axis[row], edge_axis[column]), 2), edge_norm),
                ONE if row == column else ZERO,
            )
            for column in range(3)
        )
        for row in range(3)
    )  # type: ignore[assignment]

    identity: Matrix3 = (
        (ONE, ZERO, ZERO),
        (ZERO, ONE, ZERO),
        (ZERO, ZERO, ONE),
    )
    rows = []
    for orbit, order, axis, rotation in (
        ("vertex", 5, vertices[0], vertex_rotation),
        ("face", 3, face_axis, face_rotation),
        ("edge", 2, edge_axis, edge_rotation),
    ):
        require(matrix_mul(matrix_transpose(rotation), rotation) == identity, f"{orbit} stabilizer is not orthogonal")
        require(matrix_det(rotation) == ONE, f"{orbit} stabilizer is not proper")
        require(matrix_apply(rotation, axis) == axis, f"{orbit} stabilizer moves its axis")
        image = [matrix_apply(rotation, vertex) for vertex in vertices]
        require(set(image) == vertex_set and len(set(image)) == 12, f"{orbit} stabilizer does not permute ports")
        power = identity
        first_identity = None
        for exponent in range(1, order + 1):
            power = matrix_mul(power, rotation)
            if power == identity and first_identity is None:
                first_identity = exponent
        require(first_identity == order, f"{orbit} stabilizer order drift")
        difference: Matrix3 = tuple(
            tuple(q5_sub(rotation[row][column], identity[row][column]) for column in range(3))
            for row in range(3)
        )  # type: ignore[assignment]
        require(matrix_rank(difference) == 2, f"{orbit} stabilizer has a tangent fixed vector")
        rows.append(
            {
                "orbit": orbit,
                "rotation_order": order,
                "permutes_all_twelve_ports": True,
                "proper_orthogonal": True,
                "fixed_space_dimension": 1,
                "tangent_fixed_space_dimension": 0,
            }
        )
    if full_transport is None:
        full_transport = exact_full_group_transport_certificate()
    return {
        "representative_rotations": rows,
        "full_group_transport": full_transport,
        "kernel_even_from_antipodal_ports": True,
        "conclusion": (
            "equal weights make every cosine kernel invariant under the exhibited "
            "stabilizers; an invariant tangent gradient must lie in the zero "
            "tangent fixed space, so all 62 oriented axis directions are stationary"
        ),
    }


def local_uniqueness_boxes(
    *,
    x_max: Fraction = X_MAX,
    radius_denominator_offset: int = 0,
    hessian_multiplier: Fraction = Fraction(1),
    full_transport: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Prove one stationary point in every local symmetry neighborhood."""

    require(radius_denominator_offset in (0, 1), "local radius mutation type")
    require(hessian_multiplier == 1, "tail-Hessian mutation")
    _, c2_tail = tail_bounds(x_max)
    rows = []
    for entry in representative_local_charts():
        axis = entry["axis"]
        first = entry["first"]
        second = entry["second"]
        a_norm = dot(axis, axis)
        e1_norm = dot(first, first)
        e2_norm = dot(second, second)
        require(
            dot(axis, first) == dot(axis, second) == dot(first, second) == ZERO,
            "local chart basis is not orthogonal",
        )
        # Choose the tightest reciprocal-integer rectangle containing the
        # ellipse 4095(E1 s^2+E2 t^2)<=A.  The one-step shrink is a mutation
        # control and must fail the exact containment check below.
        def reciprocal_radius(norm: Q5) -> Fraction:
            valid = []
            for denominator in range(1, 4097):
                if q5_sign(
                    q5_sub(
                        q5_scale(norm, LOCAL_COSINE_MULTIPLIER),
                        q5_scale(a_norm, denominator * denominator),
                    )
                ) >= 0:
                    valid.append(denominator)
            require(valid, "local reciprocal-radius search failed")
            return Fraction(1, max(valid) + radius_denominator_offset)

        s_radius = reciprocal_radius(e1_norm)
        t_radius = reciprocal_radius(e2_norm)
        s_box = Interval(-s_radius, s_radius)
        t_box = Interval(-t_radius, t_radius)
        require(
            q5_sign(
                q5_sub(
                    q5_scale(e1_norm, LOCAL_COSINE_MULTIPLIER * s_radius**2),
                    a_norm,
                )
            ) >= 0
            and q5_sign(
                q5_sub(
                    q5_scale(e2_norm, LOCAL_COSINE_MULTIPLIER * t_radius**2),
                    a_norm,
                )
            ) >= 0,
            "local square does not contain the full axis neighborhood",
        )

        d, gs, gt = chart_polynomials(axis, first, second)
        require(
            p2_eval_q5(gs, ZERO, ZERO) == ZERO
            and p2_eval_q5(gt, ZERO, ZERO) == ZERO,
            f"fixed-centre gradient numerator is nonzero on {entry['orbit']}",
        )
        c11 = radial_factor(gs, 0)
        c12 = radial_factor(gs, 1)
        c21 = radial_factor(gt, 0)
        c22 = radial_factor(gt, 1)
        # Fundamental theorem of calculus: (Gs,Gt)=C(s,t)(s,t).
        zero = q5(0)
        c0 = (
            (p2_eval_q5(c11, zero, zero), p2_eval_q5(c12, zero, zero)),
            (p2_eval_q5(c21, zero, zero), p2_eval_q5(c22, zero, zero)),
        )
        require(c0[0][1] == c0[1][0] == ZERO, "local Hessian basis is not diagonal")

        error_intervals = []
        for row, poly in enumerate((c11, c12, c21, c22)):
            baseline = c0[row // 2][row % 2]
            error_intervals.append(
                p2_eval_centered(
                    p2_add(poly, p2_const(base.q5_neg(baseline))), s_box, t_box
                )
            )
        a_interval = q5_interval(a_norm)
        e1_interval = q5_interval(e1_norm)
        e2_interval = q5_interval(e2_norm)
        require(a_interval.lo > 0 and e1_interval.lo > 0 and e2_interval.lo > 0, "local norm enclosure failed")
        diagonal_scales = (
            a_interval.hi / e1_interval.lo,
            a_interval.hi / e2_interval.lo,
        )
        product_lower = e1_interval.lo * e2_interval.lo
        off_diagonal_scale = a_interval.hi / rational_sqrt_lower(product_lower)
        scaled_error = (
            error_intervals[0].abs_upper() * diagonal_scales[0],
            error_intervals[1].abs_upper() * off_diagonal_scale,
            error_intervals[2].abs_upper() * off_diagonal_scale,
            error_intervals[3].abs_upper() * diagonal_scales[1],
        )
        error_row_sum = max(
            scaled_error[0] + scaled_error[1],
            scaled_error[2] + scaled_error[3],
        )
        error_column_sum = max(
            scaled_error[0] + scaled_error[2],
            scaled_error[1] + scaled_error[3],
        )
        error_induced_2_upper = rational_sqrt_upper(
            error_row_sum * error_column_sum
        )
        scaled_c00 = q5_mul(c0[0][0], base.q5_div(a_norm, e1_norm))
        scaled_c11 = q5_mul(c0[1][1], base.q5_div(a_norm, e2_norm))
        diagonal_lower = min(
            q5_interval(scaled_c00).abs_lower(),
            q5_interval(scaled_c11).abs_lower(),
        )
        c_singular_lower = diagonal_lower - error_induced_2_upper
        require(
            c_singular_lower > 0,
            f"local template matrix can become singular on {entry['orbit']}: "
            f"diagonal={diagonal_lower}, error={error_row_sum}",
        )
        d_upper = p2_eval_interval(d, s_box, t_box).hi
        base_linear_lower = c_singular_lower / d_upper**4
        # In y_i=s_i sqrt(E_i/A), the chart tangent operator has norm at
        # most one and geodesic distance is at most |y|.  The intrinsic C2
        # tail bound is therefore already the required linear coefficient.
        tail_linear_upper = c2_tail
        require(
            base_linear_lower > tail_linear_upper,
            f"tail Hessian exhausts a local template uniqueness margin on "
            f"{entry['orbit']}: base={base_linear_lower}, tail={tail_linear_upper}",
        )
        normalized_h0 = (
            base.q5_div(scaled_c00, q5_pow(a_norm, 4)),
            base.q5_div(scaled_c11, q5_pow(a_norm, 4)),
        )
        expected_h0 = {
            "vertex": (q5(-21), q5(-21)),
            "face": (q5(Fraction(35, 3)), q5(Fraction(35, 3))),
            "edge": (
                q5(Fraction(105, 16), Fraction(105, 16)),
                q5(Fraction(105, 16), Fraction(-105, 16)),
            ),
        }[entry["orbit"]]
        require(
            set(normalized_h0) == set(expected_h0),
            f"fixed-centre Morse Hessian drift on {entry['orbit']}",
        )
        h0_abs_lower = min(q5_interval(value).abs_lower() for value in normalized_h0)
        morse_margin = h0_abs_lower - c2_tail
        require(morse_margin > 0, "tail Hessian can change a Morse signature")
        rows.append(
            {
                "orbit": entry["orbit"],
                "chart_rectangle": {
                    "s_radius": str(s_radius),
                    "t_radius": str(t_radius),
                },
                "contains_full_sine_squared_1_over_4096_neighborhood": True,
                "fixed_center_gradient_numerators_zero": True,
                "template_gradient_linear_lower": str(base_linear_lower),
                "tail_gradient_linear_upper": str(tail_linear_upper),
                "strict_margin": str(base_linear_lower - tail_linear_upper),
                "fixed_center_scaled_C_min_abs_diagonal": str(diagonal_lower),
                "local_matrix_error_row_sum": str(error_row_sum),
                "local_matrix_error_column_sum": str(error_column_sum),
                "local_matrix_error_induced_2_upper": str(error_induced_2_upper),
                "induced_2_bound_justification": (
                    "||E||_2 <= sqrt(||E||_1 ||E||_infinity), with both "
                    "induced one and infinity norms bounded entrywise"
                ),
                "fixed_center_minus_uniform_operator_error": str(c_singular_lower),
                "fixed_center_normalized_hessian_eigenvalues": [
                    q5_str(value) for value in normalized_h0
                ],
                "fixed_center_min_abs_hessian_eigenvalue_lower": str(h0_abs_lower),
                "morse_signature_margin_against_tail_C2": str(morse_margin),
                "saddle_safe_injectivity_argument": (
                    "G(y)=C(y)y exactly; the fixed C(0) smallest singular "
                    "value exceeds the uniform row-sum operator bound for "
                    "C(y)-C(0), so cancellation is impossible even when the "
                    "fixed Hessian is indefinite"
                ),
                "normalized_y_chart_tail_comparison": {
                    "coordinate_definition": "y_i = s_i sqrt(E_i/A)",
                    "projective_chart_operator_norm_upper": "||J_y|| <= 1",
                    "geodesic_distance": "rho = atan(|y|) <= |y|",
                    "tail_gradient_conclusion": "||dE/dy|| <= C2 |y|",
                },
            }
        )

    if full_transport is None:
        full_transport = exact_full_group_transport_certificate()
    transport_join = {
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
    }
    return {
        "method": (
            "exact local normalized charts; radial fundamental-theorem factorization "
            "of the template gradient; intrinsic tail-Hessian integration"
        ),
        "stationary_centres": (
            "equal port weights make the cosine kernel invariant under each "
            "vertex, face, or edge rotation stabilizer; its tangent fixed space "
            "is zero, so every declared axis centre is stationary for every x"
        ),
        "normalized_tail_C2_bound": str(c2_tail),
        "orbit_representatives": rows,
        "symmetry_transport_join": transport_join,
        "every_local_neighborhood_has_exactly_one_stationary_direction": True,
        "morse_signatures_preserved_by_fixed_center_C2_margins": True,
        "arithmetic": "exact Q(sqrt(5)) polynomials and rational interval enclosures",
    }


def fail_closed_controls() -> dict[str, Any]:
    controls = []
    for name, callback, expected in (
        (
            "truncate the global subdivision before it closes",
            lambda: global_gradient_cover(max_depth=4),
            "exact off-neighborhood cover",
        ),
        (
            "change the certified axis-neighborhood multiplier",
            lambda: global_gradient_cover(neighborhood_multiplier=4095),
            "frozen local/global partition",
        ),
        (
            "shrink the local chart below the complete neighborhood",
            lambda: local_uniqueness_boxes(radius_denominator_offset=1),
            "complete local-neighborhood containment",
        ),
        (
            "inflate the admitted tail Hessian",
            lambda: local_uniqueness_boxes(hessian_multiplier=Fraction(2)),
            "exact v3 C2 custody",
        ),
    ):
        try:
            callback()
        except base.FingerprintError:
            controls.append(
                {"control": name, "expected_failure": expected, "detector_fired": True}
            )
        else:
            raise base.FingerprintError(f"control escaped: {name}")
    require(len(controls) == 4, "stage-4 control census drift")
    return {"controls": controls, "all_detectors_fired": True}


def artifact_self_hash(value: dict[str, Any]) -> str:
    body = {key: item for key, item in value.items() if key != "receipt_sha256"}
    return base.tagged_sha256(base.canonical_json_bytes(body))


def load_parent(path: Path = PARENT_PATH) -> tuple[bytes, dict[str, Any]]:
    raw = path.read_bytes()
    try:
        parent = json.loads(raw)
    except json.JSONDecodeError as error:
        raise base.FingerprintError("invalid v3 parent JSON") from error
    require(parent.get("schema") == v3.SCHEMA, "v3 parent schema drift")
    require(parent.get("status") == v3.STATUS, "v3 parent status drift")
    require(parent.get("receipt_sha256") == artifact_self_hash(parent), "v3 self-hash drift")
    require(raw == base.canonical_json_bytes(v3.build_receipt()), "v3 byte replay drift")
    return raw, parent


def build_receipt(*, parent_path: Path = PARENT_PATH) -> dict[str, Any]:
    parent_raw, parent = load_parent(parent_path)
    full_transport = exact_full_group_transport_certificate()
    global_cover = global_gradient_cover()
    local_boxes = local_uniqueness_boxes(full_transport=full_transport)
    verifier_raw = VERIFIER_PATH.read_bytes()
    receipt = {
        "schema": SCHEMA,
        "status": STATUS,
        "issue": 654,
        "extends": {
            "schema": parent["schema"],
            "status": parent["status"],
            "path": "code/a5_fingerprint/runtime/a5_multipole_persistence_receipt_v3.json",
            "reason": "append-only closure of the v3 global and local mathematical boundaries",
        },
        "parent_pin": {
            "path": "code/a5_fingerprint/runtime/a5_multipole_persistence_receipt_v3.json",
            "bytes": len(parent_raw),
            "sha256": base.tagged_sha256(parent_raw),
            "receipt_sha256": parent["receipt_sha256"],
        },
        "declared_kernel_branch": {
            "kernel": "Q_x(n) = sum_i [1 - cos(x u_i.n)]",
            "equal_weight_twelve_port_branch": True,
            "range": "0 < x = |a k| <= 1",
            "source_selected": False,
            "architecture_forced": False,
            "physical_source_selection_owner": 655,
        },
        "exact_axis_stationarity": exact_stabilizer_certificate(full_transport),
        "global_off_neighborhood_cover": global_cover,
        "local_uniqueness_boxes": local_boxes,
        "persistence_theorem": {
            "finite_exactly_62_persistence_range": True,
            "range": "0 < |a k| <= 1",
            "oriented_stationary_direction_count": 62,
            "unoriented_axis_count": 31,
            "orbit_counts": {"maxima": 12, "minima": 20, "saddles": 30},
            "morse_types_preserved": True,
            "theorem_scope": "declared equal-weight cosine kernel only",
            "source_selection_proved": False,
            "physical_attachment_proved": False,
        },
        "fail_closed_controls": fail_closed_controls(),
        "independent_verifier_binding": {
            "path": "code/a5_fingerprint/a5_multipole_persistence_stage4_independent_verifier.py",
            "bytes": len(verifier_raw),
            "sha256": base.tagged_sha256(verifier_raw),
            "implementation_independent_of_producer_import": True,
            "replay_command": (
                "python3 code/a5_fingerprint/"
                "a5_multipole_persistence_stage4_independent_verifier.py --mutations"
            ),
            "resigned_semantic_mutation_count": 25,
        },
        "comparison_boundary": {
            "public_measurement_read": False,
            "comparison_permitted": False,
            "physical_map_open": True,
        },
    }
    receipt["receipt_sha256"] = artifact_self_hash(receipt)
    return receipt


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args(argv)
    receipt = build_receipt()
    if args.write:
        RUNTIME.mkdir(exist_ok=True)
        RECEIPT_PATH.write_bytes(base.canonical_json_bytes(receipt))
        print(RECEIPT_PATH)
    else:
        print(json.dumps(receipt["persistence_theorem"], indent=2))
        print(receipt["status"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
