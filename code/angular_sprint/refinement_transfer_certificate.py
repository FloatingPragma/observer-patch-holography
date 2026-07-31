#!/usr/bin/env python3
"""Issue #643 level-one refinement identifiability certificate.

The static certificate lane of issue #643 fixes exact angular templates on
the twelve icosahedral ports together with a four-dimensional interpolation
ambiguity in the degree<=3 band. That ambiguity statement is a property of
the twelve-port truncation. This producer decides, in exact Q(sqrt5)
arithmetic, what the first refinement level of the port tower does to it:

* the twelve icosahedron vertices (cyclic permutations of (0,+-1,+-phi))
  carry squared norm (5+sqrt5)/2; its Galois norm is 5, which is not a
  rational square, so the vertex radius is not an element of Q(sqrt5) and
  every sphere-evaluation condition at a vertex splits into two
  independent Q(sqrt5)-linear rows by parity of the homogeneity gap;
* the thirty edge-midpoint sums over adjacent vertex pairs carry squared
  norm 6+2*sqrt5 = (1+sqrt5)^2, so their unit projections have exact
  Q(sqrt5) coordinates and each one contributes a single Q(sqrt5)-linear
  evaluation row;
* rational harmonic bases for degrees {0,1,2,3,6,10,12} are exact
  Laplacian nullspaces of the homogeneous monomial spaces, with the 2L+1
  dimension count, termwise harmonicity, and full coefficient rank
  recomputed for every basis;
* exact Gaussian elimination over Q(sqrt5) records the kernel dimensions
  of the degree<=3 band under vertex-only, midpoint-only, and combined
  level-one conditions, and of the comb support band
  span{H_0,H_6,H_10,H_12} under the combined level-one conditions.

Recorded outcome: the degree<=3 band keeps a positive kernel under the
twelve-vertex split conditions and has kernel zero under the full
level-one condition set, so the static base-port ambiguity is a property
of the truncation, not of the refinement tower; the comb support band
keeps a positive kernel at level one. Whether the physical readout
exposes refined-port values is an open source premise owned by the repair
and refinement law and is not decided here. No sky transfer is selected,
no repair-semigroup intertwining theorem is proved, and no comparison or
public measurement data is read. No floating point number appears in any
computation or in the receipt.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from math import isqrt
from pathlib import Path
from typing import Any, Sequence

HERE = Path(__file__).resolve().parent
RUNTIME = HERE / "runtime"
RECEIPT_PATH = RUNTIME / "refinement_transfer_receipt.json"

SCHEMA = "oph.refinement_transfer_receipt.v1"
ISSUE = 643
STATUS = (
    "STATIC_BASE_PORT_TRANSFER_NONIDENTIFIABLE__"
    "CANONICAL_BAND_LEVEL_ONE_REFINEMENT_IDENTIFIABLE__"
    "PHYSICAL_TRANSFER_OPEN"
)

CANONICAL_DEGREES = (0, 1, 2, 3)
CANONICAL_EVEN_DEGREES = (0, 2)
CANONICAL_ODD_DEGREES = (1, 3)
COMB_DEGREES = (0, 6, 10, 12)
ALL_DEGREES = (0, 1, 2, 3, 6, 10, 12)


class CertificateError(ValueError):
    """The refinement transfer certificate refused to build."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise CertificateError(message)


def canonical_json_bytes(value: Any) -> bytes:
    return (
        json.dumps(
            value,
            allow_nan=False,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        )
        + "\n"
    ).encode("ascii")


def tagged_sha256(data: bytes) -> str:
    return "sha256:" + hashlib.sha256(data).hexdigest()


# ---------------------------------------------------------------------------
# Exact arithmetic in Q(sqrt5)
# ---------------------------------------------------------------------------


class F5:
    """An element a + b*sqrt5 of Q(sqrt5) with exact Fraction coefficients."""

    __slots__ = ("a", "b")

    def __init__(self, a: Any = 0, b: Any = 0) -> None:
        self.a = Fraction(a)
        self.b = Fraction(b)

    def __add__(self, other: "F5") -> "F5":
        return F5(self.a + other.a, self.b + other.b)

    def __sub__(self, other: "F5") -> "F5":
        return F5(self.a - other.a, self.b - other.b)

    def __neg__(self) -> "F5":
        return F5(-self.a, -self.b)

    def __mul__(self, other: "F5") -> "F5":
        return F5(
            self.a * other.a + 5 * self.b * other.b,
            self.a * other.b + self.b * other.a,
        )

    def conj(self) -> "F5":
        """The Galois conjugate sqrt5 -> -sqrt5."""

        return F5(self.a, -self.b)

    def inv(self) -> "F5":
        norm = self.a * self.a - 5 * self.b * self.b
        if norm == 0:
            raise ZeroDivisionError("zero element of Q(sqrt5)")
        return F5(self.a / norm, -self.b / norm)

    def __truediv__(self, other: "F5") -> "F5":
        return self * other.inv()

    def is_zero(self) -> bool:
        return self.a == 0 and self.b == 0

    def __eq__(self, other: object) -> bool:
        return isinstance(other, F5) and self.a == other.a and self.b == other.b

    def __hash__(self) -> int:
        return hash((self.a, self.b))

    def __repr__(self) -> str:
        return f"F5({self.a},{self.b})"

    def text(self) -> str:
        return f"{self.a}+{self.b}*sqrt5"


ZERO = F5(0)
ONE = F5(1)
PHI = F5(Fraction(1, 2), Fraction(1, 2))  # golden ratio (1+sqrt5)/2
R_VERTEX = F5(Fraction(5, 2), Fraction(1, 2))  # squared vertex norm 2+phi
R_MIDPOINT = F5(6, 2)  # squared midpoint-sum norm 6+2*sqrt5
ONE_PLUS_SQRT5 = F5(1, 1)  # exact midpoint-sum radius, (1+sqrt5)^2 = 6+2*sqrt5

Vec3 = tuple[F5, F5, F5]
RowF5 = list[F5]


def f5_pow(base: F5, exponent: int) -> F5:
    require(exponent >= 0, "negative exponent in exact power")
    out = ONE
    for _ in range(exponent):
        out = out * base
    return out


def dot(u: Vec3, v: Vec3) -> F5:
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]


# ---------------------------------------------------------------------------
# Exact linear algebra over Q(sqrt5)
# ---------------------------------------------------------------------------


def rref_f5(matrix: Sequence[RowF5]) -> tuple[list[RowF5], list[int]]:
    """Exact reduced row echelon form over Q(sqrt5); returns pivot columns.

    Pivot rows are selected among the eligible rows by fewest nonzero
    entries; the selection keeps every step exact and bounds fill-in.
    """

    m = [row[:] for row in matrix]
    rows = len(m)
    cols = len(m[0]) if rows else 0
    pivots: list[int] = []
    r = 0
    for c in range(cols):
        candidates = [i for i in range(r, rows) if not m[i][c].is_zero()]
        if not candidates:
            continue
        pivot = min(
            candidates,
            key=lambda i: sum(1 for entry in m[i] if not entry.is_zero()),
        )
        m[r], m[pivot] = m[pivot], m[r]
        scale = m[r][c].inv()
        m[r] = [entry * scale for entry in m[r]]
        for i in range(rows):
            if i != r and not m[i][c].is_zero():
                factor = m[i][c]
                m[i] = [m[i][j] - factor * m[r][j] for j in range(cols)]
        pivots.append(c)
        r += 1
        if r == rows:
            break
    return m, pivots


def rank_f5(matrix: Sequence[RowF5]) -> int:
    if not matrix:
        return 0
    return len(rref_f5(matrix)[1])


def nullspace_f5(matrix: Sequence[RowF5]) -> list[RowF5]:
    """Exact right-nullspace basis of a nonempty matrix over Q(sqrt5)."""

    reduced, pivots = rref_f5(matrix)
    cols = len(matrix[0]) if matrix else 0
    free = [c for c in range(cols) if c not in pivots]
    basis: list[RowF5] = []
    for f in free:
        vec = [ZERO for _ in range(cols)]
        vec[f] = ONE
        for row_index, p in enumerate(pivots):
            vec[p] = -reduced[row_index][f]
        basis.append(vec)
    return basis


# ---------------------------------------------------------------------------
# Point model: vertices, adjacency, level-one midpoints
# ---------------------------------------------------------------------------


def icosahedron_vertices() -> list[Vec3]:
    """The twelve unnormalized icosahedron vertices: cyclic perms of (0,+-1,+-phi)."""

    verts: list[Vec3] = []
    for s1 in (ONE, -ONE):
        for s2 in (PHI, -PHI):
            verts.append((ZERO, s1, s2))
    for s1 in (ONE, -ONE):
        for s2 in (PHI, -PHI):
            verts.append((s1, s2, ZERO))
    for s1 in (ONE, -ONE):
        for s2 in (PHI, -PHI):
            verts.append((s2, ZERO, s1))
    require(len(set(verts)) == 12, "the vertex model does not list twelve distinct points")
    require(R_VERTEX == F5(2) + PHI, "the vertex radius constant is not 2+phi")
    for vertex in verts:
        require(
            dot(vertex, vertex) == R_VERTEX,
            "a vertex squared norm is not (5+sqrt5)/2",
        )
    return verts


def adjacent_pairs(verts: Sequence[Vec3]) -> list[tuple[int, int]]:
    """The thirty edges: index pairs with inner product exactly phi."""

    pairs = [
        (i, j)
        for i in range(len(verts))
        for j in range(i + 1, len(verts))
        if dot(verts[i], verts[j]) == PHI
    ]
    require(len(pairs) == 30, f"expected 30 edges, found {len(pairs)}")
    for index in range(len(verts)):
        degree = sum(1 for a, b in pairs if index in (a, b))
        require(degree == 5, "a vertex does not have exactly five neighbors")
    return pairs


def unit_midpoints(verts: Sequence[Vec3]) -> list[Vec3]:
    """The thirty level-one unit points: midpoint sums projected to the sphere.

    Every midpoint sum m = u + v over an adjacent pair has squared norm
    6+2*sqrt5 = (1+sqrt5)^2, so the projection m/(1+sqrt5) stays inside
    Q(sqrt5) and has squared norm exactly 1.
    """

    require(
        ONE_PLUS_SQRT5 * ONE_PLUS_SQRT5 == R_MIDPOINT,
        "(1+sqrt5)^2 is not 6+2*sqrt5",
    )
    points: list[Vec3] = []
    for i, j in adjacent_pairs(verts):
        m = (
            verts[i][0] + verts[j][0],
            verts[i][1] + verts[j][1],
            verts[i][2] + verts[j][2],
        )
        require(dot(m, m) == R_MIDPOINT, "a midpoint sum squared norm is not 6+2*sqrt5")
        unit = (
            m[0] / ONE_PLUS_SQRT5,
            m[1] / ONE_PLUS_SQRT5,
            m[2] / ONE_PLUS_SQRT5,
        )
        require(dot(unit, unit) == ONE, "a projected midpoint squared norm is not 1")
        points.append(unit)
    require(len(set(points)) == 30, "the projected midpoints are not thirty distinct points")
    return points


def verify_vertex_radius_nonsquare(radius: F5) -> dict[str, Any]:
    """The squared vertex norm is not a square in Q(sqrt5).

    The Galois norm of the radius is the product with its conjugate. A
    square x*x in Q(sqrt5) has Galois norm N(x)^2, the square of a
    rational. The computed norm here is the integer 5; a rational whose
    square is an integer is itself an integer (monic polynomial x^2 - 5),
    and the integer square-root check excludes an integer root. Therefore
    no element of Q(sqrt5) squares to the radius, and 1 and the vertex
    radius are linearly independent over Q(sqrt5): the two parity parts
    of a vertex evaluation condition vanish separately.
    """

    product = radius * radius.conj()
    require(product.b == 0, "the Galois norm left the rational subfield")
    galois_norm = product.a
    require(
        galois_norm == 5,
        "the vertex radius Galois norm is not 5; the two-row vertex split premise fails",
    )
    n = int(galois_norm)
    root = isqrt(n)
    require(
        root * root != n,
        "the Galois norm is an integer square; the two-row vertex split premise fails",
    )
    return {
        "galois_norm": str(galois_norm),
        "galois_norm_equals_five": True,
        "integer_square_root_check": f"isqrt({n})^2 = {root * root} != {n}",
        "square_in_q_sqrt5": False,
        "consequence": (
            "1 and the vertex radius are linearly independent over Q(sqrt5); "
            "each vertex sphere-evaluation condition splits into two "
            "independent field-linear rows by parity of the homogeneity gap"
        ),
    }


# ---------------------------------------------------------------------------
# Exact rational harmonic bases (Laplacian nullspaces)
# ---------------------------------------------------------------------------


def monomials(degree: int) -> list[tuple[int, int, int]]:
    """The homogeneous degree-L monomial exponents in a fixed order."""

    monos = [
        (i, j, degree - i - j)
        for i in range(degree, -1, -1)
        for j in range(degree - i, -1, -1)
    ]
    require(
        len(monos) == (degree + 1) * (degree + 2) // 2,
        "the monomial count is not (L+1)(L+2)/2",
    )
    return monos


def apply_laplacian(coeffs: Sequence[Fraction], degree: int) -> list[Fraction]:
    """Exact Laplacian coefficients over the degree-(L-2) monomials.

    The Laplacian of a homogeneous polynomial of degree below two is zero
    and is returned as the empty coefficient list.
    """

    if degree < 2:
        return []
    source = monomials(degree)
    target = monomials(degree - 2)
    index = {mono: t for t, mono in enumerate(target)}
    out = [Fraction(0)] * len(target)
    for coeff, (i, j, k) in zip(coeffs, source):
        if coeff == 0:
            continue
        if i >= 2:
            out[index[(i - 2, j, k)]] += coeff * i * (i - 1)
        if j >= 2:
            out[index[(i, j - 2, k)]] += coeff * j * (j - 1)
        if k >= 2:
            out[index[(i, j, k - 2)]] += coeff * k * (k - 1)
    return out


_HARMONIC_CACHE: dict[int, tuple[list[tuple[int, int, int]], list[list[Fraction]]]] = {}


def harmonic_basis(degree: int) -> tuple[list[tuple[int, int, int]], list[list[Fraction]]]:
    """An exact rational basis of the degree-L harmonic homogeneous space.

    The basis is the exact nullspace of the Laplacian as a linear map from
    the degree-L monomial space to the degree-(L-2) monomial space,
    computed by Gaussian elimination over exact rationals. Three checks
    run on every construction: the dimension equals 2L+1, the Laplacian of
    every basis vector is termwise zero, and the coefficient matrix has
    full rank 2L+1 over Q.
    """

    cached = _HARMONIC_CACHE.get(degree)
    if cached is not None:
        return cached
    monos = monomials(degree)
    if degree < 2:
        basis = [
            [Fraction(1) if t == s else Fraction(0) for t in range(len(monos))]
            for s in range(len(monos))
        ]
    else:
        target = monomials(degree - 2)
        matrix: list[RowF5] = [[ZERO] * len(monos) for _ in target]
        for s_index in range(len(monos)):
            unit = [Fraction(0)] * len(monos)
            unit[s_index] = Fraction(1)
            image = apply_laplacian(unit, degree)
            for t_index, value in enumerate(image):
                if value != 0:
                    matrix[t_index][s_index] = F5(value)
        basis = []
        for vec in nullspace_f5(matrix):
            row: list[Fraction] = []
            for entry in vec:
                require(
                    entry.b == 0,
                    "a harmonic basis coefficient left the rational subfield",
                )
                row.append(entry.a)
            basis.append(row)
    require(
        len(basis) == 2 * degree + 1,
        f"the harmonic dimension at degree {degree} is not {2 * degree + 1}",
    )
    for vec in basis:
        require(
            all(coeff == 0 for coeff in apply_laplacian(vec, degree)),
            f"a degree-{degree} basis vector is not harmonic",
        )
    require(
        rank_f5([[F5(coeff) for coeff in vec] for vec in basis]) == 2 * degree + 1,
        f"the degree-{degree} harmonic basis is not full rank over Q",
    )
    _HARMONIC_CACHE[degree] = (monos, basis)
    return monos, basis


def space_dimension(degrees: Sequence[int]) -> int:
    return sum(2 * degree + 1 for degree in degrees)


# ---------------------------------------------------------------------------
# Evaluation rows
# ---------------------------------------------------------------------------


def coordinate_powers(point: Vec3, max_degree: int) -> list[list[F5]]:
    table: list[list[F5]] = []
    for coordinate in point:
        row = [ONE]
        for _ in range(max_degree):
            row.append(row[-1] * coordinate)
        table.append(row)
    return table


def evaluate_basis(degree: int, powers: Sequence[Sequence[F5]]) -> list[F5]:
    """Exact values of the degree-L harmonic basis at a point given its powers."""

    monos, basis = harmonic_basis(degree)
    values: list[F5] = []
    for vec in basis:
        total = ZERO
        for coeff, (i, j, k) in zip(vec, monos):
            if coeff == 0:
                continue
            total = total + F5(coeff) * (powers[0][i] * powers[1][j] * powers[2][k])
        values.append(total)
    return values


def midpoint_rows(degrees: Sequence[int], unit_points: Sequence[Vec3]) -> list[RowF5]:
    """One Q(sqrt5)-linear evaluation row per unit point over the degree set.

    The coordinates lie in Q(sqrt5) with squared norm one, so the
    sphere-evaluation condition f(p) = 0 is a single field-linear row: the
    exact value of every basis polynomial at the point. Every point is
    checked to be a unit point before its row is emitted.
    """

    max_degree = max(degrees)
    rows: list[RowF5] = []
    for point in unit_points:
        require(
            dot(point, point) == ONE,
            "a level-one evaluation point is not a unit point",
        )
        powers = coordinate_powers(point, max_degree)
        row: RowF5 = []
        for degree in degrees:
            row.extend(evaluate_basis(degree, powers))
        rows.append(row)
    return rows


def vertex_split_rows(degrees: Sequence[int], verts: Sequence[Vec3]) -> list[RowF5]:
    """Two Q(sqrt5)-linear rows per vertex, split by homogeneity-gap parity.

    The condition f(p/|p|) = 0 multiplied by |p|^Lmax is
    sum_L |p|^(Lmax-L) f_L(p) = 0. Terms with even gap Lmax-L carry a
    power of the squared norm and stay in Q(sqrt5); terms with odd gap
    carry one factor |p| besides. Because the squared vertex norm is not a
    square in Q(sqrt5), the two parts vanish separately, giving one row
    per parity class per vertex. A row is identically zero when the degree
    set has no member in its parity class.
    """

    verify_vertex_radius_nonsquare(R_VERTEX)
    max_degree = max(degrees)
    rows: list[RowF5] = []
    for point in verts:
        require(
            dot(point, point) == R_VERTEX,
            "a vertex does not carry squared norm (5+sqrt5)/2",
        )
        powers = coordinate_powers(point, max_degree)
        even_row: RowF5 = []
        odd_row: RowF5 = []
        for degree in degrees:
            gap = max_degree - degree
            values = evaluate_basis(degree, powers)
            if gap % 2 == 0:
                factor = f5_pow(R_VERTEX, gap // 2)
                even_row.extend(factor * value for value in values)
                odd_row.extend(ZERO for _ in values)
            else:
                factor = f5_pow(R_VERTEX, (gap - 1) // 2)
                odd_row.extend(factor * value for value in values)
                even_row.extend(ZERO for _ in values)
        rows.append(even_row)
        rows.append(odd_row)
    return rows


# ---------------------------------------------------------------------------
# Kernel analysis
# ---------------------------------------------------------------------------


def condition_analysis(
    degrees: Sequence[int],
    verts: Sequence[Vec3],
    unit_points: Sequence[Vec3],
    *,
    vertex_conditions: bool,
    midpoint_conditions: bool,
) -> dict[str, Any]:
    """Exact rank and kernel dimension of a degree set under condition rows."""

    rows: list[RowF5] = []
    record: dict[str, Any] = {
        "space_degrees": list(degrees),
        "space_dimension": space_dimension(degrees),
        "vertex_conditions": vertex_conditions,
        "midpoint_conditions": midpoint_conditions,
    }
    if vertex_conditions:
        vrows = vertex_split_rows(degrees, verts)
        record["identically_zero_vertex_rows"] = sum(
            1 for row in vrows if all(entry.is_zero() for entry in row)
        )
        rows.extend(vrows)
    if midpoint_conditions:
        rows.extend(midpoint_rows(degrees, unit_points))
    dimension = record["space_dimension"]
    for row in rows:
        require(len(row) == dimension, "a condition row width does not match the space")
    rank = rank_f5(rows)
    record["condition_rows"] = len(rows)
    record["matrix_rank"] = rank
    record["kernel_dimension"] = dimension - rank
    return record


# ---------------------------------------------------------------------------
# Receipt
# ---------------------------------------------------------------------------


def build_receipt() -> dict[str, Any]:
    verts = icosahedron_vertices()
    pairs = adjacent_pairs(verts)
    units = unit_midpoints(verts)
    nonsquare = verify_vertex_radius_nonsquare(R_VERTEX)

    vertex_only = condition_analysis(
        CANONICAL_DEGREES, verts, units, vertex_conditions=True, midpoint_conditions=False
    )
    vertex_only_even = condition_analysis(
        CANONICAL_EVEN_DEGREES, verts, units, vertex_conditions=True, midpoint_conditions=False
    )
    vertex_only_odd = condition_analysis(
        CANONICAL_ODD_DEGREES, verts, units, vertex_conditions=True, midpoint_conditions=False
    )
    level_one = condition_analysis(
        CANONICAL_DEGREES, verts, units, vertex_conditions=True, midpoint_conditions=True
    )
    comb_level_one = condition_analysis(
        COMB_DEGREES, verts, units, vertex_conditions=True, midpoint_conditions=True
    )
    midpoints_only = condition_analysis(
        CANONICAL_DEGREES, verts, units, vertex_conditions=False, midpoint_conditions=True
    )

    # Parity split sanity: single-parity degree sets leave one identically
    # zero row per vertex; the mixed degree<=3 set leaves none.
    require(
        vertex_only["identically_zero_vertex_rows"] == 0,
        "the mixed degree<=3 set does not fill both parity rows at every vertex",
    )
    require(
        vertex_only_even["identically_zero_vertex_rows"] == 12,
        "the even subspace does not collapse to one nontrivial row per vertex",
    )
    require(
        vertex_only_odd["identically_zero_vertex_rows"] == 12,
        "the odd subspace does not collapse to one nontrivial row per vertex",
    )
    require(
        comb_level_one["identically_zero_vertex_rows"] == 12,
        "the all-even comb set does not collapse to one nontrivial vertex row",
    )
    require(
        vertex_only["kernel_dimension"]
        == vertex_only_even["kernel_dimension"] + vertex_only_odd["kernel_dimension"],
        "the vertex-condition kernel does not decompose along the parity split",
    )

    # The frozen status and boundary strings state exactly these computed
    # facts; the build refuses to emit them against a contrary computation.
    require(
        vertex_only["kernel_dimension"] > 0,
        "the twelve-vertex kernel is zero; the static nonidentifiability "
        "status string does not apply",
    )
    require(
        level_one["kernel_dimension"] == 0,
        "the level-one kernel of the degree<=3 band is not zero; the "
        "level-one identifiability status string does not apply",
    )
    require(
        comb_level_one["kernel_dimension"]
        >= comb_level_one["space_dimension"] - comb_level_one["condition_rows"],
        "the comb kernel dimension is below the row-count bound",
    )
    require(
        comb_level_one["kernel_dimension"] > 0,
        "the comb band kernel is zero at level one; the recorded frontier "
        "statement does not apply",
    )

    receipt: dict[str, Any] = {
        "schema": SCHEMA,
        "issue": ISSUE,
        "status": STATUS,
        "arithmetic": {
            "field": "Q(sqrt5)",
            "representation": "a+b*sqrt5 with exact Fraction coefficients",
            "floating_point_used": False,
        },
        "points": {
            "icosahedron_vertices": {
                "count": len(verts),
                "construction": "cyclic permutations of (0,+-1,+-phi), phi = (1+sqrt5)/2",
                "squared_norm": R_VERTEX.text(),
                "coordinates": [[c.text() for c in v] for v in verts],
            },
            "adjacency": {
                "rule": "u adjacent to v iff u.v = phi = " + PHI.text(),
                "neighbors_per_vertex": 5,
                "edge_count": len(pairs),
            },
            "level_one_midpoints": {
                "count": len(units),
                "construction": "m = u + v over the thirty adjacent vertex pairs",
                "squared_norm": R_MIDPOINT.text(),
                "radius_identity": "(1+1*sqrt5)^2 = 6+2*sqrt5",
                "unit_projection_divisor": ONE_PLUS_SQRT5.text(),
                "unit_squared_norm": "1+0*sqrt5",
                "unit_coordinates": [[c.text() for c in p] for p in units],
            },
        },
        "vertex_radius_nonsquare": nonsquare,
        "evaluation_model": {
            "sphere_condition": (
                "f(p/|p|) = 0 multiplied by |p|^Lmax into "
                "sum_L |p|^(Lmax-L) f_L(p) = 0"
            ),
            "midpoint_rows": (
                "unit coordinates lie in Q(sqrt5); one field-linear row per point"
            ),
            "vertex_rows": (
                "two field-linear rows per vertex, split by parity of Lmax-L; "
                "the split rests on the vertex radius nonsquareness fact"
            ),
        },
        "harmonic_spaces": {
            str(degree): {
                "dimension": len(harmonic_basis(degree)[1]),
                "construction": (
                    "exact rational nullspace of the Laplacian on the "
                    f"homogeneous degree-{degree} monomial space"
                ),
                "laplacian_zero_verified": True,
                "basis_rank_verified": True,
            }
            for degree in ALL_DEGREES
        },
        "kernels": {
            "canonical_band_vertex_only": vertex_only,
            "canonical_band_even_part_vertex_only": vertex_only_even,
            "canonical_band_odd_part_vertex_only": vertex_only_odd,
            "canonical_band_level_one": level_one,
            "comb_band_level_one": comb_level_one,
            "canonical_band_midpoints_only": midpoints_only,
        },
        "parity_split_sanity": {
            "canonical_band_zero_vertex_rows": vertex_only["identically_zero_vertex_rows"],
            "canonical_even_part_zero_vertex_rows": vertex_only_even[
                "identically_zero_vertex_rows"
            ],
            "canonical_odd_part_zero_vertex_rows": vertex_only_odd[
                "identically_zero_vertex_rows"
            ],
            "comb_band_zero_vertex_rows": comb_level_one["identically_zero_vertex_rows"],
            "statement": (
                "a single-parity degree set makes one of the two rows per "
                "vertex identically zero; the mixed degree<=3 set makes both "
                "rows nontrivial at every vertex"
            ),
        },
        "refinement_readout_premise": (
            "the level-one identifiability statement consumes exact function "
            "values on the refined port set (twelve vertices and thirty "
            "projected edge midpoints); whether the physical readout exposes "
            "refined-port values is an open source premise owned by the "
            "repair/refinement law (oph-physics-sim docs/CANONICAL_REPAIR_LAW.md) "
            "and is not decided here"
        ),
        "not_claimed": [
            "no sky transfer is selected",
            "no repair-semigroup intertwining theorem is proved",
            (
                "the comb support band span(H_0,H_6,H_10,H_12) stays "
                "unidentifiable at level one with kernel dimension "
                f"{comb_level_one['kernel_dimension']}"
            ),
            "no comparison or public measurement data is read",
        ],
        "counterfamily_typing": (
            "every degree<=3 completion vanishing on the full level-one point "
            "set is identically zero (kernel dimension 0), so every surviving "
            "static counterfunction lies outside the canonical degree<=3 band "
            "or is separated by level-one data"
        ),
        "comparison_boundary": {
            "public_measurement_read": False,
            "comparison_permitted": False,
        },
    }
    receipt["receipt_sha256"] = tagged_sha256(canonical_json_bytes(receipt))
    return receipt


def write_runtime() -> Path:
    RUNTIME.mkdir(exist_ok=True)
    RECEIPT_PATH.write_bytes(canonical_json_bytes(build_receipt()))
    return RECEIPT_PATH


def verify_runtime() -> None:
    if RECEIPT_PATH.read_bytes() != canonical_json_bytes(build_receipt()):
        raise SystemExit("refinement transfer receipt is stale")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.write:
        print(write_runtime())
    if args.verify:
        verify_runtime()
        print("REFINEMENT_TRANSFER_VALID")
    if not args.write and not args.verify:
        receipt = build_receipt()
        kernels = receipt["kernels"]
        print(
            json.dumps(
                {
                    "status": receipt["status"],
                    "kernel_dimensions": {
                        name: record["kernel_dimension"]
                        for name, record in kernels.items()
                    },
                    "receipt_sha256": receipt["receipt_sha256"],
                },
                indent=2,
                sort_keys=True,
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
