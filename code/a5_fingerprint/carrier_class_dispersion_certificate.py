"""Carrier-class dispersion band certificate.

Scope. The declared positive-weight scalar cosine class consists of the
symbols

    Lambda(k, n) = sum_s w_s sum_{u in O_s} [1 - cos(a k r_s (u . n))]

where a is finite and strictly positive, each O_s is one A5 orbit of unit
directions, each shell radius r_s is positive, each active shell weight w_s
is strictly positive, and the continuum
normalization fixes the k^2 coefficient to one. Invariance forces equal
weights inside each orbit, and the cosine form is even, so reciprocity is
automatic. Writing the normalized spatial-symbol expansion as

    Lambda_norm = k^2 + C4 k^4 + (B0 + B6 I6(n)) k^6
                  + (D0 + D6 I6(n)) k^8 + O(k^10),

this certificate proves the class theorem and independently reproduces both
frozen branch rays.

1. Kernel factorization. The group-summed sixth-power kernel of the sixty
   proper rotations factors on the invariant line: for every direction u,

       sum_{g in G} ((g u) . n)^6 = (60/7) + (64/35) I6(u) I6(n)

   on the unit sphere in both arguments. Consequently the level-six content
   of every orbit sum is one multiple of the same I6, with coefficient
   proportional to the I6 value of the orbit seed. The proportionality
   constant is universal across all orbits, including generic sixty-point
   orbits.

2. Sign law and isotropic floor. Every member has C4 < 0, and

       B0 / C4^2 = (10/21) * (mu2 mu6 / mu4^2) >= 10/21,

   where mu_m is the weighted radial moment sum_s w_s |O_s| r_s^m. The
   floor is the Lagrange identity: writing W_s = w_s |O_s|,
   mu2 mu6 - mu4^2 is a positive
   combination of (r_i^2 - r_j^2)^2 terms, so equality holds exactly on
   the single-radius members. One-sign nonnegativity is load-bearing: a
   signed-weight control member violates the floor, while active shells are
   taken strictly positive.

3. Rank-six band. Every member has

       B6 / B0 = (16/75) * <I6(u_seed)>  in  [-16/135, 16/75],

   the mean taken with the positive weights w |O| r^6; on single-radius
   members the weights reduce to w |O|. The endpoints are attained exactly
   by the pure face orbit and the pure vertex orbit; the pure edge orbit
   sits at -1/15; the interior zero is attained by the vertex-face mixture
   with per-direction vertex:face weight ratio 25:27. The interval bound consumes the 62-direction
   stationary census of the fixed-point packet, which places the range of
   I6 at [-5/9, 1].

4. Eighth-order confinement and the cross-order lock. The degree-eight
   kernel factors as

       sum_{g in G} ((g u).n)^8 = (20/3) + (256/75) I6(u) I6(n)

   on the unit sphere, and the harmonic angular-rank-eight multiplicity is zero,
   so the k^8 anisotropy of every member is one multiple of the same
   rotated I6 with no new shape. Writing the k^8 term as
   (D0 + D6 I6) k^8, every pure unit-radius orbit has
   D0 = -a^6/60480 and D6/D0 = (64/125) I6(seed), so the
   eighth-order band is [-64/225, 64/125]. Every single-radius member,
   including a member whose anisotropy vanishes, obeys the genuinely
   division-free cross-order identity

       5 D6 B0 = 12 B6 D0,

   equivalently D6/D0 = (12/5) (B6/B0), since the positive-weight class
   has nonzero B0 and D0. The quotient of the two anisotropy ratios is
   undefined at the zero-anisotropy member, while the polynomial identity
   remains meaningful there.

   At common radius r, D0 = -(a r)^6/60480. A generic rank-six anisotropic model
   carries independent k^6 and k^8 amplitudes; the single-radius stratum
   fixes their relative amplitude. Multi-radius members retain correlated
   radial-moment dependence. The isotropic tower alternates in sign at every order,
   because the even moments are positive and the cosine series
   alternates.

5. Frozen-branch reproduction. The map evaluated at the vertex orbit
   reproduces the complete FZ-11 coefficient ray, and evaluated at the edge
   orbit reproduces the complete FZ-12 coefficient ray, by independent
   recomputation. The face-orbit branch completes the three-orbit table
   with C4 = -a^2/20, B0 = a^4/840, B6 = -2 a^4/14175.

The class statements sharpen the falsification surface. Under the registered
physical-symbol, sector, frame, scale, readout, and exclusivity premises, a
resolved dispersion measurement with B0/C4^2 below 10/21, a rank-six-to-
isotropic ratio outside [-16/135, 16/75], an intrinsic residue at ranks one
through five, or a rank-six residue off the rotated I6 line excludes every
member of the declared positive-weight scalar cosine class at once rather
than one frozen branch.
The physical premises are not constructed here. Their target-clean producer,
physical bridge, quantitative output, and pre-comparison registration belong
to the relevant propagation and prediction lanes.

Run with --write to refresh the committed receipt.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

import a5_multipole_fixed_point_certificate as base
import spin_six_universality_certificate as universality

q5 = base.q5
Q5 = base.Q5
ZERO = base.ZERO
ONE = base.ONE
PHI = base.PHI
q5_add = base.q5_add
q5_sub = base.q5_sub
q5_mul = base.q5_mul
q5_scale = base.q5_scale
q5_div = base.q5_div
q5_pow = base.q5_pow
q5_neg = base.q5_neg
q5_str = base.q5_str
q5_sign = base.q5_sign
require = base.require
p_zero = base.p_zero
p_const = base.p_const
p_add = base.p_add
p_scale = base.p_scale
p_scale_q5 = base.p_scale_q5
p_pow = base.p_pow
p_mul = base.p_mul
p_linear = base.p_linear
p_eval = base.p_eval
p_reduce_sphere = base.p_reduce_sphere
p_is_zero = base.p_is_zero
moment_sum = base.moment_sum

SCHEMA = "oph.carrier_class_dispersion.v1"
STATUS = (
    "CARRIER_CLASS_DISPERSION_BAND_CERTIFIED__"
    "ISOTROPIC_FLOOR_10_21_WITH_SINGLE_RADIUS_SATURATION__"
    "RANK_SIX_BAND_MINUS_16_135_TO_16_75__"
    "EIGHTH_ORDER_I6_ONLY_WITH_CROSS_ORDER_LOCK_12_5__"
    "FROZEN_VERTEX_AND_EDGE_BRANCHES_REPRODUCED__"
    "PHYSICAL_BRIDGE_PREMISES_OPEN"
)
RECEIPT_PATH = Path(__file__).resolve().parent / "runtime" / (
    "carrier_class_dispersion_receipt.json"
)

HALF = Fraction(1, 2)
# Entries of the golden rotation candidates: {1/(2 phi), 1/2, phi/2} with
# 1/(2 phi) = (sqrt5 - 1)/4 and phi/2 = (1 + sqrt5)/4.
INV_TWO_PHI = q5(Fraction(-1, 4), Fraction(1, 4))
PHI_HALF = q5(Fraction(1, 4), Fraction(1, 4))
Q5_HALF = q5(HALF)


# ---------------------------------------------------------------------------
# Exact matrix helpers over Q(sqrt5)
# ---------------------------------------------------------------------------

Mat = tuple[tuple[Q5, Q5, Q5], ...]


def mat_mul(a: Mat, b: Mat) -> Mat:
    return tuple(
        tuple(
            q5_add(
                q5_add(q5_mul(a[i][0], b[0][j]), q5_mul(a[i][1], b[1][j])),
                q5_mul(a[i][2], b[2][j]),
            )
            for j in range(3)
        )
        for i in range(3)
    )


def mat_apply(a: Mat, v: tuple[Q5, Q5, Q5]) -> tuple[Q5, Q5, Q5]:
    return tuple(
        q5_add(
            q5_add(q5_mul(a[i][0], v[0]), q5_mul(a[i][1], v[1])),
            q5_mul(a[i][2], v[2]),
        )
        for i in range(3)
    )


def mat_transpose(a: Mat) -> Mat:
    return tuple(tuple(a[j][i] for j in range(3)) for i in range(3))


def mat_det(a: Mat) -> Q5:
    total = ZERO
    for c0, c1, c2, sign in (
        (0, 1, 2, 1),
        (1, 2, 0, 1),
        (2, 0, 1, 1),
        (2, 1, 0, -1),
        (0, 2, 1, -1),
        (1, 0, 2, -1),
    ):
        term = q5_mul(a[0][c0], q5_mul(a[1][c1], a[2][c2]))
        total = q5_add(total, term if sign == 1 else q5_neg(term))
    return total


IDENTITY: Mat = tuple(
    tuple(ONE if i == j else ZERO for j in range(3)) for i in range(3)
)


def mat_is_special_orthogonal(a: Mat) -> bool:
    prod = mat_mul(a, mat_transpose(a))
    return prod == IDENTITY and mat_det(a) == ONE


# ---------------------------------------------------------------------------
# The sixty proper rotations as exact matrices
# ---------------------------------------------------------------------------


def vertex_set(verts) -> frozenset:
    return frozenset(verts)


def permutes_vertices(a: Mat, verts, vset) -> bool:
    return all(mat_apply(a, v) in vset for v in verts)


def rotation_group_certificate() -> tuple[list[Mat], dict[str, Any]]:
    """Generate and verify the sixty proper rotations over Q(sqrt5)."""

    verts = base.cartesian_vertices()
    vset = vertex_set(verts)

    cyclic: Mat = (
        (ZERO, ZERO, ONE),
        (ONE, ZERO, ZERO),
        (ZERO, ONE, ZERO),
    )
    require(
        mat_is_special_orthogonal(cyclic)
        and permutes_vertices(cyclic, verts, vset),
        "cyclic coordinate rotation rejected",
    )

    magnitudes = (INV_TWO_PHI, Q5_HALF, PHI_HALF)
    import itertools

    def cross(u, v):
        return (
            q5_sub(q5_mul(u[1], v[2]), q5_mul(u[2], v[1])),
            q5_sub(q5_mul(u[2], v[0]), q5_mul(u[0], v[2])),
            q5_sub(q5_mul(u[0], v[1]), q5_mul(u[1], v[0])),
        )

    def dot(u, v):
        return q5_add(
            q5_add(q5_mul(u[0], v[0]), q5_mul(u[1], v[1])),
            q5_mul(u[2], v[2]),
        )

    row_candidates = []
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product((1, -1), repeat=3):
            row = tuple(
                magnitudes[perm[j]]
                if signs[j] == 1
                else q5_neg(magnitudes[perm[j]])
                for j in range(3)
            )
            if dot(row, row) == ONE:
                row_candidates.append(row)

    seed = None
    for row_a in row_candidates:
        if seed is not None:
            break
        for row_b in row_candidates:
            if dot(row_a, row_b) != ZERO:
                continue
            candidate: Mat = (row_a, row_b, cross(row_a, row_b))
            if mat_is_special_orthogonal(candidate) and permutes_vertices(
                candidate, verts, vset
            ):
                seed = candidate
                break
    require(seed is not None, "no golden seed rotation found")

    group: dict[Mat, None] = {IDENTITY: None, cyclic: None, seed: None}
    frontier = [cyclic, seed]
    while frontier:
        nxt = []
        for a in list(group):
            for b in frontier:
                c = mat_mul(a, b)
                if c not in group:
                    group[c] = None
                    nxt.append(c)
                d = mat_mul(b, a)
                if d not in group:
                    group[d] = None
                    nxt.append(d)
        frontier = nxt
        require(len(group) <= 60, "rotation closure exceeds sixty elements")
    rotations = list(group)
    require(len(rotations) == 60, "rotation closure is not sixty elements")
    for a in rotations:
        require(
            mat_is_special_orthogonal(a), "non-orthogonal closure element"
        )
        require(
            permutes_vertices(a, verts, vset),
            "closure element leaves the vertex set",
        )
    order_two = sum(1 for a in rotations if mat_mul(a, a) == IDENTITY) - 1
    require(order_two == 15, "involution count is not fifteen")
    return rotations, {
        "statement": (
            "sixty exact special-orthogonal matrices over Q(sqrt5) close "
            "under product, permute the twelve vertex directions, and "
            "contain exactly fifteen involutions"
        ),
        "order": 60,
        "involutions": 15,
    }


# ---------------------------------------------------------------------------
# I6 evaluation at arbitrary exact directions
# ---------------------------------------------------------------------------


def i6_at_unit(i6_raw, v: tuple[Q5, Q5, Q5]) -> Q5:
    """Evaluate I6 at the unit direction of v, exactly.

    The raw invariant polynomial has even-degree parts only, so the
    normalized value needs only integer powers of |v|^2.
    """

    norm_sq = ZERO
    for axis in range(3):
        norm_sq = q5_add(norm_sq, q5_mul(v[axis], v[axis]))
    require(q5_sign(norm_sq) > 0, "zero direction in I6 evaluation")
    by_degree: dict[int, Q5] = {}
    for (i, j, k), coeff in i6_raw.items():
        deg = i + j + k
        require(deg % 2 == 0, "odd-degree monomial in the raw invariant")
        mono = q5_mul(
            coeff,
            q5_mul(
                q5_pow(v[0], i), q5_mul(q5_pow(v[1], j), q5_pow(v[2], k))
            ),
        )
        by_degree[deg] = q5_add(by_degree.get(deg, ZERO), mono)
    total = ZERO
    for deg, value in by_degree.items():
        total = q5_add(total, q5_div(value, q5_pow(norm_sq, deg // 2)))
    return total


# ---------------------------------------------------------------------------
# Kernel factorization on the invariant line
# ---------------------------------------------------------------------------


def kernel_factorization_certificate(
    rotations: list[Mat], i6_raw, i6_reduced
) -> tuple[Q5, dict[str, Any]]:
    """Prove sum_g ((g u).n)^6 = 60/7 + d I6(u) I6(n) with d = 64/35.

    The group-summed kernel is invariant in each argument separately, and
    the invariant table pins the even-degree-six invariant space on the
    sphere to span{1, I6}. Membership is verified constructively for each
    test seed: the reduced kernel minus its isotropic part is one exact
    multiple of the reduced I6, and the multiple is proportional to the
    I6 value of the seed with one universal constant.
    """

    probe = next(m for m in sorted(i6_reduced) if sum(m) > 0)
    verts = base.cartesian_vertices()
    orbits = universality.orbit_directions()

    test_seeds: list[tuple[str, tuple[Q5, Q5, Q5]]] = [
        ("vertex", verts[0]),
        ("edge", orbits["edge_30"]["dirs"][0]),
        ("face", orbits["face_20"]["dirs"][0]),
    ]
    for idx, raw in enumerate(
        ((1, 2, 3), (2, -1, 5), (1, 1, 7), (3, 5, -2), (0, 2, 9))
    ):
        test_seeds.append(
            (f"generic_{idx}", tuple(q5(c) for c in raw))
        )

    d_constant: Q5 | None = None
    rows = []
    for name, seed in test_seeds:
        norm_sq = ZERO
        for axis in range(3):
            norm_sq = q5_add(norm_sq, q5_mul(seed[axis], seed[axis]))
        for power, iso_constant in ((2, Fraction(20)), (4, Fraction(12))):
            low = p_zero()
            for g in rotations:
                low = p_add(
                    low, p_pow(p_linear(mat_apply(g, seed)), power)
                )
            low = p_scale_q5(
                low, q5_div(ONE, q5_pow(norm_sq, power // 2))
            )
            require(
                p_reduce_sphere(low) == {(0, 0, 0): q5(iso_constant)},
                f"{name}: degree-{power} group sum is not isotropic "
                f"{iso_constant}",
            )
        kernel = p_zero()
        for g in rotations:
            kernel = p_add(kernel, p_pow(p_linear(mat_apply(g, seed)), 6))
        kernel = p_scale_q5(kernel, q5_div(ONE, q5_pow(norm_sq, 3)))
        reduced = p_reduce_sphere(kernel)
        multiple = q5_div(reduced.get(probe, ZERO), i6_reduced[probe])
        residue = p_add(reduced, p_scale_q5(i6_reduced, q5_neg(multiple)))
        iso = residue.pop((0, 0, 0), ZERO)
        require(
            iso == q5(Fraction(60, 7)),
            f"{name}: isotropic kernel part is not 60/7",
        )
        require(
            p_is_zero(residue),
            f"{name}: kernel leaves the invariant line",
        )
        seed_value = i6_at_unit(i6_raw, seed)
        if seed_value != ZERO:
            ratio = q5_div(multiple, seed_value)
            if d_constant is None:
                d_constant = ratio
            require(
                ratio == d_constant,
                f"{name}: kernel constant is not universal",
            )
        else:
            require(
                multiple == ZERO,
                f"{name}: zero seed value with nonzero multiple",
            )
        rows.append(
            {
                "seed": name,
                "i6_at_seed": q5_str(seed_value),
                "i6_multiple": q5_str(multiple),
            }
        )
    require(d_constant is not None, "kernel constant undetermined")
    require(
        d_constant == q5(Fraction(64, 35)),
        "kernel constant is not 64/35",
    )

    # Completeness of the identity for every seed direction. The kernel is
    # invariant in the seed under right translation, so as a degree-six
    # form of the seed it lies in the invariant space of that degree. The
    # recomputed invariant table pins the harmonic invariant multiplicities
    # at L = 0..6 to (1, 0, 0, 0, 0, 0, 1), so the degree-six invariant
    # space is the two-dimensional span of |u|^6 and the harmonic I6 form.
    # Residue-zero verification at two seeds with distinct I6 values
    # therefore determines the identity for all seeds; the remaining six
    # seeds are redundancy.
    invariants = universality.invariant_table_certificate()
    table = invariants["table"]
    require(
        [table[str(level)] for level in range(7)]
        == [1, 0, 0, 0, 0, 0, 1],
        "invariant multiplicities through level six drift",
    )
    distinct_values = {
        row["i6_at_seed"] for row in rows
    }
    require(
        len(distinct_values) >= 2,
        "seed set spans fewer than two I6 values",
    )
    return d_constant, {
        "statement": (
            "sum over the sixty rotations of ((g u).n)^6 equals 60/7 plus "
            "(64/35) I6(u) I6(n) on the unit sphere in both arguments, for "
            "every seed direction"
        ),
        "isotropic_part": "60/7",
        "universal_constant": q5_str(d_constant),
        "degree_two_and_four_group_sums": (
            "isotropic with constants 20 and 12 at every seed, closing the "
            "generic-orbit isotropy of the second and fourth moments"
        ),
        "completeness_argument": (
            "the kernel is right-translation invariant in the seed, the "
            "recomputed invariant multiplicities (1,0,0,0,0,0,1) through "
            "level six pin the degree-six invariant space to the span of "
            "the radial power and the I6 form, and residue-zero at seeds "
            "with distinct I6 values determines the identity on that "
            "two-dimensional space; the extra seeds are redundancy"
        ),
        "invariant_multiplicities_l0_to_l6": [1, 0, 0, 0, 0, 0, 1],
        "seeds": rows,
    }


# ---------------------------------------------------------------------------
# Per-orbit coefficient table and frozen-branch reproduction
# ---------------------------------------------------------------------------


def orbit_branch_certificate(i6_raw, i6_reduced, d_constant) -> dict[str, Any]:
    """The three pure single-orbit branches with exact coefficients."""

    probe = next(m for m in sorted(i6_reduced) if sum(m) > 0)
    orbits = universality.orbit_directions()
    order = (("vertex_12", 12), ("edge_30", 30), ("face_20", 20))
    expected = {
        "vertex_12": {
            "beta": Fraction(64, 175),
            "B6_over_a4": Fraction(2, 7875),
            "B6_over_B0": Fraction(16, 75),
            "B6_over_C4_sq": Fraction(32, 315),
            "i6_seed": Fraction(1),
        },
        "edge_30": {
            "beta": Fraction(-2, 7),
            "B6_over_a4": Fraction(-1, 12600),
            "B6_over_B0": Fraction(-1, 15),
            "B6_over_C4_sq": Fraction(-2, 63),
            "i6_seed": Fraction(-5, 16),
        },
        "face_20": {
            "beta": Fraction(-64, 189),
            "B6_over_a4": Fraction(-2, 14175),
            "B6_over_B0": Fraction(-16, 135),
            "B6_over_C4_sq": Fraction(-32, 567),
            "i6_seed": Fraction(-5, 9),
        },
    }
    rows = []
    for name, count in order:
        data = orbits[name]
        dirs, norm_sq = data["dirs"], data["norm_sq"]
        require(len(dirs) == count, f"{name} count drift")
        for k in (1, 3, 5):
            require(
                p_is_zero(moment_sum(dirs, k)),
                f"{name}: odd moment {k} survives",
            )
        m2 = p_reduce_sphere(
            p_scale_q5(moment_sum(dirs, 2), q5_div(ONE, norm_sq))
        )
        require(
            m2 == {(0, 0, 0): q5(Fraction(count, 3))},
            f"{name}: second moment is not count/3",
        )
        m4 = p_reduce_sphere(
            p_scale_q5(moment_sum(dirs, 4), q5_div(ONE, q5_pow(norm_sq, 2)))
        )
        require(
            m4 == {(0, 0, 0): q5(Fraction(count, 5))},
            f"{name}: fourth moment is not count/5",
        )
        m6 = p_reduce_sphere(
            p_scale_q5(moment_sum(dirs, 6), q5_div(ONE, q5_pow(norm_sq, 3)))
        )
        beta = q5_div(m6.get(probe, ZERO), i6_reduced[probe])
        residue = p_add(m6, p_scale_q5(i6_reduced, q5_neg(beta)))
        iso = residue.pop((0, 0, 0), ZERO)
        require(
            iso == q5(Fraction(count, 7)),
            f"{name}: sixth-moment isotropic part is not count/7",
        )
        require(p_is_zero(residue), f"{name}: sixth moment off the line")
        want = expected[name]
        require(
            beta == q5(want["beta"]),
            f"{name}: beta is not {want['beta']}",
        )
        seed_value = i6_at_unit(i6_raw, dirs[0])
        require(
            seed_value == q5(want["i6_seed"]),
            f"{name}: seed I6 value drift",
        )
        require(
            beta
            == q5_mul(
                q5(Fraction(count, 60)), q5_mul(d_constant, seed_value)
            ),
            f"{name}: beta breaks the kernel factorization",
        )
        # Normalized single-orbit branch: prefactor 6/(count a^2) fixes the
        # k^2 coefficient to one; then
        #   C4 = -a^2/20, B0 = a^4/840, B6 = a^4 beta / (120 count).
        require(
            q5_div(beta, q5(Fraction(120 * count)))
            == q5(want["B6_over_a4"]),
            f"{name}: B6 drift",
        )
        require(
            q5(want["B6_over_a4"])
            == q5_mul(q5(want["B6_over_B0"]), q5(Fraction(1, 840))),
            f"{name}: B6/B0 inconsistent",
        )
        require(
            q5(want["B6_over_a4"])
            == q5_mul(q5(want["B6_over_C4_sq"]), q5(Fraction(1, 400))),
            f"{name}: B6/C4^2 inconsistent",
        )
        require(
            q5(want["B6_over_B0"])
            == q5_mul(q5(Fraction(16, 75)), q5(want["i6_seed"])),
            f"{name}: B6/B0 is not (16/75) I6(seed)",
        )
        rows.append(
            {
                "orbit": name,
                "count": count,
                "normalization": f"6/({count} a^2)",
                "C4_over_a2": "-1/20",
                "B0_over_a4": "1/840",
                "B6_over_a4": str(want["B6_over_a4"]),
                "B6_over_B0": str(want["B6_over_B0"]),
                "B6_over_C4_squared": str(want["B6_over_C4_sq"]),
                "B0_over_C4_squared": "10/21",
                "i6_at_seed": str(want["i6_seed"]),
                "beta": str(want["beta"]),
            }
        )
    return {
        "statement": (
            "each pure single-orbit unit-radius branch has C4 = -a^2/20, "
            "B0 = a^4/840, and B6/B0 = (16/75) I6(seed); the vertex, "
            "edge, and face branches sit at 16/75, -1/15, and -16/135"
        ),
        "branches": rows,
    }


def series_route_control() -> dict[str, Any]:
    """Second route: raw cosine series coefficients for each pure orbit.

    The Taylor coefficients of 1 - cos are 1/2, -1/24, 1/720 at x^2, x^4,
    x^6. This route rebuilds the normalized k^4 and k^6 coefficients from
    those factorials and the raw moment polynomials without the moment
    ratio shortcut, and must land on the same table.
    """

    orbits = universality.orbit_directions()
    cartesian = base.build_cartesian_frame()
    i6_raw = cartesian.pop("_i6_poly_object")
    i6_reduced = p_reduce_sphere(i6_raw)
    probe = next(m for m in sorted(i6_reduced) if sum(m) > 0)
    checks = []
    for name, count in (("vertex_12", 12), ("edge_30", 30), ("face_20", 20)):
        data = orbits[name]
        dirs, norm_sq = data["dirs"], data["norm_sq"]
        pref = Fraction(6, count)
        k4 = p_reduce_sphere(
            p_scale(
                p_scale_q5(
                    moment_sum(dirs, 4), q5_div(ONE, q5_pow(norm_sq, 2))
                ),
                pref * Fraction(-1, 24),
            )
        )
        require(
            k4 == {(0, 0, 0): q5(Fraction(-1, 20))},
            f"{name}: series k^4 coefficient drift",
        )
        k6 = p_reduce_sphere(
            p_scale(
                p_scale_q5(
                    moment_sum(dirs, 6), q5_div(ONE, q5_pow(norm_sq, 3))
                ),
                pref * Fraction(1, 720),
            )
        )
        b6 = q5_div(k6.get(probe, ZERO), i6_reduced[probe])
        residue = p_add(k6, p_scale_q5(i6_reduced, q5_neg(b6)))
        iso = residue.pop((0, 0, 0), ZERO)
        require(
            iso == q5(Fraction(1, 840)),
            f"{name}: series isotropic k^6 coefficient drift",
        )
        require(p_is_zero(residue), f"{name}: series rank-six residue")
        checks.append({"orbit": name, "b6_value": b6})
    expected_b6 = (
        Fraction(2, 7875),
        Fraction(-1, 12600),
        Fraction(-2, 14175),
    )
    require(
        all(
            row["b6_value"] == q5(value)
            for row, value in zip(checks, expected_b6)
        ),
        "series-route B6 table drift",
    )
    for row, value in zip(checks, expected_b6):
        row["B6_over_a4"] = str(value)
        del row["b6_value"]
    return {
        "statement": (
            "the factorial series route reproduces C4 = -a^2/20, "
            "B0 = a^4/840, and the three-orbit B6 table independently of "
            "the moment-ratio route"
        ),
        "branches": checks,
    }


def _validate_frozen_branch_values(fz11: dict, fz12: dict) -> dict[str, Any]:
    """Validate typed coefficient fields of the two frozen branch payloads."""

    rel11 = fz11["exact_prediction"]["scale_free_relations"]
    require(
        rel11["B0_over_C4_squared"] == "10/21"
        and rel11["B6_over_B0"] == "16/75"
        and rel11["B6_over_C4_squared"] == "32/315",
        "FZ-11 frozen relations disagree with the class map",
    )
    coeffs11 = fz11["exact_prediction"]["coefficients"]
    require(
        coeffs11["C4_over_a2"] == "-1/20"
        and coeffs11["B0_over_a4"] == "1/840"
        and coeffs11["B6_over_a4"] == "2/7875",
        "FZ-11 frozen coefficients disagree with the class map",
    )
    rel12 = fz12["conditional_physical_candidate"]["scale_free_relations"]
    require(
        rel12["B0_over_C4_squared"] == "10/21"
        and rel12["B6_over_B0"] == "-1/15"
        and rel12["B6_over_C4_squared"] == "-2/63",
        "FZ-12 frozen relations disagree with the class map",
    )
    coeffs12 = fz12["conditional_physical_candidate"]["coefficients"]
    require(
        coeffs12["C4_over_a2"] == "-1/20"
        and coeffs12["B0_over_a4"] == "1/840"
        and coeffs12["B6_over_a4"] == "-1/12600",
        "FZ-12 frozen coefficients disagree with the class map",
    )
    return {
        "fz11_relations": dict(rel11),
        "fz11_coefficients": dict(coeffs11),
        "fz12_relations": dict(rel12),
        "fz12_coefficients": dict(coeffs12),
    }


def frozen_receipt_crosscheck() -> dict[str, Any]:
    """Read-only typed agreement and custody pins for both frozen branches."""

    here = Path(__file__).resolve().parent / "runtime"
    paths = {
        "FZ-11": here / "spin_six_primitive_port_prediction_receipt.json",
        "FZ-12": here / "seam_current_edge_prediction_receipt.json",
    }
    raw = {row: path.read_bytes() for row, path in paths.items()}
    payloads = {row: json.loads(data) for row, data in raw.items()}
    values = _validate_frozen_branch_values(payloads["FZ-11"], payloads["FZ-12"])

    register = json.loads(
        (Path(__file__).resolve().parents[2] / "claims" /
         "frozen_prediction_register.json").read_text()
    )
    register_by_id = {row["id"]: row for row in register["rows"]}
    for row_id in ("FZ-11", "FZ-12"):
        digest = hashlib.sha256(raw[row_id]).hexdigest()
        require(
            register_by_id[row_id]["content_sha256"] == digest,
            f"{row_id} frozen file does not match its registered custody hash",
        )

    values.update(
        {
            "statement": (
                "the frozen FZ-11 and FZ-12 receipts carry exactly the vertex "
                "and edge values of the class coefficient map; typed fields and "
                "registered byte hashes are checked without changing frozen bytes"
            ),
            "parent_pins": [
                {
                    "row": row_id,
                    "path": str(paths[row_id].relative_to(Path(__file__).parents[2])),
                    "bytes": len(raw[row_id]),
                    "sha256": base.tagged_sha256(raw[row_id]),
                }
                for row_id in ("FZ-11", "FZ-12")
            ],
        }
    )
    return values


# ---------------------------------------------------------------------------
# Radial floor: Lagrange identity and signed-weight control
# ---------------------------------------------------------------------------


def _mu(shells: list[tuple[Fraction, Fraction]], power: int) -> Fraction:
    # shells: (effective weight W = w |O|, r_squared) pairs;
    # mu_m uses r^m = (r^2)^(m/2).
    total = Fraction(0)
    for weight, r_sq in shells:
        total += weight * r_sq ** (power // 2)
    return total


def _lagrange_gap(shells: list[tuple[Fraction, Fraction]]) -> Fraction:
    gap = Fraction(0)
    for i in range(len(shells)):
        wi, ri = shells[i]
        for j in range(i + 1, len(shells)):
            wj, rj = shells[j]
            gap += wi * wj * ri * rj * (ri - rj) ** 2
    return gap


def radial_floor_certificate() -> dict[str, Any]:
    """mu2 mu6 - mu4^2 equals the positive Lagrange sum, so the isotropic
    ratio B0/C4^2 = (10/21)(mu2 mu6/mu4^2) has floor 10/21 exactly on
    single-radius members."""

    import itertools
    import random

    rng = random.Random(20260803)
    symbolic_checks = []
    for shell_count in (2, 3, 4, 5):
        for _ in range(40):
            shells = [
                (
                    Fraction(rng.randint(1, 60), rng.randint(1, 9)),
                    Fraction(rng.randint(1, 60), rng.randint(1, 9)),
                )
                for _ in range(shell_count)
            ]
            lhs = _mu(shells, 2) * _mu(shells, 6) - _mu(shells, 4) ** 2
            require(
                lhs == _lagrange_gap(shells),
                "Lagrange identity fails on a random positive instance",
            )
            require(lhs >= 0, "positive-weight gap is negative")
        symbolic_checks.append({"shells": shell_count, "instances": 40})

    # Exact symbolic two-shell identity over polynomial coefficients:
    # coefficients of (w1, w2, s, t) in mu2 mu6 - mu4^2 and in
    # w1 w2 s t (s - t)^2 agree monomial by monomial.
    # mu2 mu6 - mu4^2 = (w1 s + w2 t)(w1 s^3 + w2 t^3) - (w1 s^2 + w2 t^2)^2
    #                 = w1 w2 (s t^3 + s^3 t - 2 s^2 t^2)
    #                 = w1 w2 s t (s - t)^2.
    # The monomial bookkeeping is checked exactly:
    from collections import Counter

    lhs_monomials = Counter()
    # (w1 s + w2 t)(w1 s^3 + w2 t^3)
    for (wa, sa) in (((1, 0), (1, 0)), ((0, 1), (0, 1))):
        for (wb, sb) in (((1, 0), (3, 0)), ((0, 1), (0, 3))):
            key = (
                wa[0] + wb[0],
                wa[1] + wb[1],
                sa[0] + sb[0],
                sa[1] + sb[1],
            )
            lhs_monomials[key] += 1
    # minus (w1 s^2 + w2 t^2)^2
    for (wa, sa) in (((1, 0), (2, 0)), ((0, 1), (0, 2))):
        for (wb, sb) in (((1, 0), (2, 0)), ((0, 1), (0, 2))):
            key = (
                wa[0] + wb[0],
                wa[1] + wb[1],
                sa[0] + sb[0],
                sa[1] + sb[1],
            )
            lhs_monomials[key] -= 1
    rhs_monomials = Counter(
        {(1, 1, 3, 1): 1, (1, 1, 1, 3): 1, (1, 1, 2, 2): -2}
    )
    lhs_clean = {k: v for k, v in lhs_monomials.items() if v}
    require(
        lhs_clean == dict(rhs_monomials),
        "two-shell Lagrange identity fails symbolically",
    )

    # Equality case: gap zero forces every cross term zero, so all radii
    # with positive weight coincide.
    single = [(Fraction(3), Fraction(7, 2)), (Fraction(5), Fraction(7, 2))]
    require(_lagrange_gap(single) == 0, "single-radius gap is nonzero")
    ratio_single = Fraction(10, 21) * (
        _mu(single, 2) * _mu(single, 6) / _mu(single, 4) ** 2
    )
    require(ratio_single == Fraction(10, 21), "single-radius floor drift")

    two_shell = [(Fraction(1), Fraction(1)), (Fraction(1), Fraction(4))]
    ratio_two = Fraction(10, 21) * (
        _mu(two_shell, 2) * _mu(two_shell, 6) / _mu(two_shell, 4) ** 2
    )
    require(ratio_two > Fraction(10, 21), "two-shell member is not above the floor")

    # Signed-weight control: with one negative weight the floor fails while
    # the branch normalizes (mu2 > 0, mu4 > 0).
    signed = [(Fraction(1), Fraction(1)), (Fraction(-1, 100), Fraction(4))]
    require(
        _mu(signed, 2) > 0 and _mu(signed, 4) > 0,
        "signed control loses normalizability",
    )
    ratio_signed = Fraction(10, 21) * (
        _mu(signed, 2) * _mu(signed, 6) / _mu(signed, 4) ** 2
    )
    require(
        ratio_signed < Fraction(10, 21),
        "signed control does not violate the floor",
    )

    return {
        "statement": (
            "mu2 mu6 - mu4^2 is the positive Lagrange combination of "
            "(r_i^2 - r_j^2)^2 terms, so B0/C4^2 >= 10/21 for every "
            "positive-weight member with equality exactly at one radius; "
            "a signed-weight control violates the floor, so one-sign "
            "nonnegativity is load-bearing and active shells are taken "
            "strictly positive"
        ),
        "floor": "10/21",
        "effective_shell_weight": "W_i = w_i |O_i|",
        "identity": "mu2 mu6 - mu4^2 = sum_{i<j} W_i W_j r_i^2 r_j^2 (r_i^2 - r_j^2)^2",
        "two_shell_symbolic": "verified monomial by monomial",
        "random_instances": symbolic_checks,
        "signed_control": {
            "shells": "weights (1, -1/100), squared radii (1, 4)",
            "ratio": str(ratio_signed),
            "verdict": "below the floor, so positivity cannot be dropped",
        },
    }


# ---------------------------------------------------------------------------
# Rank-six band and interior points
# ---------------------------------------------------------------------------


def rank_six_band_certificate(i6_raw) -> dict[str, Any]:
    """B6/B0 = (16/75) <I6(seed)> lies in [-16/135, 16/75] for every
    member, with attained endpoints, the edge point at -1/15, and an
    exact interior zero.

    General members average the seed I6 values with the positive weights
    w |O| r^6, because the rank-six numerator sums w beta r^6 with
    beta = (|O|/60)(64/35) I6(seed) and the isotropic denominator sums
    (w |O| r^6)/7. Single-radius members reduce the weights to w |O|.
    """

    orbits = universality.orbit_directions()
    values = {}
    for name in ("vertex_12", "edge_30", "face_20"):
        values[name] = i6_at_unit(i6_raw, orbits[name]["dirs"][0])
    require(
        values["vertex_12"] == ONE
        and values["edge_30"] == q5(Fraction(-5, 16))
        and values["face_20"] == q5(Fraction(-5, 9)),
        "orbit I6 values drift",
    )
    # I6 at every one of the 62 stationary directions, recomputed here;
    # completeness of the census is the fixed-point packet theorem.
    census = []
    for name, expected_value in (
        ("vertex_12", ONE),
        ("face_20", q5(Fraction(-5, 9))),
        ("edge_30", q5(Fraction(-5, 16))),
    ):
        for direction in orbits[name]["dirs"]:
            census.append(i6_at_unit(i6_raw, direction))
            require(
                census[-1] == expected_value,
                f"{name}: stationary value drift",
            )
    require(len(census) == 62, "stationary census size drift")

    # Every member averages I6 seed values with the positive weights
    # w |O| r^6, so the reachable set is the convex hull of the seed
    # values; the range statement for generic seeds consumes the
    # 62-direction stationary census of the fixed-point packet.
    band_low = Fraction(16, 75) * Fraction(-5, 9)
    band_high = Fraction(16, 75) * Fraction(1)
    require(
        band_low == Fraction(-16, 135) and band_high == Fraction(16, 75),
        "band endpoints drift",
    )

    # General-member control: a three-shell mixed-orbit member evaluated
    # by two exact routes. Route one sums the orbit rank-six multiples
    # directly; route two applies the weighted-mean formula. Both give
    # B6/B0 = -383632/5682975, inside the band.
    beta_by_orbit = {
        "vertex_12": Fraction(64, 175),
        "edge_30": Fraction(-2, 7),
        "face_20": Fraction(-64, 189),
    }
    i6_by_orbit = {
        "vertex_12": Fraction(1),
        "edge_30": Fraction(-5, 16),
        "face_20": Fraction(-5, 9),
    }
    size_by_orbit = {"vertex_12": 12, "edge_30": 30, "face_20": 20}
    mixed = [
        (Fraction(1), Fraction(1), "vertex_12"),
        (Fraction(1), Fraction(4), "edge_30"),
        (Fraction(3, 7), Fraction(9, 4), "face_20"),
    ]
    t6 = sum(w * beta_by_orbit[o] * r2**3 for w, r2, o in mixed)
    mu6 = sum(w * size_by_orbit[o] * r2**3 for w, r2, o in mixed)
    direct_ratio = 7 * t6 / mu6
    mean_numerator = sum(
        w * size_by_orbit[o] * r2**3 * i6_by_orbit[o] for w, r2, o in mixed
    )
    mean_ratio = Fraction(16, 75) * mean_numerator / mu6
    require(
        direct_ratio == mean_ratio == Fraction(-383632, 5682975),
        "general-member band routes disagree",
    )
    require(
        Fraction(-16, 135) <= direct_ratio <= Fraction(16, 75),
        "general member escapes the band",
    )
    # Interior zero: vertex weight 25, face weight 27 at one radius.
    # Per-direction weights 25 on the vertex orbit and 27 on the face
    # orbit; the orbit-sum multiples 64/175 and -64/189 carry the
    # orbit sizes, and 25*(64/175) = 27*(64/189) = 64/7.
    zero_mix = 25 * Fraction(64, 175) + 27 * Fraction(-64, 189)
    require(zero_mix == 0, "vertex-face 25:27 mixture does not cancel")

    generic = i6_at_unit(i6_raw, (q5(1), q5(2), q5(3)))
    require(
        q5_sign(q5_sub(generic, q5(Fraction(-5, 9)))) > 0
        and q5_sign(q5_sub(ONE, generic)) > 0,
        "generic seed value escapes the census range",
    )
    return {
        "statement": (
            "every member has B6/B0 = (16/75) times the mean seed I6 "
            "value with positive weights w |O| r^6, so the band is "
            "[-16/135, 16/75] for the declared class, with the face and "
            "vertex branches at the endpoints, the edge branch at -1/15, "
            "and an exact zero at the 25:27 vertex-face mixture; the "
            "census completeness input is the 62-direction stationary "
            "theorem of the fixed-point packet"
        ),
        "band": ["-16/135", "16/75"],
        "band_scope": "every member of the class",
        "general_member_control": {
            "shells": "weights (1, 1, 3/7), squared radii (1, 4, 9/4), orbits (vertex, edge, face)",
            "B6_over_B0": "-383632/5682975",
            "routes": "direct rank-six sum and weighted-mean formula agree",
        },
        "edge_point": "-1/15",
        "interior_zero": "vertex weight 25, face weight 27, one radius",
        "census_recomputed": 62,
        "census_completeness_input": (
            "a5_multipole fixed-point packet, 62-direction stationary "
            "census with values 1, -5/9, -5/16"
        ),
    }


# ---------------------------------------------------------------------------
# Eighth order: I6-only anisotropy and the cross-order lock
# ---------------------------------------------------------------------------


def eighth_order_certificate(
    rotations: list[Mat], i6_raw, i6_reduced
) -> dict[str, Any]:
    """Degree-eight kernel factorization, branch table, and the 12/5 lock.

    The degree-eight kernel is invariant in each argument; the recomputed
    invariant multiplicities give m8 = 0, so the reduced kernel lies in
    span{1, I6} in each argument and the same two-seed argument closes the
    identity for every seed. The residue-zero checks below verify the
    absence of a rank-eight component at the sampled seeds.
    """

    probe = next(m for m in sorted(i6_reduced) if sum(m) > 0)
    verts = base.cartesian_vertices()
    orbits = universality.orbit_directions()
    invariants = universality.invariant_table_certificate()
    even_multiplicities = [
        invariants["table"][str(level)] for level in (0, 2, 4, 6, 8)
    ]
    require(
        even_multiplicities == [1, 0, 0, 1, 0],
        "harmonic even-rank multiplicities through rank eight drift",
    )

    seeds = [
        ("vertex", verts[0]),
        ("edge", orbits["edge_30"]["dirs"][0]),
        ("face", orbits["face_20"]["dirs"][0]),
        ("generic_0", (q5(1), q5(2), q5(3))),
        ("generic_1", (q5(2), q5(-1), q5(5))),
    ]
    s_constant: Q5 | None = None
    rows = []
    for name, seed in seeds:
        norm_sq = ZERO
        for axis in range(3):
            norm_sq = q5_add(norm_sq, q5_mul(seed[axis], seed[axis]))
        kernel = p_zero()
        for g in rotations:
            kernel = p_add(kernel, p_pow(p_linear(mat_apply(g, seed)), 8))
        kernel = p_scale_q5(kernel, q5_div(ONE, q5_pow(norm_sq, 4)))
        reduced = p_reduce_sphere(kernel)
        multiple = q5_div(reduced.get(probe, ZERO), i6_reduced[probe])
        residue = p_add(reduced, p_scale_q5(i6_reduced, q5_neg(multiple)))
        iso = residue.pop((0, 0, 0), ZERO)
        require(
            iso == q5(Fraction(20, 3)),
            f"{name}: degree-eight isotropic part is not 20/3",
        )
        require(
            p_is_zero(residue),
            f"{name}: degree-eight kernel leaves the invariant line",
        )
        seed_value = i6_at_unit(i6_raw, seed)
        if seed_value != ZERO:
            ratio = q5_div(multiple, seed_value)
            if s_constant is None:
                s_constant = ratio
            require(
                ratio == s_constant,
                f"{name}: degree-eight constant is not universal",
            )
        rows.append({"seed": name, "i6_multiple": q5_str(multiple)})
    require(
        s_constant == q5(Fraction(256, 75)),
        "degree-eight kernel constant is not 256/75",
    )

    # Branch table: normalized k^8 coefficients per pure unit-radius orbit.
    # Prefactor 6/(|O| a^2) and the -x^8/40320 series term give
    # D0 = -a^6/60480 and D6 = -a^6 m8 / (6720 |O|) with m8 the orbit-sum
    # degree-eight I6 multiple (|O|/60)(256/75) I6(seed).
    expected = {
        "vertex_12": (Fraction(256, 375), Fraction(64, 125)),
        "edge_30": (Fraction(-8, 15), Fraction(-4, 25)),
        "face_20": (Fraction(-256, 405), Fraction(-64, 225)),
    }
    branch_rows = []
    for name, count in (("vertex_12", 12), ("edge_30", 30), ("face_20", 20)):
        data = orbits[name]
        dirs, norm_sq = data["dirs"], data["norm_sq"]
        m8 = p_reduce_sphere(
            p_scale_q5(
                moment_sum(dirs, 8), q5_div(ONE, q5_pow(norm_sq, 4))
            )
        )
        multiple = q5_div(m8.get(probe, ZERO), i6_reduced[probe])
        residue = p_add(m8, p_scale_q5(i6_reduced, q5_neg(multiple)))
        iso = residue.pop((0, 0, 0), ZERO)
        require(
            iso == q5(Fraction(count, 9)),
            f"{name}: eighth-moment isotropic part is not count/9",
        )
        require(p_is_zero(residue), f"{name}: eighth moment off the line")
        want_m8, want_ratio = expected[name]
        require(multiple == q5(want_m8), f"{name}: m8 drift")
        # D6/D0 = 9 m8 / |O|
        require(
            q5_scale(multiple, Fraction(9, count)) == q5(want_ratio),
            f"{name}: eighth-order ratio drift",
        )
        branch_rows.append(
            {
                "orbit": name,
                "m8_i6_multiple": str(want_m8),
                "D0_over_a6": "-1/60480",
                "D6_over_D0": str(want_ratio),
                "radius_scope": "unit radius; at radius r, D0/a^6 = -r^6/60480",
            }
        )

    # Cross-order lock, with its genuinely division-free polynomial form:
    # D6/D0 = (12/5) (B6/B0), equivalently
    # 5 D6 B0 = 12 B6 D0.  It is support-independent on single-radius
    # members and remains meaningful when both anisotropies vanish.
    require(
        Fraction(64, 125) / Fraction(16, 75) == Fraction(12, 5),
        "cross-order lock is not 12/5",
    )
    for name, i6_seed in (
        ("vertex_12", Fraction(1)),
        ("edge_30", Fraction(-5, 16)),
        ("face_20", Fraction(-5, 9)),
    ):
        b6_ratio = Fraction(16, 75) * i6_seed
        d6_ratio = Fraction(64, 125) * i6_seed
        require(
            d6_ratio == Fraction(12, 5) * b6_ratio,
            f"{name}: lock fails",
        )

    # Exact zero-anisotropy control.  At one common radius, vertex weight
    # 25 and face weight 27 give orbit-weighted I6 numerator
    # 25*12*1 + 27*20*(-5/9) = 0.  Both ratios vanish, while the
    # division-free identity remains defined and true.
    zero_numerator = (
        Fraction(25 * 12) * Fraction(1)
        + Fraction(27 * 20) * Fraction(-5, 9)
    )
    zero_denominator = Fraction(25 * 12 + 27 * 20)
    require(zero_numerator == 0, "25:27 zero-mixture numerator drift")
    zero_mean = zero_numerator / zero_denominator
    zero_b6_ratio = Fraction(16, 75) * zero_mean
    zero_d6_ratio = Fraction(64, 125) * zero_mean
    require(
        zero_d6_ratio == Fraction(12, 5) * zero_b6_ratio == 0,
        "cross-order lock fails at the zero-anisotropy mixture",
    )

    # A two-radius member must not satisfy the single-radius lock in general.
    # This guards against accidentally widening the theorem to the complete
    # multi-radius class.  The effective weights include orbit sizes.
    multi_b6_ratio = Fraction(16, 75) * (
        Fraction(12) * 1 ** 3 * 1
        + Fraction(20) * 4 ** 3 * Fraction(-5, 9)
    ) / (Fraction(12) * 1 ** 3 + Fraction(20) * 4 ** 3)
    multi_d6_ratio = Fraction(64, 125) * (
        Fraction(12) * 1 ** 4 * 1
        + Fraction(20) * 4 ** 4 * Fraction(-5, 9)
    ) / (Fraction(12) * 1 ** 4 + Fraction(20) * 4 ** 4)
    multi_lock_residual = multi_d6_ratio - Fraction(12, 5) * multi_b6_ratio
    require(
        multi_lock_residual == Fraction(-57344, 10360225),
        "two-radius negative-control residual drift",
    )
    require(multi_lock_residual != 0, "multi-radius member incorrectly obeys lock")

    # Vertex eighth-order coefficient against the fixed-point packet
    # template A(x) = 2 x^6 (30 - x^2)/118125: the template x^8
    # coefficient -2/118125 equals twice the normalized D6, matching the
    # factor two between the template and the 1/(2 a^2) symbol.
    d6_vertex = Fraction(-1, 1) * Fraction(256, 375) / (6720 * 12)
    require(
        d6_vertex == Fraction(-1, 118125),
        "vertex eighth-order coefficient drift",
    )
    require(
        Fraction(-2, 118125) == 2 * d6_vertex,
        "fixed-point template cross-check fails",
    )

    return {
        "statement": (
            "the degree-eight kernel factors as 20/3 plus (256/75) I6(u) "
            "I6(n), the harmonic angular-rank-eight multiplicity is zero, every pure "
            "unit-radius orbit has D0 = -a^6/60480 with D6/D0 = "
            "(64/125) I6(seed) in the band [-64/225, 64/125], and every "
            "single-radius member obeys 5 D6 B0 = 12 B6 D0, "
            "equivalently D6/D0 = (12/5)(B6/B0), including the 25:27 zero-anisotropy "
            "vertex-face mixture; at common radius r the isotropic "
            "coefficient is D0 = -(a r)^6/60480; the "
            "isotropic tower alternates in sign at every order"
        ),
        "general_isotropic_coefficient": "D0 = -(a^6/60480)(mu8/mu2)",
        "general_anisotropic_ratio": (
            "D6/D0 = (64/125) <I6(seed)> with positive weights w |O| r^8"
        ),
        "kernel_isotropic_part": "20/3",
        "kernel_universal_constant": "256/75",
        "harmonic_invariant_multiplicities_l0_l2_l4_l6_l8": (
            even_multiplicities
        ),
        "eighth_order_band": ["-64/225", "64/125"],
        "cross_order_lock": "12/5",
        "cross_order_ratio_identity": "D6/D0 = (12/5)(B6/B0)",
        "cross_order_identity": "5 D6 B0 = 12 B6 D0",
        "zero_mixture_control": {
            "member": "vertex weight 25, face weight 27, one common radius",
            "B6_over_B0": "0",
            "D6_over_D0": "0",
            "division_free_identity": "5 D6 B0 = 12 B6 D0 = 0",
        },
        "multi_radius_negative_control": {
            "member": "vertex radius 1 plus face radius 2, equal per-direction weights",
            "B6_over_B0": str(multi_b6_ratio),
            "D6_over_D0": str(multi_d6_ratio),
            "lock_residual": str(multi_lock_residual),
            "verdict": "nonzero; the cross-order lock is single-radius only",
        },
        "sign_alternation": (
            "positive even moments and the alternating cosine series give "
            "isotropic coefficients of sign (-1)^(m+1) at order k^(2m) "
            "for every member"
        ),
        "fixed_point_template_crosscheck": (
            "the vertex D6 = -a^6/118125 equals half the -2/118125 x^8 "
            "coefficient of the certified through-eighth template"
        ),
        "kernel_seeds": rows,
        "branches": branch_rows,
    }


# ---------------------------------------------------------------------------
# Receipt assembly
# ---------------------------------------------------------------------------


def build_receipt() -> dict[str, Any]:
    cartesian = base.build_cartesian_frame()
    i6_raw = cartesian.pop("_i6_poly_object")
    cartesian.pop("_vertices_object", None)
    i6_reduced = p_reduce_sphere(i6_raw)

    rotations, rotation_cert = rotation_group_certificate()
    d_constant, kernel_cert = kernel_factorization_certificate(
        rotations, i6_raw, i6_reduced
    )
    branches = orbit_branch_certificate(i6_raw, i6_reduced, d_constant)
    eighth = eighth_order_certificate(rotations, i6_raw, i6_reduced)
    series = series_route_control()
    frozen = frozen_receipt_crosscheck()
    floor = radial_floor_certificate()
    band = rank_six_band_certificate(i6_raw)

    body: dict[str, Any] = {
        "schema": SCHEMA,
        "status": STATUS,
        "class_definition": (
            "positive-weight finite mixtures of A5 orbits at finitely many "
            "positive radii with finite scale a > 0, hop symbol sum w "
            "[1 - cos(a k r u.n)], "
            "continuum k^2 normalization, with this cosine sum the full "
            "spatial symbol and hence complete through k^8; invariance forces equal weights "
            "inside each orbit and the cosine form is even"
        ),
        "class_statements": {
            "sign_law": "C4 < 0 for every member",
            "isotropic_floor": (
                "B0/C4^2 = (10/21)(mu2 mu6/mu4^2) >= 10/21, equality "
                "exactly on single-radius members"
            ),
            "rank_six_band": (
                "every member has B6/B0 = (16/75) <I6(seed)> with "
                "positive weights w |O| r^6, confined to "
                "[-16/135, 16/75]; single-radius members reduce the "
                "weights to w |O|"
            ),
            "rank_purity": (
                "ranks one through five empty and rank six on the I6 line, "
                "consumed from the universality certificate"
            ),
            "eighth_order_lock": (
                "D0 = -(a^6/60480)(mu8/mu2), and the k^8 anisotropy is the "
                "same rotated I6 with "
                "D6/D0 = (64/125) <I6(seed)> using positive weights "
                "w |O| r^8, and every single-radius "
                "member obeys 5 D6 B0 = 12 B6 D0, equivalently "
                "D6/D0 = (12/5)(B6/B0), including the zero-anisotropy "
                "member; multi-radius members need not obey this lock; "
                "D0 = -(a r)^6/60480 at common radius r"
            ),
        },
        "kill_surface": (
            "under the registered physical-symbol, sector, frame, scale, "
            "readout, and exclusivity premises, a resolved intrinsic "
            "dispersion with B0/C4^2 below 10/21, or "
            "B6/B0 outside [-16/135, 16/75], or a rank one-to-five "
            "residue, or a rank-six residue off the rotated I6 line, "
            "excludes every member of the declared positive-weight scalar "
            "cosine class; "
            "exact saturation of the floor at 10/21 certifies a single-radius "
            "carrier"
        ),
        "physical_premises_unchanged": (
            "sector, frame, scale, and exclusivity premises of FZ-11 and "
            "FZ-12 are not constructed here; their target-clean producer, "
            "physical bridge, quantitative output, and pre-comparison "
            "registration belong to the relevant propagation and prediction "
            "lanes; this "
            "certificate adds no comparison and reads no data"
        ),
        "rotation_group": rotation_cert,
        "kernel_factorization": kernel_cert,
        "single_orbit_branches": branches,
        "eighth_order": eighth,
        "series_route_control": series,
        "frozen_receipt_crosscheck": frozen,
        "radial_floor": floor,
        "rank_six_band": band,
    }
    body["receipt_sha256"] = base.tagged_sha256(
        base.canonical_json_bytes(body)
    )
    return body


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args(argv)
    receipt = build_receipt()
    if args.write:
        RECEIPT_PATH.write_text(
            base.canonical_json_bytes(receipt).decode() + "\n"
        )
        print(f"wrote {RECEIPT_PATH}")
    print(receipt["status"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
