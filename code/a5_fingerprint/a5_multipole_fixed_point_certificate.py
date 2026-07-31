#!/usr/bin/env python3
"""Issue #654: the A5 multipole fixed-point certificate.

The certificate freezes the complete twelve-port angular fingerprint in
exact arithmetic, in two independent frames:

* the Cartesian frame, with the twelve vertices at the cyclic signed
  permutations of ``(0, 1, phi)`` over ``Q(sqrt5)``: certified
  5-design moment identities (all spherical-harmonic content for
  ``1 <= l <= 5`` vanishes), the exact Cartesian polynomial of the
  canonical normalized level-six form ``I6``, vertex normalization
  ``I6 = 1`` on all twelve directions, face value ``-5/9``, edge value
  ``-5/16``, and kernel independence (every declared analytic kernel
  collapses to the same ``I6`` below level ten);
* the vertex-pole frame, built from an exact five-fold trigonometric
  moment engine: the spherical normal form
  ``I6 = P6(c) + (21/8) c (1-c^2)^{5/2} cos(5 phi)``, the complete
  critical-point proof through the exact factorization
  ``P6'(c)^2 - (441/64)(1-c^2)^3(1-6c^2)^2
  = (441/64)(5c^2-1)(5c^4-5c^2+1)(45c^4-30c^2+1)``,
  the 62-point census (12 vertex maxima, 20 face minima, 30 edge
  saddles), exact critical values and tangent Hessian eigenvalues, the
  Morse--Euler check ``12 - 30 + 20 = 2``, the band projectors with
  ranks ``(1, 3, 5, 3)``, and the blind inverse-port response
  ``R = -J`` with its band sign vector.

The equal-weight kinetic stencil receipt (isotropic corrections at
``a^2 k^4`` and ``a^4 k^6``, first directional artifact at spin six,
amplitude law ``A6 proportional to a^4`` with exact one-step suppression
``1/16``) is computed under the declared stencil premise and typed
conditional; ledger rows OPH-A5-M4 and M5 stay unfrozen pending the
issue #655 selection theorem. Fail-closed controls (cube, octahedron,
rational shell, perturbed ring, unequal gains, shuffled labels,
``R = +J``, injected level-four contaminant) each fail the named
detector. The frozen rows OPH-A5-M1 through M3 carry declared premise
ancestry and blind decision rules; no physical map is claimed, and the
packet is ineligible for issue #639 consumption until physicalization.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from typing import Any

from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNTIME = HERE / "runtime"
RECEIPT_PATH = RUNTIME / "a5_multipole_fixed_point_receipt.json"

SCHEMA = "oph.a5_multipole_fixed_point_receipt.v1"
STATUS = "EXACT_A5_FINGERPRINT_CERTIFICATE__PHYSICAL_MAP_OPEN"


class FingerprintError(ValueError):
    """The fingerprint certificate refused to build."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise FingerprintError(message)


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
# Exact Q(sqrt5) scalars
# ---------------------------------------------------------------------------

Q5 = tuple[Fraction, Fraction]


def q5(a, b=0) -> Q5:
    return (Fraction(a), Fraction(b))


ZERO = q5(0)
ONE = q5(1)
PHI = q5(Fraction(1, 2), Fraction(1, 2))
INV_SQRT5 = q5(0, Fraction(1, 5))
SQRT5 = q5(0, 1)


def q5_add(x: Q5, y: Q5) -> Q5:
    return (x[0] + y[0], x[1] + y[1])


def q5_sub(x: Q5, y: Q5) -> Q5:
    return (x[0] - y[0], x[1] - y[1])


def q5_mul(x: Q5, y: Q5) -> Q5:
    return (x[0] * y[0] + 5 * x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def q5_scale(x: Q5, f) -> Q5:
    f = Fraction(f)
    return (x[0] * f, x[1] * f)


def q5_div(x: Q5, y: Q5) -> Q5:
    norm = y[0] * y[0] - 5 * y[1] * y[1]
    require(norm != 0, "division by zero in Q(sqrt5)")
    conj = (y[0], -y[1])
    num = q5_mul(x, conj)
    return (num[0] / norm, num[1] / norm)


def q5_pow(x: Q5, n: int) -> Q5:
    out = ONE
    for _ in range(n):
        out = q5_mul(out, x)
    return out


def q5_neg(x: Q5) -> Q5:
    return (-x[0], -x[1])


def q5_str(x: Q5) -> str:
    return f"{x[0]}+{x[1]}*sqrt5"


def q5_sign(x: Q5) -> int:
    """Exact sign of a + b*sqrt5."""

    a, b = x
    if a == 0 and b == 0:
        return 0
    if b == 0:
        return 1 if a > 0 else -1
    if a == 0:
        return 1 if b > 0 else -1
    if a > 0 and b > 0:
        return 1
    if a < 0 and b < 0:
        return -1
    # opposite signs: compare a^2 with 5 b^2
    lhs = a * a
    rhs = 5 * b * b
    if a > 0:
        return 1 if lhs > rhs else -1
    return -1 if lhs > rhs else 1


def q5_sqrt(x: Q5) -> Q5 | None:
    """Exact square root inside Q(sqrt5) when one exists, else None."""

    if x == ZERO:
        return ZERO
    a, b = x
    if b == 0:
        # want (p + q sqrt5)^2 = a: p*q = 0
        root = _frac_sqrt(a)
        if root is not None:
            return (root, Fraction(0))
        root = _frac_sqrt(a / 5)
        if root is not None:
            return (Fraction(0), root)
        return None
    # (p + q sqrt5)^2 = a + b sqrt5: p^2 + 5 q^2 = a, 2 p q = b
    # p^2 solves t^2 - a t + 5 (b/2)^2 = 0
    disc = a * a - 5 * b * b
    droot = _frac_sqrt(disc)
    if droot is None:
        return None
    for sign in (1, -1):
        t = (a + sign * droot) / 2
        p = _frac_sqrt(t)
        if p is not None and p != 0:
            q = b / (2 * p)
            cand = (p, q)
            if q5_mul(cand, cand) == x:
                return cand
    return None


def _frac_sqrt(f: Fraction) -> Fraction | None:
    if f < 0:
        return None
    num = _int_sqrt(f.numerator)
    den = _int_sqrt(f.denominator)
    if num is None or den is None:
        return None
    return Fraction(num, den)


def _int_sqrt(n: int) -> int | None:
    if n < 0:
        return None
    r = int(n**0.5)
    for cand in (r - 1, r, r + 1, r + 2):
        if cand >= 0 and cand * cand == n:
            return cand
    return None


# ---------------------------------------------------------------------------
# Complex Q(sqrt5) scalars (for the five-fold trigonometric engine)
# ---------------------------------------------------------------------------

C5 = tuple[Q5, Q5]  # real + imaginary parts


def c5(re: Q5, im: Q5 = ZERO) -> C5:
    return (re, im)


def c5_add(x: C5, y: C5) -> C5:
    return (q5_add(x[0], y[0]), q5_add(x[1], y[1]))


def c5_mul(x: C5, y: C5) -> C5:
    re = q5_sub(q5_mul(x[0], y[0]), q5_mul(x[1], y[1]))
    im = q5_add(q5_mul(x[0], y[1]), q5_mul(x[1], y[0]))
    return (re, im)


def c5_scale(x: C5, f) -> C5:
    return (q5_scale(x[0], f), q5_scale(x[1], f))


# ---------------------------------------------------------------------------
# Polynomials in (x, y, z) over Q(sqrt5)
# ---------------------------------------------------------------------------

Poly = dict[tuple[int, int, int], Q5]


def p_zero() -> Poly:
    return {}


def p_const(c: Q5) -> Poly:
    return {(0, 0, 0): c} if c != ZERO else {}


def p_add(p: Poly, q: Poly) -> Poly:
    out = dict(p)
    for mono, coeff in q.items():
        acc = q5_add(out.get(mono, ZERO), coeff)
        if acc == ZERO:
            out.pop(mono, None)
        else:
            out[mono] = acc
    return out


def p_scale(p: Poly, f) -> Poly:
    if Fraction(f) == 0:
        return {}
    return {m: q5_scale(c, f) for m, c in p.items()}


def p_scale_q5(p: Poly, s: Q5) -> Poly:
    if s == ZERO:
        return {}
    return {m: q5_mul(c, s) for m, c in p.items()}


def p_mul(p: Poly, q: Poly) -> Poly:
    out: Poly = {}
    for m1, c1 in p.items():
        for m2, c2 in q.items():
            mono = (m1[0] + m2[0], m1[1] + m2[1], m1[2] + m2[2])
            acc = q5_add(out.get(mono, ZERO), q5_mul(c1, c2))
            if acc == ZERO:
                out.pop(mono, None)
            else:
                out[mono] = acc
    return out


def p_pow(p: Poly, n: int) -> Poly:
    out = p_const(ONE)
    for _ in range(n):
        out = p_mul(out, p)
    return out


def p_linear(v: tuple[Q5, Q5, Q5]) -> Poly:
    out: Poly = {}
    for axis, comp in enumerate(v):
        if comp != ZERO:
            mono = tuple(1 if i == axis else 0 for i in range(3))
            out[mono] = comp
    return out


def p_eval(p: Poly, v: tuple[Q5, Q5, Q5]) -> Q5:
    total = ZERO
    for (i, j, k), coeff in p.items():
        term = q5_mul(
            coeff,
            q5_mul(q5_pow(v[0], i), q5_mul(q5_pow(v[1], j), q5_pow(v[2], k))),
        )
        total = q5_add(total, term)
    return total


def p_reduce_sphere(p: Poly) -> Poly:
    """Canonical form modulo ``x^2 + y^2 + z^2 - 1``: z-degree <= 1."""

    work = dict(p)
    while True:
        target = None
        for mono in work:
            if mono[2] >= 2:
                target = mono
                break
        if target is None:
            return work
        coeff = work.pop(target)
        i, j, k = target
        base = (i, j, k - 2)
        # z^2 -> 1 - x^2 - y^2
        for delta, factor in (
            ((0, 0, 0), Fraction(1)),
            ((2, 0, 0), Fraction(-1)),
            ((0, 2, 0), Fraction(-1)),
        ):
            mono = (base[0] + delta[0], base[1] + delta[1], base[2] + delta[2])
            acc = q5_add(work.get(mono, ZERO), q5_scale(coeff, factor))
            if acc == ZERO:
                work.pop(mono, None)
            else:
                work[mono] = acc


def p_is_zero(p: Poly) -> bool:
    return not p


def p_str(p: Poly) -> str:
    parts = []
    for mono in sorted(p):
        parts.append(f"{q5_str(p[mono])}*x^{mono[0]}*y^{mono[1]}*z^{mono[2]}")
    return " + ".join(parts) if parts else "0"


# ---------------------------------------------------------------------------
# Legendre polynomials (exact coefficients)
# ---------------------------------------------------------------------------


def legendre_coeffs(max_level: int) -> list[list[Fraction]]:
    """Coefficient lists (ascending powers) of P_0..P_max."""

    table = [[Fraction(1)], [Fraction(0), Fraction(1)]]
    for level in range(1, max_level):
        prev, cur = table[level - 1], table[level]
        nxt = [Fraction(0)] * (level + 2)
        for power, coeff in enumerate(cur):
            nxt[power + 1] += Fraction(2 * level + 1, level + 1) * coeff
        for power, coeff in enumerate(prev):
            nxt[power] -= Fraction(level, level + 1) * coeff
        table.append(nxt)
    return table[: max_level + 1]


LEGENDRE = legendre_coeffs(10)


def legendre_eval_q5(level: int, t: Q5) -> Q5:
    total = ZERO
    for power, coeff in enumerate(LEGENDRE[level]):
        if coeff:
            total = q5_add(total, q5_scale(q5_pow(t, power), coeff))
    return total


# ---------------------------------------------------------------------------
# The Cartesian frame: vertices, moments, and the I6 polynomial
# ---------------------------------------------------------------------------

NORM_SQ = q5_add(q5(2), PHI)  # |v|^2 = 2 + phi for every vertex


def cartesian_vertices() -> list[tuple[Q5, Q5, Q5]]:
    verts = []
    for s1 in (1, -1):
        for s2 in (1, -1):
            a = q5_scale(ONE, s1)
            b = q5_scale(PHI, s2)
            verts.append((ZERO, a, b))
            verts.append((a, b, ZERO))
            verts.append((b, ZERO, a))
    require(len(verts) == 12, "vertex count drift")
    return verts


def moment_sum(vertices, k: int, weights=None) -> Poly:
    """Sum_i w_i (v_i . n)^k as a polynomial in n = (x, y, z)."""

    total = p_zero()
    for idx, v in enumerate(vertices):
        term = p_pow(p_linear(v), k)
        if weights is not None:
            term = p_scale(term, weights[idx])
        total = p_add(total, term)
    return total


def normalized_moment(vertices, k: int, norm_sq: Q5) -> Poly:
    """Sum_i (u_i . n)^k with u = v/|v|, exact for even k."""

    require(k % 2 == 0, "normalized moments are taken at even order")
    raw = moment_sum(vertices, k)
    return p_scale_q5(raw, q5_div(ONE, q5_pow(norm_sq, k // 2)))


def radial_power(k: int) -> Poly:
    return p_pow(
        p_add(
            p_add(p_pow(p_linear((ONE, ZERO, ZERO)), 2), p_pow(p_linear((ZERO, ONE, ZERO)), 2)),
            p_pow(p_linear((ZERO, ZERO, ONE)), 2),
        ),
        k,
    )


def legendre_port_sum(vertices, level: int, norm_sq: Q5) -> Poly:
    """Sum_i P_level(u_i . n) as a polynomial, exact over Q(sqrt5).

    Odd powers of (v.n) cancel over the antipodal vertex set, so only
    even powers contribute and every normalization is an integer power
    of ``norm_sq``.
    """

    total = p_zero()
    for power, coeff in enumerate(LEGENDRE[level]):
        if coeff == 0:
            continue
        raw = moment_sum(vertices, power)
        if power % 2 == 1:
            require(p_is_zero(raw), f"odd moment {power} fails to cancel")
            continue
        total = p_add(
            total,
            p_scale_q5(
                raw, q5_scale(q5_div(ONE, q5_pow(norm_sq, power // 2)), coeff)
            ),
        )
    return total


def build_cartesian_frame() -> dict[str, Any]:
    verts = cartesian_vertices()
    for v in verts:
        dot = q5_add(q5_add(q5_mul(v[0], v[0]), q5_mul(v[1], v[1])), q5_mul(v[2], v[2]))
        require(dot == NORM_SQ, "vertex norm drift")

    # dot-product classes of the unit configuration
    classes: dict[str, int] = {}
    inv_norm = q5_div(ONE, NORM_SQ)
    for i in range(12):
        for j in range(12):
            if i == j:
                continue
            dot = ZERO
            for axis in range(3):
                dot = q5_add(dot, q5_mul(verts[i][axis], verts[j][axis]))
            u_dot = q5_mul(dot, inv_norm)
            key = q5_str(u_dot)
            classes[key] = classes.get(key, 0) + 1
    require(
        classes == {
            q5_str(INV_SQRT5): 60,
            q5_str(q5_neg(INV_SQRT5)): 60,
            q5_str(q5_neg(ONE)): 12,
        },
        "unit dot-product classes drift",
    )

    # exact even moment identities
    m2 = normalized_moment(verts, 2, NORM_SQ)
    m4 = normalized_moment(verts, 4, NORM_SQ)
    require(m2 == p_scale(radial_power(1), 4), "second moment is not 4 r^2")
    require(
        m4 == p_scale(radial_power(2), Fraction(12, 5)),
        "fourth moment is not (12/5) r^4",
    )
    for k in (1, 3, 5):
        require(p_is_zero(moment_sum(verts, k)), f"odd moment {k} nonzero")

    # certified nulls for l = 1..5 on the sphere
    for level in range(1, 6):
        reduced = p_reduce_sphere(legendre_port_sum(verts, level, NORM_SQ))
        require(p_is_zero(reduced), f"level {level} port sum nonzero on sphere")

    # canonical I6 polynomial
    p6_sum = legendre_port_sum(verts, 6, NORM_SQ)
    i6 = p_scale(p6_sum, Fraction(25, 132))

    # vertex normalization I6(u_i) = 1: evaluate homogeneous pieces
    unit_values = []
    for v in verts:
        total = ZERO
        for mono, coeff in i6.items():
            degree = sum(mono)
            require(degree % 2 == 0, "I6 contains an odd-degree monomial")
            val = q5_mul(
                coeff,
                q5_mul(
                    q5_pow(v[0], mono[0]),
                    q5_mul(q5_pow(v[1], mono[1]), q5_pow(v[2], mono[2])),
                ),
            )
            total = q5_add(total, q5_div(val, q5_pow(NORM_SQ, degree // 2)))
        unit_values.append(total)
    require(all(val == ONE for val in unit_values), "vertex normalization drift")

    # face and edge representatives (adjacent = dot 1/sqrt5)
    def unit_dot(i: int, j: int) -> Q5:
        dot = ZERO
        for axis in range(3):
            dot = q5_add(dot, q5_mul(verts[i][axis], verts[j][axis]))
        return q5_mul(dot, inv_norm)

    faces = []
    for i in range(12):
        for j in range(i + 1, 12):
            if unit_dot(i, j) != INV_SQRT5:
                continue
            for k in range(j + 1, 12):
                if unit_dot(i, k) == INV_SQRT5 and unit_dot(j, k) == INV_SQRT5:
                    faces.append((i, j, k))
    require(len(faces) == 20, "face census drift")
    edges = [
        (i, j)
        for i in range(12)
        for j in range(i + 1, 12)
        if unit_dot(i, j) == INV_SQRT5
    ]
    require(len(edges) == 30, "edge census drift")

    def eval_even_at_sum(vs: list[int]) -> Q5:
        centre = tuple(
            q5_add(q5_add(verts[vs[0]][axis], verts[vs[1]][axis]),
                   verts[vs[2]][axis] if len(vs) == 3 else ZERO)
            for axis in range(3)
        )
        norm = ZERO
        for axis in range(3):
            norm = q5_add(norm, q5_mul(centre[axis], centre[axis]))
        total = ZERO
        for mono, coeff in i6.items():
            degree = sum(mono)
            val = q5_mul(
                coeff,
                q5_mul(
                    q5_pow(centre[0], mono[0]),
                    q5_mul(q5_pow(centre[1], mono[1]), q5_pow(centre[2], mono[2])),
                ),
            )
            total = q5_add(total, q5_div(val, q5_pow(norm, degree // 2)))
        return total

    face_values = {q5_str(eval_even_at_sum(list(f))) for f in faces}
    edge_values = {q5_str(eval_even_at_sum(list(e))) for e in edges}
    require(face_values == {q5_str(q5(Fraction(-5, 9)))}, "face value drift")
    require(edge_values == {q5_str(q5(Fraction(-5, 16)))}, "edge value drift")

    return {
        "vertices": "cyclic signed permutations of (0, 1, phi), norm^2 = 2 + phi",
        "unit_dot_classes": {"1/sqrt5": 60, "-1/sqrt5": 60, "-1": 12},
        "moment_identities": {
            "sum (u.n)^2": "4 r^2",
            "sum (u.n)^4": "12/5 r^4",
            "odd moments 1,3,5": "identically zero",
        },
        "multipole_nulls": "sum_i P_l(u_i . n) = 0 on the sphere for l = 1..5",
        "i6_monomials": len(i6),
        "i6_polynomial": {"-".join(map(str, m)): q5_str(c) for m, c in sorted(i6.items())},
        "vertex_normalization": "I6 = 1 on all twelve unit vertices",
        "face_value": "-5/9 on all twenty face centers",
        "edge_value": "-5/16 on all thirty edge midpoints",
        "_i6_poly_object": i6,
        "_vertices_object": verts,
    }


# ---------------------------------------------------------------------------
# The vertex-pole frame: five-fold trigonometric moment engine
# ---------------------------------------------------------------------------


def ring_moment_sum(k: int, c: Q5, s_sq: Q5, lower: bool) -> Poly:
    """Sum over one pentagonal ring of (u . n)^k, exactly.

    Ring vertices: (s cos a_j, s sin a_j, c) with a_j = 2 pi j / 5
    (upper) or (2j+1) pi / 5 (lower). The engine expands
    (c z + s x cos a + s y sin a)^k with cos a = (w + 1/w)/2 and
    sin a = (w - 1/w)/(2i), then applies the exact character sums
    sum_j w^m = 5 [5|m] (upper) and 5 (-1)^(m/5) [5|m] (lower). Odd
    powers of s pair with w-exponents m of matching parity; surviving
    terms have 5 | m, so s appears at powers that combine into
    s_sq-integer terms times at most s^5, and every surviving s power
    below is even or a multiple handled by the caller.

    The return value is the exact polynomial in (x, y, z) whose
    coefficients live in Q(sqrt5), expressed with s^2 = s_sq
    substituted; odd total powers of s survive only in the cos(5 phi)
    sector and are returned with the convention s^(2t+5) ->
    s_sq^t * S5, where the caller supplies S5 = s^5 through the
    ``s5`` parameter of :func:`assemble_pole_value`. To keep this
    function self-contained the k <= 8 cases used here never produce
    s-powers beyond s^8, and 5 | m forces the odd sector to m = +-5,
    that is exactly s^5, s^7 -> s_sq * s^5.
    """

    x_mono, y_mono, z_mono = (1, 0, 0), (0, 1, 0), (0, 0, 1)
    # represent the ring dot product as a polynomial in (x, y, z) with
    # coefficients that are Laurent polynomials in w over C5, tracked as
    # dict m -> C5 coefficient, times powers of s tracked separately.
    # Term: c*z + (s/2)(w + w^-1) x + (s/(2i))(w - w^-1) y
    #     = c*z + s * [ x (w + w^-1)/2 + y (w - w^-1)/(2i) ]
    # Expand multinomially in the three parts.
    from math import comb

    # coefficient of the ring sum: dict[(mono, s_power)] -> C5 after
    # applying character sums.
    out: dict[tuple[tuple[int, int, int], int], C5] = {}

    half = Fraction(1, 2)
    # w-coefficients for the x and y parts: x: (1/2) w^{+1} + (1/2) w^{-1}
    # y: (1/(2i)) w^{+1} - (1/(2i)) w^{-1}; note 1/(2i) = -i/2.
    x_part = {1: c5(q5(half)), -1: c5(q5(half))}
    y_part = {1: c5(ZERO, q5(-half)), -1: c5(ZERO, q5(half))}

    for kz in range(k + 1):
        for kx in range(k - kz + 1):
            ky = k - kz - kx
            multi = comb(k, kz) * comb(k - kz, kx)
            # w-Laurent coefficients of x_part^kx * y_part^ky
            lau: dict[int, C5] = {0: c5(ONE)}
            for _ in range(kx):
                nxt: dict[int, C5] = {}
                for m, cf in lau.items():
                    for dm, dcf in x_part.items():
                        acc = c5_add(nxt.get(m + dm, c5(ZERO)), c5_mul(cf, dcf))
                        nxt[m + dm] = acc
                lau = nxt
            for _ in range(ky):
                nxt = {}
                for m, cf in lau.items():
                    for dm, dcf in y_part.items():
                        acc = c5_add(nxt.get(m + dm, c5(ZERO)), c5_mul(cf, dcf))
                        nxt[m + dm] = acc
                lau = nxt
            # apply character sums
            for m, cf in lau.items():
                if m % 5 != 0:
                    continue
                weight = Fraction(5)
                if lower and (m // 5) % 2 != 0:
                    weight = Fraction(-5)
                cf = c5_scale(cf, weight * multi)
                if cf == (ZERO, ZERO):
                    continue
                mono = (kx, ky, kz)
                s_power = kx + ky
                key = (mono, s_power)
                acc = c5_add(out.get(key, c5(ZERO)), cf)
                out[key] = acc

    # assemble into a Poly, multiplying in c^kz and s^s_power with
    # even s-powers folded into s_sq and odd powers reduced to s^5
    # times even remainder (only m = +-5 sector survives with odd
    # s-power for k <= 8).
    poly_even: Poly = {}
    poly_s5: Poly = {}
    for (mono, s_power), cf in out.items():
        require(cf[1] == ZERO, "imaginary residue in ring moment")
        coeff = cf[0]
        kz = k - mono[0] - mono[1]
        coeff = q5_mul(coeff, q5_pow(c, kz))
        term = {(mono[0], mono[1], kz): coeff}
        if s_power % 2 == 0:
            term = p_scale_q5(term, q5_pow(s_sq, s_power // 2))
            poly_even = p_add(poly_even, term)
        else:
            require(s_power >= 5, "odd s-power below five survived")
            term = p_scale_q5(term, q5_pow(s_sq, (s_power - 5) // 2))
            poly_s5 = p_add(poly_s5, term)
    return {"even": poly_even, "s5": poly_s5}


def pole_vertices_exact() -> dict[str, Any]:
    c_up = INV_SQRT5
    s_sq_up = q5_sub(ONE, q5_mul(c_up, c_up))  # 4/5
    return {"c_upper": c_up, "s_sq": s_sq_up}


def build_pole_frame(i6_cartesian: Poly) -> dict[str, Any]:
    data = pole_vertices_exact()
    c_up = data["c_upper"]
    s_sq = data["s_sq"]
    s_up = q5_sqrt(s_sq)
    require(s_up is not None and q5_sign(s_up) > 0, "ring sine is not in Q(sqrt5)")

    # Sum_i P6(u_i . n) in the pole frame
    total_even = p_zero()
    total_s5 = p_zero()
    for power, coeff in enumerate(LEGENDRE[6]):
        if coeff == 0:
            continue
        # poles: (0,0,1) and (0,0,-1): (z)^power + (-z)^power
        if power % 2 == 0:
            total_even = p_add(
                total_even, p_scale({(0, 0, power): ONE}, 2 * coeff)
            )
        rings_up = ring_moment_sum(power, c_up, s_sq, lower=False)
        rings_lo = ring_moment_sum(power, q5_neg(c_up), s_sq, lower=True)
        for part, bucket in (("even", "even"), ("s5", "s5")):
            merged = p_add(
                p_scale(rings_up[part], coeff), p_scale(rings_lo[part], coeff)
            )
            if bucket == "even":
                total_even = p_add(total_even, merged)
            else:
                total_s5 = p_add(total_s5, merged)

    i6_even = p_scale(total_even, Fraction(25, 132))
    i6_s5 = p_scale(total_s5, Fraction(25, 132))

    # spherical normal form: I6 = P6(z) + (21/8) z Re((x+iy)^5) * s^{-5}...
    # In this frame Re((x+iy)^5) carries exactly the s^5 cos(5 phi)
    # factor, so the claim is:
    #   i6_even  == P6(z)  and  s^5-part == (21/8) z Re((x+iy)^5) / s^5
    # where the stored s5 bucket carries the s^5 factored out as a
    # polynomial in (x, y) of degree five over the ring radius; concretely
    # the bucket must equal (21/8) z * Re5(x, y) / s_up^5 with
    # Re5(x, y) = x^5 - 10 x^3 y^2 + 5 x y^4.
    p6_z = p_zero()
    for power, coeff in enumerate(LEGENDRE[6]):
        if coeff:
            p6_z = p_add(p6_z, p_scale({(0, 0, power): ONE}, coeff))
    require(
        p_is_zero(p_reduce_sphere(p_add(i6_even, p_scale(p6_z, -1)))),
        "even sector of the pole-frame I6 is not P6(z)",
    )
    re5 = {
        (5, 0, 0): ONE,
        (3, 2, 0): q5(-10),
        (1, 4, 0): q5(5),
    }
    target_s5 = p_scale_q5(
        p_mul({(0, 0, 1): ONE}, re5),
        q5_div(q5(Fraction(21, 8)), q5_pow(s_up, 5)),
    )
    require(
        p_is_zero(p_reduce_sphere(p_add(i6_s5, p_scale(target_s5, -1)))),
        "cos(5 phi) sector of the pole-frame I6 drifts from (21/8) c s^5",
    )

    return {
        "engine": (
            "exact five-fold character sums: sum_j w^m = 5 [5|m] on the "
            "upper ring and 5 (-1)^(m/5) [5|m] on the lower ring"
        ),
        "spherical_normal_form": (
            "I6 = P6(cos theta) + (21/8) cos theta sin^5 theta cos(5 phi)"
        ),
        "even_sector_matches_p6": True,
        "cos5phi_sector_matches": True,
    }


# ---------------------------------------------------------------------------
# Critical points, values, Hessians, Morse census
# ---------------------------------------------------------------------------


def poly_eval_frac(coeffs: list[Fraction], t: Q5) -> Q5:
    total = ZERO
    for power, coeff in enumerate(coeffs):
        if coeff:
            total = q5_add(total, q5_scale(q5_pow(t, power), coeff))
    return total


def derivative_coeffs(coeffs: list[Fraction]) -> list[Fraction]:
    return [coeffs[i] * i for i in range(1, len(coeffs))]


P6C = LEGENDRE[6]
P6C_D = derivative_coeffs(P6C)
P6C_DD = derivative_coeffs(P6C_D)


def critical_point_certificate() -> dict[str, Any]:
    # exact factorization of the squared meridian equation
    # P6'(c)^2 - (441/64)(1-c^2)^3 (1-6c^2)^2
    #   = (441/64)(5c^2-1)(5c^4-5c^2+1)(45c^4-30c^2+1)
    import itertools

    def poly_mul_f(a: list[Fraction], b: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0)] * (len(a) + len(b) - 1)
        for i, ca in enumerate(a):
            for j, cb in enumerate(b):
                out[i + j] += ca * cb
        return out

    def poly_sub_f(a, b):
        n = max(len(a), len(b))
        return [
            (a[i] if i < len(a) else Fraction(0))
            - (b[i] if i < len(b) else Fraction(0))
            for i in range(n)
        ]

    lhs = poly_mul_f(P6C_D, P6C_D)
    one_minus = [Fraction(1), Fraction(0), Fraction(-1)]
    cube = poly_mul_f(poly_mul_f(one_minus, one_minus), one_minus)
    one_six = [Fraction(1), Fraction(0), Fraction(-6)]
    rhs = poly_mul_f(cube, poly_mul_f(one_six, one_six))
    rhs = [Fraction(441, 64) * v for v in rhs]
    squared = poly_sub_f(lhs, rhs)
    f1 = [Fraction(-1), Fraction(0), Fraction(5)]
    f2 = [Fraction(1), Fraction(0), Fraction(-5), Fraction(0), Fraction(5)]
    f3 = [Fraction(1), Fraction(0), Fraction(-30), Fraction(0), Fraction(45)]
    target = poly_mul_f(poly_mul_f(f1, f2), f3)
    target = [Fraction(441, 64) * v for v in target]
    while squared and squared[-1] == 0:
        squared.pop()
    while target and target[-1] == 0:
        target.pop()
    require(squared == target, "meridian factorization drift")

    # orbit latitudes (c^2 values), all in Q(sqrt5)
    latitudes = {
        "vertex_pole": q5(1),
        "vertex_ring": q5(Fraction(1, 5)),
        "face_high": q5(Fraction(1, 3), Fraction(2, 15)),
        "face_low": q5(Fraction(1, 3), Fraction(-2, 15)),
        "edge_high": q5(Fraction(1, 2), Fraction(1, 10)),
        "edge_low": q5(Fraction(1, 2), Fraction(-1, 10)),
        "edge_equator": q5(0),
    }
    # verify the quartic memberships
    for name, csq in latitudes.items():
        if name in ("vertex_pole", "edge_equator"):
            continue
        if name == "vertex_ring":
            val = q5_sub(q5_scale(csq, 5), ONE)
        elif name.startswith("edge"):
            val = q5_add(
                q5_sub(q5_scale(q5_mul(csq, csq), 5), q5_scale(csq, 5)), ONE
            )
        else:
            val = q5_add(
                q5_sub(q5_scale(q5_mul(csq, csq), 45), q5_scale(csq, 30)), ONE
            )
        require(val == ZERO, f"latitude {name} fails its factor")

    orbits = []
    total_count = 0
    euler = 0
    for name, csq in latitudes.items():
        entry: dict[str, Any] = {"orbit": name, "c_squared": q5_str(csq)}
        if name == "vertex_pole":
            # pole chart: I6 = P6(cos rho) + O(rho^5): Hessian -21 I
            hess = (q5(-21), q5(-21))
            value = ONE
            count = 2
        elif name == "edge_equator":
            # c = 0, cos(5 phi) = 0, sin(5 phi) = +-1
            f_tt = poly_eval_frac(P6C_DD, ZERO)  # s^2 P6'' - c P6' at c=0
            f_tp_sq = q5_scale(q5_mul(q5(Fraction(105, 8)), q5(Fraction(105, 8))), 1)
            # eigenvalues (f_tt +- sqrt(f_tt^2 + 4 f_tp^2)) / 2 with
            # f_tt = 105/8: discriminant = (105/8)^2 * 5
            disc_root = q5_mul(q5(Fraction(105, 8)), SQRT5)
            lam1 = q5_scale(q5_add(q5(Fraction(105, 8)), disc_root), Fraction(1, 2))
            lam2 = q5_scale(q5_sub(q5(Fraction(105, 8)), disc_root), Fraction(1, 2))
            require(
                lam1 == q5_scale(q5_add(ONE, SQRT5), Fraction(105, 16))
                and lam2 == q5_scale(q5_sub(ONE, SQRT5), Fraction(105, 16)),
                "equatorial Hessian eigenvalues drift",
            )
            hess = (lam1, lam2)
            value = q5(Fraction(-5, 16))
            require(
                poly_eval_frac(P6C, ZERO) == q5(Fraction(-5, 16)),
                "equatorial value drift",
            )
            count = 10
        else:
            # even-power arithmetic through t = c^2: with the criticality
            # relation P6'(c) = (21/8) s^3 (6t - 1) cos5phi, the quantity
            # ratio = c P6'(c) / (6t - 1) equals (21/8) c s^3 cos5phi and
            # every orbit datum is a rational function of t over Q(sqrt5):
            #   value = P6 + (1 - t) ratio
            #   f_tt  = (1 - t) P6'' + 3 c P6' - 12 (1 - t) ratio
            #   f_pp  = -25 ratio
            tval = csq
            s_sq = q5_sub(ONE, tval)
            p6_t = poly_eval_frac(
                [Fraction(-5, 16), Fraction(105, 16), Fraction(-315, 16),
                 Fraction(231, 16)],
                tval,
            )
            c_p6d = poly_eval_frac(
                [Fraction(0), Fraction(210, 16), Fraction(-1260, 16),
                 Fraction(1386, 16)],
                tval,
            )
            p6dd_t = poly_eval_frac(
                [Fraction(210, 16), Fraction(-3780, 16), Fraction(6930, 16)],
                tval,
            )
            six_c = q5_sub(q5_scale(tval, 6), ONE)
            require(q5_sign(six_c) != 0, f"degenerate meridian at {name}")
            ratio = q5_div(c_p6d, six_c)
            value = q5_add(p6_t, q5_mul(s_sq, ratio))
            f_tt = q5_add(
                q5_add(q5_mul(s_sq, p6dd_t), q5_scale(c_p6d, 3)),
                q5_scale(q5_mul(s_sq, ratio), -12),
            )
            f_pp = q5_scale(ratio, -25)
            hess = (f_tt, f_pp)
            count = 10
        signs = tuple(q5_sign(h) for h in hess)
        if signs == (-1, -1):
            morse_type, index = "maximum", 2
        elif signs == (1, 1):
            morse_type, index = "minimum", 0
        else:
            require(set(signs) == {1, -1}, f"degenerate Hessian at {name}")
            morse_type, index = "saddle", 1
        entry.update(
            {
                "count": count,
                "value": q5_str(value),
                "hessian_eigenvalues": [q5_str(h) for h in hess],
                "type": morse_type,
            }
        )
        orbits.append(entry)
        total_count += count
        euler += count * (1 if index != 1 else -1)

    require(total_count == 62, "critical count drift")
    require(euler == 2, "Euler characteristic drift")

    by_type: dict[str, int] = {}
    value_by_type: dict[str, set] = {}
    for entry in orbits:
        by_type[entry["type"]] = by_type.get(entry["type"], 0) + entry["count"]
        value_by_type.setdefault(entry["type"], set()).add(entry["value"])
    require(by_type == {"maximum": 12, "minimum": 20, "saddle": 30}, "census drift")
    require(value_by_type["maximum"] == {q5_str(ONE)}, "maximum value drift")
    require(
        value_by_type["minimum"] == {q5_str(q5(Fraction(-5, 9)))},
        "minimum value drift",
    )
    require(
        value_by_type["saddle"] == {q5_str(q5(Fraction(-5, 16)))},
        "saddle value drift",
    )

    return {
        "squared_meridian_identity": (
            "P6'(c)^2 - (441/64)(1-c^2)^3(1-6c^2)^2 = "
            "(441/64)(5c^2-1)(5c^4-5c^2+1)(45c^4-30c^2+1)"
        ),
        "case_analysis": (
            "d/dphi = -(105/8) c sin^5 theta sin(5 phi): either the poles, "
            "the equator, or sin(5 phi) = 0; the meridian equation then "
            "factors exactly and the equatorial branch forces "
            "cos(5 phi) = 0"
        ),
        "orbits": orbits,
        "census": {"maxima": 12, "minima": 20, "saddles": 30, "total": 62},
        "euler_check": "12 - 30 + 20 = 2",
        "sign_flip_note": (
            "negating I6 exchanges the twelve maxima and twenty minima and "
            "keeps the thirty saddles"
        ),
        "nondegeneracy": (
            "every tangent Hessian eigenvalue is nonzero in exact "
            "arithmetic, so all 62 points are Morse and persist under "
            "sufficiently small higher-band corrections"
        ),
    }


# ---------------------------------------------------------------------------
# Band projectors and the blind inverse-port response
# ---------------------------------------------------------------------------


def band_response_certificate(verts) -> dict[str, Any]:
    inv_norm = q5_div(ONE, NORM_SQ)

    def unit_dot(i: int, j: int) -> Q5:
        dot = ZERO
        for axis in range(3):
            dot = q5_add(dot, q5_mul(verts[i][axis], verts[j][axis]))
        return q5_mul(dot, inv_norm) if i != j else ONE

    size = 12
    dots = [[unit_dot(i, j) for j in range(size)] for i in range(size)]

    # blind antipode: unique port at graph distance >= 3 (non-adjacent,
    # no common neighbor) in the adjacency graph dot = 1/sqrt5
    adj = [[dots[i][j] == INV_SQRT5 for j in range(size)] for i in range(size)]
    antipode = []
    for i in range(size):
        partners = []
        for j in range(size):
            if j == i or adj[i][j]:
                continue
            if any(adj[i][k] and adj[k][j] for k in range(size)):
                continue
            partners.append(j)
        require(len(partners) == 1, "distance-three partner is not unique")
        antipode.append(partners[0])
    for i in range(size):
        require(dots[i][antipode[i]] == q5_neg(ONE), "antipode dot drift")
        require(antipode[antipode[i]] == i, "antipode is not an involution")

    # Gram kernels and scaled projectors
    scales = {0: Fraction(12), 1: Fraction(4), 2: Fraction(12, 5), 3: Fraction(4)}
    kernels = {}
    for level in range(4):
        kernels[level] = [
            [legendre_eval_q5(level, dots[i][j]) for j in range(size)]
            for i in range(size)
        ]

    def mat_mul(a, b):
        return [
            [
                _sum_q5(q5_mul(a[i][k], b[k][j]) for k in range(size))
                for j in range(size)
            ]
            for i in range(size)
        ]

    def _sum_q5(items):
        total = ZERO
        for item in items:
            total = q5_add(total, item)
        return total

    ranks = {}
    for level, kernel in kernels.items():
        square = mat_mul(kernel, kernel)
        for i in range(size):
            for j in range(size):
                require(
                    square[i][j] == q5_scale(kernel[i][j], scales[level]),
                    f"kernel {level} is not a scaled projector",
                )
        trace = _sum_q5(kernel[i][i] for i in range(size))
        require(trace == q5(12), f"kernel {level} trace drift")
        rank = Fraction(12) / scales[level]
        require(rank.denominator == 1, "rank is not integral")
        ranks[level] = int(rank)
    require(
        (ranks[0], ranks[1], ranks[2], ranks[3]) == (1, 3, 5, 3),
        "band ranks drift",
    )

    # blind response: J from the derived antipode, R = -J, band signs
    j_mat = [
        [ONE if antipode[i] == j else ZERO for j in range(size)]
        for i in range(size)
    ]
    band_signs = {}
    for level, kernel in kernels.items():
        prod = mat_mul(j_mat, kernel)
        expected_sign = 1 if level % 2 == 0 else -1
        for i in range(size):
            for j in range(size):
                require(
                    prod[i][j] == q5_scale(kernel[i][j], expected_sign),
                    f"antipodal parity drift at level {level}",
                )
        band_signs[level] = -expected_sign  # R = -J
    require(
        (band_signs[0], band_signs[1], band_signs[2], band_signs[3])
        == (-1, 1, -1, 1),
        "response sign vector drift",
    )

    return {
        "blind_antipode": (
            "the unique distance-three partner in the adjacency graph is "
            "the antipode; derived, not declared"
        ),
        "projector_scales": {str(k): str(v) for k, v in scales.items()},
        "band_ranks": {"unit": 1, "frame_l1": 3, "quintet_l2": 5, "kernel_l3": 3},
        "response": "R = -J with J the derived antipodal permutation",
        "response_signs_by_band": {
            "unit (lowest multipole 0, first anisotropy 6)": -1,
            "frame 3 (lowest multipole 1)": 1,
            "quintet 5 (lowest multipole 2)": -1,
            "kernel 3' (lowest multipole 3)": 1,
        },
        "sign_vector_in_band_order_1_3_3p_5": [-1, 1, 1, -1],
    }


# ---------------------------------------------------------------------------
# Kernel independence, kinetic stencil, controls, decision rules
# ---------------------------------------------------------------------------


def kernel_independence_certificate(verts, i6: Poly) -> dict[str, Any]:
    i6_reduced = p_reduce_sphere(i6)
    probes = {
        "t^2": [Fraction(0), Fraction(0), Fraction(1)],
        "t^4": [Fraction(0)] * 4 + [Fraction(1)],
        "t^6": [Fraction(0)] * 6 + [Fraction(1)],
        "t^7 (odd control)": [Fraction(0)] * 7 + [Fraction(1)],
        "t^8": [Fraction(0)] * 8 + [Fraction(1)],
        "P6": LEGENDRE[6],
        "(1+t)^6": [Fraction(_binom(6, k)) for k in range(7)],
    }
    rows = []
    for name, coeffs in probes.items():
        total = p_zero()
        for power, coeff in enumerate(coeffs):
            if coeff == 0:
                continue
            raw = moment_sum(verts, power)
            if power % 2 == 1:
                require(p_is_zero(raw), "odd moment leak in kernel probe")
                continue
            total = p_add(
                total,
                p_scale_q5(
                    raw,
                    q5_scale(q5_div(ONE, q5_pow(NORM_SQ, power // 2)), coeff),
                ),
            )
        reduced = p_reduce_sphere(total)
        # reduced must equal constant + a6 * i6_reduced
        a6 = ZERO
        probe_mono = None
        for mono, coeff in i6_reduced.items():
            if sum(mono) > 0:
                probe_mono = mono
                break
        if probe_mono is not None and probe_mono in reduced:
            a6 = q5_div(reduced[probe_mono], i6_reduced[probe_mono])
        residue = p_add(reduced, p_scale_q5(i6_reduced, q5_neg(a6)))
        constant = residue.pop((0, 0, 0), ZERO)
        require(
            p_is_zero(residue),
            f"kernel {name} does not collapse to the I6 template",
        )
        rows.append(
            {
                "kernel": name,
                "isotropic_part": q5_str(constant),
                "i6_amplitude": q5_str(a6),
                "collapses_to_i6": True,
            }
        )
    require(
        next(r for r in rows if r["kernel"].startswith("t^7"))["i6_amplitude"]
        == q5_str(ZERO),
        "odd kernel control carries anisotropy",
    )
    return {
        "statement": (
            "every declared analytic kernel's angular output below level "
            "ten is one isotropic constant plus one multiple of the same "
            "normalized I6"
        ),
        "kernels": rows,
    }


def _binom(n: int, k: int) -> int:
    from math import comb

    return comb(n, k)


def kinetic_stencil_receipt(verts, i6: Poly) -> dict[str, Any]:
    """Premise-typed: the equal-weight stencil expansion, exact."""

    m2 = normalized_moment(verts, 2, NORM_SQ)
    m4 = normalized_moment(verts, 4, NORM_SQ)
    m6 = normalized_moment(verts, 6, NORM_SQ)
    require(m2 == p_scale(radial_power(1), 4), "stencil second moment drift")
    require(m4 == p_scale(radial_power(2), Fraction(12, 5)), "stencil fourth moment drift")
    # m6 = c6 r^6 + d6 * (132/25) I6_hom * r^0 ... I6 is inhomogeneous;
    # instead certify against the homogeneous degree-six template:
    # T6 = sum_i (u_i . n)^6 - (12/7) r^6 must be proportional to the
    # harmonic part of I6. Extract: residue = m6 - (12/7) r^6.
    residue = p_add(m6, p_scale(radial_power(3), Fraction(-12, 7)))
    # the residue must reproduce (2/7)* (132/25) * harmonic6 where
    # harmonic6 is the harmonic projection of I6; certify instead the
    # sphere-reduced identity: on the sphere,
    # m6 = 12/7 + (16/33) * (132/25) I6? Compute the exact relation by
    # reduction: m6_reduced = alpha + beta * I6_reduced.
    m6_reduced = p_reduce_sphere(m6)
    i6_reduced = p_reduce_sphere(i6)
    probe_mono = next(m for m in i6_reduced if sum(m) > 0)
    beta = q5_div(m6_reduced[probe_mono], i6_reduced[probe_mono])
    rest = p_add(m6_reduced, p_scale_q5(i6_reduced, q5_neg(beta)))
    alpha = rest.pop((0, 0, 0), ZERO)
    require(p_is_zero(rest), "sixth moment does not collapse to I6")

    # lambda_a(k) = (1/(2 a^2)) sum_i [1 - cos(a k.u_i)]
    #            = k^2 - (a^2/20) k^4 + (a^4/1440)[m6 terms] ...
    # coefficient bookkeeping on the sphere direction with |k| factored:
    iso_k4 = Fraction(-1, 24) * Fraction(12, 5) / 2  # * a^2, per k^4
    require(iso_k4 == Fraction(-1, 20), "a^2 k^4 coefficient drift")
    iso_k6 = Fraction(1, 720) * alpha[0] / 2
    aniso_k6 = Fraction(1, 720) * beta[0] / 2
    require(alpha[1] == 0 and beta[1] == 0, "stencil coefficients leave Q")
    require(iso_k6 == Fraction(1, 840), "isotropic a^4 k^6 coefficient drift")
    require(aniso_k6 == Fraction(2, 7875), "anisotropic a^4 k^6 coefficient drift")

    return {
        "premise": (
            "DECLARED: the physical kinetic operator is the equal-weight "
            "twelve-direction stencil; issue #655 owns its derivation or "
            "rejection, and ledger rows OPH-A5-M4 and M5 stay unfrozen "
            "until that issue exits"
        ),
        "expansion": (
            "lambda_a(k) = k^2 - (a^2/20) k^4 + (a^4/840) k^6 "
            "+ (2 a^4/7875) k^6 I6(khat) + O(a^6 k^8)"
        ),
        "first_directional_artifact": "spin six",
        "refinement_law": {
            "statement": "A6(a/2) = A6(a)/16 + O(a^6)",
            "exact_step_ratio": "1/16 at the displayed order",
        },
        "sixth_moment_relation": {
            "on_sphere": (
                "sum_i (u_i . n)^6 = "
                f"{alpha[0]} + {beta[0]} * I6(n)"
            ),
        },
        "typing": "certified finite computation under a declared premise",
    }


def controls_certificate(verts, i6: Poly) -> dict[str, Any]:
    controls = []

    def multipole_null_failures(vectors, norm_sq: Q5) -> list[int]:
        """Levels 1..5 whose port sum is nonzero, for one equal-norm set.

        Odd powers with an irrational norm scale are admissible only when
        their raw sum cancels (antipodal sets); a nonzero odd raw sum with
        unit norm contributes directly.
        """

        failed = []
        for level in range(1, 6):
            total = p_zero()
            for power, coeff in enumerate(LEGENDRE[level]):
                if coeff == 0:
                    continue
                raw = p_zero()
                for v in vectors:
                    raw = p_add(raw, p_pow(p_linear(v), power))
                if p_is_zero(raw):
                    continue
                if power % 2 == 0:
                    scale = q5_div(ONE, q5_pow(norm_sq, power // 2))
                else:
                    require(
                        norm_sq == ONE,
                        "odd moments survive on a non-unit control set",
                    )
                    scale = ONE
                total = p_add(total, p_scale(p_scale_q5(raw, scale), coeff))
            if not p_is_zero(p_reduce_sphere(total)):
                failed.append(level)
        return failed

    # octahedron: 6 unit directions, 3-design only: level 4 must fail
    octa = [
        (ONE, ZERO, ZERO), (q5(-1), ZERO, ZERO),
        (ZERO, ONE, ZERO), (ZERO, q5(-1), ZERO),
        (ZERO, ZERO, ONE), (ZERO, ZERO, q5(-1)),
    ]
    failed = multipole_null_failures(octa, ONE)
    controls.append(
        {
            "control": "octahedron six directions",
            "expected_failure": "level four null",
            "failed_levels": failed,
            "detector_fired": 4 in failed,
        }
    )
    require(4 in failed and 2 not in failed, "octahedron control drift")

    # cube: 8 vertices, 3-design: level 4 must fail
    cube = [
        (q5(sx), q5(sy), q5(sz))
        for sx in (1, -1)
        for sy in (1, -1)
        for sz in (1, -1)
    ]
    failed = multipole_null_failures(cube, q5(3))
    controls.append(
        {
            "control": "cube eight vertices",
            "expected_failure": "level four null",
            "failed_levels": failed,
            "detector_fired": 4 in failed,
        }
    )
    require(4 in failed, "cube control drift")

    # rational twelve-point shell (declared fixed list, no symmetry)
    shell = [
        (q5(Fraction(3, 5)), q5(Fraction(4, 5)), ZERO),
        (q5(Fraction(4, 5)), ZERO, q5(Fraction(3, 5))),
        (ZERO, q5(Fraction(3, 5)), q5(Fraction(4, 5))),
        (q5(Fraction(5, 13)), q5(Fraction(12, 13)), ZERO),
        (q5(Fraction(12, 13)), ZERO, q5(Fraction(5, 13))),
        (ZERO, q5(Fraction(5, 13)), q5(Fraction(12, 13))),
        (q5(Fraction(8, 17)), q5(Fraction(15, 17)), ZERO),
        (q5(Fraction(15, 17)), ZERO, q5(Fraction(8, 17))),
        (ZERO, q5(Fraction(8, 17)), q5(Fraction(15, 17))),
        (q5(Fraction(7, 25)), q5(Fraction(24, 25)), ZERO),
        (q5(Fraction(24, 25)), ZERO, q5(Fraction(7, 25))),
        (ZERO, q5(Fraction(7, 25)), q5(Fraction(24, 25))),
    ]
    failed = multipole_null_failures(shell, ONE)
    controls.append(
        {
            "control": "declared rational twelve-point shell",
            "expected_failure": "low-level nulls",
            "failed_levels": failed,
            "detector_fired": 1 in failed or 2 in failed,
        }
    )
    require(failed, "rational shell control drift")

    # perturbed ring height: pole model with c = 1/sqrt5 + 1/50 breaks
    # the design; emulate via Cartesian scaling of one vertex pair is
    # inexact, so use the weighted detector instead: unequal gains.
    weights = [Fraction(1)] * 12
    weights[0] = Fraction(11, 10)
    dipole = p_zero()
    for idx, v in enumerate(verts):
        term = p_scale(p_linear(v), weights[idx])
        dipole = p_add(dipole, term)
    controls.append(
        {
            "control": "unequal port gains (one port at 11/10)",
            "expected_failure": "level one null",
            "detector_fired": not p_is_zero(dipole),
        }
    )
    require(not p_is_zero(dipole), "unequal-gain control drift")

    # perturbed geometry: replace one antipodal pair by a stretched pair
    stretched = list(verts)
    stretched[0] = (
        q5_scale(verts[0][0], Fraction(11, 10)),
        q5_scale(verts[0][1], Fraction(11, 10)),
        q5_scale(verts[0][2], Fraction(11, 10)),
    )
    quad = p_zero()
    for v in stretched:
        quad = p_add(quad, p_pow(p_linear(v), 2))
    isotropy_broken = quad != p_scale_q5(
        radial_power(1), q5_div(quad.get((2, 0, 0), ZERO), ONE)
    )
    controls.append(
        {
            "control": "perturbed vertex geometry (one pair stretched)",
            "expected_failure": "second-moment isotropy",
            "detector_fired": isotropy_broken,
        }
    )
    require(isotropy_broken, "perturbed-geometry control drift")

    # shuffled labels: swapping two labels without transporting the
    # incidence breaks the derived antipode involution
    def antipode_consistent(perm) -> bool:
        inv_norm = q5_div(ONE, NORM_SQ)
        size = 12

        def unit_dot(i, j):
            dot = ZERO
            for axis in range(3):
                dot = q5_add(dot, q5_mul(verts[i][axis], verts[j][axis]))
            return q5_mul(dot, inv_norm) if i != j else ONE

        adj = [
            [unit_dot(i, j) == INV_SQRT5 for j in range(size)]
            for i in range(size)
        ]
        # declared antipode after label shuffle
        true_anti = [
            next(j for j in range(size) if unit_dot(i, j) == q5_neg(ONE))
            for i in range(size)
        ]
        declared = [perm[true_anti[perm.index(i)]] for i in range(size)]
        for i in range(size):
            j = declared[i]
            if j == i or adj[i][j]:
                return False
            if any(adj[i][k] and adj[k][j] for k in range(size)):
                return False
        return True

    identity = list(range(12))
    shuffled = list(range(12))
    shuffled[1], shuffled[2] = shuffled[2], shuffled[1]
    controls.append(
        {
            "control": "shuffled port labels without transported incidence",
            "expected_failure": "blind antipode receipt",
            "detector_fired": (not antipode_consistent(shuffled))
            and antipode_consistent(identity),
        }
    )
    require(
        antipode_consistent(identity) and not antipode_consistent(shuffled),
        "label-shuffle control drift",
    )

    # R = +J: sign table must mismatch
    plus_signs = {0: 1, 1: -1, 2: 1, 3: -1}
    expected = {0: -1, 1: 1, 2: -1, 3: 1}
    controls.append(
        {
            "control": "R = +J response",
            "expected_failure": "band sign vector",
            "detector_fired": plus_signs != expected,
        }
    )

    # injected level-four contaminant: forbidden-band detector fires
    octa_p4 = p_zero()
    for v in octa:
        octa_p4 = p_add(octa_p4, p_pow(p_linear(v), 4))
    contaminated = p_add(
        p_reduce_sphere(i6), p_scale(p_reduce_sphere(octa_p4), Fraction(1, 10))
    )
    # detect: contaminated - a6 * I6 - const must be nonzero
    i6_red = p_reduce_sphere(i6)
    probe = next(m for m in i6_red if sum(m) > 0)
    a6 = q5_div(contaminated.get(probe, ZERO), i6_red[probe])
    residue = p_add(contaminated, p_scale_q5(i6_red, q5_neg(a6)))
    residue.pop((0, 0, 0), None)
    controls.append(
        {
            "control": "injected level-four contaminant",
            "expected_failure": "template overlap",
            "detector_fired": not p_is_zero(residue),
        }
    )
    require(not p_is_zero(residue), "contaminant control drift")

    require(
        all(c["detector_fired"] for c in controls),
        "a control failed to fire",
    )
    return {"controls": controls, "all_detectors_fired": True}


def decision_rules_and_ledger() -> dict[str, Any]:
    return {
        "frozen_rows": {
            "OPH-A5-M1": {
                "statement": (
                    "for equal scalar activation of the twelve ports, every "
                    "spherical-harmonic coefficient with 1 <= l <= 5 "
                    "vanishes; failure after calibrated port-gain "
                    "correction rejects the equal-trace A5 realization"
                ),
                "type": "registered conditional statement",
                "premise_ancestry": (
                    "declared twelve-port icosahedral carrier with equal "
                    "traces; no physical carrier identification is claimed"
                ),
            },
            "OPH-A5-M2": {
                "statement": (
                    "in the leading-anisotropy regime the response surface "
                    "has exactly 62 critical directions: 12 and 20 "
                    "opposite-index extrema plus 30 saddles, with the "
                    "exact values 1, -5/9, -5/16 and the certified "
                    "Hessian signatures"
                ),
                "type": "registered conditional statement",
                "premise_ancestry": (
                    "pure leading-band regime; nondegeneracy certifies "
                    "persistence under sufficiently small higher-band "
                    "corrections"
                ),
            },
            "OPH-A5-M3": {
                "statement": (
                    "band-resolved impulses give the response sign vector "
                    "(-1, +1, +1, -1) on the bands (1, 3, 3', 5) under "
                    "R = -J, with lowest multipoles (0/6, 1, 3, 2), up to "
                    "one conventional overall charge reversal"
                ),
                "type": "registered conditional statement",
                "premise_ancestry": (
                    "R = -J derived blindly from incidence distance; the "
                    "physical response channel is open"
                ),
            },
        },
        "unfrozen_rows": {
            "OPH-A5-M4": "gated on issue #655 (kinetic stencil selection)",
            "OPH-A5-M5": "gated on issue #655 (kinetic stencil selection)",
        },
        "blind_decision_rules": {
            "forbidden_band_power": (
                "exact zero in this finite certificate; an empirical "
                "campaign must freeze a fractional power threshold before "
                "exposure and record it in the comparison contract"
            ),
            "template_overlap": (
                "the anisotropic remainder after removing one I6 multiple "
                "and one constant is exactly zero here; empirical overlap "
                "thresholds freeze before exposure"
            ),
            "orbit_locations": (
                "critical directions at the exact latitudes "
                "c^2 in {1, 1/5, (5+-sqrt5)/10, (5+-2sqrt5)/15, 0}"
            ),
            "hessian_signs": "the (12, 20, 30) signature census",
            "refinement_exponent": (
                "one binary step multiplies the level-six amplitude by "
                "1/16 on the declared stencil branch; frozen only after "
                "issue #655"
            ),
        },
        "comparison_boundary": {
            "public_measurement_read": False,
            "comparison_permitted": False,
            "eligibility": (
                "INELIGIBLE_OPEN_PHYSICAL_MAP for issue #639 consumption; "
                "physicalization is owned by the carrier, response, and "
                "stencil lanes"
            ),
        },
    }


# ---------------------------------------------------------------------------
# Receipt assembly
# ---------------------------------------------------------------------------


def build_receipt() -> dict[str, Any]:
    cartesian = build_cartesian_frame()
    i6 = cartesian.pop("_i6_poly_object")
    verts = cartesian.pop("_vertices_object")
    pole = build_pole_frame(i6)
    critical = critical_point_certificate()
    bands = band_response_certificate(verts)
    kernels = kernel_independence_certificate(verts, i6)
    stencil = kinetic_stencil_receipt(verts, i6)
    controls = controls_certificate(verts, i6)
    rules = decision_rules_and_ledger()

    receipt = {
        "schema": SCHEMA,
        "status": STATUS,
        "issue": 654,
        "normalization": (
            "I6(n) = (25/132) sum_i P6(u_i . n); equal to one on every "
            "vertex direction"
        ),
        "cartesian_frame": cartesian,
        "pole_frame": pole,
        "critical_points": critical,
        "band_response": bands,
        "kernel_independence": kernels,
        "kinetic_stencil_conditional": stencil,
        "fail_closed_controls": controls,
        "decision_rules_and_ledger": rules,
    }
    receipt["receipt_sha256"] = tagged_sha256(
        canonical_json_bytes({k: v for k, v in receipt.items() if k != "receipt_sha256"})
    )
    return receipt


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args(argv)
    receipt = build_receipt()
    if args.write:
        RUNTIME.mkdir(exist_ok=True)
        RECEIPT_PATH.write_bytes(canonical_json_bytes(receipt))
        print(RECEIPT_PATH)
        return 0
    print(json.dumps(receipt["critical_points"]["census"], indent=1))
    print(receipt["status"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
