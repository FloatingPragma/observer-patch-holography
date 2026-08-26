#!/usr/bin/env python3
"""Exact global classification of the quartic-truncated W5 entropy packet.

The twelve-port quintet consists exactly of six antipodal-pair values
``a_i``, each repeated twice, with ``sum(a_i)=0`` and
``2*sum(a_i**2)=1``.  After removing a positive common factor, the
quartic-truncated entropy objective is

    G_r(a) = -S3(a) + r*S4(a)/2,
    S_k(a) = 2*sum_i a_i**k,

on the strict record domain ``1+r*a_i>0``.  Exact SymPy identities prove:

* the stationarity equation is cubic;
* every stationary point with three distinct coordinate values is a saddle;
* every possible interior minimum therefore has two coordinate values, whose
  five multiplicity classes can be compared exactly;
* with ``r_c=(32*sqrt(15)-20*sqrt(6))/27``, the strict-domain minimizer is
  the multiplicity-one C5 orbit for ``0<r<r_c``, both the multiplicity-one
  and multiplicity-two orbits at ``r_c``, and the multiplicity-two golden
  orbit for ``r_c<r<sqrt(24)``;
* for ``sqrt(24)<=r<sqrt(60)`` the strict-domain infimum is not attained.

The compact closed simplex is also classified exactly after setting
``b_i=1+r*a_i``.  For ``sqrt(24)<r<sqrt(60)`` its minimizer has four zero
weights and two nonzero weights
``3 +/- sqrt(r**2-24)/2``.  Its quadrupole spectrum is proportional to
``(-2,1-d,1+d)``, ``d=sqrt((r**2-15)/5)``, and hence supplies a continuous
simple-spectrum gap-ratio branch.  This is a viable boundary/regularization
route, not a physical prediction: the entropy Taylor expansion is least
controlled at ``q_i=0``, the source does not emit ``r``, and the
quadrupole-to-physical-log-mass attachment remains open.

Measured charged-lepton values occur only in the compare-only block.  No
finite-seed search or numerical optimizer is used as proof, no physical
promotion is allowed, and ``quartic_packet_globally_excluded`` is always
false.
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any, Sequence

import sympy as sp

ROOT = Path(__file__).resolve().parents[2]
OUT_PATH = ROOT / "particles" / "runs" / "flavor" / "entropy_w5_shape_certificate.json"

SCHEMA = "oph.entropy_w5_shape_certificate.v3"
ISSUE_CONTEXT = [546, 591]


class CertificateError(ValueError):
    def __init__(self, code: str, message: str) -> None:
        super().__init__(f"{code}: {message}")
        self.code = code


def require(condition: bool, code: str, message: str) -> None:
    if not condition:
        raise CertificateError(code, message)


class Q5:
    """Exact a + b*sqrt5 with Fraction coefficients."""

    __slots__ = ("a", "b")

    def __init__(self, a: Any = 0, b: Any = 0) -> None:
        self.a = Fraction(a)
        self.b = Fraction(b)

    def __add__(self, o: "Q5") -> "Q5":
        return Q5(self.a + o.a, self.b + o.b)

    def __sub__(self, o: "Q5") -> "Q5":
        return Q5(self.a - o.a, self.b - o.b)

    def __neg__(self) -> "Q5":
        return Q5(-self.a, -self.b)

    def __mul__(self, o: "Q5") -> "Q5":
        return Q5(self.a * o.a + 5 * self.b * o.b, self.a * o.b + self.b * o.a)

    def inv(self) -> "Q5":
        n = self.a * self.a - 5 * self.b * self.b
        require(n != 0, "Q5_DIVZERO", "inverse of zero")
        return Q5(self.a / n, -self.b / n)

    def __truediv__(self, o: "Q5") -> "Q5":
        return self * o.inv()

    def __eq__(self, o: object) -> bool:
        return isinstance(o, Q5) and self.a == o.a and self.b == o.b

    def is_zero(self) -> bool:
        return self.a == 0 and self.b == 0

    def is_positive(self) -> bool:
        if self.b == 0:
            return self.a > 0
        if self.a == 0:
            return self.b > 0
        if self.a > 0 and self.b > 0:
            return True
        if self.a < 0 and self.b < 0:
            return False
        if self.a > 0:
            return self.a * self.a > 5 * self.b * self.b
        return 5 * self.b * self.b > self.a * self.a

    def to_float(self) -> float:
        return float(self.a) + float(self.b) * math.sqrt(5.0)

    def render(self) -> str:
        return f"{self.a} + {self.b}*sqrt5"


ZERO = Q5(0)
ONE = Q5(1)
PHI = Q5(Fraction(1, 2), Fraction(1, 2))
PHI2 = PHI * PHI  # phi^2 = phi + 1
PHI4 = PHI2 * PHI2


def q5_lt(x: Q5, y: Q5) -> bool:
    return (y - x).is_positive()


# ---------------------------------------------------------------------------
# The diagonal (D2) stratum in exact arithmetic
# ---------------------------------------------------------------------------
#
# Icosahedron vertices (0, ±1, ±phi) and cyclic permutations; the three
# coordinate axes are mutually orthogonal two-fold axes, and a diagonal
# traceless Q = diag(q1, q2, q3) is exactly the D2-fixed stratum.  Each
# vertex family contributes the squared-coordinate profile (0, 1, phi^2)
# in cyclic position, so with |v|^2 = 2 + phi the twelve port values
# collapse to three values, each on four ports:
#
#   y_a = (q2 + phi^2 q3) / (2 + phi)     [vertices (0, ±1, ±phi)]
#   y_b = (q1 + phi^2 q2) / (2 + phi)     [vertices (±1, ±phi, 0)]
#   y_c = (q3 + phi^2 q1) / (2 + phi)     [vertices (±phi, 0, ±1)]
#
# The port sums are S2 = 4(y_a^2 + y_b^2 + y_c^2), S3 = 4 sum y^3,
# S4 = 4 sum y^4, and sum y = 0 exactly when tr Q = 0.


NORM = Q5(2) + PHI  # |v|^2 = 2 + phi = (5 + sqrt5)/2


def port_profile(q: Sequence[Q5]) -> tuple[Q5, Q5, Q5]:
    q1, q2, q3 = q
    inv = NORM.inv()
    return (
        (q2 + PHI2 * q3) * inv,
        (q1 + PHI2 * q2) * inv,
        (q3 + PHI2 * q1) * inv,
    )


def port_sums(q: Sequence[Q5]) -> tuple[Q5, Q5, Q5]:
    ys = port_profile(q)
    s2 = ZERO
    s3 = ZERO
    s4 = ZERO
    four = Q5(4)
    for y in ys:
        y2 = y * y
        s2 = s2 + four * y2
        s3 = s3 + four * y2 * y
        s4 = s4 + four * y2 * y2
    return s2, s3, s4


def invert_profile(y: Sequence[Q5]) -> tuple[Q5, Q5, Q5]:
    """Solve the cyclic system for diag(q) with the given port profile."""

    y_a, y_b, y_c = (value * NORM for value in y)
    # q2 + phi^2 q3 = y_a ; q1 + phi^2 q2 = y_b ; q3 + phi^2 q1 = y_c.
    # From eq2: q1 = y_b - phi^2 q2; substituted into eq3: q3 = y_c - phi^2 q1;
    # eq1 then fixes q2:  q2 (1 + phi^6) = y_a - phi^2 y_c + phi^4 y_b.
    p2 = PHI2
    q2 = (y_a - p2 * y_c + p2 * p2 * y_b) / (ONE + p2 * p2 * p2)
    q1 = y_b - p2 * q2
    q3 = y_c - p2 * q1
    for got, want in zip(port_profile((q1, q2, q3)), y):
        require((got - want).is_zero(), "PROFILE_ROUNDTRIP", "profile inversion must reproduce the profile")
    require((q1 + q2 + q3).is_zero(), "TRACELESS", "inverted diagonal must be traceless")
    return q1, q2, q3


# ---------------------------------------------------------------------------
# The two exact branch orbits
# ---------------------------------------------------------------------------


def prolate_orbit() -> dict[str, Any]:
    """The vertex-axis (C5) prolate orbit, computed on its own stratum.

    Q = a a^T / |a|^2 - I/3 for a vertex axis a = (0, 1, phi).  The twelve
    port values split as two poles at 2/3 and ten ring ports; the exact
    values follow from (v . a)^2 / (|v|^2 |a|^2) being 1 at the poles and
    exactly 1/sqrt5-quadratic on the ring.
    """

    # (v.a)^2/(|v|^2|a|^2) for v over the twelve vertices, a = (0,1,phi):
    # poles (v = ±a): 1.  Ring: the icosahedral frame gives cos^2 = 1/sqrt5
    # weight: (v.a)^2 = ((5+sqrt5)/2)^2 / 5 for the ten ring vertices.
    # Port value x = cos^2 - 1/3.
    pole = ONE - Q5(Fraction(1, 3))
    # Ring cosine: v = (1, phi, 0), a = (0, 1, phi), v.a = phi,
    # |v|^2 = |a|^2 = 2 + phi, and (2 + phi)^2 = 5 phi^2, so
    # cos^2 = phi^2 / (2 + phi)^2 = 1/5 exactly, for all ten ring vertices.
    cos2_ring = (PHI * PHI) / (NORM * NORM)
    require(cos2_ring == Q5(Fraction(1, 5)), "PROLATE_RING", "the ring cosine must be exactly one fifth")
    ring = cos2_ring - Q5(Fraction(1, 3))
    # Verify the design identity: 2*pole + 10*ring must vanish (trace zero).
    require((Q5(2) * pole + Q5(10) * ring).is_zero(), "PROLATE_TRACE", "prolate port profile must sum to zero")
    s2 = Q5(2) * pole * pole + Q5(10) * ring * ring
    s3 = Q5(2) * pole * pole * pole + Q5(10) * ring * ring * ring
    s4 = Q5(2) * pole * pole * pole * pole + Q5(10) * ring * ring * ring * ring
    return {
        "orbit": "vertex_axis_prolate_C5",
        "eigenvalue_pattern": "(2, -1, -1)/3 about the vertex axis",
        "double_eigenvalue": True,
        "port_values": {"pole": pole.render(), "ring": ring.render(), "pole_count": 2, "ring_count": 10},
        "S2": s2,
        "S3": s3,
        "S4": s4,
        "min_port_value": ring,
    }


def golden_orbit() -> dict[str, Any]:
    """The D2-stratum cubic-extremal orbit with simple spectrum.

    On the diagonal stratum the port profile is three values with sum
    zero; the cubic extremal profile at fixed square sum is the symmetric
    pattern (2, -1, -1) up to scale and permutation.  Inverting the exact
    cyclic map produces the diagonal triple, whose centered eigenvalues
    come out proportional to (-(phi^4 - 1), -1, phi^4) with sorted-gap
    ratio exactly phi.
    """

    profile = (Q5(2), Q5(-1), Q5(-1))
    q = invert_profile(profile)
    total = q[0] + q[1] + q[2]
    require(total.is_zero(), "GOLDEN_TRACE", "the golden diagonal must be traceless")
    # Sort exactly by the Q(sqrt5) order.
    triple = sorted(q, key=lambda v: (v.to_float()))
    for left, right in zip(triple, triple[1:]):
        require(q5_lt(left, right), "GOLDEN_SIMPLE", "the golden orbit must have simple spectrum")
    gap_low = triple[1] - triple[0]
    gap_high = triple[2] - triple[1]
    ratio = gap_high / gap_low
    require(ratio == PHI, "GOLDEN_RATIO", "the sorted-gap ratio must equal phi exactly")
    scaled = [value / (-triple[1]) for value in triple]
    require(
        scaled[0] == -(PHI4 - ONE) and scaled[1] == -ONE and scaled[2] == PHI4,
        "GOLDEN_PATTERN",
        "the centered triple must be proportional to (-(phi^4-1), -1, phi^4)",
    )
    s2, s3, s4 = port_sums(q)
    min_port = min(port_profile(q), key=lambda v: v.to_float())
    require(min_port == Q5(-1), "GOLDEN_MIN_PORT", "the golden profile minimum must be -1")
    return {
        "orbit": "diagonal_D2_cubic_extremal",
        "eigenvalues": [value.render() for value in triple],
        "eigenvalue_pattern": "proportional to (-(phi^4 - 1), -1, phi^4)",
        "sorted_gap_ratio": "phi = (1 + sqrt5)/2, exactly",
        "double_eigenvalue": False,
        "symmetric_criticality": (
            "the profile (2, -1, -1) is a critical point of the cubic at fixed "
            "square sum on the three-value plane, and D2-equivariance lifts "
            "stratum criticality to criticality of the full constrained "
            "functional"
        ),
        "S2": s2,
        "S3": s3,
        "S4": s4,
        "q_diagonal": [value.render() for value in q],
        "min_port_value": min_port,
    }


def normalized_invariants(branch: dict[str, Any]) -> tuple[Q5, Q5]:
    """Return (s3n, s4n) = (S3/S2^{3/2}, S4/S2^2) squared-rationalized.

    S2^{3/2} is irrational in general; the certificate stores the exact
    squared normalization s3n2 = S3^2 / S2^3 (a Q(sqrt5) number) together
    with the exact s4n = S4 / S2^2, and every comparison below is arranged
    to use only these exact quantities.
    """

    s2 = branch["S2"]
    s3 = branch["S3"]
    s4 = branch["S4"]
    s3n2 = (s3 * s3) / (s2 * s2 * s2)
    s4n = s4 / (s2 * s2)
    return s3n2, s4n


def _legacy_candidate_crossing_display(
    prolate: dict[str, Any], golden: dict[str, Any]
) -> dict[str, Any]:
    """Legacy two-candidate display, not used by the v3 proof or payload.

    With each branch normalized to unit S2, the shape energy is
    E(r) = -s3n r^3/6 + s4n r^4/12, with s3n = sqrt(s3n2) > 0 for both
    branches.  The crossing r_c solves E_p(r) = E_g(r):
    r_c = 2 (s3n_p - s3n_g) / (s4n_p - s4n_g).  The certificate stores the
    exact squared data and the crossing as a root of the recorded equation.
    Its binary-float rendering is diagnostic only.  The v3 certificate below
    proves and serializes the crossing with exact SymPy arithmetic.
    """

    p32, p4 = normalized_invariants(prolate)
    g32, g4 = normalized_invariants(golden)
    require(
        prolate["S3"].is_positive() and golden["S3"].is_positive(),
        "BRANCH_CUBIC_SIGN",
        "both branch cubics must be strictly positive in the chosen orientation",
    )
    require(not (p4 - g4).is_zero(), "BRANCH_QUARTIC_TIE", "the quartic invariants must differ")
    # Display floats, exactness carried by the stored equation.
    s3p = math.sqrt(p32.to_float())
    s3g = math.sqrt(g32.to_float())
    r_c = 2.0 * (s3p - s3g) / (p4 - g4).to_float()
    # Positivity domain: the record q_i = (1 + r w_i)/12 with unit-S2 profile
    # w stays a probability vector exactly when r < 1/|min_i w_i|, i.e.
    # r^2 < S2 / (min_i port value)^2, an exact Q(sqrt5) number per branch.
    positivity = {}
    r2_bound_float = {}
    for name, branch in (("prolate", prolate), ("golden", golden)):
        min_port = branch["min_port_value"]
        require(not min_port.is_positive() and not min_port.is_zero(), "POSITIVITY_SIGN", "the minimal port value must be negative")
        r2_bound = branch["S2"] / (min_port * min_port)
        r2_bound_float[name] = r2_bound.to_float()
        positivity[name] = {
            "r_squared_bound": r2_bound.render(),
            "r_bound_display": f"{math.sqrt(r2_bound_float[name]):.12f}",
        }
    require(
        r_c * r_c < r2_bound_float["golden"] and r_c * r_c < r2_bound_float["prolate"],
        "CROSSING_OUTSIDE_POSITIVITY",
        "the crossing must lie inside both positivity domains",
    )
    return {
        "equation": "2*(sqrt(s3n2_prolate) - sqrt(s3n2_golden)) = r_c * (s4n_prolate - s4n_golden)",
        "s3n2_prolate": p32.render(),
        "s3n2_golden": g32.render(),
        "s4n_prolate": p4.render(),
        "s4n_golden": g4.render(),
        "r_c_display": f"{r_c:.12f}",
        "positivity_domain": {
            "condition": "r < 1/|min_i w_i| for the unit-S2 profile w, i.e. r^2 < S2/(min port value)^2",
            "branches": positivity,
        },
        "selection": (
            "on the compared strata and inside each branch's positivity "
            "domain: below r_c the prolate branch has the lower shape "
            "energy; between r_c and the golden positivity bound the golden "
            "branch does; above that bound the golden record leaves the "
            "probability simplex and this comparison emits no admissible "
            "displayed simple-spectrum branch; the v3 global theorem below "
            "supersedes this legacy two-candidate display"
        ),
        "other_critical_orbits": (
            "the C3-axis orbit is critical by symmetric criticality with "
            "S3 = 0 exactly (two-value antipodal profile), S4/S2^2 = 1/12, "
            "positivity bound r^2 < 12, and shape energy r^4/144; it "
            "undercuts the golden branch only for r^2 > 96, outside every "
            "positivity bound, so it never enters this compared set; this "
            "calculation does not classify generic, C2-fixed, or V4-fixed "
            "critical orbits"
        ),
        "display_only_note": "r_c_display and r_bound_display are float renderings of the recorded exact quantities",
    }


# ---------------------------------------------------------------------------
# Exact global quartic classification
# ---------------------------------------------------------------------------


def _sympy_text(expression: sp.Expr) -> str:
    """Stable plain-text rendering after exact simplification/factorization."""

    return sp.sstr(sp.factor(sp.simplify(expression)))


def _is_exact_zero(expression: sp.Expr) -> bool:
    return sp.simplify(expression) == 0


def _two_value_rows(r: sp.Symbol) -> tuple[list[dict[str, Any]], dict[int, sp.Expr]]:
    """Verify and serialize every two-value stationary class exactly."""

    rows: list[dict[str, Any]] = []
    objectives: dict[int, sp.Expr] = {}
    for multiplicity in range(1, 6):
        other = 6 - multiplicity
        high = sp.sqrt(sp.Rational(other, 12 * multiplicity))
        low = -sp.sqrt(sp.Rational(multiplicity, 12 * other))
        s3 = sp.simplify(2 * (multiplicity * high**3 + other * low**3))
        s4 = sp.simplify(2 * (multiplicity * high**4 + other * low**4))
        s3_expected = sp.Rational(3 - multiplicity, 1) / sp.sqrt(
            3 * multiplicity * other
        )
        s4_expected = sp.Rational(3, multiplicity * other) - sp.Rational(1, 4)
        checks = {
            "zero_sum": _is_exact_zero(multiplicity * high + other * low),
            "unit_s2": _is_exact_zero(
                2 * (multiplicity * high**2 + other * low**2) - 1
            ),
            "s3_formula": _is_exact_zero(s3 - s3_expected),
            "s4_formula": _is_exact_zero(s4 - s4_expected),
        }
        require(
            all(checks.values()),
            "TWO_VALUE_IDENTITY",
            f"multiplicity {multiplicity} failed an exact identity",
        )
        objective = sp.simplify(-s3 + r * s4 / 2)
        objectives[multiplicity] = objective
        rows.append(
            {
                "positive_root_multiplicity": multiplicity,
                "negative_root_multiplicity": other,
                "positive_root": _sympy_text(high),
                "negative_root": _sympy_text(low),
                "S3": _sympy_text(s3),
                "S4": _sympy_text(s4),
                "G_r": _sympy_text(objective),
                "strict_positivity_bound_r_squared": str(
                    sp.Rational(12 * other, multiplicity)
                ),
                "exact_checks": checks,
            }
        )
    return rows, objectives


def _three_root_saddle_certificate(r: sp.Symbol) -> dict[str, Any]:
    """Exact second-variation exclusions for all three-root stationary points."""

    alpha, beta, gamma = sp.symbols("alpha beta gamma", real=True)
    gap_a, gap_b = sp.symbols("A B", positive=True)
    outer_m, outer_p = sp.symbols("M P", integer=True, positive=True)

    # If the middle root occurs at least twice, the split direction (+1,-1)
    # within that group is tangent to both constraints.  Its quadratic form
    # is twice the derivative of the stationarity cubic at beta.
    middle_split = sp.factor(
        2 * 4 * r * (beta - alpha) * (beta - gamma)
    ).subs({alpha: beta - gap_a, gamma: beta + gap_b})
    middle_split_expected = -8 * r * gap_a * gap_b

    # If beta has multiplicity one, let M and P be the multiplicities of
    # alpha and gamma.  A group-constant tangent vector has components
    # P*B, -M*P*(A+B), M*A on the three groups.  Its exact constrained
    # second variation factors as below.
    tangent_alpha = outer_p * (gamma - beta)
    tangent_beta = -outer_m * outer_p * (gamma - alpha)
    tangent_gamma = outer_m * (beta - alpha)
    lagrange_second = {
        alpha: 4 * r * (alpha - beta) * (alpha - gamma),
        beta: 4 * r * (beta - alpha) * (beta - gamma),
        gamma: 4 * r * (gamma - alpha) * (gamma - beta),
    }
    group_quadratic = sp.factor(
        outer_m * lagrange_second[alpha] * tangent_alpha**2
        + lagrange_second[beta] * tangent_beta**2
        + outer_p * lagrange_second[gamma] * tangent_gamma**2
    ).subs({alpha: beta - gap_a, gamma: beta + gap_b})
    group_expected = sp.factor(
        4
        * r
        * outer_m
        * outer_p
        * gap_a
        * gap_b
        * (gap_a + gap_b)
        * (
            outer_m * (1 - outer_p) * gap_a
            + outer_p * (1 - outer_m) * gap_b
        )
    )
    outer_sign_rows = []
    all_brackets_negative = True
    for m_value in range(1, 5):
        p_value = 5 - m_value
        bracket = sp.expand(
            group_expected
            / (4 * r * outer_m * outer_p * gap_a * gap_b * (gap_a + gap_b))
        ).subs({outer_m: m_value, outer_p: p_value})
        polynomial = sp.Poly(bracket, gap_a, gap_b)
        coefficients = polynomial.coeffs()
        negative_for_positive_gaps = (
            all(coefficient <= 0 for coefficient in coefficients)
            and any(coefficient < 0 for coefficient in coefficients)
        )
        all_brackets_negative = all_brackets_negative and negative_for_positive_gaps
        outer_sign_rows.append(
            {
                "M": m_value,
                "P": p_value,
                "bracket": _sympy_text(bracket),
                "strictly_negative_for_A_B_positive": negative_for_positive_gaps,
            }
        )

    checks = {
        "middle_split_factor_identity": _is_exact_zero(
            middle_split - middle_split_expected
        ),
        "middle_split_strictly_negative": bool(
            middle_split_expected.is_negative
        ),
        "group_constant_factor_identity": _is_exact_zero(
            group_quadratic - group_expected
        ),
        "all_outer_multiplicity_brackets_strictly_negative": (
            all_brackets_negative
        ),
    }
    require(
        all(checks.values()),
        "THREE_ROOT_SADDLE",
        "the exact three-root saddle factorization failed",
    )
    return {
        "ordered_roots": "alpha < beta < gamma",
        "stationarity_cubic": (
            "4*r*z^3 - 6*z^2 - 2*mu*z - lambda = 0; its roots are "
            "the coordinate values of a constrained stationary point"
        ),
        "middle_multiplicity_at_least_two": {
            "tangent": "split two beta coordinates by +epsilon and -epsilon",
            "second_variation_factor": _sympy_text(middle_split_expected),
            "sign": "strictly_negative for r,A,B>0",
        },
        "middle_multiplicity_one": {
            "notation": (
                "M and P are outer-root multiplicities, M+P=5; "
                "A=beta-alpha>0 and B=gamma-beta>0"
            ),
            "group_constant_tangent": "(P*B, -M*P*(A+B), M*A)",
            "second_variation_factor": (
                "4*r*M*P*A*B*(A+B)*"
                "(M*(1-P)*A+P*(1-M)*B)"
            ),
            "outer_multiplicity_sign_checks": outer_sign_rows,
            "sign": "strictly_negative for every M=1,2,3,4 and P=5-M",
        },
        "conclusion": (
            "every three-root stationary point on the strict six-pair "
            "sphere is a saddle, so an attained global minimum has at most "
            "two coordinate values"
        ),
        "checks": checks,
    }


def _closed_simplex_certificate(
    r: sp.Symbol,
    objectives: dict[int, sp.Expr],
) -> dict[str, Any]:
    """Exact compact-simplex classification, including the boundary branch."""

    K = sp.symbols("K", real=True)
    h = lambda value: value**4 - 6 * value**3
    h2 = sp.factor(K * (K - 36) / 2)

    face_gap_a, face_gap_b = sp.symbols("A_face B_face", positive=True)
    face_three_root_sign_rows = []
    face_three_root_group_negative = True
    for support in range(4, 7):
        for outer_m in range(1, support - 1):
            outer_p = support - 1 - outer_m
            bracket = sp.expand(
                outer_m * (1 - outer_p) * face_gap_a
                + outer_p * (1 - outer_m) * face_gap_b
            )
            coefficients = sp.Poly(bracket, face_gap_a, face_gap_b).coeffs()
            negative = (
                all(coefficient <= 0 for coefficient in coefficients)
                and any(coefficient < 0 for coefficient in coefficients)
            )
            face_three_root_group_negative = (
                face_three_root_group_negative and negative
            )
            face_three_root_sign_rows.append(
                {
                    "support_size": support,
                    "M": outer_m,
                    "P": outer_p,
                    "bracket": _sympy_text(bracket),
                    "strictly_negative": negative,
                }
            )
    face_three_root_checks = {
        "middle_split_negative": bool(
            (-8 * face_gap_a * face_gap_b).is_negative
        ),
        "group_tangent_negative_for_support_4_to_6": (
            face_three_root_group_negative
        ),
        "support_three_vieta_contradiction": sp.Rational(9, 2) != 6,
    }

    # On a face with n positive coordinates and 18<K<36, feasibility forces
    # exactly one high root: if q>=2 high roots, positivity would require
    # q*K<36, impossible.  The remaining candidates have one x_n and n-1
    # copies of y_n.  Their exact excess over the two-support value H2 is
    # positive by the following factorizations.
    expected_factors = {
        3: lambda t: sp.Rational(4, 9) * (t - 6) ** 2 * (t + 3),
        4: lambda t: sp.Rational(3, 64) * (t - 6) ** 2 * (t + 6) ** 2,
        5: lambda t: sp.Rational(12, 125)
        * (t - 6) ** 2
        * (t**2 + 6 * t + 18),
        6: lambda t: sp.Rational(5, 36)
        * (t - 6) ** 2
        * (t**2 + 4 * t + 12),
    }
    face_rows = []
    face_checks = []
    for support in range(3, 7):
        t = sp.symbols(f"t_{support}", positive=True)
        k_from_t = sp.simplify((36 + (support - 1) * t**2) / support)
        high = sp.simplify((6 + (support - 1) * t) / support)
        low = sp.simplify((6 - t) / support)
        h_support = sp.expand(h(high) + (support - 1) * h(low))
        excess = sp.factor(h_support - h2.subs(K, k_from_t))
        expected = sp.factor(expected_factors[support](t))
        checks = {
            "sum_is_six": _is_exact_zero(high + (support - 1) * low - 6),
            "sum_squares_is_K": _is_exact_zero(
                high**2 + (support - 1) * low**2 - k_from_t
            ),
            "H_excess_factor_identity": _is_exact_zero(excess - expected),
            "factor_strictly_positive_for_0_t_less_than_6": bool(
                (expected / (t - 6) ** 2).is_positive
            ),
        }
        face_checks.extend(checks.values())
        face_rows.append(
            {
                "support_size": support,
                "t_definition": f"sqrt(({support}*K-36)/({support}-1))",
                "high_value_x_k": _sympy_text(high),
                "low_value_y_k": _sympy_text(low),
                "H_k_minus_H_2_factor": _sympy_text(expected),
                "domain": "18 < K < 36 implies 0 < t_k < 6",
                "exact_checks": checks,
            }
        )

    # For 6<K<18, an exact one-variable KKT enumeration proves that no
    # proper-face two-root candidate satisfies positivity, both split
    # second-order conditions, and the zero-coordinate activation condition
    # lambda<=0 simultaneously.  Uniform one-root face points are checked
    # separately at their isolated K=36/n values.
    x = sp.symbols("x", real=True)
    low_kkt_rows = []
    low_kkt_checks = []
    for support in range(2, 6):
        for low_multiplicity in range(1, support):
            high_multiplicity = support - low_multiplicity
            y = sp.simplify((6 - low_multiplicity * x) / high_multiplicity)
            k_value = sp.factor(
                low_multiplicity * x**2 + high_multiplicity * y**2
            )
            low_split = sp.factor(2 * (x - y) * (4 * x + 2 * y - 9))
            high_split = sp.factor(-2 * (x - y) * (2 * x + 4 * y - 9))
            activation_lambda = sp.factor(2 * x * y * (9 - 2 * (x + y)))
            conditions = [x > 0, y > x, k_value < 18, activation_lambda <= 0]
            if low_multiplicity >= 2:
                conditions.append(low_split >= 0)
            if high_multiplicity >= 2:
                conditions.append(high_split >= 0)
            reduced = sp.reduce_inequalities(conditions, x)
            empty = reduced == sp.false
            low_kkt_checks.append(empty)
            low_kkt_rows.append(
                {
                    "support_size": support,
                    "low_multiplicity": low_multiplicity,
                    "high_multiplicity": high_multiplicity,
                    "exact_semialgebraic_result": (
                        "empty" if empty else sp.sstr(reduced)
                    ),
                }
            )

    uniform_face_rows = []
    uniform_face_checks = []
    r_c = (32 * sp.sqrt(15) - 20 * sp.sqrt(6)) / 27
    for support in range(3, 6):
        k_value = sp.Rational(36, support)
        r_value = sp.sqrt(2 * (k_value - 6))
        uniform_h = support * h(sp.Rational(6, support))
        uniform_g = sp.simplify((uniform_h + 12 * k_value - 42) / r_value**3)
        winning_multiplicity = 1 if (r_c - r_value).is_positive else 2
        excess = sp.simplify(
            uniform_g - objectives[winning_multiplicity].subs(r, r_value)
        )
        positive = bool(excess.is_positive)
        uniform_face_checks.append(positive)
        uniform_face_rows.append(
            {
                "support_size": support,
                "r_squared": _sympy_text(r_value**2),
                "G_uniform_minus_interior_global": _sympy_text(excess),
                "strictly_positive": positive,
            }
        )

    # Independent nonnegative global identity.  With pair probabilities
    # p_i=b_i/6, H=1296*(fixed(e2)+J), J=e3-4e4.  Expanding J over triples
    # gives a sum of nonnegative monomials, with equality exactly when the
    # support contains at most two indices.
    p = sp.symbols("p0:6", nonnegative=True)
    e3 = sum(p[i] * p[j] * p[k] for i in range(6) for j in range(i + 1, 6) for k in range(j + 1, 6))
    e4 = sum(
        p[i] * p[j] * p[k] * p[l]
        for i in range(6)
        for j in range(i + 1, 6)
        for k in range(j + 1, 6)
        for l in range(k + 1, 6)
    )
    triple_sum = sum(
        p[i]
        * p[j]
        * p[k]
        * (p[i] + p[j] + p[k])
        for i in range(6)
        for j in range(i + 1, 6)
        for k in range(j + 1, 6)
    )
    # Polynomial identity before imposing sum(p)=1:
    # e1*e3-4e4 = sum_T p_T*sum_{i in T}p_i.  On the probability
    # simplex e1=1, its left side is J=e3-4e4.
    j_identity = sp.expand(sum(p) * e3 - 4 * e4 - triple_sum)
    e2_symbol, e3_symbol, e4_symbol = sp.symbols("e2 e3 e4", real=True)
    power3 = 1 - 3 * e2_symbol + 3 * e3_symbol
    power4 = (
        1
        - 4 * e2_symbol
        + 2 * e2_symbol**2
        + 4 * e3_symbol
        - 4 * e4_symbol
    )
    h_elementary_identity = sp.expand(
        power4
        - power3
        - (-e2_symbol + 2 * e2_symbol**2 + e3_symbol - 4 * e4_symbol)
    )

    b_symbol = sp.symbols("b", real=True)
    coordinate_transform_identity = sp.expand(
        (b_symbol - 1) ** 4
        - 2 * (b_symbol - 1) ** 3
        - (b_symbol**4 - 6 * b_symbol**3)
        - (12 * b_symbol**2 - 10 * b_symbol + 3)
    )

    sqrt_boundary = sp.sqrt(r**2 - 24)
    b_minus = 3 - sqrt_boundary / 2
    b_plus = 3 + sqrt_boundary / 2
    x_minus = sp.simplify(b_minus - 1)
    x_plus = sp.simplify(b_plus - 1)
    g_boundary = sp.factor((r**4 - 480) / (8 * r**3))
    e_boundary = sp.factor((r**4 - 480) / 48)
    boundary_values = [x_minus / r, x_plus / r] + [-1 / r] * 4
    s3_boundary = sp.simplify(2 * sum(value**3 for value in boundary_values))
    s4_boundary = sp.simplify(2 * sum(value**4 for value in boundary_values))
    g_boundary_direct = sp.simplify(-s3_boundary + r * s4_boundary / 2)
    g1_difference = sp.factor(g_boundary - objectives[1])
    g1_expected = -(
        (r - 2 * sp.sqrt(15)) ** 2
        * (3 * r**2 + 4 * sp.sqrt(15) * r + 60)
        / (60 * r**3)
    )

    d = sp.sqrt((r**2 - 15) / 5)
    ratio = sp.simplify(2 * d / (3 - d))
    ratio_symbol = sp.symbols("R", positive=True)
    inverse_r2 = 15 + 45 * ratio_symbol**2 / (ratio_symbol + 2) ** 2
    recovered_d = sp.simplify(3 * ratio_symbol / (ratio_symbol + 2))
    spectrum_checks = {
        "spectrum_trace_zero": _is_exact_zero(-2 + (1 - d) + (1 + d)),
        "rank_two_discriminant_is_4d_squared": _is_exact_zero(
            (b_plus - b_minus) ** 2
            + 4 * b_plus * b_minus / 5
            - 4 * d**2
        ),
        "ratio_inverse_identity": _is_exact_zero(
            inverse_r2 - (15 + 5 * recovered_d**2)
        ),
        "ratio_solves_for_d": _is_exact_zero(
            2 * recovered_d / (3 - recovered_d) - ratio_symbol
        ),
        "d_upper_domain_identity": _is_exact_zero(
            9 - d**2 - (60 - r**2) / 5
        ),
        "boundary_start_ratio_is_phi": _is_exact_zero(
            ratio.subs(r, sp.sqrt(24)) - (1 + sp.sqrt(5)) / 2
        ),
    }
    checks = {
        "face_factor_identities": all(face_checks),
        "face_three_root_points_excluded": all(face_three_root_checks.values()),
        "low_amplitude_boundary_kkt_sets_empty": all(low_kkt_checks),
        "uniform_face_points_are_not_global": all(uniform_face_checks),
        "J_nonnegative_sum_identity": _is_exact_zero(j_identity),
        "H_elementary_symmetric_identity": _is_exact_zero(
            h_elementary_identity
        ),
        "coordinate_objective_transform_identity": _is_exact_zero(
            coordinate_transform_identity
        ),
        "boundary_b_sum": _is_exact_zero(b_minus + b_plus - 6),
        "boundary_b_sum_squares": _is_exact_zero(
            b_minus**2 + b_plus**2 - (6 + r**2 / 2)
        ),
        "boundary_G_identity": _is_exact_zero(g_boundary_direct - g_boundary),
        "energy_scaling_identity": _is_exact_zero(e_boundary - r**3 * g_boundary / 6),
        "boundary_beats_m1_factor_identity": _is_exact_zero(
            g1_difference - g1_expected
        ),
        "quadrupole_spectrum_and_ratio_identities": all(spectrum_checks.values()),
    }
    require(
        all(checks.values()),
        "CLOSED_SIMPLEX_CLASSIFICATION",
        "failed exact checks: "
        + ", ".join(name for name, passed in checks.items() if not passed),
    )
    return {
        "coordinate_dictionary": {
            "a_i": (
                "six normalized W5 antipodal-pair values, repeated twice; "
                "sum a_i=0 and 2*sum a_i^2=1"
            ),
            "x_i": "r*a_i, so q_port_i=(1+x_i)/12",
            "b_i": "1+x_i=1+r*a_i=12*q_port_i",
            "p_i": "b_i/6=2*q_port_i, the probability of antipodal pair i",
            "K": "sum b_i^2 = 6+r^2/2",
        },
        "domain": "b_i >= 0, sum b_i=6, sum b_i^2=K",
        "objective_transform": {
            "H": "sum_i (b_i^4-6*b_i^3)",
            "relation_to_reduced_G": "r^3*G = H+12*K-42",
            "relation_to_original_energy": "E=r^3*G/6",
        },
        "face_minimum_structure": {
            "stationarity": (
                "on each relative face, every positive b_i is a root of "
                "4*z^3-18*z^2-2*mu*z-lambda=0"
            ),
            "three_root_exclusion": (
                "the same exact split/group-tangent saddle factor excludes "
                "three roots for support at least four; at support three, "
                "Vieta gives root sum 9/2 while the face sum is 6"
            ),
            "exact_sign_checks": {
                "checks": face_three_root_checks,
                "group_tangent_cases": face_three_root_sign_rows,
            },
            "conclusion": "a face minimum has at most two positive values",
        },
        "low_amplitude_boundary_exclusion": {
            "domain": "6 < K < 18",
            "method": (
                "exact SymPy reduction of every proper-face two-root KKT "
                "case, including split and zero-activation necessary conditions"
            ),
            "cases": low_kkt_rows,
            "isolated_uniform_face_checks": uniform_face_rows,
        },
        "high_amplitude_face_comparison": {
            "domain": "18 < K < 36",
            "feasibility_reduction": (
                "a two-value positive face with q high roots requires q*K<36; "
                "K>18 therefore forces q=1"
            ),
            "H_2": _sympy_text(h2),
            "factor_table": face_rows,
        },
        "independent_nonnegative_global_certificate": {
            "definition": "J=e3(p)-4*e4(p)",
            "exact_identity": (
                "J=sum_{i<j<k} p_i*p_j*p_k*(p_i+p_j+p_k)"
            ),
            "nonnegative": _is_exact_zero(j_identity),
            "equality_condition": "support(p) has at most two indices",
            "objective_relation": (
                "H=1296*(-e2+2*e2^2+J), with e2 fixed by K"
            ),
            "scope": (
                "for K>=18 a support-two point is feasible and saturates "
                "the exact global lower bound"
            ),
        },
        "boundary_minimizer": {
            "domain": "sqrt(24) <= r < sqrt(60)",
            "b_values": {
                "four_values": "0",
                "b_minus": "3-sqrt(r^2-24)/2",
                "b_plus": "3+sqrt(r^2-24)/2",
            },
            "x_values": {
                "four_values": "-1",
                "x_minus": "2-sqrt(r^2-24)/2",
                "x_plus": "2+sqrt(r^2-24)/2",
            },
            "definition_guard": (
                "b_i=1+x_i=1+r*a_i=12*q_port_i; the b and x values "
                "must not be interchanged"
            ),
            "G_boundary": _sympy_text(g_boundary),
            "original_energy_E_boundary": _sympy_text(e_boundary),
            "G_boundary_minus_G_m1_factor": sp.sstr(g1_expected),
            "strict_sign": (
                "negative for sqrt(24)<=r<sqrt(60); it vanishes only at "
                "the excluded upper endpoint r=2*sqrt(15)=sqrt(60)"
            ),
        },
        "boundary_quadrupole": {
            "spectrum_up_to_positive_scale": ["-2", "1-d", "1+d"],
            "d": "sqrt((r^2-15)/5)",
            "domain_order": "sqrt(24)<=r<sqrt(60) gives 0<d<3",
            "sorted_gap_ratio_R": "2*d/(3-d)",
            "inverse_r_squared": "15+45*R^2/(R+2)^2",
            "derivation": (
                "the two selected icosahedral axes have squared inner "
                "product 1/5; the exact rank-two discriminant is 4*d^2"
            ),
            "exact_checks": spectrum_checks,
        },
        "closed_domain_classification": [
            {"domain": "0<r<r_c", "global_minimizer": "m=1 full-support orbit"},
            {"domain": "r=r_c", "global_minimizers": ["m=1", "m=2"]},
            {
                "domain": "r_c<r<sqrt(24)",
                "global_minimizer": "m=2 full-support golden orbit",
            },
            {
                "domain": "sqrt(24)<=r<sqrt(60)",
                "global_minimizer": "support-two boundary orbit with four zero weights",
            },
            {
                "domain": "r=sqrt(60)",
                "global_minimizer": "one pair has b=6 and the other five have b=0",
            },
        ],
        "full_global_minimizer_classification_proved": all(checks.values()),
        "checks": checks,
    }


def exact_global_quartic_certificate() -> dict[str, Any]:
    """Assemble the exact strict-domain and closed-simplex theorem."""

    r = sp.symbols("r", positive=True)
    rows, objectives = _two_value_rows(r)
    saddle = _three_root_saddle_certificate(r)
    r_c = sp.simplify((32 * sp.sqrt(15) - 20 * sp.sqrt(6)) / 27)
    diff_12 = sp.factor(objectives[1] - objectives[2])
    diff_23 = sp.factor(objectives[2] - objectives[3])
    diff_13 = sp.factor(objectives[1] - objectives[3])
    comparisons = {
        "G1_minus_G2": _sympy_text(diff_12),
        "G1_minus_G2_expected": "9*(r-r_c)/80",
        "G2_minus_G3": _sympy_text(diff_23),
        "G1_minus_G3": _sympy_text(diff_13),
    }
    exact_checks = {
        "two_value_table": all(
            all(row["exact_checks"].values()) for row in rows
        ),
        "three_root_points_are_saddles": all(saddle["checks"].values()),
        "G1_G2_crossing_identity": _is_exact_zero(
            diff_12 - sp.Rational(9, 80) * (r - r_c)
        ),
        "G2_G3_identity": _is_exact_zero(
            diff_23 - (r - 4 * sp.sqrt(6)) / 48
        ),
        "G1_G3_identity": _is_exact_zero(
            diff_13 - sp.Rational(2, 15) * (r - sp.sqrt(15))
        ),
        "crossing_positive": bool(r_c.is_positive),
        "crossing_below_sqrt24": bool((sp.sqrt(24) - r_c).is_positive),
        "m3_cannot_overtake_before_sqrt24": bool(
            (4 * sp.sqrt(6) - sp.sqrt(24)).is_positive
        ),
        "m1_beats_m3_through_first_crossing": bool(
            (sp.sqrt(15) - r_c).is_positive
        ),
        "antipodal_m4_is_strictly_above_m2": bool(
            sp.simplify(objectives[4] - objectives[2]).is_positive
        ),
        "antipodal_m5_is_strictly_above_m1": bool(
            sp.simplify(objectives[5] - objectives[1]).is_positive
        ),
    }
    require(
        all(exact_checks.values()),
        "STRICT_GLOBAL_CLASSIFICATION",
        "an exact strict-domain identity failed",
    )
    closed = _closed_simplex_certificate(r, objectives)
    all_checks_pass = all(exact_checks.values()) and all(closed["checks"].values())
    require(
        all_checks_pass,
        "QUARTIC_GLOBAL_CLASSIFICATION",
        "the global quartic certificate did not close",
    )
    return {
        "proof_kind": "exact_sympy_global_minimizer_classification",
        "six_pair_reduction": {
            "coordinates": "w=(a1,a1,...,a6,a6) after ordering antipodal pairs",
            "linear_constraint": "sum_i a_i=0",
            "unit_sphere_constraint": "2*sum_i a_i^2=1",
            "reduced_objective": "G_r=-S3+(r/2)*S4, S_k=2*sum_i a_i^k",
            "original_energy_relation": "E=(r^3/6)*G_r",
        },
        "stationarity_and_saddles": saddle,
        "two_value_stationary_table": rows,
        "crossing": {
            "r_c_exact": "(32*sqrt(15)-20*sqrt(6))/27",
            "r_c_display": f"{float(sp.N(r_c, 16)):.13f}",
            "comparisons": comparisons,
        },
        "strict_full_support_classification": {
            "domain": "0<r<sqrt(60) and 1+r*a_i>0 for every i",
            "intervals": [
                {"domain": "0<r<r_c", "global_minimizer": "m=1 C5 orbit"},
                {"domain": "r=r_c", "global_minimizers": ["m=1 C5", "m=2 golden"]},
                {"domain": "r_c<r<sqrt(24)", "global_minimizer": "m=2 golden orbit"},
                {
                    "domain": "sqrt(24)<=r<sqrt(60)",
                    "global_minimizer": None,
                    "infimum": "support-two closed-simplex boundary branch",
                    "attained": False,
                },
            ],
            "full_global_minimizer_classification_proved": all_checks_pass,
        },
        "closed_probability_simplex": closed,
        "checks": exact_checks,
        "all_exact_checks_pass": all_checks_pass,
    }


# ---------------------------------------------------------------------------
# Comparison block (compare-only; explicit ancestry)
# ---------------------------------------------------------------------------


def comparison_block() -> dict[str, Any]:
    # Measured values are deliberately scoped to this compare-only function.
    pdg_lepton_masses_gev = {
        "electron": Fraction("0.00051099895069"),
        "muon": Fraction("0.1056583755"),
        "tau": Fraction("1.77693"),
    }
    logs = {k: math.log(float(v)) for k, v in pdg_lepton_masses_gev.items()}
    ordered = sorted(logs.values())
    mean = sum(ordered) / 3.0
    centered = [value - mean for value in ordered]
    gap_low = centered[1] - centered[0]
    gap_high = centered[2] - centered[1]
    observed = gap_high / gap_low
    phi = PHI.to_float()
    mismatch_direct = abs(observed - phi) / phi
    mismatch_flipped = abs((1.0 / observed) - phi) / phi
    closest = min(mismatch_direct, mismatch_flipped)
    boundary_ratio = 1.0 / observed
    boundary_inferred_r = math.sqrt(
        15.0 + 45.0 * boundary_ratio**2 / (boundary_ratio + 2.0) ** 2
    )
    return {
        "ancestry": (
            "measured PDG charged-lepton masses, compare-only; no lepton "
            "value enters any derivation above"
        ),
        "observed_sorted_gap_ratio": f"{observed:.10f}",
        "observed_flipped": f"{1.0 / observed:.10f}",
        "golden_branch_compare": {
            "packet_output": (
                "phi = 1.6180339887 (m=2 full-support branch); "
                "the m=1 branch is degenerate"
            ),
            "relative_mismatch_direct": f"{mismatch_direct:.6f}",
            "relative_mismatch_flipped": f"{mismatch_flipped:.6f}",
            "closest_pairing_mismatch": f"{closest:.6f}",
            "verdict": "GOLDEN_BRANCH_EXCLUDED_COMPARE_ONLY",
        },
        "boundary_reentry_compare": {
            "orientation": "larger-over-smaller sorted quadrupole gap",
            "observed_R": f"{boundary_ratio:.10f}",
            "target_attached_inferred_r": f"{boundary_inferred_r:.10f}",
            "theory_map": (
                "R=2*d/(3-d), d=sqrt((r^2-15)/5), "
                "r^2=15+45*R^2/(R+2)^2"
            ),
            "domain_check": "sqrt(24)<r<sqrt(60)",
            "verdict": "BOUNDARY_BRANCH_CAN_MATCH_CENTRAL_SHAPE_COMPARE_ONLY",
            "nonclaim": (
                "r is inferred from the measured shape, not source-emitted; "
                "this is not a prediction or validation"
            ),
        },
        "verdict": "STRICT_FULL_SUPPORT_CONDITIONAL_NO_GO_BOUNDARY_REENTRY_VIABLE",
        "global_packet_verdict": "EXACTLY_CLASSIFIED_NOT_GLOBALLY_EXCLUDED",
        "verdict_statement": (
            "the exact strict full-support minima are either quadrupole-"
            "degenerate, the golden branch excluded by the compare-only "
            "shape check, or absent at high amplitude. The exact closed-"
            "simplex boundary branch continuously re-enters the observed "
            "shape at a target-attached amplitude, so the quartic packet is "
            "not globally excluded"
        ),
    }


# ---------------------------------------------------------------------------
# Controls
# ---------------------------------------------------------------------------


def control_prolate_simple_spectrum() -> dict[str, Any]:
    """The prolate branch must fail the simple-spectrum gate."""

    branch = prolate_orbit()
    try:
        require(
            branch["double_eigenvalue"] is False,
            "PROLATE_DEGENERATE",
            "the prolate branch has a double eigenvalue",
        )
    except CertificateError:
        return {
            "expected_failure": True,
            "failed": True,
            "code": "PROLATE_DEGENERATE",
        }
    return {"expected_failure": True, "failed": False}


def control_profile_mutation() -> dict[str, Any]:
    """A mutated extremal profile must break the golden-ratio identity, so
    the identity gate refuses it (the required failure)."""

    mutated = (Q5(2), Q5(-1) + Q5(Fraction(1, 7)), Q5(-1) - Q5(Fraction(1, 7)))
    q = invert_profile(mutated)
    triple = sorted(q, key=lambda v: v.to_float())
    ratio = (triple[2] - triple[1]) / (triple[1] - triple[0])
    try:
        require(ratio == PHI, "GOLDEN_RATIO", "mutated profile must not give phi")
    except CertificateError:
        return {
            "expected_failure": True,
            "failed": True,
            "code": "GOLDEN_RATIO",
            "meaning": "the golden identity consumes the exact extremal profile; a mutation is refused",
        }
    return {"expected_failure": True, "failed": False}


def control_frame_band_cubic_vanishes() -> dict[str, Any]:
    """An attempt to relocate the cubic functional to an antipode-odd
    profile must be refused: odd cubes cancel exactly."""

    values = [Q5(3), Q5(-3), Q5(Fraction(1, 2)), Q5(Fraction(-1, 2)), Q5(7), Q5(-7)]
    total = ZERO
    for value in values:
        total = total + value * value * value
    try:
        require(not total.is_zero(), "FRAME_CUBIC", "antipode-odd cubes cancel, so the cubic cannot live on the frame band")
    except CertificateError:
        return {
            "expected_failure": True,
            "failed": True,
            "code": "FRAME_CUBIC",
            "meaning": "the cubic response vanishes on any antipode-odd profile, which pins the shape functional to the quintet band",
        }
    return {"expected_failure": True, "failed": False}


def control_mutated_quartic_crossing() -> dict[str, Any]:
    """A shifted crossing must fail the exact G1-G2 identity."""

    r = sp.symbols("r", positive=True)
    _, objectives = _two_value_rows(r)
    r_c = (32 * sp.sqrt(15) - 20 * sp.sqrt(6)) / 27
    mutated = r_c + 1
    identity = sp.simplify(
        objectives[1]
        - objectives[2]
        - sp.Rational(9, 80) * (r - mutated)
    )
    try:
        require(
            identity == 0,
            "QUARTIC_CROSSING",
            "a shifted crossing cannot satisfy the exact branch identity",
        )
    except CertificateError:
        return {
            "expected_failure": True,
            "failed": True,
            "code": "QUARTIC_CROSSING",
        }
    return {"expected_failure": True, "failed": False}


def control_boundary_coordinate_mixup() -> dict[str, Any]:
    """The boundary b and x coordinates differ by one and cannot be relabeled."""

    r = sp.symbols("r", positive=True)
    b_minus = 3 - sp.sqrt(r**2 - 24) / 2
    x_minus = 2 - sp.sqrt(r**2 - 24) / 2
    try:
        require(
            _is_exact_zero(b_minus - x_minus),
            "BOUNDARY_COORDINATE_MIXUP",
            "b_i=1+x_i, so their boundary values cannot be interchanged",
        )
    except CertificateError:
        return {
            "expected_failure": True,
            "failed": True,
            "code": "BOUNDARY_COORDINATE_MIXUP",
        }
    return {"expected_failure": True, "failed": False}


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------


def build_payload() -> dict[str, Any]:
    prolate = prolate_orbit()
    golden = golden_orbit()
    exact_global = exact_global_quartic_certificate()

    controls = {
        "prolate_simple_spectrum": control_prolate_simple_spectrum(),
        "profile_mutation": control_profile_mutation(),
        "frame_band_cubic": control_frame_band_cubic_vanishes(),
        "mutated_quartic_crossing": control_mutated_quartic_crossing(),
        "boundary_coordinate_mixup": control_boundary_coordinate_mixup(),
    }
    for name, verdict in controls.items():
        require(
            verdict["expected_failure"] is True and verdict["failed"] is True,
            "CONTROL_NOT_FAILED",
            f"control {name} did not record its required failure",
        )

    def strip(branch: dict[str, Any]) -> dict[str, Any]:
        out = dict(branch)
        for key in ("S2", "S3", "S4", "min_port_value"):
            out[key] = branch[key].render()
        return out

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "issue_context": ISSUE_CONTEXT,
        "status": "EXACT_GLOBAL_QUARTIC_MINIMIZER_CLASSIFICATION",
        "row_class": "exact_symbolic_global_minimizer_theorem_with_compare_only_physical_attachment",
        "physical_claim": False,
        "promotion_allowed": False,
        "claim_boundary": (
            "Exact global-minimizer theorem for the quartic-truncated W5 "
            "entropy functional on both the strict full-support domain and "
            "its closed probability simplex. It does not establish the "
            "quadrupole-to-physical-log-mass attachment, emit the amplitude, "
            "control the entropy remainder at zero weights, or close the "
            "higher-order/full-entropy mechanism. No physical mass, ratio, "
            "or prediction is emitted."
        ),
        "exhaustiveness_boundary": {
            "minimizer_domains_classified": [
                "strict full-support probability domain for 0<r<sqrt(60)",
                "closed probability simplex for 0<r<=sqrt(60)",
            ],
            "full_critical_orbit_classification_proved": False,
            "global_minimizer_classification_proved": exact_global[
                "all_exact_checks_pass"
            ],
            "quartic_packet_globally_excluded": False,
            "strict_physical_no_go_conditional_on": [
                "strict full support q_i>0",
                "the quadrupole-to-physical-log-mass attachment",
                "the quartic truncation being the operative selector",
            ],
            "viable_routes_not_excluded": [
                "zero-weight closed-simplex boundary states",
                "strictly positive regularizations approaching the boundary",
                "higher entropy orders and the full entropy functional",
                "a source-derived amplitude law",
                "a different source-derived W5 effective action",
            ],
            "remainder_warning": (
                "the Taylor expansion about the uniform record is least "
                "controlled at boundary points with q_i=0"
            ),
        },
        "packet_derivation": {
            "expansion": "D(u(1+x)||u) = u(S2/2 - S3/6 + S4/12) + O(x^5), u = 1/12",
            "band_selection": (
                "S3 vanishes identically on the antipode-odd frame band, so "
                "the quintet band is the only nontrivial band with a cubic "
                "entropy response"
            ),
            "premise": (
                "the shape functional is the quartic-truncated expansion "
                "restricted to the quintet band at fixed band amplitude; the "
                "amplitude is a declared parameter of the premise, ranging "
                "over the probability domain; no source-derived amplitude "
                "emitter is supplied"
            ),
        },
        "exact_global_quartic_certificate": exact_global,
        "exact_candidate_geometry": {
            "role": (
                "redundant exact Q(sqrt(5)) identification of the m=1 C5 "
                "and m=2 golden quadrupole orbits; not the global proof"
            ),
            "prolate": strip(prolate),
            "golden": strip(golden),
        },
        "comparison": comparison_block(),
        "controls": controls,
        "boundary_extension": (
            "the exact closed-simplex minimizer has a continuous boundary "
            "gap-ratio branch, so the packet is not globally excluded; zero-"
            "weight, regularized, higher-order, full-entropy, and source-"
            "amplitude routes remain viable"
        ),
        "finite_seed_or_numerical_optimization_used_as_proof": False,
    }
    return payload


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=OUT_PATH)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    payload = build_payload()
    if args.verify:
        stored = json.loads(args.out.read_text(encoding="utf-8"))
        require(stored == payload, "DRIFT", "stored certificate does not match a rebuild")
        print(json.dumps({"status": "PASS"}))
        return 0
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "WROTE", "out": str(args.out), "verdict": payload["comparison"]["verdict"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
