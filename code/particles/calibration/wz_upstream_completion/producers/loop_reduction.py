#!/usr/bin/env python3
"""One-loop two-point reduction core for the direct FJ engine.

Exact symbolic reduction of one-loop two-point integrands to the
denominator-power basis

    I(a, b) = integral of 1 / (D1^a D2^b),   0 <= a, b <= 2,

with D1 = k^2 - m1sq (line V1 -> V2, momentum k) and
D2 = (k - p)^2 - m2sq (line V2 -> V1, momentum k - p), and the
space-time dimension d = 4 - 2 eps kept symbolic everywhere.  Nothing
is hand-entered: numerators arrive as polynomials in the invariants
k2 = k.k and kp = k.p from the tensor projector, and the reduction
uses only the exact identities

    k2 = D1 + m1sq,
    kp = (D1 - D2 + p2 + m1sq - m2sq) / 2,

term-by-term cancellation of D powers, and the dimensional
regularization facts that pure polynomial integrals vanish and that
odd single-denominator moments vanish after the shift q = k - p.

Squared denominators appear only through the exact partial-fraction
split of the R_xi vector propagator on a massless line; massive lines
split into single-power pieces with masses m^2 and xi m^2 before the
reduction is called.

Basis names: I(1,1) = B0, I(1,0) = A0(m1sq), I(0,1) = A0(m2sq),
I(2,0) = A0p(m1sq) = dA0/dm1sq, I(0,2) = A0p(m2sq), I(2,1), I(1,2),
I(2,2) kept as finite symbols.  Single UV poles: B0 -> Delta,
A0(m) -> m Delta, A0p -> Delta, I(2,1) = I(1,2) = I(2,2) -> 0.

A two-point tensor block is fully described by its two projections
g_{mu nu} Pi^{mu nu} and p_mu p_nu Pi^{mu nu}; the transverse and
longitudinal parts are Pi_T = (g.Pi - p.Pi/p2)/(d - 1) and
Pi_L = p.Pi/p2.
"""

from __future__ import annotations

import sympy as sp

d_sym = sp.Symbol("d")
eps = sp.Symbol("eps")
p2 = sp.Symbol("p2")
k2 = sp.Symbol("k2")
kp = sp.Symbol("kp")
m1sq = sp.Symbol("m1sq")
m2sq = sp.Symbol("m2sq")
D1 = sp.Symbol("D1")
D2 = sp.Symbol("D2")
Delta = sp.Symbol("Delta")

BASIS = {
    (1, 1): sp.Symbol("B0"),
    (1, 0): sp.Symbol("A0m1"),
    (0, 1): sp.Symbol("A0m2"),
    (2, 0): sp.Symbol("A0pm1"),
    (0, 2): sp.Symbol("A0pm2"),
    (2, 1): sp.Symbol("C21"),
    (1, 2): sp.Symbol("C12"),
    (2, 2): sp.Symbol("C22"),
}

POLE_OF = {
    sp.Symbol("B0"): Delta,
    sp.Symbol("A0m1"): m1sq * Delta,
    sp.Symbol("A0m2"): m2sq * Delta,
    sp.Symbol("A0pm1"): Delta,
    sp.Symbol("A0pm2"): Delta,
    sp.Symbol("C21"): sp.Integer(0),
    sp.Symbol("C12"): sp.Integer(0),
    sp.Symbol("C22"): sp.Integer(0),
}


def dot(a: str, b: str) -> sp.Expr:
    pair = tuple(sorted((a, b)))
    if pair == ("k", "k"):
        return k2
    if pair == ("k", "p"):
        return kp
    if pair == ("p", "p"):
        return p2
    raise ValueError(f"unknown scalar product {pair}")


def project(terms: list[tuple[sp.Expr, str | None, str | None]], projector: str) -> sp.Expr:
    """Contract a rank-two structure with g_{mu nu} or p_mu p_nu.

    Each term is (coefficient, carrier_mu, carrier_nu); a carrier is
    "k", "p", or None, and None in both slots is an explicit metric."""

    total = sp.Integer(0)
    for coefficient, car_mu, car_nu in terms:
        if car_mu is None and car_nu is None:
            factor = d_sym if projector == "g" else p2
        elif car_mu is not None and car_nu is not None:
            if projector == "g":
                factor = dot(car_mu, car_nu)
            else:
                factor = dot(car_mu, "p") * dot(car_nu, "p")
        else:
            raise ValueError("mixed open/metric term")
        total += coefficient * factor
    return sp.expand(total)


def reduce_two_point(numerator: sp.Expr, a_power: int, b_power: int) -> sp.Expr:
    """Reduce integral of numerator / (D1^a D2^b) to the basis.

    After the k2/kp substitution the integrand is a Laurent monomial
    sum D1^alpha D2^beta / (D1^a D2^b); positive residual powers on
    both slots make a polynomial (zero), residual (a', b') with
    0 <= a', b' <= 2 lands on a basis element, and residual powers
    beyond the basis are reduced by the shifted-tadpole identities:
    every D1^n / D2 tadpole reduces through D1 = Dq + 2 q.p +
    (p2 + m2sq - m1sq) with odd q moments dropped and even q moments
    reduced recursively in the single denominator Dq."""

    substituted = sp.expand(numerator.subs({
        k2: D1 + m1sq,
        kp: (D1 - D2 + p2 + m1sq - m2sq) / 2,
    }))
    poly = sp.Poly(substituted, D1, D2)
    total = sp.Integer(0)
    for (alpha, beta), coefficient in poly.terms():
        res_a = a_power - alpha
        res_b = b_power - beta
        if res_a <= 0 and res_b <= 0:
            continue  # polynomial integrand
        if res_a < 0:
            total += coefficient * excess_tadpole(-res_a, res_b, first_excess=True)
            continue
        if res_b < 0:
            total += coefficient * excess_tadpole(-res_b, res_a, first_excess=False)
            continue
        total += coefficient * BASIS[(res_a, res_b)]
    return sp.expand(total)


def excess_tadpole(excess: int, other_power: int, first_excess: bool) -> sp.Expr:
    """Integral of D1^excess / D2^other (or mirrored).

    Shift to the D2 line momentum q = k - p:
    D1 = Dq + 2 q.p + s with s = p2 + m2sq - m1sq, Dq = q^2 - m2sq.
    Odd q moments vanish; (q.p)^2 -> q^2 p2 / d = (Dq + m2sq) p2 / d.
    The engine needs excess <= 2 against other_power <= 2."""

    if other_power == 0:
        return sp.Integer(0)  # polynomial
    if first_excess:
        s = p2 + m2sq - m1sq
        base = {(1,): BASIS[(0, 1)], (2,): BASIS[(0, 2)]}[(other_power,)]
        msq = m2sq
    else:
        s = p2 + m1sq - m2sq
        base = {(1,): BASIS[(1, 0)], (2,): BASIS[(2, 0)]}[(other_power,)]
        msq = m1sq
    if excess == 1:
        # (Dq + s) / Dq^n: the Dq/Dq^n piece lowers the power.
        lowered = tadpole_symbol(other_power - 1, first_excess)
        return sp.expand(s * base + lowered)
    if excess == 2:
        # (Dq + 2 q.p + s)^2 with odd terms dropped:
        # Dq^2 + 2 Dq s + s^2 + 4 (q.p)^2
        qp2 = (p2 / d_sym) * (tadpole_symbol(other_power - 1, first_excess) + msq * base)
        return sp.expand(
            tadpole_symbol(other_power - 2, first_excess)
            + 2 * s * tadpole_symbol(other_power - 1, first_excess)
            + s ** 2 * base
            + 4 * qp2
        )
    raise ValueError("excess above the two-point rank")


def tadpole_symbol(power: int, on_second_line: bool) -> sp.Expr:
    if power <= 0:
        return sp.Integer(0)  # 1/Dq^0: polynomial
    key = (0, power) if on_second_line else (power, 0)
    return BASIS[key]


def reduce_tadpole(numerator: sp.Expr) -> sp.Expr:
    """Reduce integral of numerator/D1 to the canonical A0 basis; the
    numerator is polynomial in k2 and kp (odd kp moments vanish and
    kp^2 -> k2 p2 / d by Lorentz averaging)."""

    substituted = sp.expand(numerator.subs(kp ** 2, k2 * p2 / d_sym).subs(kp, 0))
    substituted = sp.expand(substituted.subs(k2, D1 + m1sq))
    poly = sp.Poly(substituted, D1)
    total = sp.Integer(0)
    for (a,), coefficient in poly.terms():
        if a == 0:
            total += coefficient * BASIS[(1, 0)]
    return sp.expand(total)


def transverse_longitudinal(g_projection: sp.Expr, p_projection: sp.Expr) -> tuple[sp.Expr, sp.Expr]:
    pi_l = sp.expand(p_projection / p2)
    pi_t = sp.expand((g_projection - pi_l) / (d_sym - 1))
    return pi_t, pi_l


def uv_pole(expression: sp.Expr, mass_map: dict | None = None) -> sp.Expr:
    """Single-pole coefficient (unit Delta) with d = 4 - 2 eps and the
    eps x Delta cross terms routed to the finite part exactly.  The
    optional mass_map substitutes the canonical line masses AFTER the
    divergent parts are inserted, so reduced expressions can stay in
    canonical (m1sq, m2sq) form until pole extraction."""

    with_delta = expression.subs(POLE_OF)
    if mass_map:
        with_delta = with_delta.subs(mass_map)
    with_delta = sp.expand(with_delta.subs(d_sym, 4 - 2 * eps))
    if with_delta.has(sp.Symbol("B0"), sp.Symbol("A0m1"), sp.Symbol("A0m2")):
        raise ValueError("unsubstituted basis symbol in pole extraction")
    coefficient = with_delta.coeff(Delta)
    series = sp.series(sp.together(coefficient), eps, 0, 1).removeO()
    return sp.simplify(series.subs(eps, 0))
