#!/usr/bin/env python3
"""Interval contraction certificate for the OPH P/alpha closure map.

This is the stage-2 (basin-then-contract) contraction certificate for the P
coordinate. It upgrades the sampled finite-difference surfaces emitted by
``fixed_point_witness.py`` and ``fixed_point_certificate.py`` to a rigorous
interval enclosure:

- the declared closure chain ``alpha -> P -> M_U -> alpha_U -> alpha_i(m_Z)
  -> a0(P) -> Delta_Th(P) -> alpha`` is re-implemented, formula by formula,
  from ``paper_math.py`` in mpmath binary interval arithmetic (``mpmath.iv``,
  outward rounding on every elementary operation);
- the two implicit solves in the chain (the tree-level ``m_Z`` closure and the
  D10 pixel closure for ``alpha_U``) are enclosed by verified sign-change
  brackets with epsilon-inflation, with the interval derivative of each
  residual proven sign-definite on the bracket (existence, local uniqueness,
  and C^1 dependence by the implicit function theorem);
- the derivative of the full map is enclosed by forward-mode dual-number
  automatic differentiation over intervals, with the implicit nodes handled by
  the implicit function theorem;
- the SU(2)/SU(3) edge-sum truncation tails are bounded by explicit geometric
  majorants and added to the interval enclosures, so the certificate covers
  both the declared finite-cutoff sums and their infinite-cutoff limits;
- the certificate is direct Banach: an interval ``I`` in alpha space with
  ``g(I)`` contained in the interior of ``I`` and an interval bound
  ``|g'| <= L < 1`` on ``I``, which proves existence and uniqueness of the
  fixed point in ``I``.

Two readout maps are certified:

- ``thomson_structured_running``: the declared solver mode; the fixed point is
  the source/root value near ``alpha_inv = 136.9948``;
- ``thomson_structured_running_plus_gauge_width``: the closure-ledger CL-2
  mixed map that adds the unified gauge width ``alpha_U(P)`` to the
  inverse-alpha readout; its fixed point is near ``alpha_inv = 137.03596``.

Claim boundary: this certifies the declared numerical map at the declared
representation cutoffs and declared one-loop RG/matching, tree-level m_Z,
Stage-5 continuation-mass, and exact one-loop kernel conventions. It is not an
exact fine-structure derivation and does not change the stage-3 landing
verdict (closure ledger row CL-1).

mpmath contexts are private instances (``MPContext`` / ``MPIntervalContext``);
the global ``mpmath.mp`` and ``mpmath.iv`` precision settings are never
touched.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any, Callable

from mpmath.ctx_iv import MPIntervalContext
from mpmath.ctx_mp import MPContext


ARTIFACT_DATE = "2026-07-14"
DEFAULT_OUT = (
    Path(__file__).resolve().parent
    / "runtime"
    / f"p_interval_contraction_certificate_{ARTIFACT_DATE}.json"
)

MODE_SOURCE = "thomson_structured_running"
MODE_GAUGE_WIDTH = "thomson_structured_running_plus_gauge_width"
MODES = (MODE_SOURCE, MODE_GAUGE_WIDTH)

N_C = 3
N_G = 3
BETA_EW = 4  # N_c + 1


class CertificateError(RuntimeError):
    """Raised when a rigorous verification step cannot be completed."""


# ---------------------------------------------------------------------------
# Backends: plain mp point arithmetic, iv interval arithmetic, and dual-number
# forward AD lifted over either.
# ---------------------------------------------------------------------------


class Dual:
    """Forward-mode dual number over any mpmath-like scalar type."""

    __slots__ = ("x", "d")

    def __init__(self, x: Any, d: Any) -> None:
        self.x = x
        self.d = d

    def __add__(self, other: "Dual") -> "Dual":
        return Dual(self.x + other.x, self.d + other.d)

    def __sub__(self, other: "Dual") -> "Dual":
        return Dual(self.x - other.x, self.d - other.d)

    def __neg__(self) -> "Dual":
        return Dual(-self.x, -self.d)

    def __mul__(self, other: "Dual") -> "Dual":
        return Dual(self.x * other.x, self.x * other.d + self.d * other.x)

    def __truediv__(self, other: "Dual") -> "Dual":
        q = self.x / other.x
        return Dual(q, (self.d - q * other.d) / other.x)


class _BaseBackend:
    """Shared constant table for a scalar context."""

    def _init_constants(self) -> None:
        self.zero = self.num(0)
        self.one = self.num(1)
        self.two = self.num(2)
        self.sqrt_pi = self.sqrt(self.pi)
        self.two_pi = self.two * self.pi
        self.four_pi = self.num(4) * self.pi
        self.four_pi_sq = self.four_pi * self.pi
        self.phi = (self.one + self.sqrt(self.num(5))) / self.two
        self.beta_ew = self.num(BETA_EW)
        self.b_mssm = (
            self.num(33) / self.num(5),
            self.one,
            -self.num(3),
        )

    def pow_int(self, x: Any, n: int) -> Any:
        if n < 0:
            return self.one / self.pow_int(x, -n)
        result = self.one
        base = x
        while n:
            if n & 1:
                result = result * base
            base = base * base
            n >>= 1
        return result


class MpBackend(_BaseBackend):
    kind = "mp"

    def __init__(self, dps: int) -> None:
        self.ctx = MPContext()
        self.ctx.dps = dps
        self.dps = dps
        self.pi = self.ctx.pi
        self._init_constants()

    def num(self, v: Any) -> Any:
        return self.ctx.mpf(v)

    def exp(self, x: Any) -> Any:
        return self.ctx.exp(x)

    def log(self, x: Any) -> Any:
        return self.ctx.log(x)

    def sqrt(self, x: Any) -> Any:
        return self.ctx.sqrt(x)

    def cos(self, x: Any) -> Any:
        return self.ctx.cos(x)


class IvBackend(_BaseBackend):
    kind = "iv"

    def __init__(self, dps: int) -> None:
        self.ctx = MPIntervalContext()
        self.ctx.dps = dps
        self.dps = dps
        self.pi = self.ctx.pi
        self._init_constants()

    def num(self, v: Any) -> Any:
        return self.ctx.mpf(v)

    def exp(self, x: Any) -> Any:
        return self.ctx.exp(x)

    def log(self, x: Any) -> Any:
        return self.ctx.log(x)

    def sqrt(self, x: Any) -> Any:
        return self.ctx.sqrt(x)

    def cos(self, x: Any) -> Any:
        return self.ctx.cos(x)

    def hull(self, lo: Any, hi: Any) -> Any:
        return self.ctx.mpf([lo, hi])

    def thin(self, x: Any) -> Any:
        return self.ctx.mpf(x)

    def contains_zero(self, x: Any) -> bool:
        return self.ctx.mpf(0) in x

    def sup_abs(self, x: Any) -> Any:
        neg_lo = -x.a
        hi = x.b
        return neg_lo if neg_lo > hi else hi


class DualBackend(_BaseBackend):
    """Lift a scalar backend to dual numbers (value, d/dalpha)."""

    kind = "dual"

    def __init__(self, base: _BaseBackend) -> None:
        self.base = base
        self.pi = Dual(base.pi, base.zero)
        self._init_constants()

    def num(self, v: Any) -> Dual:
        return Dual(self.base.num(v), self.base.zero)

    def constant(self, x: Any) -> Dual:
        return Dual(x, self.base.zero)

    def lift(self, x: Any, d: Any) -> Dual:
        return Dual(x, d)

    def exp(self, u: Dual) -> Dual:
        e = self.base.exp(u.x)
        return Dual(e, e * u.d)

    def log(self, u: Dual) -> Dual:
        return Dual(self.base.log(u.x), u.d / u.x)

    def sqrt(self, u: Dual) -> Dual:
        s = self.base.sqrt(u.x)
        return Dual(s, u.d / (self.two.x * s))

    def cos(self, u: Dual) -> Dual:
        c = self.base.cos(u.x)
        # sin(x) = cos(x - pi/2); avoids requiring a separate sin primitive.
        s = self.base.cos(u.x - self.base.pi / self.base.two)
        return Dual(c, -s * u.d)


# ---------------------------------------------------------------------------
# Declared chain, written once over a generic backend.
# The formulas mirror paper_math.py exactly.
# ---------------------------------------------------------------------------


def build_su2_terms(b: _BaseBackend, cutoff: int) -> list[tuple[Any, Any, Any]]:
    terms = []
    for n in range(cutoff + 1):
        dim = n + 1
        # c2 = j(j+1) with j = n/2, i.e. n(n+2)/4.
        terms.append((b.num(dim), b.num(n * (n + 2)) / b.num(4), b.log(b.num(dim))))
    return terms


def build_su3_terms(b: _BaseBackend, cutoff: int) -> list[tuple[Any, Any, Any]]:
    terms = []
    for p in range(cutoff + 1):
        for q in range(cutoff + 1):
            dim = ((p + 1) * (q + 1) * (p + q + 2)) // 2
            c2_num = p * p + q * q + p * q + 3 * p + 3 * q
            terms.append((b.num(dim), b.num(c2_num) / b.num(3), b.log(b.num(dim))))
    return terms


def ellbar(
    b: _BaseBackend,
    terms: list[tuple[Any, Any, Any]],
    t: Any,
    tail: tuple[Any, Any] | None = None,
) -> Any:
    """Weighted mean of log-dimensions: sum(dim e^{-t c2} ln dim)/sum(dim e^{-t c2}).

    ``tail`` optionally supplies one-sided enclosures ``(z_tail, s_tail)`` for
    the omitted terms beyond the cutoff, so the returned interval encloses both
    the finite-cutoff sum and the infinite sum.
    """
    z = b.zero
    s = b.zero
    for dim, c2, ln_dim in terms:
        w = dim * b.exp(-(t * c2))
        z = z + w
        s = s + w * ln_dim
    if tail is not None:
        z = z + tail[0]
        s = s + tail[1]
    return s / z


def mu_u_of_p(b: _BaseBackend, p: Any) -> Any:
    return b.exp(-b.two_pi) * b.exp(b.log(p) / b.num(6))


def v_transmutation(b: _BaseBackend, alpha_u: Any, p: Any) -> Any:
    e_cell = b.one / b.sqrt(p)
    return e_cell * b.exp(-(b.two_pi / (b.beta_ew * alpha_u)))


def alpha_run(b: _BaseBackend, alpha_u: Any, bcoef: Any, mu: Any, mu_u: Any) -> Any:
    inv = b.one / alpha_u + (bcoef / b.two_pi) * b.log(mu_u / mu)
    return b.one / inv


def mz_tree(b: _BaseBackend, v: Any, alpha1: Any, alpha2: Any) -> Any:
    alpha_y = b.num(3) * alpha1 / b.num(5)
    return (v / b.two) * b.sqrt(b.four_pi * alpha2 + b.four_pi * alpha_y)


def h_mz(b: _BaseBackend, mu: Any, alpha_u: Any, mu_u: Any, v: Any) -> Any:
    a1 = alpha_run(b, alpha_u, b.b_mssm[0], mu, mu_u)
    a2 = alpha_run(b, alpha_u, b.b_mssm[1], mu, mu_u)
    return mz_tree(b, v, a1, a2) - mu


def alpha_em_inv(b: _BaseBackend, alpha1: Any, alpha2: Any) -> Any:
    alpha_y = b.num(3) * alpha1 / b.num(5)
    return b.one / alpha2 + b.one / alpha_y


def koide_roots(b: _BaseBackend) -> tuple[Any, Any, Any]:
    """Stage-5 Koide roots, ascending (P-independent constants)."""
    delta = b.beta_ew / (b.two * b.num(N_C) * b.num(N_G))
    sqrt_two = b.sqrt(b.two)
    roots = []
    for k in range(3):
        angle = delta + b.two_pi * b.num(k) / b.num(3)
        roots.append(b.one + sqrt_two * b.cos(angle))
    return roots


def stage5_integer_vectors(mpb: MpBackend) -> dict[str, tuple[int, int, int]]:
    """Mirror paper_math._derive_stage5_integer_vectors (integer outputs only)."""
    b = mpb
    epsilon = b.one / b.num(6)
    roots = sorted(koide_roots(b))
    roots_sq = [r * r for r in roots]
    n_tau = N_G
    n_mu = N_G + 1
    best_n_e = None
    best_residual = None
    for n_e in range(n_mu + 1, n_mu + 8):
        exponents = (n_e, n_mu, n_tau)
        errors = []
        for i in range(3):
            for j in range(i + 1, 3):
                lhs = roots_sq[i] / roots_sq[j]
                rhs = b.pow_int(epsilon, exponents[i] - exponents[j])
                errors.append(abs(b.log(lhs / rhs)))
        residual = max(errors)
        if best_residual is None or residual < best_residual:
            best_residual = residual
            best_n_e = n_e
    if best_n_e is None:
        raise CertificateError("could not derive Stage-5 charged-lepton exponents")
    return {
        "n_u": (2 * N_C, N_C, 0),
        "n_d": (2 * N_C, N_C + 1, N_C - 1),
        "n_e": (best_n_e, n_mu, n_tau),
    }


def sorted_roots(b: _BaseBackend, mp_order: list[int]) -> list[Any]:
    roots = koide_roots(b)
    return [roots[i] for i in mp_order]


def quark_masses(
    b: _BaseBackend, v: Any, vectors: dict[str, tuple[int, int, int]]
) -> dict[str, Any]:
    pref = v / b.sqrt(b.two)
    epsilon = b.one / b.num(6)
    n_u = vectors["n_u"]
    n_d = vectors["n_d"]
    return {
        "u": pref * b.pow_int(epsilon, n_u[0]),
        "c": pref * b.pow_int(epsilon, n_u[1]),
        "t": pref * b.pow_int(epsilon, n_u[2]),
        "d": pref * b.pow_int(epsilon, n_d[0]),
        "s": pref * b.pow_int(epsilon, n_d[1]),
        "b": pref * b.pow_int(epsilon, n_d[2]),
    }


def lepton_masses(
    b: _BaseBackend,
    v: Any,
    vectors: dict[str, tuple[int, int, int]],
    roots_ascending: list[Any],
) -> dict[str, Any]:
    roots_sq = [r * r for r in roots_ascending]
    n_e = vectors["n_e"]
    sqrt_two = b.sqrt(b.two)
    six = b.num(6)
    log_gm_c = b.zero
    for index, exponent in enumerate(n_e):
        numerator = roots_sq[index] * sqrt_two * b.pow_int(six, exponent)
        log_gm_c = log_gm_c + b.log(numerator / v)
    log_gm_c = log_gm_c / b.num(3)
    s0 = b.exp(-log_gm_c)
    scale = s0 * b.exp(b.log(b.two) / six)
    return {
        "e": scale * roots_sq[0],
        "mu": scale * roots_sq[1],
        "tau": scale * roots_sq[2],
    }


def transport_kernel_exact(
    b: _BaseBackend, q_scale: Any, mass: Any, charge_squared: Any, multiplicity: Any
) -> Any:
    z = (q_scale * q_scale) / (mass * mass)
    a = z / b.num(4)
    sqrt_a = b.sqrt(a)
    sqrt_one_plus_a = b.sqrt(b.one + a)
    asinh_sqrt_a = b.log(sqrt_a + sqrt_one_plus_a)
    integral = (
        -(b.num(5) / b.num(18))
        + b.one / (b.num(6) * a)
        + ((b.two * a - b.one) * sqrt_one_plus_a * asinh_sqrt_a)
        / (b.num(6) * a * sqrt_a)
    )
    return b.two * multiplicity * charge_squared / b.pi * integral


def structured_delta_alpha_inv(
    b: _BaseBackend,
    mz: Any,
    v: Any,
    alpha3: Any,
    vectors: dict[str, tuple[int, int, int]],
    roots_ascending: list[Any],
) -> Any:
    quarks = quark_masses(b, v, vectors)
    leptons = lepton_masses(b, v, vectors, roots_ascending)
    screening = b.one - b.num(N_C) * alpha3 / b.pi
    one = b.one
    three = b.num(3)
    lepton_total = (
        transport_kernel_exact(b, mz, leptons["e"], one, one)
        + transport_kernel_exact(b, mz, leptons["mu"], one, one)
        + transport_kernel_exact(b, mz, leptons["tau"], one, one)
    )
    four_ninths = b.num(4) / b.num(9)
    one_ninth = b.one / b.num(9)
    quark_naive = (
        transport_kernel_exact(b, mz, quarks["u"], four_ninths, three)
        + transport_kernel_exact(b, mz, quarks["d"], one_ninth, three)
        + transport_kernel_exact(b, mz, quarks["s"], one_ninth, three)
        + transport_kernel_exact(b, mz, quarks["c"], four_ninths, three)
        + transport_kernel_exact(b, mz, quarks["b"], one_ninth, three)
    )
    return lepton_total + screening * quark_naive


# ---------------------------------------------------------------------------
# Edge-sum tail bounds (interval backend only).
#
# For the SU(2) sum over n > N: term(n) = (n+1) e^{-t c2(n)} [ln(n+1)] with
# c2(n) = n(n+2)/4 and c2(n+1) - c2(n) = (2n+3)/4 increasing. For the SU(3)
# double sum outside the [0,N]^2 square, shells m = max(p,q) > N carry at most
# 2m+1 terms, dim <= (m+1)^3, ln dim <= 3 ln(m+1), c2 >= (m^2+3m)/3, and
# c2 <= m^2 + 2m within the shell. In both cases the successive-term (or
# successive-shell) ratio is bounded by a small combinatorial constant times
# e^{-t_lo * gap}; the blanket constants below (16 and 64) dominate every
# coefficient ratio that appears, including the extra c2 factor of the
# t-derivative series, so a single geometric ratio bounds the value tails and
# the derivative tails simultaneously.
# ---------------------------------------------------------------------------


def su2_tail_bounds(ivb: IvBackend, t: Any, cutoff: int) -> dict[str, Any]:
    if not (t.a > 0):
        raise CertificateError("su2 tail bound requires t > 0")
    n0 = cutoff + 1
    c2_first = ivb.num(n0 * (n0 + 2)) / ivb.num(4)
    gap = ivb.num(2 * n0 + 3) / ivb.num(4)
    q = ivb.num(16) * ivb.exp(-(ivb.thin(t.a) * gap))
    if not (q.b < ivb.num("0.5").a):
        raise CertificateError("su2 geometric tail ratio is not < 1/2")
    geom = ivb.one / (ivb.one - q)
    z1 = ivb.num(n0 + 1) * ivb.exp(-(t * c2_first))
    s1 = z1 * ivb.log(ivb.num(n0 + 1))
    return {
        "z": z1 * geom,
        "s": s1 * geom,
        "dz": c2_first * z1 * geom,
        "ds": c2_first * s1 * geom,
        "ratio_bound": q,
        "first_tail_index": n0,
    }


def su3_tail_bounds(ivb: IvBackend, t: Any, cutoff: int) -> dict[str, Any]:
    if not (t.a > 0):
        raise CertificateError("su3 tail bound requires t > 0")
    m0 = cutoff + 1
    c2_min = ivb.num(m0 * m0 + 3 * m0) / ivb.num(3)
    c2_max_shell = ivb.num(m0 * m0 + 2 * m0)
    gap = ivb.num(2 * m0 + 4) / ivb.num(3)
    q = ivb.num(64) * ivb.exp(-(ivb.thin(t.a) * gap))
    if not (q.b < ivb.num("0.5").a):
        raise CertificateError("su3 geometric tail ratio is not < 1/2")
    geom = ivb.one / (ivb.one - q)
    count = ivb.num(2 * m0 + 1)
    dim_bound = ivb.num((m0 + 1) ** 3)
    z1 = count * dim_bound * ivb.exp(-(t * c2_min))
    s1 = z1 * ivb.num(3) * ivb.log(ivb.num(m0 + 1))
    return {
        "z": z1 * geom,
        "s": s1 * geom,
        "dz": c2_max_shell * z1 * geom,
        "ds": c2_max_shell * s1 * geom,
        "ratio_bound": q,
        "first_tail_index": m0,
    }


def dual_tail(ivb: IvBackend, bounds: dict[str, Any], value_key: str, deriv_key: str, t_deriv: Any) -> tuple[Any, Any]:
    """One-sided dual tail contribution for an ellbar accumulator.

    Value tail is in [0, eps]; its d/dalpha is -(sum c2 w ...) * dt/dalpha,
    enclosed by [-eps_d, eps_d] * |t'| via the symmetric hull below.
    """
    eps = bounds[value_key].b
    eps_d = bounds[deriv_key].b
    value = ivb.hull(ivb.num(0).a, eps)
    deriv = ivb.hull(-eps_d, eps_d) * t_deriv
    return value, deriv


# ---------------------------------------------------------------------------
# Point solver (private mp context): mirrors the declared solver's scan-and-
# bisect root selection, then polishes with a bracketed Illinois iteration.
# Used only to produce candidates and inflation scales; all rigor comes from
# the interval verification.
# ---------------------------------------------------------------------------


class PointSolver:
    def __init__(self, mpb: MpBackend, su2_cutoff: int, su3_cutoff: int) -> None:
        self.b = mpb
        self.su2_cutoff = su2_cutoff
        self.su3_cutoff = su3_cutoff
        self.su2_terms = build_su2_terms(mpb, su2_cutoff)
        self.su3_terms = build_su3_terms(mpb, su3_cutoff)
        self.vectors = stage5_integer_vectors(mpb)
        roots = koide_roots(mpb)
        self.root_order = sorted(range(3), key=lambda i: roots[i])
        self.roots_ascending = [roots[i] for i in self.root_order]
        self._mz_seed: Any | None = None
        self._u_seed: Any | None = None
        self.rel_tol = self.b.num(10) ** (-(mpb.dps - 10))

    # -- m_Z closure ---------------------------------------------------------

    def _mz_bracket(self, f: Callable[[Any], Any], mu_u: Any) -> tuple[Any, Any, Any, Any]:
        b = self.b
        seed = self._mz_seed
        if seed is not None:
            widen = b.num("1e-6")
            while widen < b.num("0.5"):
                lo = seed * (b.one - widen)
                hi = seed * (b.one + widen)
                flo, fhi = f(lo), f(hi)
                if flo * fhi < 0:
                    return lo, hi, flo, fhi
                widen = widen * b.num(32)
        log_lo = b.log(mu_u) - b.num(50)
        log_hi = b.log(mu_u)
        prev_mu = None
        prev_val = None
        for index in range(260):
            frac = b.num(index) / b.num(259)
            mu = b.exp(log_lo + frac * (log_hi - log_lo))
            val = f(mu)
            if prev_val is not None and val * prev_val < 0:
                return prev_mu, mu, prev_val, val
            prev_mu, prev_val = mu, val
        raise CertificateError("could not bracket the m_Z fixed point (point solve)")

    def solve_mz(self, alpha_u: Any, p: Any, mu_u: Any) -> Any:
        b = self.b
        v = v_transmutation(b, alpha_u, p)

        def f(mu: Any) -> Any:
            return h_mz(b, mu, alpha_u, mu_u, v)

        lo, hi, flo, fhi = self._mz_bracket(f, mu_u)
        root = self._illinois(f, lo, hi, flo, fhi)
        self._mz_seed = root
        return root

    # -- alpha_U pixel closure ------------------------------------------------

    def pixel_residual(self, alpha_u: Any, p: Any, mu_u: Any) -> Any:
        b = self.b
        mz = self.solve_mz(alpha_u, p, mu_u)
        a2 = alpha_run(b, alpha_u, b.b_mssm[1], mz, mu_u)
        a3 = alpha_run(b, alpha_u, b.b_mssm[2], mz, mu_u)
        return (
            ellbar(b, self.su2_terms, b.four_pi_sq * a2)
            + ellbar(b, self.su3_terms, b.four_pi_sq * a3)
            - p / b.num(4)
        )

    def solve_alpha_u(self, p: Any) -> tuple[Any, Any]:
        b = self.b
        mu_u = mu_u_of_p(b, p)

        def f(u: Any) -> Any:
            return self.pixel_residual(u, p, mu_u)

        seed = self._u_seed
        bracket = None
        if seed is not None:
            widen = b.num("1e-7")
            while widen < b.num("0.5"):
                lo = seed * (b.one - widen)
                hi = seed * (b.one + widen)
                try:
                    flo, fhi = f(lo), f(hi)
                except CertificateError:
                    break
                if flo * fhi < 0:
                    bracket = (lo, hi, flo, fhi)
                    break
                widen = widen * b.num(32)
        if bracket is None:
            lo_scan = b.num("0.02")
            hi_scan = b.num("0.08")
            prev_u = None
            prev_val = None
            for index in range(41):
                u = lo_scan + (hi_scan - lo_scan) * b.num(index) / b.num(40)
                try:
                    val = f(u)
                except CertificateError:
                    continue
                if prev_val is not None and val * prev_val < 0:
                    bracket = (prev_u, u, prev_val, val)
                    break
                prev_u, prev_val = u, val
        if bracket is None:
            raise CertificateError("could not bracket alpha_U (point solve)")
        root = self._illinois(f, *bracket)
        self._u_seed = root
        return root, mu_u

    def _illinois(
        self, f: Callable[[Any], Any], lo: Any, hi: Any, flo: Any, fhi: Any
    ) -> Any:
        b = self.b
        for _ in range(240):
            if abs(hi - lo) < self.rel_tol * abs(hi):
                break
            mid = hi - fhi * (hi - lo) / (fhi - flo)
            if not (lo < mid < hi or hi < mid < lo):
                mid = (lo + hi) / b.two
            fmid = f(mid)
            if fmid == 0:
                return mid
            if fmid * fhi < 0:
                lo, flo = hi, fhi
                hi, fhi = mid, fmid
            else:
                hi, fhi = mid, fmid
                flo = flo / b.two
        return (lo + hi) / b.two

    # -- full inner map --------------------------------------------------------

    def inner(self, alpha: Any, mode: str) -> dict[str, Any]:
        b = self.b
        p = b.phi + alpha * b.sqrt_pi
        alpha_u, mu_u = self.solve_alpha_u(p)
        mz = self.solve_mz(alpha_u, p, mu_u)
        a1 = alpha_run(b, alpha_u, b.b_mssm[0], mz, mu_u)
        a2 = alpha_run(b, alpha_u, b.b_mssm[1], mz, mu_u)
        a3 = alpha_run(b, alpha_u, b.b_mssm[2], mz, mu_u)
        a0_inv = alpha_em_inv(b, a1, a2)
        v = v_transmutation(b, alpha_u, p)
        delta = structured_delta_alpha_inv(
            b, mz, v, a3, self.vectors, self.roots_ascending
        )
        alpha_inv = a0_inv + delta
        if mode == MODE_GAUGE_WIDTH:
            alpha_inv = alpha_inv + alpha_u
        elif mode != MODE_SOURCE:
            raise ValueError(f"unsupported mode: {mode}")
        return {
            "p": p,
            "mu_u": mu_u,
            "alpha_u": alpha_u,
            "mz": mz,
            "v": v,
            "a0_inv": a0_inv,
            "delta": delta,
            "alpha_inv": alpha_inv,
            "g": b.one / alpha_inv,
        }

    def fixed_point(self, mode: str, seed: str = "0.0073") -> Any:
        b = self.b
        alpha = b.num(seed)
        for _ in range(60):
            new_alpha = self.inner(alpha, mode)["g"]
            if abs(new_alpha - alpha) < self.rel_tol * abs(new_alpha):
                return new_alpha
            alpha = new_alpha
        raise CertificateError("point fixed-point iteration did not converge")


# ---------------------------------------------------------------------------
# Interval chain with verified implicit nodes.
# ---------------------------------------------------------------------------


class IntervalChain:
    def __init__(
        self,
        ivb: IvBackend,
        point: PointSolver,
        su2_cutoff: int,
        su3_cutoff: int,
    ) -> None:
        self.ivb = ivb
        self.dual = DualBackend(ivb)
        self.point = point
        self.su2_cutoff = su2_cutoff
        self.su3_cutoff = su3_cutoff
        self.su2_terms = build_su2_terms(ivb, su2_cutoff)
        self.su3_terms = build_su3_terms(ivb, su3_cutoff)
        self.vectors = point.vectors
        roots = koide_roots(ivb)
        self.roots_ascending = [roots[i] for i in point.root_order]
        for i in range(2):
            if not (self.roots_ascending[i].b < self.roots_ascending[i + 1].a):
                raise CertificateError("Koide root intervals are not disjoint")
        zero = ivb.zero
        self.su2_terms_dual = [
            (Dual(a, zero), Dual(b_, zero), Dual(c, zero)) for a, b_, c in self.su2_terms
        ]
        self.su3_terms_dual = [
            (Dual(a, zero), Dual(b_, zero), Dual(c, zero)) for a, b_, c in self.su3_terms
        ]
        # Finite-difference sensitivity scales (mp numbers); set by certify_mode
        # before any interval evaluation. Guidance for bracket inflation only:
        # every bracket is verified rigorously afterwards.
        self.scales: dict[str, Any] = {}

    # -- helpers ----------------------------------------------------------------

    def _to_iv(self, x_mp: Any) -> Any:
        return self.ivb.thin(x_mp)

    def _dual_const(self, x_iv: Any) -> Dual:
        return Dual(x_iv, self.ivb.zero)

    # -- verified m_Z ------------------------------------------------------------

    def verified_mz(
        self,
        alpha_u_iv: Any,
        p_iv: Any,
        mu_u_iv: Any,
        candidate_mp: Any,
    ) -> tuple[Any, Any, dict[str, Any]]:
        """Return (M, h_m enclosure, checks) with a verified sign-change bracket.

        The initial bracket pad covers the root variation over the actual widths
        of the alpha_U and P interval arguments (finite-difference estimates,
        times a safety factor); the sign checks below are the rigorous step.
        """
        ivb = self.ivb
        mpb = self.point.b
        v_iv = v_transmutation(ivb, alpha_u_iv, p_iv)
        w_u = _width_mp(mpb, alpha_u_iv)
        w_p = _width_mp(mpb, p_iv)
        pad = (
            self.scales["dmz_rel_du"] * w_u + self.scales["dmz_rel_dp"] * w_p
        ) * mpb.num(8) + self.scales["mz_floor"]
        for attempt in range(10):
            lo_mp = candidate_mp * (mpb.one - pad)
            hi_mp = candidate_mp * (mpb.one + pad)
            lo = self._to_iv(lo_mp)
            hi = self._to_iv(hi_mp)
            h_lo = h_mz(ivb, lo, alpha_u_iv, mu_u_iv, v_iv)
            h_hi = h_mz(ivb, hi, alpha_u_iv, mu_u_iv, v_iv)
            zero = ivb.num(0)
            if h_lo.a > zero.b and h_hi.b < zero.a:
                orientation = "decreasing"
            elif h_lo.b < zero.a and h_hi.a > zero.b:
                orientation = "increasing"
            else:
                pad = pad * mpb.num(8)
                continue
            m_box = ivb.hull(lo.a, hi.b)
            # h_m over the bracket (dual seeded in mu).
            d = self.dual
            h_m = h_mz(
                d,
                Dual(m_box, ivb.one),
                self._dual_const(alpha_u_iv),
                self._dual_const(mu_u_iv),
                self._dual_const(v_iv),
            ).d
            if ivb.contains_zero(h_m):
                pad = pad * mpb.num(8)
                continue
            checks = {
                "bracket_rel_half_width": str(pad),
                "orientation": orientation,
                "endpoint_signs_verified": True,
                "h_m_sign_definite": True,
                "attempts": attempt + 1,
            }
            return m_box, h_m, checks
        raise CertificateError("verified m_Z bracket failed after inflation attempts")

    def mz_dual(
        self,
        alpha_u_dual: Dual,
        p_dual: Dual,
        mu_u_dual: Dual,
        candidate_mp: Any,
    ) -> tuple[Dual, dict[str, Any]]:
        """Verified m_Z enclosure with IFT derivative for arbitrary dual seeds."""
        ivb = self.ivb
        m_box, h_m, checks = self.verified_mz(
            alpha_u_dual.x, p_dual.x, mu_u_dual.x, candidate_mp
        )
        d = self.dual
        v_dual = v_transmutation(d, alpha_u_dual, p_dual)
        slot = h_mz(d, Dual(m_box, ivb.zero), alpha_u_dual, mu_u_dual, v_dual).d
        m_prime = -(slot / h_m)
        return Dual(m_box, m_prime), checks

    # -- pixel closure residual over duals ---------------------------------------

    def residual_dual(
        self,
        u_dual: Dual,
        p_dual: Dual,
        mu_u_dual: Dual,
        mz_candidate_mp: Any,
    ) -> tuple[Dual, dict[str, Any]]:
        d = self.dual
        ivb = self.ivb
        mz, mz_checks = self.mz_dual(u_dual, p_dual, mu_u_dual, mz_candidate_mp)
        a2 = alpha_run(d, u_dual, d.b_mssm[1], mz, mu_u_dual)
        a3 = alpha_run(d, u_dual, d.b_mssm[2], mz, mu_u_dual)
        t2 = d.four_pi_sq * a2
        t3 = d.four_pi_sq * a3
        su2_tails = su2_tail_bounds(ivb, t2.x, self.su2_cutoff)
        su3_tails = su3_tail_bounds(ivb, t3.x, self.su3_cutoff)
        z2_tail, dz2_tail = dual_tail(ivb, su2_tails, "z", "dz", t2.d)
        s2_tail, ds2_tail = dual_tail(ivb, su2_tails, "s", "ds", t2.d)
        z3_tail, dz3_tail = dual_tail(ivb, su3_tails, "z", "dz", t3.d)
        s3_tail, ds3_tail = dual_tail(ivb, su3_tails, "s", "ds", t3.d)
        e2 = ellbar(
            d, self.su2_terms_dual, t2, (Dual(z2_tail, dz2_tail), Dual(s2_tail, ds2_tail))
        )
        e3 = ellbar(
            d, self.su3_terms_dual, t3, (Dual(z3_tail, dz3_tail), Dual(s3_tail, ds3_tail))
        )
        residual = e2 + e3 - p_dual / d.num(4)
        info = {
            "mz_checks": mz_checks,
            "su2_tail_hi": _iv_upper_str(su2_tails["s"]),
            "su3_tail_hi": _iv_upper_str(su3_tails["s"]),
        }
        return residual, info

    def residual_iv(
        self,
        u_iv: Any,
        p_iv: Any,
        mu_u_iv: Any,
        mz_candidate_mp: Any,
    ) -> Any:
        zero = self.ivb.zero
        res, _ = self.residual_dual(
            Dual(u_iv, zero),
            Dual(p_iv, zero),
            Dual(mu_u_iv, zero),
            mz_candidate_mp,
        )
        return res.x

    # -- verified alpha_U ---------------------------------------------------------

    def verified_alpha_u(
        self,
        p_dual: Dual,
        mu_u_dual: Dual,
        u_candidate_mp: Any,
        pad_abs: Any,
    ) -> tuple[Dual, dict[str, Any]]:
        ivb = self.ivb
        mpb = self.point.b
        p_iv = p_dual.x
        mu_u_iv = mu_u_dual.x
        p_mid = _mid_mp(mpb, p_iv)
        zero = ivb.num(0)
        pad = pad_abs
        for attempt in range(10):
            lo_mp = u_candidate_mp - pad
            hi_mp = u_candidate_mp + pad
            mz_lo = self.point.solve_mz(lo_mp, p_mid, mu_u_of_p(mpb, p_mid))
            mz_hi = self.point.solve_mz(hi_mp, p_mid, mu_u_of_p(mpb, p_mid))
            r_lo = self.residual_iv(self._to_iv(lo_mp), p_iv, mu_u_iv, mz_lo)
            r_hi = self.residual_iv(self._to_iv(hi_mp), p_iv, mu_u_iv, mz_hi)
            if r_lo.a > zero.b and r_hi.b < zero.a:
                orientation = "decreasing"
            elif r_lo.b < zero.a and r_hi.a > zero.b:
                orientation = "increasing"
            else:
                pad = pad * mpb.num(8)
                continue
            u_box = ivb.hull(self._to_iv(lo_mp).a, self._to_iv(hi_mp).b)
            mz_c = self.point.solve_mz(u_candidate_mp, p_mid, mu_u_of_p(mpb, p_mid))
            # dR/du over the box (u seeded 1, alpha-independent inputs seeded 0).
            r_u_dual, _ = self.residual_dual(
                Dual(u_box, ivb.one),
                Dual(p_iv, ivb.zero),
                Dual(mu_u_iv, ivb.zero),
                mz_c,
            )
            r_u = r_u_dual.d
            if ivb.contains_zero(r_u):
                pad = pad * mpb.num(8)
                continue
            # dR/dalpha at fixed u over the box (u seeded 0, P carrying d/dalpha).
            r_alpha_dual, info = self.residual_dual(
                Dual(u_box, ivb.zero), p_dual, mu_u_dual, mz_c
            )
            u_prime = -(r_alpha_dual.d / r_u)
            checks = {
                "bracket_abs_half_width": str(pad),
                "orientation": orientation,
                "endpoint_signs_verified": True,
                "R_u_sign_definite": True,
                "attempts": attempt + 1,
                "su2_tail_hi": info["su2_tail_hi"],
                "su3_tail_hi": info["su3_tail_hi"],
            }
            return Dual(u_box, u_prime), checks
        raise CertificateError("verified alpha_U bracket failed after inflation attempts")

    # -- full dual map -------------------------------------------------------------

    def g_dual(
        self,
        alpha_iv: Any,
        mode: str,
    ) -> tuple[Dual, dict[str, Any]]:
        ivb = self.ivb
        d = self.dual
        mpb = self.point.b

        alpha_dual = Dual(alpha_iv, ivb.one)
        p_dual = d.phi + alpha_dual * d.sqrt_pi
        mu_u_dual = mu_u_of_p(d, p_dual)

        w_p = _width_mp(mpb, p_dual.x)
        u_pad = self.scales["dudp_abs"] * w_p * mpb.num(8) + self.scales["u_floor"]

        p_mid = _mid_mp(mpb, p_dual.x)
        u_candidate, _ = self.point.solve_alpha_u(p_mid)

        u_dual, u_checks = self.verified_alpha_u(p_dual, mu_u_dual, u_candidate, u_pad)

        mz_c = self.point.solve_mz(u_candidate, p_mid, mu_u_of_p(mpb, p_mid))
        mz_dual, mz_checks = self.mz_dual(u_dual, p_dual, mu_u_dual, mz_c)

        a1 = alpha_run(d, u_dual, d.b_mssm[0], mz_dual, mu_u_dual)
        a2 = alpha_run(d, u_dual, d.b_mssm[1], mz_dual, mu_u_dual)
        a3 = alpha_run(d, u_dual, d.b_mssm[2], mz_dual, mu_u_dual)
        a0_inv = alpha_em_inv(d, a1, a2)
        v_dual = v_transmutation(d, u_dual, p_dual)
        roots_d = [self._dual_const(r) for r in self.roots_ascending]
        delta = structured_delta_alpha_inv(
            d, mz_dual, v_dual, a3, self.vectors, roots_d
        )
        alpha_inv_dual = a0_inv + delta
        if mode == MODE_GAUGE_WIDTH:
            alpha_inv_dual = alpha_inv_dual + u_dual
        elif mode != MODE_SOURCE:
            raise ValueError(f"unsupported mode: {mode}")
        g = d.one / alpha_inv_dual
        diagnostics = {
            "alpha_u_verification": u_checks,
            "mz_verification": mz_checks,
            "alpha_u_enclosure": _iv_pair(u_dual.x),
            "mz_enclosure": _iv_pair(mz_dual.x),
            "a0_inv_enclosure": _iv_pair(a0_inv.x),
            "delta_enclosure": _iv_pair(delta.x),
            "alpha_inv_enclosure": _iv_pair(alpha_inv_dual.x),
        }
        return g, diagnostics


# ---------------------------------------------------------------------------
# Serialization helpers.
# ---------------------------------------------------------------------------


def _iv_pair(x: Any) -> dict[str, str]:
    text = str(x).strip("[]")
    lo, hi = text.split(",")
    return {"lo": lo.strip(), "hi": hi.strip()}


def _iv_upper_str(x: Any) -> str:
    return _iv_pair(x)["hi"]


def _width_mp(mpb: MpBackend, x: Any) -> Any:
    return mpb.ctx.mpf(x.delta.b._mpi_[1])


def _mid_mp(mpb: MpBackend, x: Any) -> Any:
    return mpb.ctx.mpf(x.mid.a._mpi_[0])


def _mp_str(mpb: MpBackend, x: Any, dps: int | None = None) -> str:
    return mpb.ctx.nstr(x, dps or mpb.dps, strip_zeros=False)


# ---------------------------------------------------------------------------
# Certificate construction.
# ---------------------------------------------------------------------------


MAP_DEFINITIONS = {
    MODE_SOURCE: (
        "alpha -> 1/(alpha_em^-1(m_Z^2;P) + Delta_Th(P)) with "
        "P = phi + alpha*sqrt(pi); Delta_Th is the internal Stage-5 "
        "structured Thomson continuation with the exact one-loop fermion "
        "kernel and quark screening factor 1 - N_c*alpha_3(m_Z;P)/pi"
    ),
    MODE_GAUGE_WIDTH: (
        "alpha -> 1/(alpha_em^-1(m_Z^2;P) + Delta_Th(P) + alpha_U(P)) with "
        "P = phi + alpha*sqrt(pi); the closure-ledger CL-2 mixed map that "
        "adds the finite-screen unified gauge width alpha_U(P) to the "
        "inverse-alpha readout"
    ),
}


def certify_mode(
    mode: str,
    mp_dps: int,
    iv_dps: int,
    su2_cutoff: int,
    su3_cutoff: int,
    half_width: str,
    refine_passes: int,
) -> dict[str, Any]:
    mpb = MpBackend(mp_dps)
    point = PointSolver(mpb, su2_cutoff, su3_cutoff)
    ivb = IvBackend(iv_dps)
    chain = IntervalChain(ivb, point, su2_cutoff, su3_cutoff)

    alpha_star = point.fixed_point(mode)
    inner_star = point.inner(alpha_star, mode)

    # Finite-difference scale estimates (guidance only; not part of the proof).
    hp = mpb.num(10) ** (-14)
    p_star = inner_star["p"]
    u_plus, _ = point.solve_alpha_u(p_star + hp)
    u_minus, _ = point.solve_alpha_u(p_star - hp)
    dudp = abs(u_plus - u_minus) / (2 * hp)
    mu_u_star = inner_star["mu_u"]
    u_star = inner_star["alpha_u"]
    mz_star = inner_star["mz"]
    mzp = point.solve_mz(u_star, p_star + hp, mu_u_of_p(mpb, p_star + hp))
    mzm = point.solve_mz(u_star, p_star - hp, mu_u_of_p(mpb, p_star - hp))
    dmz_rel_dp = abs(mzp - mzm) / (2 * hp) / mz_star
    hu = mpb.num(10) ** (-14)
    mzpu = point.solve_mz(u_star + hu, p_star, mu_u_star)
    mzmu = point.solve_mz(u_star - hu, p_star, mu_u_star)
    dmz_rel_du = abs(mzpu - mzmu) / (2 * hu) / mz_star
    floor = mpb.num(10) ** (-(min(mp_dps, iv_dps) - 14))
    chain.scales = {
        "dudp_abs": dudp + mpb.num(10) ** (-6),
        "dmz_rel_dp": dmz_rel_dp + mpb.num(10) ** (-6),
        "dmz_rel_du": dmz_rel_du + mpb.num(10) ** (-6),
        "u_floor": floor,
        "mz_floor": floor,
    }

    r = mpb.num(half_width)
    lo = alpha_star - r
    hi = alpha_star + r
    interval = ivb.hull(ivb.thin(lo).a, ivb.thin(hi).b)

    # Derivative enclosure over the whole interval (dual chain on the box).
    g_box, diagnostics = chain.g_dual(interval, mode)
    g_prime = g_box.d
    lipschitz = ivb.sup_abs(g_prime)
    contraction = lipschitz < ivb.one.a
    if not contraction:
        raise CertificateError(f"{mode}: interval derivative bound is not < 1")

    # Mean-value (centered) form: g(I) subset g(mid) + g'(I) * (I - mid).
    # g is C^1 on I because the implicit-function denominators are verified
    # sign-definite along the box evaluation above.
    mid_iv = ivb.thin(alpha_star)
    g_mid, mid_diag = chain.g_dual(mid_iv, mode)
    g_val = g_mid.x + g_prime * (interval - mid_iv)

    inside = (interval.a < g_val.a) and (g_val.b < interval.b)
    if not inside:
        raise CertificateError(
            f"{mode}: centered-form g(I) is not contained in the interior of I "
            f"at half-width {half_width}"
        )

    # The fixed point lies in g(I) and in I; refine the enclosure by iterating
    # the centered form (contraction with factor <= L per pass).
    passes = [
        {
            "half_width": half_width,
            "g_interval": _iv_pair(g_val),
            "enclosure_width": _iv_upper_str(g_val.delta),
        }
    ]
    enclosure = g_val
    refine_diag = mid_diag
    for _ in range(refine_passes):
        candidate = enclosure
        mid_mp = _mid_mp(mpb, candidate)
        mid2 = ivb.thin(mid_mp)
        try:
            gp2, _ = chain.g_dual(candidate, mode)
            gm2, refine_diag_try = chain.g_dual(mid2, mode)
        except CertificateError:
            break
        l2 = ivb.sup_abs(gp2.d)
        if not (l2 < ivb.one.a):
            break
        g2 = gm2.x + gp2.d * (candidate - mid2)
        if not ((candidate.a < g2.a) and (g2.b < candidate.b)):
            break
        enclosure = g2
        refine_diag = refine_diag_try
        passes.append(
            {
                "half_width": _iv_upper_str(candidate.delta / ivb.two),
                "g_interval": _iv_pair(g2),
                "enclosure_width": _iv_upper_str(g2.delta),
            }
        )

    alpha_inv_enclosure = ivb.one / enclosure
    p_enclosure = ivb.phi + ivb.sqrt_pi * enclosure

    return {
        "mode": mode,
        "map_definition": MAP_DEFINITIONS[mode],
        "outer_equation": "P = phi + alpha*sqrt(pi)",
        "space": "alpha (fine-structure coupling); P and alpha_inv enclosures are derived images",
        "backend": (
            "mpmath.iv binary interval arithmetic with outward rounding on every "
            "elementary operation; private MPIntervalContext/MPContext instances; "
            "no Decimal-pad fallback used"
        ),
        "iv_dps": iv_dps,
        "point_dps": mp_dps,
        "su2_cutoff": su2_cutoff,
        "su3_cutoff": su3_cutoff,
        "padding_policy": (
            "no manual per-operation pads; implicit roots are enclosed by verified "
            "sign-change brackets with epsilon-inflation (factor-8 growth on failure), "
            "and the residual derivative is proven sign-definite on each bracket"
        ),
        "edge_sum_tail_bounds": {
            "included": True,
            "method": (
                "geometric majorant on the terms beyond the cutoff (SU(2)) and on the "
                "max(p,q) shells beyond the cutoff square (SU(3)); the tail enclosures "
                "are added one-sidedly to the heat-kernel sums, so every interval "
                "encloses both the declared finite-cutoff sums and their "
                "infinite-cutoff limits"
            ),
            "su2_ellbar_numerator_tail_upper_bound": diagnostics["alpha_u_verification"]["su2_tail_hi"],
            "su3_ellbar_numerator_tail_upper_bound": diagnostics["alpha_u_verification"]["su3_tail_hi"],
        },
        "uniqueness_interval_alpha": {
            **_iv_pair(interval),
            "half_width": half_width,
        },
        "banach": {
            "g_of_interval": _iv_pair(g_val),
            "g_maps_interval_into_interior": True,
            "gprime_enclosure": _iv_pair(g_prime),
            "lipschitz_bound": _iv_pair(lipschitz)["hi"],
            "contraction": True,
            "existence": True,
            "uniqueness_in_interval": True,
            "argument": (
                "direct Banach with a mean-value (centered) form: the interval "
                "evaluation proves sup |g'| <= L < 1 on I and "
                "g(I) subset g(mid(I)) + g'(I)*(I - mid(I)) subset interior(I); "
                "g is C^1 on I because the implicit-function denominators h_m and "
                "R_u are interval-verified sign-definite over the box evaluation"
            ),
        },
        "refinement_passes": passes,
        "certified_enclosure": {
            "alpha": {**_iv_pair(enclosure), "width": _iv_upper_str(enclosure.delta)},
            "alpha_inv": {
                **_iv_pair(alpha_inv_enclosure),
                "width": _iv_upper_str(alpha_inv_enclosure.delta),
            },
            "P": {**_iv_pair(p_enclosure), "width": _iv_upper_str(p_enclosure.delta)},
        },
        "inner_root_certificates": {
            "alpha_U_pixel_closure": refine_diag["alpha_u_verification"],
            "m_Z_tree_closure": refine_diag["mz_verification"],
            "root_selection_scope": (
                "the certified brackets contain the roots selected by the declared "
                "solver's scan-and-bisect procedure (the point solver mirrors the "
                "declared scan and its roots lie inside the verified brackets); "
                "uniqueness inside each bracket follows from the sign-definite "
                "residual derivative; global root uniqueness over the full scan "
                "windows is not claimed"
            ),
        },
        "interval_diagnostics": {
            "alpha_u": refine_diag["alpha_u_enclosure"],
            "m_z": refine_diag["mz_enclosure"],
            "a0_inv": refine_diag["a0_inv_enclosure"],
            "delta_th": refine_diag["delta_enclosure"],
            "alpha_inv_map_value": refine_diag["alpha_inv_enclosure"],
        },
        "fixed_point_point_estimate_display_only": {
            "alpha": _mp_str(mpb, alpha_star),
            "alpha_inv": _mp_str(mpb, mpb.one / alpha_star),
            "P": _mp_str(mpb, mpb.phi + alpha_star * mpb.sqrt_pi),
        },
    }


def build_certificate(
    mp_dps: int = 60,
    iv_dps: int = 60,
    su2_cutoff: int = 120,
    su3_cutoff: int = 90,
    half_width: str = "0.000004",
    refine_passes: int = 4,
    modes: tuple[str, ...] = MODES,
    decimal_crosscheck_precision: int = 0,
) -> dict[str, Any]:
    mode_blocks = {}
    for mode in modes:
        mode_blocks[mode] = certify_mode(
            mode, mp_dps, iv_dps, su2_cutoff, su3_cutoff, half_width, refine_passes
        )

    artifact: dict[str, Any] = {
        "artifact": "oph_p_interval_contraction_certificate",
        "date": ARTIFACT_DATE,
        "claim_status": "interval_contraction_certificate_for_declared_closure_map",
        "claim_boundary": (
            "Existence and uniqueness of the closure-map fixed point are certified by "
            "interval arithmetic on the stated alpha interval, for the declared "
            "numerical map at the stated representation cutoffs (with the edge-sum "
            "truncation tails bounded, the certificate also covers the infinite-cutoff "
            "edge sums). The declared one-loop RG/matching conventions, the tree-level "
            "m_Z closure, the Stage-5 continuation masses, and the exact one-loop "
            "kernel are certified as declared, not as physical endpoint theorems. "
            "This is not an exact fine-structure derivation; the stage-3 landing "
            "verdict of closure row CL-1 is unchanged."
        ),
        "protocol_stage": (
            "stage 2 (contraction certificate) of the basin-then-contract protocol "
            "for the P coordinate; see the P-closure issue (#545)"
        ),
        "scope": (
            "certificate applies to the declared numerical closure map at the stated "
            "su2_cutoff/su3_cutoff (tail bounds extend it to the infinite edge sums) "
            "and to the stated readout mode; it does not certify any relation to the "
            "measured fine-structure constant"
        ),
        "promotion_allowed": False,
        "exact_alpha_promoted": False,
        "consumer_policy": {
            "may_feed_live_particle_predictions": False,
            "may_feed_compare_or_audit_surfaces": True,
            "hidden_external_alpha_allowed": False,
            "default_thomson_endpoint_allowed": False,
        },
        "modes": mode_blocks,
    }

    if decimal_crosscheck_precision > 0:
        artifact["decimal_chain_crosscheck"] = decimal_crosscheck(
            mp_dps=mp_dps,
            su2_cutoff=su2_cutoff,
            su3_cutoff=su3_cutoff,
            decimal_precision=decimal_crosscheck_precision,
        )
    return artifact


def decimal_crosscheck(
    mp_dps: int,
    su2_cutoff: int,
    su3_cutoff: int,
    decimal_precision: int,
    alpha_probe: str = "0.0073",
) -> dict[str, Any]:
    """Compare the mp midpoint chain against paper_math's Decimal chain.

    The Decimal chain resolves its internal bisections to roughly
    ``2^-(precision+8)`` of the initial brackets, so ``decimal_precision``
    controls the number of matching digits available.
    """
    from decimal import Decimal

    from paper_math import PaperMathContext

    mpb = MpBackend(mp_dps)
    point = PointSolver(mpb, su2_cutoff, su3_cutoff)
    mine = point.inner(mpb.num(alpha_probe), MODE_SOURCE)

    ctx = PaperMathContext(
        precision=decimal_precision, su2_cutoff=su2_cutoff, su3_cutoff=su3_cutoff
    )
    probe = ctx.evaluate_alpha_fixed_point(Decimal(alpha_probe), MODE_SOURCE)
    theirs = mpb.num(str(probe.inner_alpha))
    diff = abs(mine["g"] - theirs)
    rel = diff / abs(theirs)
    if rel > 0:
        digits = int(-mpb.ctx.log(rel, 10))
    else:
        digits = mp_dps
    return {
        "alpha_probe": alpha_probe,
        "mode": MODE_SOURCE,
        "decimal_precision": decimal_precision,
        "mp_inner_alpha": _mp_str(mpb, mine["g"]),
        "decimal_inner_alpha": str(probe.inner_alpha),
        "abs_difference": _mp_str(mpb, diff, 5),
        "agreement_digits": digits,
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Emit the interval-grade Banach contraction certificate for the OPH "
            "P/alpha closure map (stage-2 contraction certificate)."
        )
    )
    parser.add_argument("--mp-dps", type=int, default=60, help="Point-solver working digits.")
    parser.add_argument("--iv-dps", type=int, default=60, help="Interval arithmetic working digits.")
    parser.add_argument("--su2-cutoff", type=int, default=120, help="SU(2) representation cutoff.")
    parser.add_argument("--su3-cutoff", type=int, default=90, help="SU(3) representation cutoff.")
    parser.add_argument(
        "--half-width",
        default="0.000004",
        help="Half-width of the certified alpha uniqueness interval.",
    )
    parser.add_argument(
        "--refine-passes",
        type=int,
        default=4,
        help="Enclosure refinement passes after the primary certificate.",
    )
    parser.add_argument(
        "--mode",
        choices=MODES + ("both",),
        default="both",
        help="Which readout map(s) to certify.",
    )
    parser.add_argument(
        "--decimal-crosscheck-precision",
        type=int,
        default=120,
        help=(
            "Decimal precision for the paper_math midpoint cross-check; "
            "0 disables the cross-check."
        ),
    )
    parser.add_argument("--output", default=str(DEFAULT_OUT), help="Output JSON path.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    modes = MODES if args.mode == "both" else (args.mode,)
    t0 = time.time()
    artifact = build_certificate(
        mp_dps=args.mp_dps,
        iv_dps=args.iv_dps,
        su2_cutoff=args.su2_cutoff,
        su3_cutoff=args.su3_cutoff,
        half_width=args.half_width,
        refine_passes=args.refine_passes,
        modes=modes,
        decimal_crosscheck_precision=args.decimal_crosscheck_precision,
    )
    wall_seconds = round(time.time() - t0, 2)
    text = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text, encoding="utf-8")
    for mode in modes:
        block = artifact["modes"][mode]
        enclosure = block["certified_enclosure"]["alpha_inv"]
        print(f"{mode}:")
        print(f"  L (Lipschitz bound)   = {block['banach']['lipschitz_bound']}")
        print(f"  alpha_inv enclosure   = [{enclosure['lo']}, {enclosure['hi']}]")
        print(f"  enclosure width       = {enclosure['width']}")
    print(f"wall seconds: {wall_seconds} (not recorded in the artifact)")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
