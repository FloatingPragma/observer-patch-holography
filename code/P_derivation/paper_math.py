#!/usr/bin/env python3
"""Paper-math implementation of the OPH P/alpha closure experiment.

This module encodes the D10 paper equations directly, without importing the
public-facing ``code/particles`` package. The main entry point is
``build_report()``, which:

1. Solves the paper D10 pixel-closure equations for a given ``P``.
2. Extracts the source-locked electromagnetic anchor
   ``a0(P) = alpha_em^-1(m_Z^2; P)``.
3. Closes Alex's equation ``P = phi + alpha * sqrt(pi)`` using a documented
   low-energy readout choice.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from decimal import Decimal, localcontext
from functools import lru_cache
from typing import Any


def _dec(value: Decimal | int | str | float) -> Decimal:
    if isinstance(value, Decimal):
        return value
    if isinstance(value, int):
        return Decimal(value)
    return Decimal(str(value))


@lru_cache(maxsize=None)
def decimal_pi(precision: int) -> Decimal:
    """Compute pi with the Chudnovsky series."""
    with localcontext() as ctx:
        ctx.prec = precision + 8
        total = Decimal(0)
        m = 1
        l = 13591409
        x = 1
        k = 6
        c = Decimal(426880) * Decimal(10005).sqrt()
        threshold = Decimal(10) ** (-(precision + 4))
        for _ in range(max(3, precision // 13 + 2)):
            term = Decimal(m * l) / Decimal(x)
            total += term
            if abs(term) < threshold:
                break
            m = (m * (k**3 - 16 * k)) // ((_ + 1) ** 3)
            l += 545140134
            x *= -262537412640768000
            k += 12
        return +(c / total)


@dataclass(frozen=True)
class D10Point:
    p: Decimal
    n_c: int
    mu_u: Decimal
    alpha_u: Decimal
    mz_run: Decimal
    v: Decimal
    alpha1_mz: Decimal
    alpha2_mz: Decimal
    alpha3_mz: Decimal
    alpha_y_mz: Decimal
    alpha_em_mz: Decimal
    alpha_em_inv_mz: Decimal
    sin2w_mz: Decimal


@dataclass(frozen=True)
class ClosureStep:
    iteration: int
    alpha_probe: Decimal
    alpha_probe_inv: Decimal
    outer_p: Decimal
    inner_alpha: Decimal
    inner_alpha_inv: Decimal
    source_anchor_alpha_inv: Decimal
    residual_alpha: Decimal


@dataclass(frozen=True)
class ClosureResult:
    mode: str
    precision: int
    su2_cutoff: int
    su3_cutoff: int
    phi: Decimal
    sqrt_pi: Decimal
    alpha: Decimal
    alpha_inv: Decimal
    p: Decimal
    god_equation_residual: Decimal
    alpha_fixed_point_residual: Decimal
    source_anchor_alpha_inv: Decimal
    structured_running: dict[str, Any] | None
    d10: D10Point
    iterations: list[ClosureStep]


@dataclass(frozen=True)
class FixedPointProbe:
    alpha: Decimal
    alpha_inv: Decimal
    p: Decimal
    inner_alpha: Decimal
    inner_alpha_inv: Decimal
    residual_alpha: Decimal


CODATA_2022_ALPHA_INV = Decimal("137.035999177")
CODATA_2022_ALPHA_INV_STD_UNCERTAINTY = Decimal("0.000000021")


class PaperMathContext:
    def __init__(
        self,
        precision: int = 40,
        su2_cutoff: int = 120,
        su3_cutoff: int = 90,
        n_c: int = 3,
        n_g: int = 3,
    ):
        self.precision = precision
        self.work_precision = precision + 12
        self.su2_cutoff = su2_cutoff
        self.su3_cutoff = su3_cutoff
        self.n_c = n_c
        self.n_g = n_g

        with localcontext() as ctx:
            ctx.prec = self.work_precision
            self.pi = decimal_pi(self.work_precision)
            self.sqrt_pi = self.pi.sqrt()
            self.sqrt_two = Decimal(2).sqrt()
            self.two_pi = Decimal(2) * self.pi
            self.four_pi = Decimal(4) * self.pi
            self.four_pi_sq = Decimal(4) * self.pi * self.pi
            self.ln10 = Decimal(10).ln()
            self.one = Decimal(1)
            self.two = Decimal(2)
            self.five = Decimal(5)
            self.six = Decimal(6)
            self.planck_energy = self.one
            self.b_mssm = (Decimal(33) / Decimal(5), self.one, -Decimal(3))
            self.beta_ew = Decimal(n_c + 1)
            self.phi = (self.one + self.five.sqrt()) / self.two
            self._su2_terms = self._build_su2_terms()
            self._su3_terms = self._build_su3_terms()
            self.stage5_vectors = self._derive_stage5_integer_vectors()

    def _build_su2_terms(self) -> tuple[tuple[Decimal, Decimal, Decimal], ...]:
        terms: list[tuple[Decimal, Decimal, Decimal]] = []
        for n in range(self.su2_cutoff + 1):
            j = Decimal(n) / Decimal(2)
            dim = Decimal(n + 1)
            c2 = j * (j + self.one)
            terms.append((dim, c2, dim.ln()))
        return tuple(terms)

    def _build_su3_terms(self) -> tuple[tuple[Decimal, Decimal, Decimal], ...]:
        terms: list[tuple[Decimal, Decimal, Decimal]] = []
        for p in range(self.su3_cutoff + 1):
            for q in range(self.su3_cutoff + 1):
                dim_int = ((p + 1) * (q + 1) * (p + q + 2)) // 2
                dim = Decimal(dim_int)
                c2 = Decimal(p * p + q * q + p * q + 3 * p + 3 * q) / Decimal(3)
                terms.append((dim, c2, dim.ln()))
        return tuple(terms)

    def ellbar_su2(self, t: Decimal) -> Decimal:
        with localcontext() as ctx:
            ctx.prec = self.work_precision
            weights: list[tuple[Decimal, Decimal]] = []
            z = Decimal(0)
            for dim, c2, ln_dim in self._su2_terms:
                w = dim * (-(t * c2)).exp()
                weights.append((w, ln_dim))
                z += w
            total = Decimal(0)
            for w, ln_dim in weights:
                total += (w / z) * ln_dim
            return +total

    def ellbar_su3(self, t: Decimal) -> Decimal:
        with localcontext() as ctx:
            ctx.prec = self.work_precision
            weights: list[tuple[Decimal, Decimal]] = []
            z = Decimal(0)
            for dim, c2, ln_dim in self._su3_terms:
                w = dim * (-(t * c2)).exp()
                weights.append((w, ln_dim))
                z += w
            total = Decimal(0)
            for w, ln_dim in weights:
                total += (w / z) * ln_dim
            return +total

    def unification_scale_gev(self, p: Decimal) -> Decimal:
        with localcontext() as ctx:
            ctx.prec = self.work_precision
            return +(self.planck_energy * (-self.two_pi).exp() * (p.ln() / Decimal(6)).exp())

    def e_cell(self, p: Decimal) -> Decimal:
        with localcontext() as ctx:
            ctx.prec = self.work_precision
            return +(self.planck_energy / p.sqrt())

    def v_from_transmutation(self, alpha_u: Decimal, p: Decimal) -> Decimal:
        with localcontext() as ctx:
            ctx.prec = self.work_precision
            exponent = -(self.two_pi / (self.beta_ew * alpha_u))
            return +(self.e_cell(p) * exponent.exp())

    def alpha_run_1loop(self, alpha_u: Decimal, b: Decimal, mu: Decimal, mu_u: Decimal) -> Decimal:
        with localcontext() as ctx:
            ctx.prec = self.work_precision
            inv = (self.one / alpha_u) + (b / self.two_pi) * (mu_u / mu).ln()
            return +(self.one / inv)

    def run_alphas_from_unification(self, alpha_u: Decimal, mu: Decimal, mu_u: Decimal) -> tuple[Decimal, Decimal, Decimal]:
        return tuple(self.alpha_run_1loop(alpha_u, b_i, mu, mu_u) for b_i in self.b_mssm)  # type: ignore[return-value]

    def alpha_em_from_alpha1_alpha2(self, alpha1: Decimal, alpha2: Decimal) -> tuple[Decimal, Decimal, Decimal]:
        with localcontext() as ctx:
            ctx.prec = self.work_precision
            alpha_y = Decimal(3) * alpha1 / Decimal(5)
            alpha_em = self.one / ((self.one / alpha2) + (self.one / alpha_y))
            sin2w = alpha_em / alpha2
            return +alpha_y, +alpha_em, +sin2w

    def mz_tree_from_v_and_couplings(self, v: Decimal, alpha1: Decimal, alpha2: Decimal) -> Decimal:
        with localcontext() as ctx:
            ctx.prec = self.work_precision
            alpha_y = Decimal(3) * alpha1 / Decimal(5)
            g2_sq = self.four_pi * alpha2
            gy_sq = self.four_pi * alpha_y
            return +(v / self.two * (g2_sq + gy_sq).sqrt())

    def decimal_cos(self, x: Decimal) -> Decimal:
        with localcontext() as ctx:
            ctx.prec = self.work_precision
            x = x % self.two_pi
            if x > self.pi:
                x -= self.two_pi
            term = self.one
            total = self.one
            n = 0
            threshold = Decimal(10) ** (-(self.work_precision - 2))
            xx = x * x
            while True:
                n += 1
                term *= -(xx) / Decimal((2 * n - 1) * (2 * n))
                total += term
                if abs(term) < threshold:
                    return +total

    def _koide_roots(self, delta: Decimal) -> tuple[Decimal, Decimal, Decimal]:
        with localcontext() as ctx:
            ctx.prec = self.work_precision
            roots = []
            for k in range(3):
                angle = delta + self.two_pi * Decimal(k) / Decimal(3)
                roots.append(self.one + self.sqrt_two * self.decimal_cos(angle))
            roots.sort()
            return tuple(+root for root in roots)

    def _derive_stage5_integer_vectors(self) -> dict[str, Any]:
        with localcontext() as ctx:
            ctx.prec = self.work_precision
            epsilon = self.one / self.six
            delta = self.beta_ew / (self.two * Decimal(self.n_c) * Decimal(self.n_g))
            roots = self._koide_roots(delta)
            roots_sq = [root * root for root in roots]
            n_tau = self.n_g
            n_mu = self.n_g + 1
            best_n_e: int | None = None
            best_residual: Decimal | None = None
            for n_e in range(n_mu + 1, n_mu + 8):
                exponents = (n_e, n_mu, n_tau)
                errors: list[Decimal] = []
                for i in range(3):
                    for j in range(i + 1, 3):
                        lhs = roots_sq[i] / roots_sq[j]
                        rhs = epsilon ** (exponents[i] - exponents[j])
                        errors.append(abs((lhs / rhs).ln()))
                residual = max(errors)
                if best_residual is None or residual < best_residual:
                    best_residual = residual
                    best_n_e = n_e
            if best_n_e is None:
                raise RuntimeError("Could not derive Stage-5 charged-lepton exponents.")
            return {
                "epsilon": +epsilon,
                "delta": +delta,
                "roots": tuple(+root for root in roots),
                "n_u": (2 * self.n_c, self.n_c, 0),
                "n_d": (2 * self.n_c, self.n_c + 1, self.n_c - 1),
                "n_e": (best_n_e, n_mu, n_tau),
            }

    def diagonal_quark_masses(self, v: Decimal) -> dict[str, Decimal]:
        with localcontext() as ctx:
            ctx.prec = self.work_precision
            pref = v / self.sqrt_two
            epsilon = self.stage5_vectors["epsilon"]
            n_u = self.stage5_vectors["n_u"]
            n_d = self.stage5_vectors["n_d"]
            return {
                "u": +(pref * (epsilon ** n_u[0])),
                "c": +(pref * (epsilon ** n_u[1])),
                "t": +(pref * (epsilon ** n_u[2])),
                "d": +(pref * (epsilon ** n_d[0])),
                "s": +(pref * (epsilon ** n_d[1])),
                "b": +(pref * (epsilon ** n_d[2])),
            }

    def charged_lepton_masses(self, v: Decimal) -> dict[str, Decimal]:
        with localcontext() as ctx:
            ctx.prec = self.work_precision
            roots_sq = [root * root for root in self.stage5_vectors["roots"]]
            n_e = self.stage5_vectors["n_e"]
            log_gm_c = Decimal(0)
            for index, exponent in enumerate(n_e):
                numerator = roots_sq[index] * self.sqrt_two * (self.six ** exponent)
                log_gm_c += (numerator / v).ln()
            log_gm_c /= Decimal(3)
            s0 = (-log_gm_c).exp()
            scale = s0 * (self.two.ln() / Decimal(6)).exp()
            return {
                "e": +(scale * roots_sq[0]),
                "mu": +(scale * roots_sq[1]),
                "tau": +(scale * roots_sq[2]),
            }

    def fermion_transport_kernel_asymptotic(
        self,
        q_scale: Decimal,
        mass: Decimal,
        charge_squared: Decimal,
        multiplicity: Decimal,
    ) -> Decimal:
        with localcontext() as ctx:
            ctx.prec = self.work_precision
            return +(
                multiplicity
                * charge_squared
                / (Decimal(3) * self.pi)
                * (((q_scale * q_scale) / (mass * mass)).ln() - Decimal(5) / Decimal(3))
            )

    def fermion_transport_kernel_exact(
        self,
        q_scale: Decimal,
        mass: Decimal,
        charge_squared: Decimal,
        multiplicity: Decimal,
    ) -> Decimal:
        with localcontext() as ctx:
            ctx.prec = self.work_precision
            z = (q_scale * q_scale) / (mass * mass)
            threshold = Decimal(10) ** (-(self.precision + 4))

            def simpson(n_intervals: int) -> Decimal:
                h = self.one / Decimal(n_intervals)
                total = Decimal(0)
                for index in range(n_intervals + 1):
                    x = Decimal(index) * h
                    t = x * (self.one - x)
                    value = t * (self.one + z * t).ln()
                    if index == 0 or index == n_intervals:
                        weight = self.one
                    elif index % 2 == 1:
                        weight = Decimal(4)
                    else:
                        weight = self.two
                    total += weight * value
                integral = total * h / Decimal(3)
                return +(Decimal(2) * multiplicity * charge_squared / self.pi * integral)

            n_intervals = 32
            previous = simpson(n_intervals)
            while n_intervals < 4096:
                n_intervals *= 2
                current = simpson(n_intervals)
                if abs(current - previous) < threshold:
                    return +current
                previous = current
            return +previous

    def _structured_thomson_running_from_masses(
        self,
        d10: D10Point,
        *,
        quarks: dict[str, Decimal],
        leptons: dict[str, Decimal],
        mass_source: str,
        kernel: str,
    ) -> dict[str, Any]:
        with localcontext() as ctx:
            ctx.prec = self.work_precision
            screening = self.one - Decimal(self.n_c) * d10.alpha3_mz / self.pi

            fermion_term = (
                self.fermion_transport_kernel_exact
                if kernel == "exact_1loop"
                else self.fermion_transport_kernel_asymptotic
            )

            lepton_total = (
                fermion_term(d10.mz_run, leptons["e"], self.one, self.one)
                + fermion_term(d10.mz_run, leptons["mu"], self.one, self.one)
                + fermion_term(d10.mz_run, leptons["tau"], self.one, self.one)
            )
            quark_naive = (
                fermion_term(d10.mz_run, quarks["u"], Decimal(4) / Decimal(9), Decimal(3))
                + fermion_term(d10.mz_run, quarks["d"], Decimal(1) / Decimal(9), Decimal(3))
                + fermion_term(d10.mz_run, quarks["s"], Decimal(1) / Decimal(9), Decimal(3))
                + fermion_term(d10.mz_run, quarks["c"], Decimal(4) / Decimal(9), Decimal(3))
                + fermion_term(d10.mz_run, quarks["b"], Decimal(1) / Decimal(9), Decimal(3))
            )
            quark_screened = screening * quark_naive
            total = lepton_total + quark_screened
            return {
                "transport_kernel": kernel,
                "epsilon_z6": self.stage5_vectors["epsilon"],
                "koide_delta": self.stage5_vectors["delta"],
                "integer_vectors": {
                    "n_u": self.stage5_vectors["n_u"],
                    "n_d": self.stage5_vectors["n_d"],
                    "n_e": self.stage5_vectors["n_e"],
                },
                "masses_gev": {
                    "u": quarks["u"],
                    "d": quarks["d"],
                    "s": quarks["s"],
                    "c": quarks["c"],
                    "b": quarks["b"],
                    "e": leptons["e"],
                    "mu": leptons["mu"],
                    "tau": leptons["tau"],
                },
                "quark_screening_factor": +screening,
                "mass_source": mass_source,
                "lepton_delta_alpha_inv": +lepton_total,
                "quark_delta_alpha_inv_naive": +quark_naive,
                "quark_delta_alpha_inv_screened": +quark_screened,
                "total_delta_alpha_inv": +total,
            }

    def structured_thomson_running(self, d10: D10Point) -> dict[str, Any]:
        quarks = self.diagonal_quark_masses(d10.v)
        leptons = self.charged_lepton_masses(d10.v)
        exact = self._structured_thomson_running_from_masses(
            d10,
            quarks=quarks,
            leptons=leptons,
            mass_source="internal_stage5_continuation",
            kernel="exact_1loop",
        )
        asymptotic = self._structured_thomson_running_from_masses(
            d10,
            quarks=quarks,
            leptons=leptons,
            mass_source="internal_stage5_continuation",
            kernel="asymptotic_log",
        )
        exact["legacy_asymptotic_comparison"] = {
            "lepton_delta_alpha_inv": asymptotic["lepton_delta_alpha_inv"],
            "quark_delta_alpha_inv_naive": asymptotic["quark_delta_alpha_inv_naive"],
            "quark_delta_alpha_inv_screened": asymptotic["quark_delta_alpha_inv_screened"],
            "total_delta_alpha_inv": asymptotic["total_delta_alpha_inv"],
        }
        exact["kernel_upgrade_delta_alpha_inv"] = +(
            exact["total_delta_alpha_inv"] - asymptotic["total_delta_alpha_inv"]
        )
        return exact

    def structured_thomson_running_asymptotic(self, d10: D10Point) -> dict[str, Any]:
        quarks = self.diagonal_quark_masses(d10.v)
        leptons = self.charged_lepton_masses(d10.v)
        return self._structured_thomson_running_from_masses(
            d10,
            quarks=quarks,
            leptons=leptons,
            mass_source="internal_stage5_continuation",
            kernel="asymptotic_log",
        )

    def solve_mz_fixed_point_tree(self, alpha_u: Decimal, p: Decimal, mu_u: Decimal) -> tuple[Decimal, Decimal, Decimal, Decimal, Decimal]:
        with localcontext() as ctx:
            ctx.prec = self.work_precision
            v_ev = self.v_from_transmutation(alpha_u, p)
            bisection_steps = max(32, self.precision + 8)

            def f(mu: Decimal) -> Decimal:
                alpha1, alpha2, alpha3 = self.run_alphas_from_unification(alpha_u, mu, mu_u)
                return self.mz_tree_from_v_and_couplings(v_ev, alpha1, alpha2) - mu

            previous_mu: Decimal | None = None
            previous_value: Decimal | None = None
            log_lo = mu_u.ln() - Decimal(50)
            log_hi = mu_u.ln()
            for index in range(260):
                frac = Decimal(index) / Decimal(259)
                mu = +(log_lo + frac * (log_hi - log_lo)).exp()
                value = f(mu)
                if previous_value is not None and value * previous_value < 0:
                    lo = previous_mu
                    hi = mu
                    lo_value = previous_value
                    for _ in range(bisection_steps):
                        mid = (lo * hi).sqrt()
                        mid_value = f(mid)
                        if lo_value * mid_value > 0:
                            lo = mid
                            lo_value = mid_value
                        else:
                            hi = mid
                    mz_run = +(lo + hi) / self.two
                    alpha1, alpha2, alpha3 = self.run_alphas_from_unification(alpha_u, mz_run, mu_u)
                    return +mz_run, +v_ev, +alpha1, +alpha2, +alpha3
                previous_mu = mu
                previous_value = value
        raise RuntimeError("Could not bracket the m_Z fixed point.")

    def pixel_residual(self, alpha_u: Decimal, p: Decimal, mu_u: Decimal) -> tuple[Decimal, tuple[Decimal, Decimal, Decimal, Decimal, Decimal]]:
        with localcontext() as ctx:
            ctx.prec = self.work_precision
            mz_run, v_ev, alpha1, alpha2, alpha3 = self.solve_mz_fixed_point_tree(alpha_u, p, mu_u)
            t2 = self.four_pi_sq * alpha2
            t3 = self.four_pi_sq * alpha3
            residual = self.ellbar_su2(t2) + self.ellbar_su3(t3) - p / Decimal(4)
            return +residual, (+mz_run, +v_ev, +alpha1, +alpha2, +alpha3)

    def solve_alpha_u_from_p(self, p: Decimal) -> tuple[Decimal, tuple[Decimal, Decimal, Decimal, Decimal, Decimal], Decimal]:
        with localcontext() as ctx:
            ctx.prec = self.work_precision
            mu_u = self.unification_scale_gev(p)
            lo = Decimal("0.02")
            hi = Decimal("0.08")
            scan_points = 40
            bisection_steps = max(32, self.precision + 8)
            previous_alpha: Decimal | None = None
            previous_residual: Decimal | None = None
            previous_rep: tuple[Decimal, Decimal, Decimal, Decimal, Decimal] | None = None
            bracket: tuple[Decimal, Decimal] | None = None
            bracket_data: tuple[
                tuple[Decimal, tuple[Decimal, Decimal, Decimal, Decimal, Decimal]],
                tuple[Decimal, tuple[Decimal, Decimal, Decimal, Decimal, Decimal]],
            ] | None = None
            for index in range(scan_points + 1):
                frac = Decimal(index) / Decimal(scan_points)
                alpha_u = lo + (hi - lo) * frac
                try:
                    residual, rep = self.pixel_residual(alpha_u, p, mu_u)
                except RuntimeError:
                    continue
                if previous_residual is not None and residual * previous_residual < 0:
                    bracket = (previous_alpha, alpha_u)
                    bracket_data = ((previous_residual, previous_rep), (residual, rep))
                    break
                previous_alpha = alpha_u
                previous_residual = residual
                previous_rep = rep
            if bracket is None:
                raise RuntimeError("Could not bracket alpha_U in the D10 pixel solve.")

            lo_alpha, hi_alpha = bracket
            if bracket_data is None:
                raise RuntimeError("Missing bracket data for alpha_U.")
            (lo_residual, lo_rep), (hi_residual, hi_rep) = bracket_data
            final_rep = hi_rep
            for _ in range(bisection_steps):
                mid = (lo_alpha + hi_alpha) / self.two
                mid_residual, mid_rep = self.pixel_residual(mid, p, mu_u)
                final_rep = mid_rep
                if lo_residual * mid_residual <= 0:
                    hi_alpha = mid
                    hi_residual = mid_residual
                else:
                    lo_alpha = mid
                    lo_residual = mid_residual
            alpha_u = +(lo_alpha + hi_alpha) / self.two
            return +alpha_u, final_rep, +mu_u

    def build_d10_from_p(self, p: Decimal | int | str | float) -> D10Point:
        with localcontext() as ctx:
            ctx.prec = self.work_precision
            p_dec = _dec(p)
            alpha_u, rep, mu_u = self.solve_alpha_u_from_p(p_dec)
            mz_run, v_ev, alpha1, alpha2, alpha3 = rep
            alpha_y, alpha_em, sin2w = self.alpha_em_from_alpha1_alpha2(alpha1, alpha2)
            return D10Point(
                p=+p_dec,
                n_c=self.n_c,
                mu_u=+mu_u,
                alpha_u=+alpha_u,
                mz_run=+mz_run,
                v=+v_ev,
                alpha1_mz=+alpha1,
                alpha2_mz=+alpha2,
                alpha3_mz=+alpha3,
                alpha_y_mz=+alpha_y,
                alpha_em_mz=+alpha_em,
                alpha_em_inv_mz=+(self.one / alpha_em),
                sin2w_mz=+sin2w,
            )

    def outer_p_from_alpha(self, alpha: Decimal) -> Decimal:
        with localcontext() as ctx:
            ctx.prec = self.work_precision
            return +(self.phi + alpha * self.sqrt_pi)

    def outer_alpha_from_p(self, p: Decimal) -> Decimal:
        with localcontext() as ctx:
            ctx.prec = self.work_precision
            return +((p - self.phi) / self.sqrt_pi)

    def observed_p_from_alpha_inv(self, alpha_inv: Decimal | int | str | float = CODATA_2022_ALPHA_INV) -> Decimal:
        with localcontext() as ctx:
            ctx.prec = self.work_precision
            return +self.outer_p_from_alpha(self.one / _dec(alpha_inv))

    def alpha_external_from_d10(
        self,
        d10: D10Point,
        mode: str,
    ) -> tuple[Decimal, Decimal, dict[str, Any] | None]:
        with localcontext() as ctx:
            ctx.prec = self.work_precision
            if mode == "mz_anchor":
                alpha_inv = d10.alpha_em_inv_mz
                return +(self.one / alpha_inv), +alpha_inv, None
            if mode == "thomson_structured_running":
                running = self.structured_thomson_running(d10)
                alpha_inv = d10.alpha_em_inv_mz + running["total_delta_alpha_inv"]
                return +(self.one / alpha_inv), +alpha_inv, running
            if mode == "thomson_structured_running_asymptotic":
                running = self.structured_thomson_running_asymptotic(d10)
                alpha_inv = d10.alpha_em_inv_mz + running["total_delta_alpha_inv"]
                return +(self.one / alpha_inv), +alpha_inv, running
        raise ValueError(f"Unsupported mode: {mode}")

    def evaluate_alpha_fixed_point(self, alpha_probe: Decimal | int | str | float, mode: str) -> FixedPointProbe:
        with localcontext() as ctx:
            ctx.prec = self.work_precision
            alpha = _dec(alpha_probe)
            p_outer = self.outer_p_from_alpha(alpha)
            d10 = self.build_d10_from_p(p_outer)
            inner_alpha, inner_alpha_inv, _running = self.alpha_external_from_d10(d10, mode)
            return FixedPointProbe(
                alpha=+alpha,
                alpha_inv=+(self.one / alpha),
                p=+p_outer,
                inner_alpha=+inner_alpha,
                inner_alpha_inv=+inner_alpha_inv,
                residual_alpha=+(inner_alpha - alpha),
            )

    def solve_closure(
        self,
        mode: str = "thomson_structured_running",
        max_iterations: int = 80,
        tolerance: Decimal | None = None,
        scan_points: int = 60,
    ) -> ClosureResult:
        with localcontext() as ctx:
            ctx.prec = self.work_precision
            tol = tolerance or Decimal(10) ** (-(self.precision - 4))
            alpha_min = Decimal("0.005")
            alpha_max = Decimal("0.01")
            if scan_points < 2:
                raise ValueError("scan_points must be at least 2")

            def evaluate(alpha_probe: Decimal) -> tuple[Decimal, Decimal, D10Point, Decimal, dict[str, Any] | None]:
                p_outer = self.outer_p_from_alpha(alpha_probe)
                d10 = self.build_d10_from_p(p_outer)
                inner_alpha, inner_alpha_inv, running = self.alpha_external_from_d10(d10, mode)
                residual = inner_alpha - alpha_probe
                return +residual, +inner_alpha, d10, +inner_alpha_inv, running

            lo_alpha: Decimal | None = None
            hi_alpha: Decimal | None = None
            lo_eval: tuple[Decimal, Decimal, D10Point, Decimal, dict[str, Any] | None] | None = None
            hi_eval: tuple[Decimal, Decimal, D10Point, Decimal, dict[str, Any] | None] | None = None
            previous_alpha: Decimal | None = None
            previous_eval: tuple[Decimal, Decimal, D10Point, Decimal, dict[str, Any] | None] | None = None

            for index in range(scan_points + 1):
                frac = Decimal(index) / Decimal(scan_points)
                alpha_probe = alpha_min + (alpha_max - alpha_min) * frac
                current_eval = evaluate(alpha_probe)
                if previous_eval is not None and current_eval[0] * previous_eval[0] <= 0:
                    lo_alpha = previous_alpha
                    hi_alpha = alpha_probe
                    lo_eval = previous_eval
                    hi_eval = current_eval
                    break
                previous_alpha = alpha_probe
                previous_eval = current_eval

            if lo_alpha is None or hi_alpha is None or lo_eval is None or hi_eval is None:
                raise RuntimeError("Could not bracket the alpha closure root.")

            steps: list[ClosureStep] = []
            final_alpha = hi_alpha
            final_d10 = hi_eval[2]
            final_alpha_inv = hi_eval[3]
            final_running = hi_eval[4]
            for iteration in range(max_iterations):
                alpha_mid = (lo_alpha + hi_alpha) / self.two
                residual_mid, inner_alpha_mid, d10_mid, inner_alpha_inv_mid, running_mid = evaluate(alpha_mid)
                steps.append(
                    ClosureStep(
                        iteration=iteration,
                        alpha_probe=+alpha_mid,
                        alpha_probe_inv=+(self.one / alpha_mid),
                        outer_p=+self.outer_p_from_alpha(alpha_mid),
                        inner_alpha=+inner_alpha_mid,
                        inner_alpha_inv=+inner_alpha_inv_mid,
                        source_anchor_alpha_inv=+d10_mid.alpha_em_inv_mz,
                        residual_alpha=+residual_mid,
                    )
                )
                final_alpha = alpha_mid
                final_d10 = d10_mid
                final_alpha_inv = inner_alpha_inv_mid
                final_running = running_mid
                if abs(residual_mid) < tol:
                    break
                if lo_eval[0] * residual_mid <= 0:
                    hi_alpha = alpha_mid
                    hi_eval = (residual_mid, inner_alpha_mid, d10_mid, inner_alpha_inv_mid, running_mid)
                else:
                    lo_alpha = alpha_mid
                    lo_eval = (residual_mid, inner_alpha_mid, d10_mid, inner_alpha_inv_mid, running_mid)

            p = self.outer_p_from_alpha(final_alpha)
            alpha_check, _, _ = self.alpha_external_from_d10(final_d10, mode)
            return ClosureResult(
                mode=mode,
                precision=self.precision,
                su2_cutoff=self.su2_cutoff,
                su3_cutoff=self.su3_cutoff,
                phi=+self.phi,
                sqrt_pi=+self.sqrt_pi,
                alpha=+final_alpha,
                alpha_inv=+final_alpha_inv,
                p=+p,
                god_equation_residual=+(p - self.outer_p_from_alpha(final_alpha)),
                alpha_fixed_point_residual=+(alpha_check - final_alpha),
                source_anchor_alpha_inv=+final_d10.alpha_em_inv_mz,
                structured_running=final_running,
                d10=final_d10,
                iterations=steps,
            )


def to_serializable(obj: Any) -> Any:
    if isinstance(obj, Decimal):
        return format(obj, "f")
    if isinstance(obj, list):
        return [to_serializable(item) for item in obj]
    if isinstance(obj, tuple):
        return [to_serializable(item) for item in obj]
    if hasattr(obj, "__dataclass_fields__"):
        return {key: to_serializable(value) for key, value in asdict(obj).items()}
    if isinstance(obj, dict):
        return {key: to_serializable(value) for key, value in obj.items()}
    return obj


def build_report(
    precision: int = 40,
    mode: str = "thomson_structured_running",
    su2_cutoff: int = 120,
    su3_cutoff: int = 90,
    max_iterations: int = 20,
    scan_points: int = 60,
) -> dict[str, Any]:
    ctx = PaperMathContext(precision=precision, su2_cutoff=su2_cutoff, su3_cutoff=su3_cutoff)
    result = ctx.solve_closure(mode=mode, max_iterations=max_iterations, scan_points=scan_points)
    return to_serializable(result)


def build_fixed_point_witness(
    precision: int = 40,
    mode: str = "thomson_structured_running",
    su2_cutoff: int = 120,
    su3_cutoff: int = 90,
    max_iterations: int = 20,
    scan_points: int = 60,
    derivative_step: str = "0.000001",
    sample_points: int = 5,
) -> dict[str, Any]:
    """Build a numerical witness for the alpha -> alpha fixed-point map.

    This is intentionally not an interval-arithmetic proof. It records a
    reproducible numerical fixed-point witness and finite-difference diagnostics
    so public claims can distinguish "candidate emitted by the code" from
    "uniqueness certified on an interval".
    """
    ctx = PaperMathContext(precision=precision, su2_cutoff=su2_cutoff, su3_cutoff=su3_cutoff)
    result = ctx.solve_closure(mode=mode, max_iterations=max_iterations, scan_points=scan_points)
    with localcontext() as dec_ctx:
        dec_ctx.prec = ctx.work_precision
        h = _dec(derivative_step)
        if h <= 0:
            raise ValueError("derivative_step must be positive")
        if sample_points < 1:
            raise ValueError("sample_points must be at least 1")

        center = result.alpha
        start = center - h * Decimal(sample_points // 2)
        probes: list[dict[str, Any]] = []
        slopes: list[Decimal] = []
        for index in range(sample_points):
            alpha = start + h * Decimal(index)
            probe = ctx.evaluate_alpha_fixed_point(alpha, mode)
            probes.append(to_serializable(probe))
            left = ctx.evaluate_alpha_fixed_point(alpha - h, mode)
            right = ctx.evaluate_alpha_fixed_point(alpha + h, mode)
            slope = (right.inner_alpha - left.inner_alpha) / (Decimal(2) * h)
            slopes.append(+slope)

        observed_p = ctx.observed_p_from_alpha_inv(CODATA_2022_ALPHA_INV)
        codata_alpha = ctx.one / CODATA_2022_ALPHA_INV
        codata_delta_alpha_inv = result.alpha_inv - CODATA_2022_ALPHA_INV
        codata_delta_p = result.p - observed_p

        max_abs_sample_slope = max(abs(slope) for slope in slopes)
        witness = {
            "claim_status": "numerical_witness_not_interval_certificate",
            "claim_boundary": (
                "This report samples the declared closure map and estimates local finite-difference "
                "slopes. It does not prove interval-wide uniqueness or a Banach contraction bound."
            ),
            "mode": mode,
            "precision": precision,
            "su2_cutoff": su2_cutoff,
            "su3_cutoff": su3_cutoff,
            "scan_points": scan_points,
            "derivative_step": h,
            "sample_points": sample_points,
            "fixed_point": result,
            "finite_difference": {
                "max_abs_sample_slope": +max_abs_sample_slope,
                "slopes": slopes,
                "probes": probes,
            },
            "codata_2022_compare_only": {
                "alpha_inv": CODATA_2022_ALPHA_INV,
                "alpha_inv_standard_uncertainty": CODATA_2022_ALPHA_INV_STD_UNCERTAINTY,
                "alpha": +codata_alpha,
                "p_from_alpha_inv": +observed_p,
                "fixed_point_minus_codata_alpha_inv": +codata_delta_alpha_inv,
                "fixed_point_p_minus_codata_p": +codata_delta_p,
            },
        }
        return to_serializable(witness)


def compute_alpha(
    precision: int = 40,
    mode: str = "thomson_structured_running",
    su2_cutoff: int = 120,
    su3_cutoff: int = 90,
    max_iterations: int = 20,
    scan_points: int = 60,
) -> tuple[Decimal, Decimal]:
    ctx = PaperMathContext(precision=precision, su2_cutoff=su2_cutoff, su3_cutoff=su3_cutoff)
    result = ctx.solve_closure(mode=mode, max_iterations=max_iterations, scan_points=scan_points)
    return result.alpha, result.alpha_inv
