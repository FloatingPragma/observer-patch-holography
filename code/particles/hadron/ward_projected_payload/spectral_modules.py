#!/usr/bin/env python3
"""First-principles spectral-measure modules for the payload harness.

LABEL (mandatory): development bracket, non-blind environment; the protocol
pass requires an isolated re-run.

All modules are computable on this machine from the Stage-5 internal chain
plus standard perturbative QCD structure. None of them is the missing
nonperturbative Ward-projected source measure. Inputs:

- Quark masses: internal Stage-5 continuation (paper_math
  ``diagonal_quark_masses``), never PDG.
- alpha_s: the D10 chain's own alpha_3(m_Z;P*) run with the four-loop MS-bar
  beta function and three-loop decoupling, ported from the source-only lane
  ``code/particles/qcd/derive_lambda_qcd_source_transmutation.py``. Flavor
  thresholds are placed at the internal Stage-5 m_c, m_b, m_t.
- Lambda_QCD scale: the source-only transmutation lane artifact
  ``code/particles/runs/qcd/lambda_qcd_source_transmutation.json``. Only the
  dimensionless ratio Lambda/mZ is consumed, so the unclosed clock candidate
  cancels. The declared interval covers the lane's threshold-bracket sweep.

Declared branches (stated, not tuned):

- parton_free: free-quark dispersion density, thresholds at 4 m_q^2 with
  Stage-5 masses. Exactly reproduces the chain's naive one-loop quark sum.
- pqcd: parton density multiplied by the massless MS-bar R-ratio series
  1 + a + c2 a^2 + c3 a^3 (a = alpha_s(s)/pi) above a declared IR cutoff
  s0 = (k * Lambda3)^2. Below s0 the density is either the free parton
  density ("free") or zero ("zero", confinement-motivated support cut).
- constituent: free-quark dispersion with dressed masses
  M_q = sqrt(m_q^2 + (kappa * Lambda3)^2), kappa = 1 declared.
"""

from __future__ import annotations

import math
from bisect import bisect_left
from functools import lru_cache
from typing import Any, Callable

from payload_harness import (
    N_C,
    QUARK_CHARGES_SQ,
    QUARK_ORDER,
    EvaluationPoint,
    beta_of_y,
    parton_moment,
)

ZETA3 = 1.2020569031595942

# Source-only Lambda_QCD transmutation lane (dimensionless ratios to mZ_run;
# the lane's display-GeV clock candidate cancels in the ratio).
LAMBDA_LANE_MZ_DISPLAY_GEV = 91.58833750192969
LAMBDA3_CENTRAL_OVER_MZ = 0.33481516768895897 / LAMBDA_LANE_MZ_DISPLAY_GEV
LAMBDA3_LO_OVER_MZ = 0.31949328977395597 / LAMBDA_LANE_MZ_DISPLAY_GEV
LAMBDA3_HI_OVER_MZ = 0.34974968299682957 / LAMBDA_LANE_MZ_DISPLAY_GEV
LAMBDA3_GRID = {
    "lane_lo": LAMBDA3_LO_OVER_MZ,
    "lane_central": LAMBDA3_CENTRAL_OVER_MZ,
    "lane_hi": LAMBDA3_HI_OVER_MZ,
}

Y_TAIL_MAX = 1.0e12  # dispersion tail handled analytically above this y


# ----------------------------------------------------------------------------
# Four-loop alpha_s engine (ported form of the source-only lane, scipy-free).
# ----------------------------------------------------------------------------


def beta_coeffs(nf: int) -> tuple[float, float, float, float]:
    b0 = 11.0 - 2.0 * nf / 3.0
    b1 = 102.0 - 38.0 * nf / 3.0
    b2 = 2857.0 / 2.0 - 5033.0 * nf / 18.0 + 325.0 * nf**2 / 54.0
    b3 = (
        149753.0 / 6.0
        + 3564.0 * ZETA3
        - (1078361.0 / 162.0 + 6508.0 * ZETA3 / 27.0) * nf
        + (50065.0 / 162.0 + 6472.0 * ZETA3 / 81.0) * nf**2
        + 1093.0 * nf**3 / 729.0
    )
    return b0, b1, b2, b3


def _decouple_down(alpha_nf: float, nl: int) -> float:
    a = alpha_nf / math.pi
    c2 = -11.0 / 72.0
    c3 = -564731.0 / 124416.0 + 82043.0 * ZETA3 / 27648.0 + 2633.0 * nl / 31104.0
    return alpha_nf * (1.0 + c2 * a * a + c3 * a**3)


def _decouple_up(alpha_nl: float, nl: int) -> float:
    """Invert the downward matching: find alpha_nf with decouple_down = alpha_nl."""
    alpha_nf = alpha_nl
    for _ in range(6):
        ratio = _decouple_down(alpha_nf, nl) / alpha_nf
        alpha_nf = alpha_nl / ratio
    return alpha_nf


class AlphaSTable:
    """Deterministic four-loop running table in t = ln(mu/mZ)."""

    def __init__(
        self,
        alpha_mz: float,
        t_c: float,
        t_b: float,
        t_t: float,
        t_min: float = -5.2,
        t_max: float = 15.0,
        h: float = 0.001,
        alpha_cap: float = 6.0,
    ):
        # t_min sits just below the smallest declared IR cutoff mu0 = 2*Lambda3.
        # The walk clamps at alpha_cap; clamped entries lie below every
        # declared cutoff and are never consumed by a correction integral.
        self.alpha_cap = alpha_cap
        self.h = h
        self.t_min = t_min
        self.t_max = t_max
        self.thresholds = {"c": t_c, "b": t_b, "t": t_t}
        self.ts: list[float] = []
        self.alphas: list[float] = []
        self._walk(alpha_mz, t_c, t_b, t_t)

    @staticmethod
    def _deriv(a: float, coeffs: tuple[float, float, float, float]) -> float:
        b0, b1, b2, b3 = coeffs
        return -2.0 * a * a * (b0 + b1 * a + b2 * a * a + b3 * a**3)

    def _rk4_leg(
        self, alpha0: float, t0: float, t1: float, nf: int, out_t: list[float], out_a: list[float]
    ) -> float:
        coeffs = beta_coeffs(nf)
        n_steps = max(1, int(round(abs(t1 - t0) / self.h)))
        dt = (t1 - t0) / n_steps
        a_cap = self.alpha_cap / (4.0 * math.pi)
        a = alpha0 / (4.0 * math.pi)
        clamped = False
        for i in range(n_steps):
            out_t.append(t0 + i * dt)
            out_a.append(4.0 * math.pi * a)
            if clamped:
                continue
            k1 = self._deriv(a, coeffs)
            k2 = self._deriv(a + 0.5 * dt * k1, coeffs)
            k3 = self._deriv(a + 0.5 * dt * k2, coeffs)
            k4 = self._deriv(a + dt * k3, coeffs)
            a_next = a + dt / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
            if not math.isfinite(a_next) or a_next > a_cap or a_next <= 0.0:
                a = a_cap
                clamped = True
            else:
                a = a_next
        out_t.append(t1)
        out_a.append(4.0 * math.pi * a)
        return 4.0 * math.pi * a

    def _walk(self, alpha_mz: float, t_c: float, t_b: float, t_t: float) -> None:
        # Downward: mZ (nf=5) -> m_b -> (nf=4) -> m_c -> (nf=3) -> t_min.
        down_t: list[float] = []
        down_a: list[float] = []
        a_mb = self._rk4_leg(alpha_mz, 0.0, t_b, 5, down_t, down_a)
        a_mb4 = _decouple_down(a_mb, 4)
        a_mc = self._rk4_leg(a_mb4, t_b, t_c, 4, down_t, down_a)
        a_mc3 = _decouple_down(a_mc, 3)
        self._rk4_leg(a_mc3, t_c, self.t_min, 3, down_t, down_a)
        # Upward: mZ (nf=5) -> m_t -> (nf=6) -> t_max.
        up_t: list[float] = []
        up_a: list[float] = []
        a_mt = self._rk4_leg(alpha_mz, 0.0, t_t, 5, up_t, up_a)
        a_mt6 = _decouple_up(a_mt, 5)
        self._rk4_leg(a_mt6, t_t, self.t_max, 6, up_t, up_a)

        pairs = sorted(zip(down_t, down_a)) + sorted(zip(up_t, up_a))
        pairs.sort(key=lambda item: item[0])
        self.ts = [t for t, _ in pairs]
        self.alphas = [a for _, a in pairs]

    def alpha_s(self, t: float) -> float:
        if t <= self.ts[0]:
            return self.alphas[0]
        if t >= self.ts[-1]:
            return self.alphas[-1]
        idx = bisect_left(self.ts, t)
        t0, t1 = self.ts[idx - 1], self.ts[idx]
        a0, a1 = self.alphas[idx - 1], self.alphas[idx]
        if t1 == t0:
            return a1
        w = (t - t0) / (t1 - t0)
        return a0 + w * (a1 - a0)

    def nf_at(self, t: float) -> int:
        if t < self.thresholds["c"]:
            return 3
        if t < self.thresholds["b"]:
            return 4
        if t < self.thresholds["t"]:
            return 5
        return 6


@lru_cache(maxsize=8)
def _alpha_table_cached(
    alpha_mz: float, t_c: float, t_b: float, t_t: float
) -> AlphaSTable:
    return AlphaSTable(alpha_mz, t_c, t_b, t_t)


def alpha_table_for(ep: EvaluationPoint) -> AlphaSTable:
    t_c = math.log(ep.quark_masses["c"] / ep.mz_run)
    t_b = math.log(ep.quark_masses["b"] / ep.mz_run)
    t_t = math.log(ep.quark_masses["t"] / ep.mz_run)
    return _alpha_table_cached(ep.alpha3_mz, t_c, t_b, t_t)


# ----------------------------------------------------------------------------
# Massless MS-bar R-ratio correction series.
# ----------------------------------------------------------------------------

ACTIVE_CHARGES = {
    3: ("u", "d", "s"),
    4: ("u", "d", "s", "c"),
    5: ("u", "d", "s", "c", "b"),
    6: ("u", "d", "s", "c", "b"),  # top excluded from the transport sum
}
QUARK_CHARGES = {"u": 2.0 / 3.0, "d": -1.0 / 3.0, "s": -1.0 / 3.0, "c": 2.0 / 3.0, "b": -1.0 / 3.0}


def r_ratio_coefficients(nf: int) -> tuple[float, float, float]:
    c1 = 1.0
    c2 = 365.0 / 24.0 - 11.0 * ZETA3 - nf * (11.0 / 12.0 - 2.0 * ZETA3 / 3.0)
    flavors = ACTIVE_CHARGES[min(nf, 6)]
    sum_q = sum(QUARK_CHARGES[q] for q in flavors)
    sum_q2 = sum(QUARK_CHARGES[q] ** 2 for q in flavors)
    singlet = -1.2395 * (sum_q * sum_q) / (3.0 * sum_q2)
    c3 = -6.63694 - 1.20013 * nf - 0.00518 * nf * nf + singlet
    return c1, c2, c3


def qcd_correction_factor(a_pi: float, nf: int, order: int) -> float:
    c1, c2, c3 = r_ratio_coefficients(nf)
    factor = 1.0
    if order >= 1:
        factor += c1 * a_pi
    if order >= 2:
        factor += c2 * a_pi * a_pi
    if order >= 3:
        factor += c3 * a_pi**3
    return factor


# ----------------------------------------------------------------------------
# Modules.
# ----------------------------------------------------------------------------


class SpectralModule:
    def __init__(self, module_id: str, declared_branch: dict[str, Any], builder: Callable):
        self.module_id = module_id
        self.declared_branch = declared_branch
        self._builder = builder

    def segments(self, ep: EvaluationPoint) -> list[dict[str, Any]]:
        return self._builder(ep)


def _parton_segments(masses_of: Callable[[EvaluationPoint], dict[str, float]]):
    def build(ep: EvaluationPoint) -> list[dict[str, Any]]:
        masses = masses_of(ep)
        segments = []
        for name in QUARK_ORDER:
            y_thr = 4.0 * (masses[name] / ep.mz_run) ** 2
            segments.append(
                {
                    "type": "parton",
                    "label": f"free_parton_{name}",
                    "y_thr": y_thr,
                    "ncq2": N_C * QUARK_CHARGES_SQ[name],
                }
            )
        return segments

    return build


def make_parton_free() -> SpectralModule:
    return SpectralModule(
        "parton_free",
        {
            "variant": "free_quark_parton_one_loop",
            "masses": "internal_stage5_continuation",
            "qcd_correction": "none",
        },
        _parton_segments(lambda ep: ep.quark_masses),
    )


def make_constituent(lambda3_key: str, kappa: float = 1.0) -> SpectralModule:
    lambda3 = LAMBDA3_GRID[lambda3_key]

    def masses_of(ep: EvaluationPoint) -> dict[str, float]:
        dress = kappa * lambda3 * ep.mz_run
        return {
            name: math.sqrt(ep.quark_masses[name] ** 2 + dress * dress)
            for name in QUARK_ORDER
        }

    return SpectralModule(
        f"constituent_kappa{kappa:g}_{lambda3_key}",
        {
            "variant": "constituent_mass_one_loop",
            "masses": "internal_stage5_continuation_dressed",
            "dressing_rule": "M_q = sqrt(m_q^2 + (kappa*Lambda3)^2), declared branch",
            "kappa": kappa,
            "lambda3_over_mz": lambda3,
            "lambda3_source": "oph_lambda_qcd_source_transmutation (" + lambda3_key + ")",
            "qcd_correction": "none",
        },
        _parton_segments(masses_of),
    )


def make_pqcd(
    lambda3_key: str,
    k_cut: float,
    below: str,
    order: int,
) -> SpectralModule:
    """Parton density with massless MS-bar QCD correction above s0=(k*Lambda3)^2."""
    lambda3 = LAMBDA3_GRID[lambda3_key]
    if below not in ("free", "zero"):
        raise ValueError("below must be 'free' or 'zero'")

    def build(ep: EvaluationPoint) -> list[dict[str, Any]]:
        table = alpha_table_for(ep)
        y0 = (k_cut * lambda3) ** 2  # (k*Lambda3/mZ)^2 in y = s/mZ^2 units
        segments: list[dict[str, Any]] = []
        thresholds = {
            name: 4.0 * (ep.quark_masses[name] / ep.mz_run) ** 2 for name in QUARK_ORDER
        }

        # Free-parton base, full support (closed form).
        for name in QUARK_ORDER:
            segments.append(
                {
                    "type": "parton",
                    "label": f"free_parton_{name}",
                    "y_thr": thresholds[name],
                    "ncq2": N_C * QUARK_CHARGES_SQ[name],
                }
            )

        # Support cut: remove the free density below y0 (closed form).
        if below == "zero":
            for name in QUARK_ORDER:
                y_thr = thresholds[name]
                if y0 > y_thr:
                    segments.append(
                        {
                            "type": "parton",
                            "label": f"support_cut_{name}",
                            "y_thr": y_thr,
                            "ncq2": -N_C * QUARK_CHARGES_SQ[name],
                            "y_hi": y0,
                        }
                    )

        # QCD correction density above y0 (numeric).
        def rho_corr(y: float) -> float:
            t = 0.5 * math.log(y)
            a_pi = table.alpha_s(t) / math.pi
            nf = table.nf_at(t)
            factor = qcd_correction_factor(a_pi, nf, order) - 1.0
            base = 0.0
            for name in QUARK_ORDER:
                b = beta_of_y(thresholds[name], y)
                if b > 0.0:
                    base += N_C * QUARK_CHARGES_SQ[name] * b * (3.0 - b * b) / 2.0
            return base * factor

        # Piece boundaries at quark thresholds inside [y0, Y_TAIL_MAX].
        edges = [y0]
        for name in ("c", "b"):
            y_thr = thresholds[name]
            if y0 < y_thr < Y_TAIL_MAX:
                edges.append(y_thr)
        edges.append(Y_TAIL_MAX)
        edges = sorted(set(edges))
        for lo, hi in zip(edges[:-1], edges[1:]):
            sqrt_left = any(abs(lo - thresholds[name]) < 1e-30 for name in ("c", "b"))
            segments.append(
                {
                    "type": "density",
                    "label": f"pqcd_correction_[{lo:.6e},{hi:.6e}]",
                    "y_lo": lo,
                    "y_hi": hi,
                    "rho": rho_corr,
                    "sqrt_left": sqrt_left,
                    "signed": True,
                }
            )

        # Analytic tail with the correction frozen at Y_TAIL_MAX.
        t_tail = 0.5 * math.log(Y_TAIL_MAX)
        a_tail = table.alpha_s(t_tail) / math.pi
        r_inf = sum(N_C * QUARK_CHARGES_SQ[name] for name in QUARK_ORDER) * (
            qcd_correction_factor(a_tail, table.nf_at(t_tail), order) - 1.0
        )
        segments.append(
            {"type": "tail", "label": "pqcd_correction_tail", "y_lo": Y_TAIL_MAX, "r_inf": r_inf}
        )
        return segments

    return SpectralModule(
        f"pqcd_o{order}_k{k_cut:g}_{below}_{lambda3_key}",
        {
            "variant": "parton_plus_massless_msbar_r_ratio_series",
            "masses": "internal_stage5_continuation",
            "alpha_s": "four_loop_from_internal_alpha3_mz_thresholds_stage5",
            "order_in_a_pi": order,
            "ir_cutoff_rule": "s0 = (k*Lambda3)^2, declared sweep grid",
            "k_cut": k_cut,
            "below_cutoff": below,
            "lambda3_over_mz": lambda3,
            "lambda3_source": "oph_lambda_qcd_source_transmutation (" + lambda3_key + ")",
        },
        build,
    )


class SyntheticAtomModule:
    """Synthetic module for round-trip tests: finite list of positive atoms."""

    def __init__(self, atoms: list[tuple[float, float]], module_id: str = "synthetic_atoms"):
        self.module_id = module_id
        self.declared_branch = {"variant": "synthetic_atoms_test_only"}
        self._atoms = list(atoms)

    def segments(self, ep: EvaluationPoint) -> list[dict[str, Any]]:
        return [
            {"type": "atom", "label": f"atom_{i}", "y": y, "weight": w}
            for i, (y, w) in enumerate(self._atoms)
        ]


def get_module(name: str) -> Any:
    if name == "parton_free":
        return make_parton_free()
    if name.startswith("constituent_"):
        key = name.split("constituent_", 1)[1]
        return make_constituent(key)
    raise ValueError(f"unknown module name: {name}")
