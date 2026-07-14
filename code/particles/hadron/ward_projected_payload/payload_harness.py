#!/usr/bin/env python3
"""Ward-projected hadronic transport payload harness (generator G1).

LABEL (mandatory): development bracket, non-blind environment; the protocol
pass requires an isolated re-run.

Declared contract implemented here, per
``code/P_derivation/THOMSON_TRANSPORT_THEOREMS.md`` and
``code/P_derivation/SOURCE_SPECTRAL_THEOREM.md``:

- Object emitted: the source transport packet in inverse-alpha units,

      Delta_source(P) = Delta_lep(P) + Delta_had(P) + Delta_EW(P),

  where Delta_had is the once-subtracted dispersion moment of a positive
  Ward-projected U(1)_Q spectral density rho_Q(s;P),

      Delta_had(P) = mZ(P)^2/(3*pi) * Integral rho_Q(s;P)/(s*(s+mZ(P)^2)) ds.

- Scheme: same subtraction as a0(P) = alpha_em^-1(m_Z^2;P), Ward-projected
  U(1)_Q lane, D10 source family ``d10_running_tree``.

- Evaluation point: the implemented internal fixed-point root P* of the
  Stage-5 chain (mode ``thomson_structured_running``), rebuilt locally from
  ``paper_math.PaperMathContext`` at declared precision and cutoffs. No
  CODATA/NIST value, no measured hadronic data, and no PDG hadron data enter
  any computation in this module.

- Screening coordinate: x = N_c * alpha_3(m_Z;P)/pi. The chain's implemented
  screen is S_impl = 1 - x. The declared residual form is
  S = 1 - x + c_Q x^2. This module reports, for every payload,

      S_eff  = Delta_had / Delta_quark_naive_one_loop,
      c_Q    = (S_eff - (1 - x)) / x^2.

Spectral modules feed the harness through ``emit_delta_source``. A module is
any object with attributes ``module_id`` (str), ``declared_branch`` (dict),
and a method ``segments(ep)`` returning a list of segment dicts:

  {"type": "atom",    "y": <s/mZ^2>, "weight": <positive weight>}
  {"type": "parton",  "y_thr": <4m^2/mZ^2>, "ncq2": <N_c*Q^2>,
                      "y_lo": None, "y_hi": None}              (closed form)
  {"type": "density", "y_lo": ..., "y_hi": ..., "rho": f(y)->float,
                      "sqrt_left": bool}                        (numeric)
  {"type": "tail",    "y_lo": ..., "r_inf": <constant R value>} (closed form)

Densities use the R-ratio normalization: a free massless parton species of
color factor N_c and charge Q contributes rho = N_c * Q^2.
"""

from __future__ import annotations

import hashlib
import json
import math
import sys
from dataclasses import dataclass
from decimal import Decimal
from functools import lru_cache
from pathlib import Path
from typing import Any, Callable

HERE = Path(__file__).resolve().parent
P_DERIVATION = HERE.parents[2] / "P_derivation"
if str(P_DERIVATION) not in sys.path:
    sys.path.insert(0, str(P_DERIVATION))

NON_BLIND_LABEL = (
    "development bracket, non-blind environment; the protocol pass requires "
    "an isolated re-run"
)

# Internal fixed-point root of the implemented Stage-5 chain, mode
# thomson_structured_running, precision 100, from
# code/P_derivation/runtime/full_p_alpha_report_current.json.
# This is a Stage-5 internal chain quantity; it is not the CODATA-mapped P_C.
P_STAR_INTERNAL = (
    "1.6309720958588973769645139035069556298539032554195104407930688323138732094230"
)

DEFAULT_PRECISION = 40
DEFAULT_SU2_CUTOFF = 60
DEFAULT_SU3_CUTOFF = 45

QUARK_CHARGES_SQ = {"u": 4.0 / 9.0, "d": 1.0 / 9.0, "s": 1.0 / 9.0, "c": 4.0 / 9.0, "b": 1.0 / 9.0}
QUARK_ORDER = ("u", "d", "s", "c", "b")
N_C = 3


# ----------------------------------------------------------------------------
# Deterministic Gauss-Legendre nodes.
# ----------------------------------------------------------------------------


@lru_cache(maxsize=None)
def gauss_legendre(n: int) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Nodes and weights on [-1, 1] by Newton iteration on P_n."""
    nodes: list[float] = []
    weights: list[float] = []
    for i in range(1, n + 1):
        x = math.cos(math.pi * (i - 0.25) / (n + 0.5))
        for _ in range(100):
            p0, p1 = 1.0, x
            for k in range(2, n + 1):
                p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
            dp = n * (x * p1 - p0) / (x * x - 1.0)
            dx = p1 / dp
            x -= dx
            if abs(dx) < 1e-15:
                break
        p0, p1 = 1.0, x
        for k in range(2, n + 1):
            p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
        dp = n * (x * p1 - p0) / (x * x - 1.0)
        nodes.append(x)
        weights.append(2.0 / ((1.0 - x * x) * dp * dp))
    return tuple(nodes), tuple(weights)


def _gl_integrate(f: Callable[[float], float], lo: float, hi: float, n: int) -> float:
    nodes, weights = gauss_legendre(n)
    mid = 0.5 * (lo + hi)
    half = 0.5 * (hi - lo)
    return half * math.fsum(w * f(mid + half * x) for x, w in zip(nodes, weights))


# ----------------------------------------------------------------------------
# Evaluation point (Stage-5 internal chain rebuild).
# ----------------------------------------------------------------------------


@dataclass(frozen=True)
class EvaluationPoint:
    p: str
    precision: int
    su2_cutoff: int
    su3_cutoff: int
    mz_run: float
    v: float
    alpha3_mz: float
    a0_alpha_em_inv_mz: float
    quark_masses: dict[str, float]
    lepton_masses: dict[str, float]

    @property
    def x_screen(self) -> float:
        return N_C * self.alpha3_mz / math.pi

    def to_json(self) -> dict[str, Any]:
        return {
            "p": self.p,
            "p_provenance": (
                "internal fixed-point root of the implemented Stage-5 chain "
                "(mode thomson_structured_running); rebuilt D10 point via "
                "paper_math.PaperMathContext.build_d10_from_p"
            ),
            "precision": self.precision,
            "su2_cutoff": self.su2_cutoff,
            "su3_cutoff": self.su3_cutoff,
            "mz_run_planck_units": self.mz_run,
            "v_planck_units": self.v,
            "alpha3_mz": self.alpha3_mz,
            "a0_alpha_em_inv_mz": self.a0_alpha_em_inv_mz,
            "x_screen": self.x_screen,
            "mass_source": "internal_stage5_continuation",
            "quark_masses_planck_units": dict(self.quark_masses),
            "lepton_masses_planck_units": dict(self.lepton_masses),
        }


def build_evaluation_point(
    p_value: str = P_STAR_INTERNAL,
    precision: int = DEFAULT_PRECISION,
    su2_cutoff: int = DEFAULT_SU2_CUTOFF,
    su3_cutoff: int = DEFAULT_SU3_CUTOFF,
) -> EvaluationPoint:
    from paper_math import PaperMathContext  # Stage-5 internal chain

    ctx = PaperMathContext(precision=precision, su2_cutoff=su2_cutoff, su3_cutoff=su3_cutoff)
    d10 = ctx.build_d10_from_p(Decimal(p_value))
    quarks = {k: float(v) for k, v in ctx.diagonal_quark_masses(d10.v).items()}
    leptons = {k: float(v) for k, v in ctx.charged_lepton_masses(d10.v).items()}
    return EvaluationPoint(
        p=p_value,
        precision=precision,
        su2_cutoff=su2_cutoff,
        su3_cutoff=su3_cutoff,
        mz_run=float(d10.mz_run),
        v=float(d10.v),
        alpha3_mz=float(d10.alpha3_mz),
        a0_alpha_em_inv_mz=float(d10.alpha_em_inv_mz),
        quark_masses=quarks,
        lepton_masses=leptons,
    )


# ----------------------------------------------------------------------------
# Closed-form parton dispersion moments.
#
# For one species (color factor N_c, squared charge Q^2, threshold
# y_thr = 4 m^2 / mZ^2), with beta = sqrt(1 - y_thr/y) and c = 1 + y_thr:
#
#   Delta = ncq2/(3*pi) * [F(y_hi) - F(y_lo)],
#   F(y)  = b^3/3 + (c - 3) b + c(3 - c)/(2 sqrt(c)) * ln((sqrt(c)+b)/(sqrt(c)-b)),
#
# evaluated at b = beta(y). The log term uses sqrt(c)-b = (c-b^2)/(sqrt(c)+b)
# with c - b^2 = y_thr + y_thr/y composed from exact ratios, so the light-mass
# limit stays stable in double precision. This is the exact once-subtracted
# dispersion moment of the free-fermion spectral density
# rho = ncq2 * beta (3 - beta^2)/2 against the kernel
# mZ^2/(3 pi s (s + mZ^2)), identical to the closed one-loop kernel in
# paper_math.fermion_transport_kernel_exact when y runs over [y_thr, inf).
# ----------------------------------------------------------------------------


def parton_moment(
    y_thr: float,
    ncq2: float,
    y_lo: float | None = None,
    y_hi: float | None = None,
) -> float:
    """Dispersion moment over [y_lo, y_hi]; None means [y_thr, infinity)."""
    c = 1.0 + y_thr
    sc = math.sqrt(c)

    def antiderivative(w: float) -> float:
        # w = y_thr / y in (0, 1]; w -> 0 is the y -> infinity limit.
        b = math.sqrt(1.0 - w)
        c_minus_b2 = y_thr + w  # equals c - b^2 without cancellation
        log_term = math.log((sc + b) * (sc + b) / c_minus_b2)
        return b * b * b / 3.0 + (c - 3.0) * b + c * (3.0 - c) / (2.0 * sc) * log_term

    def at(y: float | None, infinity_default: bool) -> float:
        if y is None:
            return antiderivative(0.0) if infinity_default else 0.0
        if y <= y_thr:
            return 0.0
        return antiderivative(y_thr / y)

    upper = at(y_hi, infinity_default=True)
    lower = at(y_lo, infinity_default=False)
    return ncq2 / (3.0 * math.pi) * (upper - lower)


def beta_of_y(y_thr: float, y: float) -> float:
    if y <= y_thr:
        return 0.0
    return math.sqrt(1.0 - y_thr / y)


def kernel_moment_atom(y: float) -> float:
    """Kernel value for a unit-weight atom at y = s/mZ^2, inverse-alpha units."""
    return 1.0 / (3.0 * math.pi * y * (1.0 + y))


def tail_moment(y_lo: float, r_inf: float) -> float:
    """Exact tail integral of a constant density r_inf above y_lo."""
    return r_inf / (3.0 * math.pi) * math.log1p(1.0 / y_lo)


# ----------------------------------------------------------------------------
# Numeric density integration (used for QCD correction densities and
# synthetic modules). Deterministic piecewise Gauss-Legendre in t = ln y,
# with a sqrt-removing substitution at declared threshold-side edges.
# ----------------------------------------------------------------------------


def density_moment(
    rho: Callable[[float], float],
    y_lo: float,
    y_hi: float,
    *,
    sqrt_left: bool = False,
    gauss_n: int = 48,
    splits_per_decade: int = 2,
) -> float:
    if y_hi <= y_lo:
        return 0.0
    t_lo = math.log(y_lo)
    t_hi = math.log(y_hi)

    def f_t(t: float) -> float:
        y = math.exp(t)
        return rho(y) / (1.0 + y)

    total = 0.0
    if sqrt_left:
        # First slice in tau = sqrt(t - t_lo) to remove the sqrt kink at y_lo.
        t_cut = min(t_hi, t_lo + min(1.0, t_hi - t_lo))
        tau_hi = math.sqrt(t_cut - t_lo)

        def f_tau(tau: float) -> float:
            return 2.0 * tau * f_t(t_lo + tau * tau)

        total += _gl_integrate(f_tau, 0.0, tau_hi, gauss_n)
        t_lo = t_cut
    if t_hi > t_lo:
        n_sub = max(1, math.ceil((t_hi - t_lo) / math.log(10.0) * splits_per_decade))
        edges = [t_lo + (t_hi - t_lo) * i / n_sub for i in range(n_sub + 1)]
        for a, b in zip(edges[:-1], edges[1:]):
            total += _gl_integrate(f_t, a, b, gauss_n)
    return total / (3.0 * math.pi)


# ----------------------------------------------------------------------------
# Baselines from the Stage-5 chain form.
# ----------------------------------------------------------------------------


def lepton_transport(ep: EvaluationPoint) -> float:
    total = 0.0
    for name in ("e", "mu", "tau"):
        y_thr = 4.0 * (ep.lepton_masses[name] / ep.mz_run) ** 2
        total += parton_moment(y_thr, 1.0)
    return total


def quark_naive_transport(ep: EvaluationPoint, masses: dict[str, float] | None = None) -> float:
    masses = masses if masses is not None else ep.quark_masses
    total = 0.0
    for name in QUARK_ORDER:
        y_thr = 4.0 * (masses[name] / ep.mz_run) ** 2
        total += parton_moment(y_thr, N_C * QUARK_CHARGES_SQ[name])
    return total


def implemented_screen(ep: EvaluationPoint) -> float:
    return 1.0 - ep.x_screen


# ----------------------------------------------------------------------------
# Harness emission.
# ----------------------------------------------------------------------------


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def emit_delta_source(
    module: Any,
    ep: EvaluationPoint,
    *,
    gauss_n: int = 48,
    splits_per_decade: int = 4,
) -> dict[str, Any]:
    """Evaluate one spectral-measure module against the declared contract."""
    delta_had = 0.0
    segment_report: list[dict[str, Any]] = []
    positivity_ok = True
    for seg in module.segments(ep):
        kind = seg["type"]
        if kind == "atom":
            value = seg["weight"] * kernel_moment_atom(seg["y"])
            if seg["weight"] < 0.0:
                positivity_ok = False
        elif kind == "parton":
            value = parton_moment(
                seg["y_thr"], seg["ncq2"], seg.get("y_lo"), seg.get("y_hi")
            )
        elif kind == "density":
            value = density_moment(
                seg["rho"],
                seg["y_lo"],
                seg["y_hi"],
                sqrt_left=bool(seg.get("sqrt_left", False)),
                gauss_n=gauss_n,
                splits_per_decade=splits_per_decade,
            )
            if bool(seg.get("signed", False)):
                positivity_ok = positivity_ok and value >= 0.0
        elif kind == "tail":
            value = tail_moment(seg["y_lo"], seg["r_inf"])
        else:
            raise ValueError(f"unknown segment type: {kind}")
        delta_had += value
        segment_report.append({"type": kind, "label": seg.get("label", ""), "value": value})

    delta_lep = lepton_transport(ep)
    naive = quark_naive_transport(ep)
    x = ep.x_screen
    s_eff = delta_had / naive
    c_q_implied = (s_eff - (1.0 - x)) / (x * x)
    delta_ew = 0.0
    delta_source_total = delta_lep + delta_had + delta_ew
    quark_screened_impl = implemented_screen(ep) * naive
    r_q_residual = delta_had + delta_ew - quark_screened_impl

    payload = {
        "artifact": "oph_ward_projected_payload_harness_delta_source",
        "label": NON_BLIND_LABEL,
        "source_family_id": "d10_running_tree",
        "current": "U1_Q",
        "scheme": {
            "same_subtraction_as_a0": True,
            "scheme_id": "d10_ward_projected_once_subtracted_at_mZ2",
            "normalization_convention": "R_ratio_massless_parton_NcQ2",
            "kernel": "mZ^2/(3*pi*s*(s+mZ^2))",
        },
        "evaluation_point": ep.to_json(),
        "module": {
            "module_id": module.module_id,
            "declared_branch": module.declared_branch,
        },
        "components_alpha_inv": {
            "delta_lep": delta_lep,
            "delta_had": delta_had,
            "delta_EW": {
                "value": delta_ew,
                "branch": "declared_zero_branch_unproven",
                "status": "theorem_4_open_zero_identity_not_established",
            },
        },
        "delta_source_alpha_inv": delta_source_total,
        "diagnostics": {
            "quark_delta_alpha_inv_naive_one_loop": naive,
            "quark_delta_alpha_inv_screened_impl": quark_screened_impl,
            "implemented_screen_1_minus_x": implemented_screen(ep),
            "x_screen": x,
            "s_effective": s_eff,
            "c_q_implied": c_q_implied,
            "r_q_residual_vs_implemented_screen": r_q_residual,
            "positivity_ok": positivity_ok,
        },
        "integration": {
            "gauss_n": gauss_n,
            "splits_per_decade": splits_per_decade,
            "segments": segment_report,
        },
        "promotion_allowed": False,
        "external_inputs_used": False,
    }
    digest_source = dict(payload)
    payload["content_sha256"] = hashlib.sha256(
        _canonical_json(digest_source).encode("utf-8")
    ).hexdigest()
    return payload


def main() -> int:
    import argparse

    import spectral_modules

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--module", default="parton_free")
    parser.add_argument("--precision", type=int, default=DEFAULT_PRECISION)
    parser.add_argument("--output", default="")
    args = parser.parse_args()

    ep = build_evaluation_point(precision=args.precision)
    module = spectral_modules.get_module(args.module)
    payload = emit_delta_source(module, ep)
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output:
        out = Path(args.output)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(text, encoding="utf-8")
        print(f"saved: {out}")
    else:
        print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
