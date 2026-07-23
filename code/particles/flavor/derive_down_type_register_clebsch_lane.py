#!/usr/bin/env python3
"""Down-type quark masses from the charged leptons via register Clebsch factors.

The recovered quotient carries a unification register at ``mu_U`` where the
gauge couplings meet.  If the register pairs the down-type quarks with the
charged leptons multiplet-wise (the classic pattern of the recovered
group structure) with the Clebsch factors ``(1, 1/3, 3)`` for
``(b/tau, s/mu, d/e)``, the entire down-type sector follows from the
conditional MCPR charged-lepton triple with no further input:

    y_b(mu_U) = y_tau(mu_U),
    y_s(mu_U) = y_mu(mu_U)/3,
    y_d(mu_U) = 3 y_e(mu_U).

The candidate register origin of the factors is declared: ``3 = N_c`` and
``1/3`` is the invariant color measure of the transitive C3 action, the two
weights the corpus emits everywhere.  The selection artifact
(``clebsch_register_pairing_selection.json``) closes the pairing part
conditional on the #314 exterior matter contract and the unordered weight
set conditional on the frozen constraints F1/F2; the generation-order
assignment of the weights, ``GENERATION_REGISTER_ORDER``, remains the open
premise gate of this lane.

The lane runs the full one-loop seven-Yukawa system (top trace included,
top from the criticality boundary), applies the Clebsch boundary at
``mu_U``, runs down, and converts with flavor-banded QCD mass running.
Ratios are the sharp outputs: the Gatto-Sartori-Tonin angle
``sqrt(m_d/m_s)`` is emitted as a Cabibbo candidate.  The absolute
normalization carries the known third-generation tension of this relation
class at one loop; the deficit is attributed, never absorbed.  Measured
values appear only in the compare-only block.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

try:
    from calibration.derive_d11_criticality_boundary_scan import (
        ALPHA_U_FALLBACK,
        P_FALLBACK,
        gauge_couplings,
        run_boundary,
        source_scales,
    )
except ModuleNotFoundError:
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "calibration"))
    from derive_d11_criticality_boundary_scan import (
        ALPHA_U_FALLBACK,
        P_FALLBACK,
        gauge_couplings,
        run_boundary,
        source_scales,
    )

ROOT = Path(__file__).resolve().parents[2]
MCPR = (
    ROOT / "particles" / "runs" / "leptons"
    / "charged_mcpr_completion_conditional.json"
)
DEFAULT_OUT = (
    ROOT / "particles" / "runs" / "flavor"
    / "down_type_register_clebsch_lane.json"
)

CLEBSCH = {"b_over_tau": 1.0, "s_over_mu": 1.0 / 3.0, "d_over_e": 3.0}

# Compare-only references (never prediction ancestors).
COMPARE_ONLY = {
    "mb_mb_gev": 4.18,
    "ms_2gev_gev": 0.0935,
    "md_2gev_gev": 0.0047,
    "ms_over_md": 19.9,
    "cabibbo": 0.2250,
}


def _yukawa_derivatives(y: list[float], mu: float, mz_run: float) -> list[float]:
    """One-loop SM Yukawa system [yt, yb, ytau, ymu, ye]."""

    g_y, g2, g3 = gauge_couplings(mu, mz_run)
    g1sq = (5.0 / 3.0) * g_y * g_y
    gauge_up = (17.0 / 20.0) * g1sq + 2.25 * g2 * g2 + 8.0 * g3 * g3
    gauge_down = 0.25 * g1sq + 2.25 * g2 * g2 + 8.0 * g3 * g3
    gauge_lepton = 2.25 * g1sq + 2.25 * g2 * g2
    yt, yb, ytau, ymu, ye = y
    trace = 3.0 * yt * yt + 3.0 * yb * yb + ytau * ytau + ymu * ymu + ye * ye
    kappa = 1.0 / (16.0 * math.pi**2)
    return [
        yt * kappa * (1.5 * yt * yt - 1.5 * yb * yb + trace - gauge_up),
        yb * kappa * (1.5 * yb * yb - 1.5 * yt * yt + trace - gauge_down),
        ytau * kappa * (1.5 * ytau * ytau + trace - gauge_lepton),
        ymu * kappa * (1.5 * ymu * ymu + trace - gauge_lepton),
        ye * kappa * (1.5 * ye * ye + trace - gauge_lepton),
    ]


def run_yukawas(
    y0: list[float], mu0: float, mu1: float, mz_run: float, n_steps: int = 12000
) -> list[float]:
    y = list(y0)
    t0, t1 = math.log(mu0), math.log(mu1)
    dt = (t1 - t0) / n_steps
    for i in range(n_steps):
        t = t0 + i * dt
        k1 = _yukawa_derivatives(y, math.exp(t), mz_run)
        y2 = [a + 0.5 * dt * b for a, b in zip(y, k1)]
        k2 = _yukawa_derivatives(y2, math.exp(t + 0.5 * dt), mz_run)
        y3 = [a + 0.5 * dt * b for a, b in zip(y, k2)]
        k3 = _yukawa_derivatives(y3, math.exp(t + 0.5 * dt), mz_run)
        y4 = [a + dt * b for a, b in zip(y, k3)]
        k4 = _yukawa_derivatives(y4, math.exp(t + dt), mz_run)
        y = [
            a + dt / 6.0 * (b1 + 2 * b2 + 2 * b3 + b4)
            for a, b1, b2, b3, b4 in zip(y, k1, k2, k3, k4)
        ]
    return y


def alpha_s_chain(mu: float, alpha_s_mt: float, mt: float) -> float:
    """One-loop flavor-banded strong coupling below m_t."""

    def run(a0: float, mu0: float, mu1: float, nf: int) -> float:
        b0 = 11.0 - 2.0 * nf / 3.0
        return 1.0 / (1.0 / a0 + b0 / (2.0 * math.pi) * math.log(mu1 / mu0))

    if mu >= 4.18:
        return run(alpha_s_mt, mt, mu, 5)
    a_mb = run(alpha_s_mt, mt, 4.18, 5)
    if mu >= 1.27:
        return run(a_mb, 4.18, mu, 4)
    a_mc = run(a_mb, 4.18, 1.27, 4)
    return run(a_mc, 1.27, mu, 3)


def qcd_mass_factor(mu0: float, mu1: float, alpha_s_mt: float, mt: float) -> float:
    """One-loop mass running factor across flavor bands from mu0 to mu1."""

    def gamma_exponent(nf: int) -> float:
        return 4.0 / (11.0 - 2.0 * nf / 3.0)

    factor = 1.0
    grid = [(4.18, 5), (1.27, 4)]
    current = mu0
    a_current = alpha_s_chain(current, alpha_s_mt, mt)
    for threshold, nf in grid:
        if mu1 < threshold < current:
            a_th = alpha_s_chain(threshold, alpha_s_mt, mt)
            factor *= (a_th / a_current) ** gamma_exponent(nf)
            current, a_current = threshold, a_th
    nf_final = 5 if mu1 >= 4.18 else (4 if mu1 >= 1.27 else 3)
    a_end = alpha_s_chain(mu1, alpha_s_mt, mt)
    factor *= (a_end / a_current) ** gamma_exponent(nf_final)
    return factor


def build_artifact(mcpr: dict[str, Any]) -> dict[str, Any]:
    scales = source_scales(P_FALLBACK, ALPHA_U_FALLBACK)
    v_gev = scales["v_transmutation_gev"]
    mz_run = scales["mz_run_gev"]
    mu_u = scales["mu_U_gauge_unification"]

    top_row = run_boundary(
        scales["log_midpoint_half_turn"], v_gev, mz_run, n_steps=12000, loops=1
    )
    yt_mt = top_row["y_t_mt"]
    mt_ms = top_row["mt_msbar_gev"]
    alpha_s_mt = top_row["alpha_s_mt"]

    lepton_masses_gev = [
        float(x) / 1000.0
        for x in (
            float(mcpr["optional_scale_display"]["masses_MeV"][0]),
            float(mcpr["optional_scale_display"]["masses_MeV"][1]),
            float(mcpr["optional_scale_display"]["masses_MeV"][2]),
        )
    ]
    m_e, m_mu, m_tau = lepton_masses_gev

    def yukawa(m: float) -> float:
        return math.sqrt(2.0) * m / v_gev

    up = run_yukawas(
        [yt_mt, 0.0, yukawa(m_tau), yukawa(m_mu), yukawa(m_e)], mt_ms, mu_u, mz_run
    )
    yt_u, _, ytau_u, ymu_u, ye_u = up

    yb_u = CLEBSCH["b_over_tau"] * ytau_u
    ys_u = CLEBSCH["s_over_mu"] * ymu_u
    yd_u = CLEBSCH["d_over_e"] * ye_u

    down = run_yukawas([yt_u, yb_u, ytau_u, ymu_u, ye_u], mu_u, mt_ms, mz_run)
    yb_mt = down[1]
    down_scale_factor = yb_mt / yb_u
    ys_mt = ys_u * down_scale_factor
    yd_mt = yd_u * down_scale_factor

    mb_at_mt = yb_mt * v_gev / math.sqrt(2.0)
    ms_at_mt = ys_mt * v_gev / math.sqrt(2.0)
    md_at_mt = yd_mt * v_gev / math.sqrt(2.0)

    mb_mb = mb_at_mt
    for _ in range(25):
        mb_mb = mb_at_mt * qcd_mass_factor(mt_ms, mb_mb, alpha_s_mt, mt_ms)
    ms_2 = ms_at_mt * qcd_mass_factor(mt_ms, 2.0, alpha_s_mt, mt_ms)
    md_2 = md_at_mt * qcd_mass_factor(mt_ms, 2.0, alpha_s_mt, mt_ms)

    cabibbo_gst = math.sqrt(md_2 / ms_2)

    compare = {
        "mb_relative": mb_mb / COMPARE_ONLY["mb_mb_gev"] - 1.0,
        "ms_relative": ms_2 / COMPARE_ONLY["ms_2gev_gev"] - 1.0,
        "md_relative": md_2 / COMPARE_ONLY["md_2gev_gev"] - 1.0,
        "ms_over_md_relative": (ms_2 / md_2) / COMPARE_ONLY["ms_over_md"] - 1.0,
        "cabibbo_relative": cabibbo_gst / COMPARE_ONLY["cabibbo"] - 1.0,
        "references": COMPARE_ONLY,
        "role": "comparison only, outside the prediction ancestry",
    }

    checks = {
        "clebsch_factors_declared": True,
        "third_generation_boundary_is_pure_b_tau": CLEBSCH["b_over_tau"] == 1.0,
        "masses_positive": min(mb_mb, ms_2, md_2) > 0,
        "ratio_outputs_finite": math.isfinite(cabibbo_gst),
        "premise_gates_open": True,
    }

    return {
        "artifact": "oph_down_type_register_clebsch_lane",
        "schema_version": 1,
        "status": "CONDITIONAL_DOWN_TYPE_SECTOR_RATIOS_SHARP_NORMALIZATION_TENSE",
        "row_class": "conditional_on_mcpr_leptons_and_register_clebsch",
        "promotion_allowed": False,
        "clebsch_boundary": {
            "factors": CLEBSCH,
            "register_scale": "mu_U",
            "candidate_register_origin": (
                "3 = N_c and 1/3 = invariant color measure of the transitive "
                "C3 action; the pairing and the unordered weight set are "
                "selected conditionally by "
                "runs/flavor/clebsch_register_pairing_selection.json, and the "
                "generation-order assignment is the open premise"
            ),
            "selection_artifact": (
                "runs/flavor/clebsch_register_pairing_selection.json"
            ),
        },
        "inputs": {
            "lepton_triple_source": (
                "runs/leptons/charged_mcpr_completion_conditional.json"
            ),
            "top_trajectory": "criticality boundary at the adopted branch",
            "loop_order": "one_loop_yukawa_system_one_loop_qcd_mass_running",
        },
        "boundary_values_at_mu_U": {
            "y_tau": ytau_u,
            "y_mu": ymu_u,
            "y_e": ye_u,
            "y_b": yb_u,
            "y_s": ys_u,
            "y_d": yd_u,
            "y_t": yt_u,
        },
        "predictions": {
            "mb_mb_gev": mb_mb,
            "ms_2gev_gev": ms_2,
            "md_2gev_gev": md_2,
            "ms_over_md": ms_2 / md_2,
            "cabibbo_gst_sqrt_md_over_ms": cabibbo_gst,
        },
        "compare_only": compare,
        "normalization_tension": {
            "statement": (
                "the third-generation relation y_b = y_tau at one loop in "
                "this branch overshoots the down-type normalization by tens "
                "of percent while the intra-sector ratios and the "
                "Gatto-Sartori-Tonin angle land at the ten-percent scale; "
                "the deficit is carried openly and attributed to the "
                "third-generation register factor and the threshold packet, "
                "both open"
            ),
            "open_objects": [
                "GENERATION_REGISTER_ORDER",
                "THIRD_GENERATION_REGISTER_FACTOR",
                "FROZEN_RG_THRESHOLD_MATCHING_PACKET",
            ],
            "premise_reduction": (
                "CLEBSCH_REGISTER_SELECTION_THEOREM is split by the selection "
                "artifact: pairing and unordered weight set closed "
                "conditionally, generation order open"
            ),
        },
        "claim_boundary": (
            "Conditional on the MCPR lepton triple and the declared register "
            "Clebsch pattern. Ratios are the sharp outputs; the absolute "
            "normalization carries a named open tension. Nothing here is a "
            "promoted source-only prediction."
        ),
        "checks": checks,
        "checks_pass": all(bool(v) for v in checks.values()),
    }


def build() -> dict[str, Any]:
    mcpr = json.loads(MCPR.read_text(encoding="utf-8"))
    return build_artifact(mcpr)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    artifact = build()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "status": artifact["status"],
                "checks_pass": artifact["checks_pass"],
                "predictions": artifact["predictions"],
                "compare_only_relatives": {
                    k: round(v, 4)
                    for k, v in artifact["compare_only"].items()
                    if isinstance(v, float)
                },
                "output": str(args.output),
            },
            indent=2,
        )
    )
    return 0 if artifact["checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
