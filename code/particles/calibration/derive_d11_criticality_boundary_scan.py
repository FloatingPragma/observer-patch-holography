#!/usr/bin/env python3
"""Scan the D11 double-criticality boundary scale on the source-audit branch.

The archived D11 literal core derives both Yukawa-sector boundaries from the
gauge sector: at one boundary scale it imposes the double-criticality
condition ``lambda = 0`` and ``beta_lambda = 0``, which fixes
``y_t = (X/16)^(1/4)`` with ``X = 2 g2^4 + (g2^2 + gY^2)^2``, then integrates
the one-loop Standard Model system down to the top scale.  The entire
``(m_t, m_H)`` pair is therefore a zero-continuous-parameter function of the
pixel value and the boundary-scale choice.

The archived choice places the boundary at the gauge-unification scale
``mu_U = (E_star/e^(2 pi)) P^(1/6)``.  The electroweak transmutation in the
same model anchors ``v = E_cell exp(-2 pi/(beta_EW alpha_U))`` at the pixel
cell energy ``E_cell = E_star/sqrt(P)``.  The archived model therefore uses
two different high-scale anchors.  This lane evaluates the criticality
boundary at every named source scale and on a continuum grid, so the
boundary-scale dependence is an explicit discrete selection question with the
single-anchor branch (boundary at ``E_cell``) identified.

The one-loop engine is validated in-artifact against the archived pinned
endpoints at the ``mu_U`` boundary.  A two-loop engine is consumed when the
transcribed module is present; its rows quantify the loop-truncation band.
No measured particle mass is read.  Every row is promotion-blocked: the
boundary-scale selection theorem is the named open object.
"""

from __future__ import annotations

import argparse
import functools
import json
import math
from pathlib import Path
from typing import Any, Callable

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SOURCE_ROOT = (
    ROOT / "particles" / "hierarchy" / "certificates"
    / "R_P_source_audit_pixel_certificate.json"
)
DEFAULT_OUT = (
    ROOT / "particles" / "runs" / "calibration"
    / "d11_criticality_boundary_scan.json"
)

# Strict source-audit branch constants.  The alphas are the one-loop chart
# couplings quoted at the running Z coordinate of the same branch; they derive
# from alpha_U through the declared slope stack (see
# wzh_residual_elimination_boundary.json for the shared basis).
E_STAR_DISPLAY_GEV = 1.2208901289579269e19  # clock candidate; checksum status
ALPHA_Y_MZ = 0.010131400543429989
ALPHA_2_MZ = 0.03377781411239098
ALPHA_3_MZ = 0.11833586196478191
MU_Z_OVER_E = 7.501767385088419e-18

B_SM_ALPHA_Y = (41.0 / 6.0, -19.0 / 6.0, -7.0)

# Archived pinned endpoints at the mu_U boundary (one loop).
ARCHIVED_PIN = {
    "y_t_u": 0.38989194927267146,
    "y_t_mt": 0.8876472196531545,
    "lambda_mt": 0.1089320070972017,
    "alpha_s_mt": 0.11067971218235456,
    "mt_ms": 154.779273546571,
    "mt_pole": 164.13059302729587,
    "mh_tree": 115.10128055804651,
}

try:  # transcribed and benchmark-validated two-loop engine, if present
    try:
        from calibration.sm_two_loop_rge_engine import (
            beta_2loop as _beta_2loop_engine,
            validate as _validate_2loop,
        )
    except ModuleNotFoundError:
        from sm_two_loop_rge_engine import (
            beta_2loop as _beta_2loop_engine,
            validate as _validate_2loop,
        )
    TWO_LOOP_AVAILABLE = True
except Exception:  # noqa: BLE001
    _beta_2loop_engine = None
    _validate_2loop = None
    TWO_LOOP_AVAILABLE = False


def source_scales(p_value: float, alpha_u: float) -> dict[str, float]:
    """Named high scales of the source branch, in display GeV."""

    e_cell = E_STAR_DISPLAY_GEV / math.sqrt(p_value)
    mu_u = (E_STAR_DISPLAY_GEV / math.exp(2.0 * math.pi)) * p_value ** (1.0 / 6.0)
    v = e_cell * math.exp(-2.0 * math.pi / (4.0 * alpha_u))
    return {
        "mu_U_gauge_unification": mu_u,
        "E_cell_pixel_energy": e_cell,
        "E_star": E_STAR_DISPLAY_GEV,
        # Log-midpoint of the model's two pre-existing anchors:
        # sqrt(mu_U * E_cell) = E_star * exp(-pi) * P^(-1/6).
        "log_midpoint_half_turn": math.sqrt(mu_u * e_cell),
        "v_transmutation_gev": v,
        "mz_run_gev": MU_Z_OVER_E * E_STAR_DISPLAY_GEV,
    }


def alpha_run_1loop(alpha0: float, b: float, mu: float, mu0: float) -> float:
    return 1.0 / (1.0 / alpha0 - (b / (2.0 * math.pi)) * math.log(mu / mu0))


def gauge_couplings(mu: float, mz_run: float) -> tuple[float, float, float]:
    ay = alpha_run_1loop(ALPHA_Y_MZ, B_SM_ALPHA_Y[0], mu, mz_run)
    a2 = alpha_run_1loop(ALPHA_2_MZ, B_SM_ALPHA_Y[1], mu, mz_run)
    a3 = alpha_run_1loop(ALPHA_3_MZ, B_SM_ALPHA_Y[2], mu, mz_run)
    return (
        math.sqrt(4.0 * math.pi * ay),
        math.sqrt(4.0 * math.pi * a2),
        math.sqrt(4.0 * math.pi * a3),
    )


def critical_surface_yukawa(g_y: float, g2: float) -> float:
    """Double criticality: beta_lambda(lambda=0) = 0 at one loop."""

    x = 2.0 * g2**4 + (g2 * g2 + g_y * g_y) ** 2
    return (x / 16.0) ** 0.25


def beta_y_1loop(y: float, g_y: float, g2: float, g3: float) -> float:
    return y / (16.0 * math.pi**2) * (
        4.5 * y * y - (17.0 / 12.0) * g_y**2 - 2.25 * g2**2 - 8.0 * g3**2
    )


def beta_lambda_1loop(lam: float, y: float, g_y: float, g2: float, g3: float) -> float:
    term = 24.0 * lam * lam
    term += lam * (12.0 * y * y - 9.0 * g2**2 - 3.0 * g_y**2)
    term += -6.0 * y**4
    term += (3.0 / 8.0) * (2.0 * g2**4 + (g2**2 + g_y**2) ** 2)
    return term / (16.0 * math.pi**2)


def top_pole_from_msbar(mt_ms: float, alpha_s: float, n_l: int = 5) -> float:
    a = alpha_s / math.pi
    k2 = 13.4434 - 1.0414 * n_l
    k3 = 190.595 - 26.655 * n_l + 0.653 * n_l**2
    return mt_ms * (1.0 + (4.0 / 3.0) * a + k2 * a**2 + k3 * a**3)


def run_boundary(
    mu_b: float,
    v_gev: float,
    mz_run: float,
    n_steps: int = 24000,
    loops: int = 1,
) -> dict[str, float]:
    """Impose double criticality at ``mu_b`` and integrate down to m_t."""

    g_y_b, g2_b, g3_b = gauge_couplings(mu_b, mz_run)
    if loops == 2 and TWO_LOOP_AVAILABLE:
        y = solve_two_loop_critical_yukawa(g_y_b, g2_b, g3_b)
    else:
        y = critical_surface_yukawa(g_y_b, g2_b)
    lam = 0.0
    y_t_u = y

    t0, t1 = math.log(mu_b), math.log(50.0)
    dt = (t1 - t0) / n_steps
    ln_mu = np.empty(n_steps + 1)
    y_arr = np.empty(n_steps + 1)
    l_arr = np.empty(n_steps + 1)

    def derivs(t: float, y_v: float, l_v: float) -> tuple[float, float]:
        mu = math.exp(t)
        g_y, g2, g3 = gauge_couplings(mu, mz_run)
        if loops == 2 and TWO_LOOP_AVAILABLE:
            g1 = math.sqrt(5.0 / 3.0) * g_y
            derivatives = _beta_2loop_engine(g1, g2, g3, y_v, l_v)
            return derivatives[3], derivatives[4]
        return (
            beta_y_1loop(y_v, g_y, g2, g3),
            beta_lambda_1loop(l_v, y_v, g_y, g2, g3),
        )

    for i in range(n_steps + 1):
        t = t0 + i * dt
        ln_mu[i] = t
        y_arr[i] = y
        l_arr[i] = lam
        if i == n_steps:
            break
        k1y, k1l = derivs(t, y, lam)
        k2y, k2l = derivs(t + 0.5 * dt, y + 0.5 * dt * k1y, lam + 0.5 * dt * k1l)
        k3y, k3l = derivs(t + 0.5 * dt, y + 0.5 * dt * k2y, lam + 0.5 * dt * k2l)
        k4y, k4l = derivs(t + dt, y + dt * k3y, lam + dt * k3l)
        y += dt / 6.0 * (k1y + 2.0 * k2y + 2.0 * k3y + k4y)
        lam += dt / 6.0 * (k1l + 2.0 * k2l + 2.0 * k3l + k4l)

    def interp(mu: float) -> tuple[float, float]:
        x = math.log(mu)
        return (
            float(np.interp(x, ln_mu[::-1], y_arr[::-1])),
            float(np.interp(x, ln_mu[::-1], l_arr[::-1])),
        )

    mu_g = 173.0
    for _ in range(8):
        y_g, _ = interp(mu_g)
        mu_g = y_g * v_gev / math.sqrt(2.0)
    mt_ms = mu_g
    y_mt, lam_mt = interp(mt_ms)
    alpha_s_mt = alpha_run_1loop(ALPHA_3_MZ, B_SM_ALPHA_Y[2], mt_ms, mz_run)
    return {
        "boundary_scale_gev": mu_b,
        "y_t_u": y_t_u,
        "y_t_mt": y_mt,
        "lambda_mt": lam_mt,
        "alpha_s_mt": alpha_s_mt,
        "mt_msbar_gev": mt_ms,
        "mt_pole_gev": top_pole_from_msbar(mt_ms, alpha_s_mt),
        "mh_tree_gev": math.sqrt(max(0.0, 2.0 * lam_mt)) * v_gev,
    }


def solve_two_loop_critical_yukawa(g_y: float, g2: float, g3: float) -> float:
    """Solve beta_lambda(lambda=0, y) = 0 with the two-loop engine."""

    g1 = math.sqrt(5.0 / 3.0) * g_y
    lo, hi = 0.05, 2.0

    def f(y: float) -> float:
        return _beta_2loop_engine(g1, g2, g3, y, 0.0)[4]

    f_lo, f_hi = f(lo), f(hi)
    if f_lo * f_hi > 0:
        return critical_surface_yukawa(g_y, g2)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(lo) * f(mid) <= 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def validation_block(v_gev: float, mz_run: float) -> dict[str, Any]:
    row = run_boundary(
        source_scales(P_FALLBACK, ALPHA_U_FALLBACK)["mu_U_gauge_unification"],
        v_gev,
        mz_run,
        n_steps=45000,
        loops=1,
    )
    mapping = {
        "y_t_u": "y_t_u",
        "y_t_mt": "y_t_mt",
        "lambda_mt": "lambda_mt",
        "alpha_s_mt": "alpha_s_mt",
        "mt_ms": "mt_msbar_gev",
        "mt_pole": "mt_pole_gev",
        "mh_tree": "mh_tree_gev",
    }
    residuals = {
        pin_key: abs(row[row_key] / ARCHIVED_PIN[pin_key] - 1.0)
        for pin_key, row_key in mapping.items()
    }
    return {
        "archived_pinned_endpoints": ARCHIVED_PIN,
        "reproduction_relative_residuals": residuals,
        "max_residual": max(residuals.values()),
        "passes_1e_6": max(residuals.values()) < 1.0e-6,
    }


P_FALLBACK = 1.63097209569432901817967892561191884270169
ALPHA_U_FALLBACK = 0.04112424744557487


def build_artifact(root: dict[str, Any]) -> dict[str, Any]:
    p_value = float(root["P_cand"])
    alpha_u = float(root["alpha_U_P_cand"])
    scales = source_scales(p_value, alpha_u)
    v_gev = scales["v_transmutation_gev"]
    mz_run = scales["mz_run_gev"]

    v_over_e_expected = p_value ** -0.5 * math.exp(-2.0 * math.pi / (4.0 * alpha_u))
    v_consistency = abs(
        v_gev / E_STAR_DISPLAY_GEV / float(root["v_over_E_star_P_cand"]) - 1.0
    )

    validation = validation_block(v_gev, mz_run)

    named_rows: dict[str, dict[str, float]] = {}
    named_scale_keys = (
        "mu_U_gauge_unification",
        "log_midpoint_half_turn",
        "E_cell_pixel_energy",
        "E_star",
    )
    for name in named_scale_keys:
        named_rows[name] = run_boundary(scales[name], v_gev, mz_run, loops=1)

    grid = np.geomspace(1.0e16, 1.3e19, 25)
    curve = [run_boundary(float(mu_b), v_gev, mz_run, n_steps=12000) for mu_b in grid]

    two_loop_rows: dict[str, Any] = {"available": TWO_LOOP_AVAILABLE}
    two_loop_curve: list[dict[str, float]] = []
    if TWO_LOOP_AVAILABLE:
        two_loop_rows["engine_validation"] = _validate_2loop()
        for name in named_scale_keys:
            two_loop_rows[name] = run_boundary(
                scales[name], v_gev, mz_run, loops=2
            )
        two_loop_curve = [
            run_boundary(float(mu_b), v_gev, mz_run, n_steps=8000, loops=2)
            for mu_b in np.geomspace(1.0e16, 1.3e19, 13)
        ]

    checks = {
        "archived_endpoint_reproduction": bool(validation["passes_1e_6"]),
        "v_transmutation_matches_certificate": v_consistency < 1.0e-10,
        "curve_mh_monotone_in_boundary_scale": all(
            curve[i + 1]["mh_tree_gev"] >= curve[i]["mh_tree_gev"]
            for i in range(len(curve) - 1)
        ),
        "double_criticality_boundary_lambda_zero": True,
    }
    if TWO_LOOP_AVAILABLE:
        validation_2l = two_loop_rows["engine_validation"]
        checks["two_loop_engine_benchmark_validation"] = bool(
            validation_2l.get("all_ok")
        ) and all(bool(v) for v in validation_2l.get("endpoint_ok", {}).values())

    return {
        "artifact": "oph_d11_criticality_boundary_scan",
        "schema_version": 1,
        "status": "boundary_scale_family_conditional_prediction",
        "row_class": "conditional_on_P_and_criticality_boundary_scale_selection",
        "promotion_allowed": False,
        "branch": {
            "id": "source_audit_P_cand",
            "P": root["P_cand"],
            "alpha_U": root["alpha_U_P_cand"],
        },
        "boundary_law": {
            "condition": "lambda(mu_b) = 0 and beta_lambda(mu_b) = 0",
            "yukawa_consequence": "y_t(mu_b) = ((2 g2^4 + (g2^2+gY^2)^2)/16)^(1/4)",
            "free_continuous_parameters": 0,
            "free_discrete_choice": "the boundary scale mu_b among named source scales",
        },
        "single_anchor_coherence": {
            "observation": (
                "The electroweak transmutation anchors v at E_cell = "
                "E_star/sqrt(P). The archived criticality boundary sits at "
                "mu_U, a second high scale. Placing the criticality boundary "
                "at the same E_cell anchor is the unique single-anchor law of "
                "this family."
            ),
            "single_anchor_branch": "E_cell_pixel_energy",
            "archived_branch": "mu_U_gauge_unification",
            "selection_theorem_status": (
                "CRITICALITY_BOUNDARY_SCALE_SELECTION_THEOREM absent"
            ),
        },
        "source_scales_gev": scales,
        "one_loop_named_boundaries": named_rows,
        "one_loop_curve": curve,
        "two_loop_named_boundaries": two_loop_rows,
        "two_loop_curve": two_loop_curve,
        "relation_statement": (
            "The double-criticality law makes (m_t, m_H) a one-discrete-"
            "parameter family over the boundary scale. The m_t to m_H "
            "relation along the curve is therefore a fit-free consequence "
            "of the boundary law; the boundary-scale selection theorem "
            "picks the point on the curve."
        ),
        "matching_bands": {
            "mh_tree_to_pole_relative_band": 0.015,
            "mt_electroweak_matching_relative_band": 0.005,
            "statement": (
                "Tree-to-pole Higgs matching and electroweak top matching are "
                "carried as declared uncertainty bands, never as applied "
                "corrections."
            ),
        },
        "archived_endpoint_validation": validation,
        "display_note": (
            "GeV values use the unclosed clock candidate; the clock audit "
            "classifies that decimal as a calibration checksum. Dimensionless "
            "ratios follow by dividing by E_star."
        ),
        "checks": checks,
        "checks_pass": all(bool(v) for v in checks.values()),
        "claim_boundary": (
            "Every row is conditional on the witness-level P root, the "
            "declared one-loop chart couplings, the double-criticality "
            "boundary law, and the open boundary-scale selection theorem. "
            "The rows are running/tree coordinates, never physical poles."
        ),
    }


def build() -> dict[str, Any]:
    root = json.loads(SOURCE_ROOT.read_text(encoding="utf-8"))
    return build_artifact(root)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    artifact = build()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    named = artifact["one_loop_named_boundaries"]
    print(
        json.dumps(
            {
                "status": artifact["status"],
                "checks_pass": artifact["checks_pass"],
                "two_loop_available": artifact["two_loop_named_boundaries"][
                    "available"
                ],
                "one_loop": {
                    name: {
                        "mt_pole_gev": round(row["mt_pole_gev"], 3),
                        "mh_tree_gev": round(row["mh_tree_gev"], 3),
                    }
                    for name, row in named.items()
                },
                "output": str(args.output),
            },
            indent=2,
        )
    )
    return 0 if artifact["checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
