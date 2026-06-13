#!/usr/bin/env python3
"""Joint (P, N_CRC) fixed-point certificate helper for OPH issue #338."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mpmath as mp

mp.mp.dps = 100

P_PUBLIC = mp.mpf(
    "1.6309682094039593248792798477826489413359828516279250606661507533907793398933432"
)
P_CAND_SOURCE = mp.mpf("1.63097209569432901817967892561191884270169")
ALPHA_CAND_INV = mp.mpf("136.994835164621649457949994585787193262029")
SOURCE_V_OVER_ESTAR = mp.mpf("2.0199803239725553e-17")
N_DISPLAY = mp.mpf("3.31e122")


def s(x: mp.mpf, digits: int = 50) -> str:
    return mp.nstr(x, digits)


def alpha_from_p(P: mp.mpf) -> mp.mpf:
    phi = (1 + mp.sqrt(5)) / 2
    return (P - phi) / mp.sqrt(mp.pi)


def tick_factor(N: mp.mpf) -> mp.mpf:
    return (N / mp.pi) ** (-mp.mpf(1) / 48)


def v_over_ecell_resonance(P: mp.mpf, N: mp.mpf) -> mp.mpf:
    return (N / mp.pi) ** (-P / 12)


def v_over_estar_resonance(P: mp.mpf, N: mp.mpf) -> mp.mpf:
    return P ** (-mp.mpf(1) / 2) * v_over_ecell_resonance(P, N)


def n_backsolved_from_source_v(P: mp.mpf, v_over_estar: mp.mpf) -> mp.mpf:
    v_over_ecell = v_over_estar * mp.sqrt(P)
    return mp.pi * (1 / v_over_ecell) ** (12 / P)


def product_contraction(qP: mp.mpf, qN: mp.mpf) -> dict[str, object]:
    q = max(qP, qN)
    return {"q_product": s(q, 30), "passes": bool(q < 1), "qP": s(qP, 30), "qN": s(qN, 30)}


def coupled_contraction(a: mp.mpf, b: mp.mpf, c: mp.mpf, d: mp.mpf, r: mp.mpf) -> dict[str, object]:
    q = max(a + b / r, d + r * c)
    return {
        "a": s(a, 30),
        "b": s(b, 30),
        "c": s(c, 30),
        "d": s(d, 30),
        "r": s(r, 30),
        "q_weighted_sup": s(q, 30),
        "passes": bool(q < 1),
        "condition": "max(a+b/r, d+r*c) < 1",
    }


def numeric_report(N: mp.mpf = N_DISPLAY) -> dict[str, object]:
    alpha = alpha_from_p(P_PUBLIC)
    return {
        "issue": 338,
        "status": "closed_product_branch_theorem_with_explicit_coupled_branch_boundary",
        "theorem": "joint fixed-point and stability theorem for the product-separated (P,N_CRC) source branch",
        "P_public": s(P_PUBLIC, 90),
        "alpha_from_P_public": s(alpha, 90),
        "alpha_inverse_from_P_public": s(1 / alpha, 90),
        "P_source_audit": s(P_CAND_SOURCE, 70),
        "alpha_source_audit_inverse": s(ALPHA_CAND_INV, 70),
        "N_display_used_for_readout_demo": s(N, 30),
        "N_display_warning": "3.31e122 is a rounded capacity-scale label, not a high-precision N_CRC certificate.",
        "global_tick_gprime_at_N_display": s(tick_factor(N), 50),
        "v_over_Ecell_resonance_at_display_N": s(v_over_ecell_resonance(P_PUBLIC, N), 50),
        "v_over_Estar_resonance_at_display_N": s(v_over_estar_resonance(P_PUBLIC, N), 50),
        "log_sensitivity_tick_to_logN": s(-mp.mpf(1) / 48, 30),
        "log_sensitivity_vEcell_to_logN": s(-P_PUBLIC / 12, 30),
        "log_sensitivity_vEcell_to_P_at_display_N": s(-mp.log(N / mp.pi) / 12, 30),
        "log_sensitivity_vEstar_to_P_at_display_N": s(-1 / (2 * P_PUBLIC) - mp.log(N / mp.pi) / 12, 30),
        "N_backsolved_from_source_v_over_Estar": s(n_backsolved_from_source_v(P_PUBLIC, SOURCE_V_OVER_ESTAR), 40),
        "N_backsolved_warning": "CIRCULAR_DIAGNOSTIC_ONLY: this solves N from a weak-scale readout and is not allowed as an OPH source certificate.",
        "claim_boundary": {
            "product_branch": "The source map is J(P,x)=(Gamma(P),C_hat(x)) on I_P x I_x. Component contraction certificates imply a unique joint fixed point with q=max(qP,qN)<1.",
            "coupled_branch": "If a future source map contains cross-feedback, supply derivative bounds a,b,c,d and r>0 with max(a+b/r,d+r*c)<1, or keep the coupled branch as residual freedom.",
            "forbidden_input": "The weak scale, Higgs mass, W/Z masses, G, Planck area, and measured Lambda are not allowed as fixed-point inputs.",
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--N", default=str(N_DISPLAY), help="N value for readout demo; default rounded 3.31e122")
    parser.add_argument("--qP", default=None, help="optional local P contraction constant")
    parser.add_argument("--qN", default=None, help="optional global log-N contraction constant")
    parser.add_argument("--a", default=None, help="coupled map bound |dP'|_P")
    parser.add_argument("--b", default=None, help="coupled map bound |dP'|_x")
    parser.add_argument("--c", default=None, help="coupled map bound |dx'|_P")
    parser.add_argument("--d", default=None, help="coupled map bound |dx'|_x")
    parser.add_argument("--r", default="1", help="weighted metric parameter for coupled map")
    parser.add_argument("--output", default=None, help="write JSON report")
    args = parser.parse_args()

    N = mp.mpf(args.N)
    out = numeric_report(N)

    if args.qP is not None and args.qN is not None:
        out["product_contraction_certificate"] = product_contraction(mp.mpf(args.qP), mp.mpf(args.qN))
    else:
        out["product_contraction_certificate"] = {
            "status": "conditional_on_component_contractions",
            "required": "qP<1 and qN<1 from the local P and global N source certificates; then q=max(qP,qN)<1.",
        }

    if all(v is not None for v in [args.a, args.b, args.c, args.d]):
        out["coupled_contraction_certificate"] = coupled_contraction(
            mp.mpf(args.a), mp.mpf(args.b), mp.mpf(args.c), mp.mpf(args.d), mp.mpf(args.r)
        )
    else:
        out["coupled_contraction_certificate"] = {
            "status": "residual_coupled_branch_boundary",
            "required": "find r>0 such that max(a+b/r, d+r*c)<1, or equivalently a weighted-sup contraction envelope.",
        }

    text = json.dumps(out, indent=2) + "\n"
    if args.output:
        Path(args.output).write_text(text, encoding="utf-8")
    print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
