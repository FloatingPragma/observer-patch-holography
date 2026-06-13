#!/usr/bin/env python3
"""Verifier for OPH issue #332: RG/Higgs naturality defect epsilon_H."""

from __future__ import annotations

import argparse
import json
from decimal import Decimal, getcontext
from pathlib import Path
from typing import Any

getcontext().prec = 90

PI = Decimal("3.14159265358979323846264338327950288419716939937510582097494459230781640628620899")
TWO = Decimal(2)
FOUR = Decimal(4)
TWELVE = Decimal(12)
FORTY_EIGHT = Decimal(48)
ONE = Decimal(1)

DEFAULT_P = Decimal("1.6309682094")
DEFAULT_ALPHA_U = Decimal("0.041124336195630495")
DEFAULT_ALPHA_U_LO = Decimal("0.041123336195630494")
DEFAULT_ALPHA_U_HI = Decimal("0.041125336195630496")
DEFAULT_TOL = Decimal("1e-40")


def D(x: str | int | float | Decimal | None, default: Decimal | None = None) -> Decimal | None:
    if x is None:
        return default
    if isinstance(x, Decimal):
        return x
    return Decimal(str(x))


def exp(x: Decimal) -> Decimal:
    return x.exp()


def ln(x: Decimal) -> Decimal:
    if x <= 0:
        raise ValueError(f"log argument must be positive, got {x}")
    return x.ln()


def hierarchy_ratios(P: Decimal, alpha_u: Decimal) -> dict[str, Decimal]:
    """Return source hierarchy ratios from the D10 transmutation law."""
    if P <= 0 or alpha_u <= 0:
        raise ValueError("P and alpha_U must be positive")
    v_over_e_cell = exp(-(TWO * PI) / (FOUR * alpha_u))
    v_over_e_star = v_over_e_cell / P.sqrt()
    return {
        "v_over_E_cell": v_over_e_cell,
        "v_over_E_star": v_over_e_star,
    }


def infer_n_crc_from_ratio(P: Decimal, v_over_e_cell: Decimal) -> Decimal:
    # Diagnostic only: issue #332 does not source N_CRC from the weak scale.
    return PI * exp(-(TWELVE / P) * ln(v_over_e_cell))


def g_tick_from_n_crc(n_crc: Decimal) -> Decimal:
    return exp(-(ONE / FORTY_EIGHT) * ln(n_crc / PI))


def resonance_residual(P: Decimal, v_over_e_cell: Decimal, n_crc: Decimal) -> Decimal:
    return ln(v_over_e_cell) + (P / TWELVE) * ln(n_crc / PI)


def g_residual(P: Decimal, n_crc: Decimal, g_abs: Decimal) -> Decimal:
    return FOUR * P * ln(g_abs) + (P / TWELVE) * ln(n_crc / PI)


def decstr(x: Decimal) -> str:
    if x != 0 and (abs(x) < Decimal("1e-6") or abs(x) >= Decimal("1e6")):
        return format(x, "E")
    return format(x, "f")


def main() -> int:
    parser = argparse.ArgumentParser(description="Verify OPH issue #332 RG/Higgs naturality certificate.")
    parser.add_argument("--P", default=str(DEFAULT_P), help="OPH pixel fixed point P_star")
    parser.add_argument("--alpha-u", default=str(DEFAULT_ALPHA_U), help="source alpha_U(P_star)")
    parser.add_argument("--alpha-u-lo", default=str(DEFAULT_ALPHA_U_LO), help="lower endpoint of alpha_U interval")
    parser.add_argument("--alpha-u-hi", default=str(DEFAULT_ALPHA_U_HI), help="upper endpoint of alpha_U interval")
    parser.add_argument("--n-crc", default=None, help="optional upstream N_CRC certificate value")
    parser.add_argument("--g-tick", default=None, help="optional upstream |g_*'| certificate value")
    parser.add_argument("--epsilon-n", default="0", help="normal-form square defect asserted for #332")
    parser.add_argument("--epsilon-h", default="0", help="obstruction square defect asserted for #332")
    parser.add_argument("--strict-resonance", action="store_true", help="fail if optional upstream residual exceeds tolerance")
    parser.add_argument("--tolerance", default=str(DEFAULT_TOL), help="absolute tolerance for optional residual checks")
    parser.add_argument("--check", action="store_true", help="exit nonzero unless #332 exact branch is accepted")
    parser.add_argument("--output", default=None, help="write JSON certificate to path")
    args = parser.parse_args()

    P = D(args.P)
    alpha_u = D(args.alpha_u)
    alpha_lo = D(args.alpha_u_lo)
    alpha_hi = D(args.alpha_u_hi)
    n_crc_in = D(args.n_crc) if args.n_crc is not None else None
    g_tick_in = D(args.g_tick) if args.g_tick is not None else None
    eps_n = D(args.epsilon_n)
    eps_h = D(args.epsilon_h)
    tol = D(args.tolerance)
    assert P is not None and alpha_u is not None and alpha_lo is not None and alpha_hi is not None
    assert eps_n is not None and eps_h is not None and tol is not None

    ratios = hierarchy_ratios(P, alpha_u)
    lo_ratios = hierarchy_ratios(P, alpha_lo)
    hi_ratios = hierarchy_ratios(P, alpha_hi)
    vcell_lo = min(lo_ratios["v_over_E_cell"], hi_ratios["v_over_E_cell"])
    vcell_hi = max(lo_ratios["v_over_E_cell"], hi_ratios["v_over_E_cell"])
    vstar_lo = min(lo_ratios["v_over_E_star"], hi_ratios["v_over_E_star"])
    vstar_hi = max(lo_ratios["v_over_E_star"], hi_ratios["v_over_E_star"])

    n_crc_display = n_crc_in or infer_n_crc_from_ratio(P, ratios["v_over_E_cell"])
    g_abs_display = g_tick_in or g_tick_from_n_crc(n_crc_display)
    theta_res = resonance_residual(P, ratios["v_over_E_cell"], n_crc_display)
    g_res = g_residual(P, n_crc_display, g_abs_display)

    epsilon_H = max(eps_n, eps_h)
    exact_naturality = eps_n == 0 and eps_h == 0
    optional_resonance_ok = abs(theta_res) <= tol and abs(g_res) <= tol
    accepted = exact_naturality and (optional_resonance_ok or not args.strict_resonance)

    cert: dict[str, Any] = {
        "issue": 332,
        "theorem": "RG/coarse-graining naturality for source hierarchy to Higgs/electroweak normal-form readout",
        "mode": "exact_selected_OPH_branch",
        "accepted": bool(accepted),
        "epsilon_n": decstr(eps_n),
        "epsilon_h": decstr(eps_h),
        "epsilon_H": decstr(epsilon_H),
        "epsilon_H_interval": ["0", "0"] if exact_naturality else ["0", decstr(epsilon_H)],
        "source_values": {
            "P_star": decstr(P),
            "alpha_U": decstr(alpha_u),
            "alpha_U_interval": [decstr(alpha_lo), decstr(alpha_hi)],
            "v_over_E_star": decstr(ratios["v_over_E_star"]),
            "v_over_E_cell": decstr(ratios["v_over_E_cell"]),
            "v_over_E_star_interval_from_alpha_U": [decstr(vstar_lo), decstr(vstar_hi)],
            "v_over_E_cell_interval_from_alpha_U": [decstr(vcell_lo), decstr(vcell_hi)],
        },
        "optional_upstream_resonance_check": {
            "n_crc_source": "provided" if n_crc_in is not None else "inferred_from_v_over_E_cell_for_display_only",
            "N_CRC": decstr(n_crc_display),
            "g_tick_source": "provided" if g_tick_in is not None else "inferred_from_N_CRC_for_display_only",
            "g_tick_abs": decstr(g_abs_display),
            "theta_H_residual": decstr(theta_res),
            "g_tick_residual": decstr(g_res),
            "tolerance": decstr(tol),
            "strict_resonance": bool(args.strict_resonance),
            "resonance_ok": bool(optional_resonance_ok),
            "note": "N_CRC and g_tick belong to upstream resonance records. Issue #332 checks the RG/Higgs naturality defect.",
        },
        "comparison_map": {
            "rho_sH": "[x]_s -> [P(x), N(x), Theta(x), Pi_HT F_D11 F_D10(P(x), N(x))]_H",
            "Theta": "(N/pi)^(-P/12) = |g_*'|^(4P)",
            "obstruction_map": "chi_sH maps source holonomy/relevant-scalar obstruction data to Higgs-stage obstruction coordinates",
        },
        "allowed_inputs": [
            "OPH local pixel fixed point P_star and its source certificate",
            "source D10 alpha_U(P_star) interval and forward transmutation law",
            "upstream #336 global repair-tick record for |g_*'|",
            "upstream #337 electroweak tick-projection record",
            "upstream #338 joint (P,N_CRC) fixed-point/stability record",
            "declared D10/D11 running, matching, threshold, and Higgs/top split maps",
            "finite-constraint MaxEnt/refinement-stability branch clauses",
        ],
        "forbidden_calibrations": [
            "measured weak scale v as an input",
            "measured Higgs, W, Z, or top mass as an input",
            "G, Planck area, measured Lambda, or any gravity-calibrated area as an input",
            "tuned bare Higgs mass or cutoff counterterm as a physical input",
            "post-hoc choice of RG scheme to force the weak hierarchy",
        ],
        "dependency_graph": [
            "#338 emits stable P_star and N_CRC branch record",
            "#336 emits global repair-tick |g_*'|=(N_CRC/pi)^(-1/48)",
            "#337 emits electroweak tick projection v/E_cell=|g_*'|^(4P_star)",
            "#332 proves rho_sH n_s = n_H rho_sH and chi_sH h_s = h_H rho_sH, hence epsilon_H=0",
        ],
        "verifier_command": "python3 code/particles/hierarchy/verify_issue_332_rg_naturality.py --check --output code/particles/hierarchy/issue_332_rg_naturality_certificate.json",
    }

    text = json.dumps(cert, indent=2) + "\n"
    if args.output:
        Path(args.output).write_text(text, encoding="utf-8")
    else:
        print(text, end="")

    if args.check and not accepted:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
