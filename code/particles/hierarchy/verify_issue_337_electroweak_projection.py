#!/usr/bin/env python3
"""Verifier for OPH issue #337: electroweak tick-projection bridge."""

from __future__ import annotations

import argparse
import json
from decimal import Decimal, getcontext
from pathlib import Path
from typing import Any

getcontext().prec = 110

PI = Decimal("3.14159265358979323846264338327950288419716939937510582097494459230781640628620899")
TWO = Decimal(2)
FOUR = Decimal(4)
SIX = Decimal(6)
TWELVE = Decimal(12)
TWENTY_FOUR = Decimal(24)
FORTY_EIGHT = Decimal(48)
ZERO = Decimal(0)

DEFAULT_P = Decimal("1.6309682094039593248792798477826489413359828516279250606661507533907793398933432")
DEFAULT_ALPHA_U = Decimal("0.041124336195630495")
DEFAULT_ALPHA_U_LO = Decimal("0.041123336195630494")
DEFAULT_ALPHA_U_HI = Decimal("0.041125336195630496")
DEFAULT_ROUNDED_N = Decimal("3.31e122")
DEFAULT_TOL = Decimal("1e-40")


def D(value: str | int | Decimal | None, default: Decimal | None = None) -> Decimal | None:
    if value is None:
        return default
    if isinstance(value, Decimal):
        return value
    return Decimal(str(value))


def ln(value: Decimal) -> Decimal:
    if value <= 0:
        raise ValueError(f"log argument must be positive, got {value}")
    return value.ln()


def exp(value: Decimal) -> Decimal:
    return value.exp()


def decstr(value: Decimal) -> str:
    if value != 0 and (abs(value) < Decimal("1e-6") or abs(value) >= Decimal("1e6")):
        return format(value, "E")
    return format(value, "f")


def hierarchy_ratios(p_star: Decimal, alpha_u: Decimal) -> dict[str, Decimal]:
    v_over_e_cell = exp(-(TWO * PI) / (FOUR * alpha_u))
    return {
        "v_over_E_cell": v_over_e_cell,
        "v_over_E_star": v_over_e_cell / p_star.sqrt(),
        "t_tr": -ln(v_over_e_cell),
    }


def bridge_log_n_over_pi(p_star: Decimal, alpha_u: Decimal) -> Decimal:
    return SIX * PI / (p_star * alpha_u)


def n_from_log_over_pi(log_n_over_pi: Decimal) -> Decimal:
    return PI * exp(log_n_over_pi)


def g_tick_from_log(log_n_over_pi: Decimal) -> Decimal:
    return exp(-log_n_over_pi / FORTY_EIGHT)


def projection_exponent(alpha_u: Decimal, log_n_over_pi: Decimal) -> Decimal:
    return TWENTY_FOUR * PI / (alpha_u * log_n_over_pi)


def bridge_residual(alpha_u: Decimal, p_star: Decimal, log_n_over_pi: Decimal) -> Decimal:
    return alpha_u * log_n_over_pi - SIX * PI / p_star


def exact_branch(p_star: Decimal, alpha_u: Decimal) -> dict[str, Decimal]:
    log_n = bridge_log_n_over_pi(p_star, alpha_u)
    ratios = hierarchy_ratios(p_star, alpha_u)
    n_crc = n_from_log_over_pi(log_n)
    g_abs = g_tick_from_log(log_n)
    target_exponent = FOUR * p_star
    projected = projection_exponent(alpha_u, log_n)
    v_from_n = exp(-(p_star / TWELVE) * log_n)
    return {
        "log_N_over_pi": log_n,
        "N_EW": n_crc,
        "g_tick_abs": g_abs,
        "projection_exponent": projected,
        "target_exponent": target_exponent,
        "projection_exponent_error": projected - target_exponent,
        "bridge_residual": bridge_residual(alpha_u, p_star, log_n),
        "v_over_E_cell": ratios["v_over_E_cell"],
        "v_over_E_star": ratios["v_over_E_star"],
        "v_from_bridge": v_from_n,
        "v_bridge_error": v_from_n - ratios["v_over_E_cell"],
        "t_tr": ratios["t_tr"],
    }


def n_diagnostic(p_star: Decimal, alpha_u: Decimal, n_crc: Decimal) -> dict[str, Decimal]:
    log_n = ln(n_crc / PI)
    ratios = hierarchy_ratios(p_star, alpha_u)
    g_abs = g_tick_from_log(log_n)
    projected = projection_exponent(alpha_u, log_n)
    target_exponent = FOUR * p_star
    return {
        "N": n_crc,
        "log_N_over_pi": log_n,
        "g_tick_abs": g_abs,
        "projection_exponent": projected,
        "target_exponent": target_exponent,
        "projection_exponent_error": projected - target_exponent,
        "bridge_residual": bridge_residual(alpha_u, p_star, log_n),
        "v_from_N": exp(-(p_star / TWELVE) * log_n),
        "v_over_E_cell_source": ratios["v_over_E_cell"],
        "v_error": exp(-(p_star / TWELVE) * log_n) - ratios["v_over_E_cell"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Verify OPH issue #337 electroweak tick projection.")
    parser.add_argument("--P", default=str(DEFAULT_P), help="OPH pixel fixed point P_star")
    parser.add_argument("--alpha-u", default=str(DEFAULT_ALPHA_U), help="source alpha_U(P_star)")
    parser.add_argument("--alpha-u-lo", default=str(DEFAULT_ALPHA_U_LO), help="lower endpoint of alpha_U interval")
    parser.add_argument("--alpha-u-hi", default=str(DEFAULT_ALPHA_U_HI), help="upper endpoint of alpha_U interval")
    parser.add_argument("--rounded-n", default=str(DEFAULT_ROUNDED_N), help="rounded capacity display for diagnostic guard")
    parser.add_argument("--tolerance", default=str(DEFAULT_TOL), help="absolute tolerance for decimal residual checks")
    parser.add_argument("--check", action="store_true", help="exit nonzero unless the exact bridge certificate passes")
    parser.add_argument("--output", default=None, help="write JSON certificate to path")
    args = parser.parse_args()

    p_star = D(args.P)
    alpha_u = D(args.alpha_u)
    alpha_lo = D(args.alpha_u_lo)
    alpha_hi = D(args.alpha_u_hi)
    rounded_n = D(args.rounded_n)
    tol = D(args.tolerance)
    assert p_star is not None and alpha_u is not None and alpha_lo is not None and alpha_hi is not None
    assert rounded_n is not None and tol is not None

    exact = exact_branch(p_star, alpha_u)
    rounded = n_diagnostic(p_star, alpha_u, rounded_n)
    accepted = abs(exact["projection_exponent_error"]) <= tol and abs(exact["v_bridge_error"]) <= tol

    cert: dict[str, Any] = {
        "issue": 337,
        "artifact": "R_EW_tick_projection_certificate",
        "status": "closed_projection_map_with_exact_bridge_condition",
        "accepted": bool(accepted),
        "theorem": "electroweak tick-projection bridge for the local/global hierarchy resonance",
        "definitions": {
            "v_over_E_cell": "exp[-2*pi/(4*alpha_U(P))]",
            "g_tick_abs": "(N/pi)^(-1/48)",
            "Pi_EW": "-log(v/E_cell)/(-log|g_*'|) = 24*pi/(alpha_U(P)*log(N/pi))",
            "target_projection": "Pi_EW(P_star,N_CRC)=4*P_star",
            "bridge_residual": "B_EW(P,N)=alpha_U(P)*log(N/pi)-6*pi/P",
            "exact_bridge_capacity": "N_EW(P)=pi*exp[6*pi/(P*alpha_U(P))]",
        },
        "source_values": {
            "P_star": decstr(p_star),
            "alpha_U": decstr(alpha_u),
            "alpha_U_interval": [decstr(alpha_lo), decstr(alpha_hi)],
            "beta_EW": "4",
            "tick_denominator": "48",
            "v_over_E_cell": decstr(exact["v_over_E_cell"]),
            "v_over_E_star": decstr(exact["v_over_E_star"]),
            "t_tr": decstr(exact["t_tr"]),
        },
        "exact_bridge": {
            "N_EW": decstr(exact["N_EW"]),
            "log_N_EW_over_pi": decstr(exact["log_N_over_pi"]),
            "g_tick_abs": decstr(exact["g_tick_abs"]),
            "projection_exponent": decstr(exact["projection_exponent"]),
            "target_exponent_4P": decstr(exact["target_exponent"]),
            "projection_exponent_error": decstr(exact["projection_exponent_error"]),
            "bridge_residual": decstr(exact["bridge_residual"]),
            "v_from_bridge": decstr(exact["v_from_bridge"]),
            "v_bridge_error": decstr(exact["v_bridge_error"]),
        },
        "rounded_capacity_diagnostic": {
            "N_display": decstr(rounded["N"]),
            "status": "diagnostic_only_not_exact_bridge_certificate",
            "g_tick_abs": decstr(rounded["g_tick_abs"]),
            "projection_exponent": decstr(rounded["projection_exponent"]),
            "target_exponent_4P": decstr(rounded["target_exponent"]),
            "projection_exponent_error": decstr(rounded["projection_exponent_error"]),
            "bridge_residual": decstr(rounded["bridge_residual"]),
            "v_from_N": decstr(rounded["v_from_N"]),
            "v_over_E_cell_source": decstr(rounded["v_over_E_cell_source"]),
            "v_error": decstr(rounded["v_error"]),
        },
        "allowed_inputs": [
            "OPH local pixel fixed point P_star",
            "source D10 alpha_U(P_star) interval and forward transmutation law",
            "electroweak transmutation multiplicity beta_EW=N_c+1=4",
            "global repair-tick record |g_*'|=(N_CRC/pi)^(-1/48)",
            "joint (P,N_CRC) same-pair stability record",
        ],
        "forbidden_calibrations": [
            "measured weak scale v as an input",
            "measured W, Z, Higgs, or top mass as an input",
            "measured G or Planck area as an input",
            "measured Lambda as an input",
            "rounded N_CRC display as a high-precision bridge certificate",
        ],
        "claim_boundary": {
            "closed_here": "unique projection map, origin of 4P, exact bridge residual, and exact bridge target N_EW(P_star)",
            "source_certificate_required": "The issue-344 exact global capacity certificate supplies B_EW(P_star,N_CRC^EW)=0 for the bridge-refined capacity object.",
            "rounded_display": "3.31e122 is a capacity-scale label and fails the exact bridge diagnostic.",
        },
        "verifier_command": "python3 code/particles/hierarchy/verify_issue_337_electroweak_projection.py --check --output code/particles/hierarchy/certificates/R_EW_tick_projection_certificate.json",
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
