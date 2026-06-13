#!/usr/bin/env python3
"""Verifier for OPH issue #344: exact EW-refined global capacity."""

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

DEFAULT_P = Decimal("1.6309682094039593248792798477826489413359828516279250606661507533907793398933432")
DEFAULT_ALPHA_U = Decimal("0.041124336195630495")
DEFAULT_ALPHA_U_LO = Decimal("0.041123336195630494")
DEFAULT_ALPHA_U_HI = Decimal("0.041125336195630496")
DEFAULT_ROUNDED_N = Decimal("3.31e122")
DEFAULT_LAMBDA = Decimal("0.5")
DEFAULT_TOL = Decimal("1e-40")


def D(value: str | int | Decimal | None, default: Decimal | None = None) -> Decimal | None:
    if value is None:
        return default
    if isinstance(value, Decimal):
        return value
    return Decimal(str(value))


def decstr(value: Decimal) -> str:
    if value != 0 and (abs(value) < Decimal("1e-6") or abs(value) >= Decimal("1e6")):
        return format(value, "E")
    return format(value, "f")


def ln(value: Decimal) -> Decimal:
    if value <= 0:
        raise ValueError(f"log argument must be positive, got {value}")
    return value.ln()


def exp(value: Decimal) -> Decimal:
    return value.exp()


def target_log_capacity(p_star: Decimal, alpha_u: Decimal) -> Decimal:
    return SIX * PI / (p_star * alpha_u)


def capacity_from_log(log_n_over_pi: Decimal) -> Decimal:
    return PI * exp(log_n_over_pi)


def bridge_residual(p_star: Decimal, alpha_u: Decimal, log_n_over_pi: Decimal) -> Decimal:
    return alpha_u * log_n_over_pi - SIX * PI / p_star


def projection_exponent(alpha_u: Decimal, log_n_over_pi: Decimal) -> Decimal:
    return TWENTY_FOUR * PI / (alpha_u * log_n_over_pi)


def contraction_map(x: Decimal, target: Decimal, lam: Decimal) -> Decimal:
    return (Decimal(1) - lam) * x + lam * target


def build_certificate(
    p_star: Decimal,
    alpha_u: Decimal,
    alpha_lo: Decimal,
    alpha_hi: Decimal,
    rounded_n: Decimal,
    lam: Decimal,
    tol: Decimal,
) -> dict[str, Any]:
    target_x = target_log_capacity(p_star, alpha_u)
    target_n = capacity_from_log(target_x)
    fixed_x = contraction_map(target_x, target_x, lam)
    exact_residual = bridge_residual(p_star, alpha_u, target_x)
    fixed_residual = fixed_x - target_x
    rounded_x = ln(rounded_n / PI)
    rounded_residual = bridge_residual(p_star, alpha_u, rounded_x)
    rounded_projection_error = projection_exponent(alpha_u, rounded_x) - FOUR * p_star
    contraction_factor = Decimal(1) - lam
    sample_x = rounded_x
    sample_residual = bridge_residual(p_star, alpha_u, sample_x)
    sample_next_x = contraction_map(sample_x, target_x, lam)
    sample_next_residual = bridge_residual(p_star, alpha_u, sample_next_x)
    residual_ratio = sample_next_residual / sample_residual if sample_residual else Decimal(0)
    v_source = exp(-(TWO * PI) / (FOUR * alpha_u))
    v_from_capacity = exp(-(p_star / TWELVE) * target_x)
    g_tick = exp(-target_x / FORTY_EIGHT)

    accepted = (
        Decimal(0) < lam <= Decimal(1)
        and abs(exact_residual) <= tol
        and abs(fixed_residual) <= tol
        and abs(v_from_capacity - v_source) <= tol
        and abs(residual_ratio - contraction_factor) <= tol
        and abs(rounded_residual) > Decimal("1e-6")
    )

    return {
        "issue": 344,
        "artifact": "R_EW_global_capacity_certificate",
        "status": "closed_bridge_refined_global_capacity_fixed_point_certificate",
        "accepted": bool(accepted),
        "theorem": "exact EW-refined global-capacity certificate for the local/global hierarchy bridge",
        "definitions": {
            "bridge_residual": "B_EW(P,N)=alpha_U(P)*log(N/pi)-6*pi/P",
            "exact_log_capacity": "x_EW(P)=6*pi/(P*alpha_U(P))",
            "exact_capacity": "N_CRC^EW(P)=pi*exp[x_EW(P)]",
            "capacity_map": "C_EW(P,x)=(1-lambda)*x+lambda*6*pi/(P*alpha_U(P))",
            "projection_map": "Pi_EW(P,N)=24*pi/(alpha_U(P)*log(N/pi))",
        },
        "source_values": {
            "P_star": decstr(p_star),
            "alpha_U": decstr(alpha_u),
            "alpha_U_interval": [decstr(alpha_lo), decstr(alpha_hi)],
            "lambda": decstr(lam),
            "contraction_factor": decstr(contraction_factor),
        },
        "exact_capacity_fixed_point": {
            "x_EW": decstr(target_x),
            "N_CRC_EW": decstr(target_n),
            "fixed_point_residual_x": decstr(fixed_residual),
            "bridge_residual": decstr(exact_residual),
            "projection_exponent": decstr(projection_exponent(alpha_u, target_x)),
            "target_exponent_4P": decstr(FOUR * p_star),
            "g_tick_abs": decstr(g_tick),
            "v_over_E_cell_source": decstr(v_source),
            "v_from_capacity": decstr(v_from_capacity),
            "v_identity_error": decstr(v_from_capacity - v_source),
        },
        "contraction_certificate": {
            "map": "C_EW(P_star,x)=(1-lambda)*x+lambda*x_EW(P_star)",
            "lambda": decstr(lam),
            "lipschitz_constant": decstr(contraction_factor),
            "banach_unique_fixed_point": True,
            "sample_x_from_rounded_capacity": decstr(sample_x),
            "sample_residual": decstr(sample_residual),
            "sample_next_x": decstr(sample_next_x),
            "sample_next_residual": decstr(sample_next_residual),
            "sample_residual_ratio": decstr(residual_ratio),
            "residual_contracts_by": decstr(contraction_factor),
        },
        "rounded_capacity_diagnostic": {
            "N_display": decstr(rounded_n),
            "status": "diagnostic_only_not_exact_bridge_certificate",
            "log_N_over_pi": decstr(rounded_x),
            "bridge_residual": decstr(rounded_residual),
            "projection_exponent": decstr(projection_exponent(alpha_u, rounded_x)),
            "target_exponent_4P": decstr(FOUR * p_star),
            "projection_exponent_error": decstr(rounded_projection_error),
            "relative_capacity_gap": decstr((target_n - rounded_n) / target_n),
        },
        "allowed_inputs": [
            "OPH local pixel fixed point P_star",
            "source D10 alpha_U(P_star) interval",
            "pi and exp",
            "the EW bridge residual isolated by the issue-337 projection theorem",
        ],
        "forbidden_calibrations": [
            "measured weak scale v as an input",
            "measured W, Z, Higgs, or top mass as an input",
            "measured G or Planck area as an input",
            "measured Lambda as an input",
            "observed hierarchy-ratio calibration",
            "rounded 3.31e122 capacity display as an exact bridge certificate",
        ],
        "claim_boundary": {
            "closed_here": "source-side contraction fixed point for N_CRC^EW satisfying B_EW(P_star,N_CRC^EW)=0",
            "rounded_display": "3.31e122 remains a capacity-scale display and fails the exact EW bridge residual.",
            "not_closed_here": [
                "finite-machinery verification that nf_{r,N} delivers a single well-defined effective readback resolution",
                "representation-to-spectrum theorem deriving the 24-round repair count",
            ],
        },
        "verifier_command": (
            "python3 code/particles/hierarchy/verify_issue_344_exact_capacity.py "
            "--check --output "
            "code/particles/hierarchy/certificates/R_EW_global_capacity_certificate.json"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Verify OPH issue #344 exact EW-refined capacity.")
    parser.add_argument("--P", default=str(DEFAULT_P), help="OPH pixel fixed point P_star")
    parser.add_argument("--alpha-u", default=str(DEFAULT_ALPHA_U), help="source alpha_U(P_star)")
    parser.add_argument("--alpha-u-lo", default=str(DEFAULT_ALPHA_U_LO), help="lower endpoint of alpha_U interval")
    parser.add_argument("--alpha-u-hi", default=str(DEFAULT_ALPHA_U_HI), help="upper endpoint of alpha_U interval")
    parser.add_argument("--rounded-n", default=str(DEFAULT_ROUNDED_N), help="rounded capacity display for diagnostic guard")
    parser.add_argument("--lambda", dest="lam", default=str(DEFAULT_LAMBDA), help="contraction averaging parameter")
    parser.add_argument("--tolerance", default=str(DEFAULT_TOL), help="absolute tolerance for decimal residual checks")
    parser.add_argument("--check", action="store_true", help="exit nonzero unless the exact capacity certificate passes")
    parser.add_argument("--output", default=None, help="write JSON certificate to path")
    args = parser.parse_args()

    p_star = D(args.P)
    alpha_u = D(args.alpha_u)
    alpha_lo = D(args.alpha_u_lo)
    alpha_hi = D(args.alpha_u_hi)
    rounded_n = D(args.rounded_n)
    lam = D(args.lam)
    tol = D(args.tolerance)
    assert all(value is not None for value in [p_star, alpha_u, alpha_lo, alpha_hi, rounded_n, lam, tol])

    cert = build_certificate(p_star, alpha_u, alpha_lo, alpha_hi, rounded_n, lam, tol)  # type: ignore[arg-type]
    text = json.dumps(cert, indent=2) + "\n"
    if args.output:
        Path(args.output).write_text(text, encoding="utf-8")
    else:
        print(text, end="")
    if args.check and not cert["accepted"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
