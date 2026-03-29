#!/usr/bin/env python3
"""Emit quark gap pairs, centered even-log profiles, and leading readbacks from sigma_u and sigma_d."""

from __future__ import annotations

import argparse
import json
import math


RHO_ORD = 1.2942849363777058
G_SHARED = 0.9231656602589083

U_PROFILE_RAY = (
    -(2 * RHO_ORD + 1) / (3 * (1 + RHO_ORD)),
    (RHO_ORD - 1) / (3 * (1 + RHO_ORD)),
    (RHO_ORD + 2) / (3 * (1 + RHO_ORD)),
)
D_PROFILE_RAY = (
    -(RHO_ORD + 2) / (3 * (1 + RHO_ORD)),
    (1 - RHO_ORD) / (3 * (1 + RHO_ORD)),
    (2 * RHO_ORD + 1) / (3 * (1 + RHO_ORD)),
)


def gap_pair(sig: float, sector: str) -> dict[str, float]:
    if sector == "u":
        return {
            "gamma21_log_per_side": RHO_ORD * sig / (1 + RHO_ORD),
            "gamma32_log_per_side": sig / (1 + RHO_ORD),
        }
    if sector == "d":
        return {
            "gamma21_log_per_side": sig / (1 + RHO_ORD),
            "gamma32_log_per_side": RHO_ORD * sig / (1 + RHO_ORD),
        }
    raise ValueError(f"unknown sector: {sector}")


def centered_profile(sig: float, sector: str) -> list[float]:
    ray = U_PROFILE_RAY if sector == "u" else D_PROFILE_RAY
    return [sig * entry for entry in ray]


def leading_singular_values(sig: float, sector: str) -> list[float]:
    return [G_SHARED * math.exp(2.0 * entry) for entry in centered_profile(sig, sector)]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sigma-u", type=float, required=True, help="sigma_u total log spread per side")
    parser.add_argument("--sigma-d", type=float, required=True, help="sigma_d total log spread per side")
    args = parser.parse_args()

    payload = {
        "artifact": "oph_quark_even_profile_rays",
        "rho_ord": RHO_ORD,
        "shared_norm_value": G_SHARED,
        "sigma_u_total_log_per_side": args.sigma_u,
        "sigma_d_total_log_per_side": args.sigma_d,
        "u_profile_ray": U_PROFILE_RAY,
        "d_profile_ray": D_PROFILE_RAY,
        "gap_pair_u": gap_pair(args.sigma_u, "u"),
        "gap_pair_d": gap_pair(args.sigma_d, "d"),
        "E_u_log": centered_profile(args.sigma_u, "u"),
        "E_d_log": centered_profile(args.sigma_d, "d"),
        "leading_singular_values_u": leading_singular_values(args.sigma_u, "u"),
        "leading_singular_values_d": leading_singular_values(args.sigma_d, "d"),
        "sigma_seed_ud_candidate": (args.sigma_u + args.sigma_d) / 2.0,
        "eta_ud_candidate": (args.sigma_u - args.sigma_d) / 2.0,
    }
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
