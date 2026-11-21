#!/usr/bin/env python3
"""Self-contained correction of the OPH no-EFE wide-binary observables.

The shipped handoff stores sqrt(nu) under a generic ``boost`` key and then
applies a second square root when reporting the Chae-style logarithmic
parameter.  This script keeps acceleration, velocity, and logarithmic
quantities distinct.  It deliberately does not implement a scalar MOND-EFE
comparator: AQUAL/QUMOND EFE predictions require a vector/tensor forward model.
"""
from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass

A_EFF = 1.1847373881e-10  # m s^-2, archived OPH effective scale
ARCHIVED_G_N = 8.895394264866785e-11
DEEP_G_N = 1.0e-11


@dataclass(frozen=True)
class Point:
    g_n_m_s2: float
    acceleration_boost_nu: float
    velocity_boost_sqrt_nu: float
    gamma_log10_velocity_boost: float
    acceleration_enhancement_fraction: float
    velocity_enhancement_fraction: float


def nu_oph(g_n: float, a_eff: float = A_EFF) -> float:
    if g_n <= 0 or a_eff <= 0:
        raise ValueError("g_n and a_eff must be positive")
    y = math.sqrt(g_n / a_eff)
    return 1.0 / (-math.expm1(-y))


def evaluate(g_n: float, a_eff: float = A_EFF) -> Point:
    nu = nu_oph(g_n, a_eff)
    v = math.sqrt(nu)
    return Point(
        g_n_m_s2=g_n,
        acceleration_boost_nu=nu,
        velocity_boost_sqrt_nu=v,
        gamma_log10_velocity_boost=math.log10(v),
        acceleration_enhancement_fraction=nu - 1.0,
        velocity_enhancement_fraction=v - 1.0,
    )


def payload() -> dict:
    archived = evaluate(ARCHIVED_G_N)
    deep = evaluate(DEEP_G_N)
    shipped_gamma = 0.5 * math.log10(archived.velocity_boost_sqrt_nu)
    return {
        "schema": "OPH gravity unit correction v2",
        "a_eff_m_s2": A_EFF,
        "points": {
            "archived_sample": asdict(archived),
            "deep_tail_1e-11": asdict(deep),
        },
        "shipped_nested_sqrt_error": {
            "stored_value_1_3135_is": "velocity boost sqrt(nu), not G_eff/G",
            "shipped_Gamma": shipped_gamma,
            "correct_Gamma": archived.gamma_log10_velocity_boost,
            "factor_in_log_parameter": archived.gamma_log10_velocity_boost / shipped_gamma,
        },
        "comparator_policy": {
            "newton": "acceleration boost = 1",
            "mond_efe": "must be generated with a published AQUAL/QUMOND forward model; scalar g_N+g_ext substitution is not accepted",
            "oph_no_efe": "conditional on the closed-cut activation premise",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    data = payload()
    if args.json:
        print(json.dumps(data, indent=2))
        return
    for name, row in data["points"].items():
        print(name)
        print(f"  g_N                      = {row['g_n_m_s2']:.12e} m s^-2")
        print(f"  acceleration boost nu    = {row['acceleration_boost_nu']:.12f}")
        print(f"  velocity boost sqrt(nu)  = {row['velocity_boost_sqrt_nu']:.12f}")
        print(f"  Gamma=log10(sqrt(nu))    = {row['gamma_log10_velocity_boost']:.12f}")
    print("\nDo not use nu(g_N+g_ext) as an AQUAL/QUMOND EFE likelihood.")


if __name__ == "__main__":
    main()
