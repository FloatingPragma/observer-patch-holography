#!/usr/bin/env python3
"""Dimensional guard for the proposed Jacobi-clock hysteresis discriminator."""
from __future__ import annotations

import json
import math

G = 6.67430e-11
M_SUN = 1.98847e30
KPC = 3.0856775814913673e19


def jacobi_rate(mass_kg: float, distance_m: float) -> float:
    return math.sqrt(G * mass_kg / distance_m**3)


def payload() -> dict:
    mass = 1.0e12 * M_SUN
    d1 = 50.0 * KPC
    d2 = 100.0 * KPC
    g1 = jacobi_rate(mass, d1)
    g2 = jacobi_rate(mass, d2)
    return {
        "schema": "Jacobi scaling dimensional guard v1",
        "rate_definition": "Gamma_J=sqrt(GM/d^3)",
        "rate_scaling": "d^(-3/2)",
        "squared_rate_or_tidal_tensor_scaling": "d^(-3)",
        "dimensionless_population_coordinate": "eta=(M_host/M_sat)*(r_sat/d)^3, equivalently a Jacobi-radius ratio",
        "worked_ratio_50_to_100_kpc": {
            "Gamma_50_over_Gamma_100": g1 / g2,
            "expected": (d2 / d1) ** 1.5,
            "Gamma2_50_over_Gamma2_100": (g1 / g2) ** 2,
            "expected_squared": (d2 / d1) ** 3,
        },
        "correction": "A claim that the relaxation rate itself scales as d^-3 is dimensionally wrong. Only Gamma_J^2 or the tidal tensor does.",
    }


if __name__ == "__main__":
    print(json.dumps(payload(), indent=2))
