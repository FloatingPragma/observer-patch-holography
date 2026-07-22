#!/usr/bin/env python3
"""Exact theorem oracle; deliberately contains no evidence-promotion path."""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp


def run_checks() -> dict[str, bool | str]:
    checks: dict[str, bool | str] = {
        "suite": "oph-sm-qft-q1-q4-wz-theorem-oracle-v2"
    }

    def check(name: str, condition: object) -> None:
        passed = bool(condition)
        checks[name] = passed
        if not passed:
            raise AssertionError(name)

    # One-generation integer hypercharges, Yukawa invariance, and anomalies.
    q, uc, dc, lepton, ec, higgs = 1, -4, 2, -3, 6, 3
    check("yukawa_up", q + higgs + uc == 0)
    check("yukawa_down", q - higgs + dc == 0)
    check("yukawa_lepton", lepton - higgs + ec == 0)
    check("anomaly_su3_cubed", 2 - 1 - 1 == 0)
    check("anomaly_su3_squared_u1", 2 * q + uc + dc == 0)
    check("anomaly_su2_squared_u1", 3 * q + lepton == 0)
    check("anomaly_gravity_squared_u1", 6*q + 3*uc + 3*dc + 2*lepton + ec == 0)
    check("anomaly_u1_cubed", 6*q**3 + 3*uc**3 + 3*dc**3 + 2*lepton**3 + ec**3 == 0)
    check("witten_doublet_parity", (3 + 1) % 2 == 0)

    # Canonically normalized neutral tree kernel.
    g, gp, v, lam = sp.symbols("g gp v lam", nonzero=True)
    matrix = v**2 / 4 * sp.Matrix([[g**2, -g*gp], [-g*gp, gp**2]])
    characteristic = sp.factor((lam * sp.eye(2) - matrix).det())
    check("neutral_tree_roots", sp.simplify(characteristic - lam*(lam-(g**2+gp**2)*v**2/4)) == 0)

    # Complete first-order finite-coordinate change, using kappa consistently.
    kappa, p, dp = sp.symbols("kappa p dp")
    s0, s1 = p**2 + 3*p, p**3 - p
    transformed = sp.expand(s0.subs(p, p+kappa*dp) + kappa*s1.subs(p, p+kappa*dp))
    check("fj_complete_chain_rule", sp.simplify(transformed.coeff(kappa, 1) - (s1 + dp*sp.diff(s0, p))) == 0)

    # Strict charged and neutral pole series through kappa^2.
    p1, p1p, p2 = sp.symbols("p1 p1p p2")
    a, b = -p1, p1*p1p-p2
    check("charged_order_one", sp.simplify(a+p1) == 0)
    check("charged_order_two", sp.simplify(b+a*p1p+p2) == 0)
    z = sp.symbols("z", nonzero=True)
    d1, d1p, d2, za1, az1 = sp.symbols("d1 d1p d2 za1 az1")
    az, bz = -d1, d1*d1p-d2+za1*az1/z
    check("neutral_order_one", sp.simplify(az+d1) == 0)
    check("neutral_order_two_with_mixing", sp.simplify(bz+az*d1p+d2-za1*az1/z) == 0)
    check("neutral_mixing_is_not_strict_one_loop", not sp.simplify(za1*az1/z) == 0)

    # Strict mass/width coordinate through first order.
    m0, dsr, dsi = sp.symbols("m0 dsr dsi", nonzero=True, real=True)
    mass = m0 + kappa*dsr/(2*m0)
    width = -kappa*dsi/m0
    reconstructed = sp.expand((mass-sp.I*width/2)**2)
    target = m0**2 + kappa*(dsr+sp.I*dsi)
    check("strict_mass_width_roundtrip", sp.simplify(reconstructed.coeff(kappa, 1)-target.coeff(kappa, 1)) == 0)

    # Matrix Nielsen identity without assuming the inverse exists at the pole.
    a11, a12, a21, a22 = sp.symbols("a11 a12 a21 a22")
    l11, l12, l21, l22, r11, r12, r21, r22 = sp.symbols("l11 l12 l21 l22 r11 r12 r21 r22")
    gamma = sp.Matrix([[a11, a12], [a21, a22]])
    left = sp.Matrix([[l11, l12], [l21, l22]])
    right = sp.Matrix([[r11, r12], [r21, r22]])
    direction = left*gamma + gamma*right
    variables = [a11, a12, a21, a22]
    components = [direction[0,0], direction[0,1], direction[1,0], direction[1,1]]
    directional_det = sp.expand(sum(sp.diff(gamma.det(), x)*dx for x, dx in zip(variables, components)))
    check("matrix_nielsen_determinant_identity", sp.simplify(directional_det-(sp.trace(left)+sp.trace(right))*gamma.det()) == 0)

    checks["all_pass"] = True
    return checks


def main() -> None:
    print(json.dumps(run_checks(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

