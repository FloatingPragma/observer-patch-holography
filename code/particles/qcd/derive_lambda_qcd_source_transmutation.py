#!/usr/bin/env python3
"""Derive Lambda_QCD from the source strong coupling by dimensional transmutation.

The strict source-audit branch emits the strong coupling
``alpha_3(M_Z_run) = 0.11833586196478191`` through the same slope stack that
produces the electroweak chart.  Four-loop MS-bar running plus standard
decoupling then determines the QCD scale ``Lambda_MSbar^(nf)`` with no
hadronic input.  This is the perturbative half of every light-hadron mass:
``m_h/E_star = (m_h/Lambda) * (Lambda/E_star)``, where the first factor is
nonperturbative and lives outside this lane.

Quark threshold locations enter only logarithmically; they are declared
external scheme inputs and are swept over wide brackets so the output is an
interval.  The four-loop machinery is validated against the published world
values as a code certification block; that comparison never enters the source
rows.  No measured hadron mass is read.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = (
    ROOT / "particles" / "runs" / "qcd" / "lambda_qcd_source_transmutation.json"
)

# Source inputs (strict source-audit branch, shared with the D10/D11 lanes).
ALPHA_3_MZ = 0.11833586196478191
MU_Z_OVER_E = 7.501767385088419e-18
E_STAR_DISPLAY_GEV = 1.2208901289579269e19  # clock candidate; checksum status
MZ_RUN_GEV = MU_Z_OVER_E * E_STAR_DISPLAY_GEV

# Declared external threshold scheme inputs (log-suppressed sensitivity).
MB_MB_GEV = 4.18
MC_MC_GEV = 1.27
MB_BRACKET = (3.5, 5.0)
MC_BRACKET = (1.0, 1.6)

ZETA3 = 1.2020569031595942


def beta_coeffs(nf: int) -> tuple[float, float, float, float]:
    b0 = 11.0 - 2.0 * nf / 3.0
    b1 = 102.0 - 38.0 * nf / 3.0
    b2 = 2857.0 / 2.0 - 5033.0 * nf / 18.0 + 325.0 * nf**2 / 54.0
    b3 = (
        149753.0 / 6.0
        + 3564.0 * ZETA3
        - (1078361.0 / 162.0 + 6508.0 * ZETA3 / 27.0) * nf
        + (50065.0 / 162.0 + 6472.0 * ZETA3 / 81.0) * nf**2
        + 1093.0 * nf**3 / 729.0
    )
    return b0, b1, b2, b3


def run_alpha_s(alpha0: float, mu0: float, mu1: float, nf: int, n_steps: int = 4000) -> float:
    """Four-loop RK4 running of a = alpha_s/(4 pi) in ln mu."""

    b0, b1, b2, b3 = beta_coeffs(nf)

    def deriv(a: float) -> float:
        return -2.0 * a * a * (b0 + b1 * a + b2 * a * a + b3 * a**3)

    a = alpha0 / (4.0 * math.pi)
    t0, t1 = math.log(mu0), math.log(mu1)
    dt = (t1 - t0) / n_steps
    for _ in range(n_steps):
        k1 = deriv(a)
        k2 = deriv(a + 0.5 * dt * k1)
        k3 = deriv(a + 0.5 * dt * k2)
        k4 = deriv(a + dt * k3)
        a += dt / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    return 4.0 * math.pi * a


def decouple_down(alpha_nf: float, nl: int) -> float:
    """Match alpha_s^(nf) to alpha_s^(nl=nf-1) at mu = m_q(m_q)."""

    a = alpha_nf / math.pi
    c2 = -11.0 / 72.0
    c3 = -564731.0 / 124416.0 + 82043.0 * ZETA3 / 27648.0 + 2633.0 * nl / 31104.0
    return alpha_nf * (1.0 + c2 * a * a + c3 * a**3)


def lambda_msbar(alpha_s: float, mu: float, nf: int) -> float:
    """Exact four-loop MS-bar Lambda by the standard asymptotic convention.

    t(a) = 1/(b0 a) + (b1/b0^2) ln(b0 a)
           + int_0^a [ -1/(x^2 B(x)) + 1/(b0 x^2) - b1/(b0^2 x) ] dx,
    with B(x) = b0 + b1 x + b2 x^2 + b3 x^3 and t = ln(mu^2/Lambda^2).
    """

    b0, b1, b2, b3 = beta_coeffs(nf)
    a = alpha_s / (4.0 * math.pi)

    def integrand(x: float) -> float:
        big_b = b0 + b1 * x + b2 * x * x + b3 * x**3
        return -1.0 / (x * x * big_b) + 1.0 / (b0 * x * x) - b1 / (b0 * b0 * x)

    tail, _ = quad(integrand, 0.0, a, limit=200)
    t = 1.0 / (b0 * a) + (b1 / (b0 * b0)) * math.log(b0 * a) + tail
    return mu * math.exp(-0.5 * t)


def chain(alpha_mz: float, mz: float, mb: float, mc: float) -> dict[str, float]:
    """Run 5 -> 4 -> 3 flavors and extract Lambda per flavor number."""

    lam5 = lambda_msbar(alpha_mz, mz, 5)
    alpha_mb = run_alpha_s(alpha_mz, mz, mb, 5)
    alpha_mb_4 = decouple_down(alpha_mb, 4)
    lam4 = lambda_msbar(alpha_mb_4, mb, 4)
    alpha_mc = run_alpha_s(alpha_mb_4, mb, mc, 4)
    alpha_mc_3 = decouple_down(alpha_mc, 3)
    lam3 = lambda_msbar(alpha_mc_3, mc, 3)
    # Stability probe: extract Lambda^(3) again after running to 1.5 GeV.
    alpha_15 = run_alpha_s(alpha_mc_3, mc, 1.5, 3)
    lam3_probe = lambda_msbar(alpha_15, 1.5, 3)
    return {
        "lambda5_gev": lam5,
        "lambda4_gev": lam4,
        "lambda3_gev": lam3,
        "lambda3_probe_gev": lam3_probe,
        "alpha_s_mb_5f": alpha_mb,
        "alpha_s_mc_3f": alpha_mc_3,
    }


def machinery_validation() -> dict[str, Any]:
    """Certify the running/extraction code against published world values.

    This block consumes the published world coupling as a code benchmark
    only; nothing here enters the source rows.
    """

    world = chain(0.1179, 91.1876, MB_MB_GEV, MC_MC_GEV)
    checks = {
        "lambda5_near_210_mev": abs(world["lambda5_gev"] - 0.210) < 0.012,
        "lambda3_near_338_mev": abs(world["lambda3_gev"] - 0.338) < 0.020,
        "lambda3_extraction_stable": abs(
            world["lambda3_probe_gev"] / world["lambda3_gev"] - 1.0
        )
        < 0.02,
    }
    return {
        "role": "code_certification_only_never_a_source_row",
        "world_input_alpha_s_mz": 0.1179,
        "world_chain": world,
        "published_compare": {"lambda5_gev": 0.210, "lambda3_gev": 0.338},
        "checks": checks,
    }


def build_artifact() -> dict[str, Any]:
    central = chain(ALPHA_3_MZ, MZ_RUN_GEV, MB_MB_GEV, MC_MC_GEV)

    corners = []
    for mb in MB_BRACKET:
        for mc in MC_BRACKET:
            corners.append(chain(ALPHA_3_MZ, MZ_RUN_GEV, mb, mc))
    lam3_values = [c["lambda3_gev"] for c in corners] + [central["lambda3_gev"]]
    lam5_values = [c["lambda5_gev"] for c in corners] + [central["lambda5_gev"]]
    lam3_interval = (min(lam3_values), max(lam3_values))
    lam5_interval = (min(lam5_values), max(lam5_values))

    validation = machinery_validation()

    checks = {
        "machinery_certified": all(validation["checks"].values()),
        "lambda3_positive": lam3_interval[0] > 0,
        "extraction_scale_stability": abs(
            central["lambda3_probe_gev"] / central["lambda3_gev"] - 1.0
        )
        < 0.02,
        "threshold_brackets_swept": len(corners) == 4,
    }

    return {
        "artifact": "oph_lambda_qcd_source_transmutation",
        "schema_version": 1,
        "status": "perturbative_transmutation_from_source_coupling",
        "row_class": "conditional_on_P_and_declared_threshold_inputs",
        "promotion_allowed": False,
        "source_input": {
            "alpha_3_mz_run": ALPHA_3_MZ,
            "mz_run_gev_display": MZ_RUN_GEV,
            "provenance": (
                "strict source-audit branch slope stack; shared basis with "
                "wzh_residual_elimination_boundary.json"
            ),
        },
        "declared_external_inputs": {
            "mb_mb_gev": MB_MB_GEV,
            "mc_mc_gev": MC_MC_GEV,
            "brackets": {"mb": MB_BRACKET, "mc": MC_BRACKET},
            "ancestry_note": (
                "threshold locations are declared external quark scheme "
                "masses; sensitivity is logarithmic and covered by the "
                "bracket sweep"
            ),
        },
        "loop_order": "four_loop_beta_with_three_loop_decoupling",
        "central": central,
        "lambda3_interval_gev": lam3_interval,
        "lambda5_interval_gev": lam5_interval,
        "over_E_star": {
            "lambda3": central["lambda3_gev"] / E_STAR_DISPLAY_GEV,
            "lambda5": central["lambda5_gev"] / E_STAR_DISPLAY_GEV,
        },
        "machinery_validation": validation,
        "display_note": (
            "GeV values use the unclosed clock candidate; the clock audit "
            "classifies that decimal as a calibration checksum."
        ),
        "checks": checks,
        "checks_pass": all(bool(v) for v in checks.values()),
        "claim_boundary": (
            "Lambda is the perturbative transmutation scale of the source "
            "coupling. Hadron masses additionally require the nonperturbative "
            "ratio m_h/Lambda, which this lane does not supply."
        ),
    }


def build() -> dict[str, Any]:
    return build_artifact()


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
                "lambda3_gev": artifact["central"]["lambda3_gev"],
                "lambda3_interval_gev": artifact["lambda3_interval_gev"],
                "lambda5_gev": artifact["central"]["lambda5_gev"],
                "world_compare_lambda3": artifact["machinery_validation"][
                    "world_chain"
                ]["lambda3_gev"],
                "output": str(args.output),
            },
            indent=2,
        )
    )
    return 0 if artifact["checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
