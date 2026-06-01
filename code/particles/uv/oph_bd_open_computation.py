#!/usr/bin/env python3
"""Low-energy target equations for the OPH-selected BD n=1 branch.

This script reproduces the numerical threshold targets reported in the
2026-06-01 OPH/BD open-computation handoff. It does not compute the raw
Bouchard-Donagi algebraic geometry, bundle cohomology, or string thresholds.
It only turns the declared OPH electroweak target vector into the corresponding
MSSM tree-level, stop-threshold, and one-loop unification proxy targets.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class OPHElectroweakTargets:
    alpha2_mz: float = 0.03377843630219015
    alphaY_mz: float = 0.010131601067241624
    alpha3_reference_mz: float = 0.1184
    v_gev: float = 246.76711732749683
    mH_gev: float = 125.1995304097179
    mt_gev: float = 172.35235532883115
    # Inferred from the unification proxy values in the 2026-06-01 handoff.
    mZ_reference_gev: float = 91.187978085123


def top_yukawa(targets: OPHElectroweakTargets) -> float:
    return math.sqrt(2.0) * targets.mt_gev / targets.v_gev


def higgs_quartic(targets: OPHElectroweakTargets) -> float:
    return targets.mH_gev**2 / (2.0 * targets.v_gev**2)


def mssm_tree_quartic_max(targets: OPHElectroweakTargets) -> float:
    # lambda = (g_2^2 + g_Y^2) / 8, with g_i^2 = 4*pi*alpha_i.
    return (math.pi / 2.0) * (targets.alpha2_mz + targets.alphaY_mz)


def stop_threshold_proxy(
    targets: OPHElectroweakTargets, tan_beta: float, x_t_over_m_s: float
) -> float:
    """Return the required M_S in GeV in the one-loop top/stop proxy."""

    lambda_tree_max = mssm_tree_quartic_max(targets)
    m_tree_max = math.sqrt(2.0 * lambda_tree_max * targets.v_gev**2)
    cos_2beta = (1.0 - tan_beta**2) / (1.0 + tan_beta**2)
    tree_m2 = (m_tree_max * cos_2beta) ** 2
    delta_m2 = targets.mH_gev**2 - tree_m2
    prefactor = 3.0 * targets.mt_gev**4 / (2.0 * math.pi**2 * targets.v_gev**2)
    bracket = delta_m2 / prefactor
    mixing = x_t_over_m_s**2 * (1.0 - x_t_over_m_s**2 / 12.0)
    return targets.mt_gev * math.exp((bracket - mixing) / 2.0)


def gauge_unification_proxy(
    targets: OPHElectroweakTargets, m_susy_gev: float = 1000.0
) -> dict[str, float]:
    """One-loop SM-to-MSSM unification proxy with SU(5)-normalized alpha_1."""

    alpha1_mz = (5.0 / 3.0) * targets.alphaY_mz
    inv_mz = [1.0 / alpha1_mz, 1.0 / targets.alpha2_mz, 1.0 / targets.alpha3_reference_mz]
    b_sm = [41.0 / 10.0, -19.0 / 6.0, -7.0]
    b_mssm = [33.0 / 5.0, 1.0, -3.0]
    log_susy_over_mz = math.log(m_susy_gev / targets.mZ_reference_gev)

    inv_susy = [
        inv_mz[i] - b_sm[i] / (2.0 * math.pi) * log_susy_over_mz
        for i in range(3)
    ]
    log_mu_over_susy = (
        (inv_susy[0] - inv_susy[1])
        * 2.0
        * math.pi
        / (b_mssm[0] - b_mssm[1])
    )
    m_unification_gev = m_susy_gev * math.exp(log_mu_over_susy)
    alpha_u_inv = inv_susy[0] - b_mssm[0] / (2.0 * math.pi) * log_mu_over_susy

    inv_alpha3_susy_pred = alpha_u_inv + b_mssm[2] / (2.0 * math.pi) * log_mu_over_susy
    inv_alpha3_mz_pred = inv_alpha3_susy_pred + b_sm[2] / (2.0 * math.pi) * log_susy_over_mz
    alpha3_pred_mz = 1.0 / inv_alpha3_mz_pred
    delta_inv_alpha3_needed = (1.0 / targets.alpha3_reference_mz) - inv_alpha3_mz_pred

    return {
        "M_SUSY_GeV": m_susy_gev,
        "log10_MU_GeV": math.log10(m_unification_gev),
        "M_U_GeV": m_unification_gev,
        "alpha_U_inverse": alpha_u_inv,
        "alpha3_pred_mZ": alpha3_pred_mz,
        "delta_alpha3_inverse_needed": delta_inv_alpha3_needed,
    }


def operator_audit() -> dict[str, Any]:
    hypercharge = {
        "Q": 1.0 / 6.0,
        "Hu": 1.0 / 2.0,
        "Hd": -1.0 / 2.0,
        "uc": -2.0 / 3.0,
        "dc": 1.0 / 3.0,
        "L": -1.0 / 2.0,
        "ec": 1.0,
    }
    operators = {
        "Q Hu uc": ["Q", "Hu", "uc"],
        "Q Hd dc": ["Q", "Hd", "dc"],
        "L Hd ec": ["L", "Hd", "ec"],
        "L L ec": ["L", "L", "ec"],
        "L Q dc": ["L", "Q", "dc"],
        "uc dc dc": ["uc", "dc", "dc"],
        "Q Q Q L": ["Q", "Q", "Q", "L"],
        "uc uc dc ec": ["uc", "uc", "dc", "ec"],
    }
    return {
        name: {
            "hypercharge_sum": sum(hypercharge[field] for field in fields),
            "hypercharge_neutral": abs(sum(hypercharge[field] for field in fields)) < 1e-12,
        }
        for name, fields in operators.items()
    }


def compute_report() -> dict[str, Any]:
    targets = OPHElectroweakTargets()
    lambda_h = higgs_quartic(targets)
    lambda_tree_max = mssm_tree_quartic_max(targets)
    stop_rows = []
    for tan_beta in (2.0, 5.0, 10.0, 50.0):
        for x_t in (0.0, math.sqrt(6.0)):
            stop_rows.append(
                {
                    "tan_beta": tan_beta,
                    "X_t_over_M_S": x_t,
                    "required_M_S_GeV": stop_threshold_proxy(targets, tan_beta, x_t),
                }
            )

    return {
        "status": "target_equation_reduction_not_full_BD_geometry_solve",
        "selected_representative": "BD_n=1_OPH_Bouchard-Donagi_E8xE8_heterotic_SU5_one_Higgs_stratum",
        "inputs": asdict(targets),
        "low_energy_targets": {
            "y_t_OPH": top_yukawa(targets),
            "lambda_H_OPH": lambda_h,
            "lambda_MSSM_tree_max": lambda_tree_max,
            "Delta_lambda_min": lambda_h - lambda_tree_max,
            "m_h_tree_max_GeV": math.sqrt(2.0 * lambda_tree_max * targets.v_gev**2),
        },
        "stop_threshold_proxy": stop_rows,
        "gauge_unification_proxy": gauge_unification_proxy(targets),
        "operator_audit": operator_audit(),
        "open_gates": [
            "raw_BD_geometry_bundle_cohomology",
            "Yukawa_superpotential_cup_products_and_instantons",
            "heavy_spectrum_and_string_GUT_thresholds",
            "hidden_sector_SUSY_breaking_soft_terms",
            "full_moduli_lock_map_F_BD_n1",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, help="Optional JSON output path.")
    args = parser.parse_args()
    report = compute_report()
    encoded = json.dumps(report, indent=2, sort_keys=True)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(encoded + "\n", encoding="utf-8")
    else:
        print(encoded)


if __name__ == "__main__":
    main()
