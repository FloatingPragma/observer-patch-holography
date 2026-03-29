#!/usr/bin/env python3
"""Export the current sector-even quark excitation-evaluator candidate."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_GENERATOR = ROOT / "particles" / "runs" / "flavor" / "generation_bundle_branch_generator.json"
DEFAULT_DESCENT = ROOT / "particles" / "runs" / "flavor" / "quark_sector_descent.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "family_excitation_evaluator.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _ctr(values: np.ndarray) -> list[float]:
    centered = values - np.mean(values)
    return [float(x) for x in centered.tolist()]


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the family excitation-evaluator candidate artifact.")
    parser.add_argument("--generator", default=str(DEFAULT_GENERATOR))
    parser.add_argument("--descent", default=str(DEFAULT_DESCENT))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    generator = json.loads(Path(args.generator).read_text(encoding="utf-8"))
    descent = json.loads(Path(args.descent).read_text(encoding="utf-8"))
    centered = generator["centered_compressed_branch_generator"]
    matrix = np.asarray(centered["real"], dtype=float) + 1j * np.asarray(centered["imag"], dtype=float)
    eigenvalues = np.linalg.eigvalsh(matrix)
    lam1, lam2, lam3 = [float(x) for x in eigenvalues.tolist()]
    gap21 = lam2 - lam1
    gap32 = lam3 - lam2
    x = np.asarray(
        [
            2.0 * (lam1 - lam1) / (lam3 - lam1) - 1.0,
            2.0 * (lam2 - lam1) / (lam3 - lam1) - 1.0,
            2.0 * (lam3 - lam1) / (lam3 - lam1) - 1.0,
        ],
        dtype=float,
    )
    x2 = float(x[1])
    rho_ord = float(3.0 * gap32 / (2.0 * gap32 + gap21))
    rho_ord_reciprocal = float(1.0 / rho_ord)
    diagnostic_sigma_u = 5.5905
    diagnostic_sigma_d = 3.3049
    diagnostic_gap_pair_u = {
        "gamma21_log_per_side": float(rho_ord * diagnostic_sigma_u / (1.0 + rho_ord)),
        "gamma32_log_per_side": float(diagnostic_sigma_u / (1.0 + rho_ord)),
    }
    diagnostic_gap_pair_d = {
        "gamma21_log_per_side": float(diagnostic_sigma_d / (1.0 + rho_ord)),
        "gamma32_log_per_side": float(rho_ord * diagnostic_sigma_d / (1.0 + rho_ord)),
    }
    diagnostic_coeffs_u = {
        "linear": (diagnostic_gap_pair_u["gamma21_log_per_side"] + diagnostic_gap_pair_u["gamma32_log_per_side"]) / 2.0,
        "quadratic": ((1.0 + x2) * diagnostic_gap_pair_u["gamma32_log_per_side"] - (1.0 - x2) * diagnostic_gap_pair_u["gamma21_log_per_side"]) / (2.0 * (1.0 - x2**2)),
    }
    diagnostic_coeffs_d = {
        "linear": (diagnostic_gap_pair_d["gamma21_log_per_side"] + diagnostic_gap_pair_d["gamma32_log_per_side"]) / 2.0,
        "quadratic": ((1.0 + x2) * diagnostic_gap_pair_d["gamma32_log_per_side"] - (1.0 - x2) * diagnostic_gap_pair_d["gamma21_log_per_side"]) / (2.0 * (1.0 - x2**2)),
    }
    e_u = _ctr(diagnostic_coeffs_u["linear"] * x + diagnostic_coeffs_u["quadratic"] * x**2)
    e_d = _ctr(diagnostic_coeffs_d["linear"] * x + diagnostic_coeffs_d["quadratic"] * x**2)

    artifact = {
        "artifact": "oph_family_excitation_evaluator",
        "generated_utc": _timestamp(),
        "proof_status": "exact_reduced_family_on_simple_ordered_three_point_spectrum",
        "theorem_candidate": "every self-adjoint trace-zero sector-even evaluator diagonal in the ordered projector algebra on a simple three-point spectrum is quadratic in X_ord",
        "input_kind": "ordered_branch_generator_spectral_package",
        "log_insertion_convention": "add_E_q_log_to_both_left_and_right_diagonal_logs",
        "family_coordinate_kind": "affine_normalized_ordered_branch_coordinate",
        "parameterization_kind": "sector_total_spread_primary",
        "smallest_constructive_missing_object": "oph_family_excitation_spread_map",
        "next_live_mover": "ordered_spectrum_to_two_spread_emitter",
        "quadratic_even_family_closed": True,
        "live_blocker": "quark_spread_map_missing",
        "ordered_family_eigenvalues": [lam1, lam2, lam3],
        "ordered_spectral_gaps": [gap21, gap32],
        "family_coordinate_x": [float(v) for v in x.tolist()],
        "centering_mode": "trace_zero",
        "quadratic_family_formula": "E_q_log = ctr(a_q * X_ord + b_q * X_ord^2)",
        "rho_ord": rho_ord,
        "rho_ord_reciprocal": rho_ord_reciprocal,
        "ratio_law_status": "closed_from_ordered_branch_generator_spectral_package",
        "reduced_primary_object": "oph_family_excitation_spread_map",
        "spread_map": {
            "sigma_u_total_log_per_side": None,
            "sigma_d_total_log_per_side": None,
            "gap_pair_u": {
                "gamma21_log_per_side": "rho_ord * sigma_u / (1 + rho_ord)",
                "gamma32_log_per_side": "sigma_u / (1 + rho_ord)",
            },
            "gap_pair_d": {
                "gamma21_log_per_side": "sigma_d / (1 + rho_ord)",
                "gamma32_log_per_side": "rho_ord * sigma_d / (1 + rho_ord)",
            },
        },
        "gap_pair_u": {"gamma21_log_per_side": None, "gamma32_log_per_side": None},
        "gap_pair_d": {"gamma21_log_per_side": None, "gamma32_log_per_side": None},
        "gap_pair_reduction": {
            "gamma21_q": "E_q,2 - E_q,1",
            "gamma32_q": "E_q,3 - E_q,2",
            "eigenvalue_recovery": {
                "e1": "-(2*gamma21_q + gamma32_q)/3",
                "e2": "(gamma21_q - gamma32_q)/3",
                "e3": "(gamma21_q + 2*gamma32_q)/3",
            },
            "coefficient_recovery": {
                "a_q": "(gamma21_q + gamma32_q)/2",
                "b_q": "((1 + x2)*gamma32_q - (1 - x2)*gamma21_q) / (2*(1 - x2^2))",
            },
        },
        "quadratic_polynomial_coeffs_u": {"linear": None, "quadratic": None},
        "quadratic_polynomial_coeffs_d": {"linear": None, "quadratic": None},
        "diagnostic_witness_sigma_u": diagnostic_sigma_u,
        "diagnostic_witness_sigma_d": diagnostic_sigma_d,
        "diagnostic_witness_gap_pair_u": diagnostic_gap_pair_u,
        "diagnostic_witness_gap_pair_d": diagnostic_gap_pair_d,
        "diagnostic_witness_coeffs_u": diagnostic_coeffs_u,
        "diagnostic_witness_coeffs_d": diagnostic_coeffs_d,
        "diagnostic_witness_E_u_log": e_u,
        "diagnostic_witness_E_d_log": e_d,
        "diagnostic_forward_consequence": {
            "singular_values_u": [0.00270846, 1.48784982, 194.714124],
            "singular_values_d": [0.03884097, 0.69477844, 28.971755],
        },
        "diagnostic_required_gap_ratio_u": rho_ord,
        "diagnostic_required_gap_ratio_d": rho_ord_reciprocal,
        "shared_norm_origin": descent.get("shared_norm_origin"),
        "shared_norm_value": descent.get("g_ch"),
        "promotion_blocker_cleared": "quark_even_excitation_evaluator_missing",
        "gap_pair_value_count_required": 2,
        "notes": [
            "The smallest reduced live family is not a free three-weight projector ansatz but the exact quadratic family forced by a simple ordered three-point spectrum plus trace-zero centering.",
            "The ordered package already fixes the canonical ratio law rho_ord, so the four-gap map reduces constructively to a two-spread emitter sigma_u and sigma_d.",
            "Under the current factorized forward builder, actual quark mass movement now depends only on emitting those two sector-total spreads; the gap pairs, quadratic coefficients, and centered E_q_log vectors are then deterministic readback.",
            "The diagnostic witness gap pairs and coefficients are reference-facing only; they show the hierarchy the missing spread map would generate without being promoted to theorem status.",
        ],
    }

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
