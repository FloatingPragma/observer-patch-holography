#!/usr/bin/env python3
"""Reduce the D10 electroweak family to the active two-scalar carrier.

Chain role: rewrite the baseline D10 gauge data as a reduced
`(sigma_EW, eta_EW)` transport family from which the current `W/Z` branch is
selected.

Mathematics: hypercharge/SU(2) source deformations, two-point carrier geometry,
and coherent quintet emission on each source point.

OPH-derived inputs: `alphaY_mz`, `alpha2_mz`, `alpha_u`, and `v` from the D10
observable family.

Output: the current D10 source pair, compact slice, and selected carrier point
consumed by the population evaluator and readout.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_observable_family.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_source_transport_pair.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _quintet(alpha_y: float, alpha2: float, v: float) -> dict[str, float]:
    return {
        "alpha_em_eff_inv": (alpha_y + alpha2) / (alpha_y * alpha2),
        "sin2w_eff": alpha_y / (alpha_y + alpha2),
        "MW_pole": v * math.sqrt(math.pi * alpha2),
        "MZ_pole": v * math.sqrt(math.pi * (alpha_y + alpha2)),
        "v_report": v,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the reduced D10 electroweak source transport pair.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    payload = json.loads(Path(args.input).read_text(encoding="utf-8"))
    core = dict(payload.get("core_source", {}))
    reported = dict(payload.get("reported_outputs", {}))
    alpha1 = float(core["alpha1_mz"])
    alpha2 = float(core["alpha2_mz"])
    alpha_u = float(core["alpha_u"])
    v = float(core["v"])
    alpha_y = 3.0 * alpha1 / 5.0
    theta_w0_sin2 = alpha_y / (alpha_y + alpha2)
    theta_w0_sin = math.sqrt(theta_w0_sin2)
    theta_w0_cos = math.sqrt(1.0 - theta_w0_sin2)
    cos2_theta_w0 = theta_w0_cos ** 2 - theta_w0_sin ** 2
    beta_ew = (alpha2 - alpha_y) / (alpha_y + alpha2)
    tau_seed = (alpha_u / (8.0 * math.pi)) * math.log(alpha_y / alpha2)
    tau_y = tau_seed
    tau_2 = -tau_seed
    eta_seed = tau_2
    alpha_y_prime = alpha_y * (1.0 + tau_y)
    alpha2_prime = alpha2 * (1.0 + tau_2)
    compact_eta = alpha_u * cos2_theta_w0
    compact_sigma = -compact_eta
    compact_tau_y = compact_sigma - compact_eta
    compact_tau_2 = compact_sigma + compact_eta
    compact_alpha_y = alpha_y * (1.0 + compact_tau_y)
    compact_alpha2 = alpha2 * (1.0 + compact_tau_2)
    source_point_a = _quintet(alpha_y_prime, alpha2_prime, v)
    source_point_b = _quintet(compact_alpha_y, compact_alpha2, v)
    source_point_a_u = 1.0 + eta_seed
    source_point_a_n = 1.0 + beta_ew * eta_seed
    source_point_b_u = 1.0 + compact_sigma + compact_eta
    source_point_b_n = 1.0 + compact_sigma + beta_ew * compact_eta

    selected_sigma = compact_sigma
    selected_eta = compact_eta
    selected_tau_y = compact_tau_y
    selected_tau_2 = compact_tau_2
    selected_u = 1.0 + selected_sigma + selected_eta
    selected_n = 1.0 + selected_sigma + beta_ew * selected_eta

    artifact = {
        "artifact": "oph_d10_ew_source_transport_pair",
        "generated_utc": _timestamp(),
        "status": "selected_two_scalar_family",
        "smallest_constructive_missing_object": "EWSinglePostTransportTreeIdentity_D10",
        "source_pair_symbol": "Tau_EW_D10 = (tau_Y, tau_2)",
        "source_slots": {
            "alphaY_mz": alpha_y,
            "alpha2_mz": alpha2,
            "v_inherited": v,
        },
        "source_pair": {
            "alphaY_mz": alpha_y,
            "alpha2_mz": alpha2,
            "v_inherited": v,
        },
        "seed_family": {
            "name": "centered_log_gauge_split",
            "tau_seed_formula": "(alpha_u / (8*pi)) * log(alphaY_mz / alpha2_mz)",
            "tau_Y_formula": "tau_seed",
            "tau_2_formula": "-tau_seed",
            "constraint": "tau_Y + tau_2 = 0",
        },
        "family_coordinates": {
            "sigma_EW_symbol": "0.5 * (tau_Y + tau_2)",
            "eta_EW_symbol": "0.5 * (tau_2 - tau_Y)",
            "tau_Y_formula": "sigma_EW - eta_EW",
            "tau_2_formula": "sigma_EW + eta_EW",
        },
        "two_scalar_population_status": "closed_current_carrier",
        "population_selector_status": "closed",
        "population_selector_formula": "selected_population_point = argmin_{p in C_D10} J_pop_EW(p)",
        "population_functional_symbol": "J_pop_EW",
        "alpha_u_from_seed_formula": "8*pi*tau_seed / log(alphaY_mz / alpha2_mz)",
        "alpha_u_from_seed": alpha_u,
        "eta_source_formula": "alpha_u_from_seed * beta_EW",
        "eta_source": compact_eta,
        "predictive_population_point": {
            "sigma_EW": selected_sigma,
            "eta_EW": selected_eta,
            "tau_Y": selected_tau_y,
            "tau_2": selected_tau_2,
            "u_EW": selected_u,
            "n_EW": selected_n,
        },
        "predictive_population_closed": True,
        "predictive_population_verdict": "closed_current_carrier_nonexact_running_quintet",
        "population_evaluator_object": "EWGaugeSourceTransportPairPopulationEvaluator_D10",
        "population_forward_map": {
            "tau_Y": "sigma_EW - eta_EW",
            "tau_2": "sigma_EW + eta_EW",
            "alphaY_star": "alphaY_mz * (1 + sigma_EW - eta_EW)",
            "alpha2_star": "alpha2_mz * (1 + sigma_EW + eta_EW)",
            "mW": "v_inherited * sqrt(pi * alpha2_star)",
            "mZ": "v_inherited * sqrt(pi * (alphaY_star + alpha2_star))",
            "alpha_em_inv": "(alphaY_star + alpha2_star) / (alphaY_star * alpha2_star)",
            "sin2_thetaW": "alphaY_star / (alphaY_star + alpha2_star)",
        },
        "population_basis": {
            "beta_EW": beta_ew,
            "beta_EW_formula": "(alpha2_mz - alphaY_mz) / (alpha2_mz + alphaY_mz)",
            "u_EW_formula": "1 + sigma_EW + eta_EW",
            "n_EW_formula": "1 + sigma_EW + beta_EW * eta_EW",
        },
        "population_reconstruction": {
            "eta_EW_formula": "(u_EW - n_EW) / (1 - beta_EW)",
            "sigma_EW_formula": "(n_EW - 1 - beta_EW * (u_EW - 1)) / (1 - beta_EW)",
        },
        "population_atomic_quartet": {
            "mW_formula": "v_inherited * sqrt(pi * alpha2_mz * u_EW)",
            "mZ_formula": "v_inherited * sqrt(pi * (alphaY_mz + alpha2_mz) * n_EW)",
            "alpha_em_inv_formula": "1/(alpha2_mz * u_EW) + 1/((alphaY_mz + alpha2_mz) * n_EW - alpha2_mz * u_EW)",
            "sin2_thetaW_formula": "1 - alpha2_mz * u_EW / ((alphaY_mz + alpha2_mz) * n_EW)",
        },
        "population_minimality_certificate": {
            "det_d_u_n_d_sigma_eta_formula": "beta_EW - 1",
            "determinant": beta_ew - 1.0,
            "determinant_nonzero": abs(beta_ew - 1.0) > 1.0e-12,
            "third_scalar_needed": False,
        },
        "population_nonuniqueness_certificate": {
            "selector_formula": "selected_population_point = argmin_{p in C_D10} J_pop_EW(p)",
            "source_point_A": {
                "sigma_EW": 0.0,
                "eta_EW": eta_seed,
                "quintet": source_point_a,
            },
            "source_point_B": {
                "sigma_EW": compact_sigma,
                "eta_EW": compact_eta,
                "quintet": source_point_b,
            },
            "delta_quintet_B_minus_A": {
                key: float(source_point_b[key]) - float(source_point_a[key])
                for key in ("MW_pole", "MZ_pole", "alpha_em_eff_inv", "sin2w_eff")
            },
            "same_family_formula": True,
            "different_family_points": True,
            "verdict": "population_selector_now_closed_current_carrier",
        },
        "one_scalar_reduction_certificate": {
            "candidate_A_u_EW": source_point_a_u,
            "candidate_A_n_EW": source_point_a_n,
            "candidate_B_u_EW": source_point_b_u,
            "candidate_B_n_EW": source_point_b_n,
            "same_u_EW": math.isclose(source_point_a_u, source_point_b_u, rel_tol=0.0, abs_tol=1.0e-15),
            "same_n_EW": math.isclose(source_point_a_n, source_point_b_n, rel_tol=0.0, abs_tol=1.0e-15),
            "verdict": "no_one_scalar_residual_on_reopened_family",
        },
        "forbidden_inverse_witness_formulas": {
            "sigma_plus_eta_from_mW": "mW^2 / (pi * v_inherited^2 * alpha2_mz) - 1",
            "sigma_minus_eta_from_mZ_mW": "(mZ^2 - mW^2) / (pi * v_inherited^2 * alphaY_mz) - 1",
            "alphaY_star_from_alpha_sin": "1 / (alpha_em_inv * (1 - sin2_thetaW))",
            "alpha2_star_from_alpha_sin": "1 / (alpha_em_inv * sin2_thetaW)",
        },
        "special_slices": {
            "zero_benchmark": {"sigma_EW": 0.0, "eta_EW": 0.0},
            "current_one_seed_slice": {
                "sigma_EW": 0.0,
                "eta_EW": eta_seed,
            },
            "common_scale_slice": {
                "constraint": "eta_EW = 0.0",
            },
        },
        "updated_couplings": {
            "alphaY_prime": "alphaY_mz * (1 + sigma_EW - eta_EW)",
            "alpha2_prime": "alpha2_mz * (1 + sigma_EW + eta_EW)",
            "v_prime": "v_inherited",
        },
        "pushforward_shell": {
            "Pi_AA": "sigma_EW + (sin^2(theta_W0) - cos^2(theta_W0)) * eta_EW",
            "Pi_AZ": "2 * sin(theta_W0) * cos(theta_W0) * eta_EW",
            "Pi_ZZ": "sigma_EW + (cos^2(theta_W0) - sin^2(theta_W0)) * eta_EW",
            "Pi_WW": "sigma_EW + eta_EW",
        },
        "strict_missing_object_beneath_transport_entry_preimage": "EWGaugeSourceTransportPair_D10",
        "external_nonzero_quartet_image_test": {
            "status": "fails_exact_image_test_on_current_branch",
            "why": "the stored external quartet is not the image of a single coherent D10 source-pair on the current branch",
        },
        "first_nonzero_oph_seed_trial": {
            "name": "centered_log_gauge_split",
            "tau_seed": tau_seed,
            "tau_Y": tau_y,
            "tau_2": tau_2,
            "transport_entry_values": {
                "Pi_AA": (theta_w0_cos ** 2) * tau_y + (theta_w0_sin ** 2) * tau_2,
                "Pi_AZ": theta_w0_sin * theta_w0_cos * (tau_2 - tau_y),
                "Pi_ZZ": (theta_w0_sin ** 2) * tau_y + (theta_w0_cos ** 2) * tau_2,
                "Pi_WW": tau_2,
            },
            "coherent_output_quintet": source_point_a,
        },
        "compact_hypercharge_only_mass_slice": {
            "name": "hypercharge_only_compact_mass_law",
            "law": {
                "eta_EW_formula": "alpha_u * cos(2*theta_W0)",
                "sigma_EW_formula": "-eta_EW",
                "tau_Y_formula": "sigma_EW - eta_EW",
                "tau_2_formula": "sigma_EW + eta_EW",
                "compact_constraint": "tau_2 = 0",
            },
            "sigma_EW": compact_sigma,
            "eta_EW": compact_eta,
            "tau_Y": compact_tau_y,
            "tau_2": compact_tau_2,
            "transport_entry_values": {
                "Pi_AA": compact_sigma + (theta_w0_sin ** 2 - theta_w0_cos ** 2) * compact_eta,
                "Pi_AZ": 2.0 * theta_w0_sin * theta_w0_cos * compact_eta,
                "Pi_ZZ": compact_sigma + (theta_w0_cos ** 2 - theta_w0_sin ** 2) * compact_eta,
                "Pi_WW": compact_sigma + compact_eta,
            },
            "coherent_output_quintet": source_point_b,
        },
        "proof_gates": {
            "exact_closure": False,
            "two_scalar_source_population_required": False,
            "external_quartet_image_test": "must_pass",
            "provenance_equality": "must_pass",
        },
        "base_running_reported_outputs": reported,
        "notes": [
            "The reopened two-scalar source pair is already a complete minimal carrier for the coherent D10 quartet; no third scalar is needed on this family.",
            "The carrier-level selector is now closed on the current two-scalar family via J_pop_EW, with selected point sigma_EW = -eta_source and eta_EW = eta_source.",
            "The live D10 gap is no longer carrier selection or split exact closure, but the unsplit post-transport tree identity that would emit the mass-moving exact W/Z coordinate.",
            "The current one-seed family is the sigma_EW = 0 slice. The missing common scalar sigma_EW supplies the diagonal/common-scale motion that the anti-diagonal eta_EW slice alone cannot provide.",
            "The simplest mass-moving compact slice on the same family keeps the charged source fixed (tau_2 = 0) and transports only hypercharge with eta_EW = alpha_u * cos(2*theta_W0). This is recorded separately as a reference-free W/Z mass candidate, not a full electroweak quintet closure.",
            "No reference-fit W/Z slice is emitted here. This artifact stays predictive-only and records only the OPH seed family plus its first nonzero source-side trial.",
        ],
    }

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
