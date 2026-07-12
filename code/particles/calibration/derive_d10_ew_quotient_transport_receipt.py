#!/usr/bin/env python3
"""Build the conditional QT1--QT5 D10 quotient-transport receipt.

The canonical numerical evaluator is imported from
``derive_d10_ew_target_free_repair_value_law``.  This receipt adds domain,
fibre-identity, chart-Jacobian, and provenance checks.  It does not manufacture
the finite quotient-path certificate: QT1--QT5 remain explicitly assumed and
not enumerated by the current repository.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path

try:
    from calibration.derive_d10_ew_target_free_repair_value_law import (
        evaluate_candidate_from_source_basis,
    )
except ModuleNotFoundError:  # direct execution from calibration/
    from derive_d10_ew_target_free_repair_value_law import (
        evaluate_candidate_from_source_basis,
    )


ROOT = Path(__file__).resolve().parents[2]
SOURCE_PAIR_JSON = ROOT / "particles" / "runs" / "calibration" / "d10_ew_source_transport_pair.json"
VALUE_LAW_JSON = (
    ROOT / "particles" / "runs" / "calibration" / "d10_ew_target_free_repair_value_law.json"
)
DEFAULT_OUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_quotient_transport_receipt.json"
CERTIFICATE_SCHEMA = (
    ROOT / "particles" / "calibration" / "d10_ew_quotient_path_certificate.schema.json"
)
EXPECTED_CERTIFICATE = (
    ROOT / "particles" / "runs" / "calibration" / "d10_ew_quotient_path_certificate.json"
)


def _now_utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _max_abs(values: list[float]) -> float:
    return max((abs(value) for value in values), default=0.0)


def build_artifact(source_pair: dict, value_law: dict) -> dict:
    pair = source_pair["source_pair"]
    compact = source_pair["compact_hypercharge_only_mass_slice"]
    quintet = compact["coherent_output_quintet"]
    alpha_2 = float(pair["alpha2_mz"])
    alpha_y = float(pair["alphaY_mz"])
    eta = float(compact["eta_EW"])
    v_value = float(quintet["v_report"])

    evaluated = evaluate_candidate_from_source_basis(
        alpha_2=alpha_2,
        alpha_y=alpha_y,
        eta_source=eta,
        v_value=v_value,
    )
    basis = evaluated["basis"]
    chart = evaluated["repair_chart"]
    couplings = evaluated["repaired_couplings"]
    outputs = evaluated["coherent_emitted_quintet"]

    tau_2 = float(chart["tau2_tree_exact"])
    delta_n = float(chart["delta_n_tree_exact"])
    tau_y = float(chart["tauY_fiber"])
    alpha_sum = alpha_2 + alpha_y
    n_fib = 1.0 + (alpha_2 * tau_2 + alpha_y * tau_y) / alpha_sum
    jacobian_det = float(outputs["MW_pole"]) * float(outputs["MZ_pole"]) / (
        4.0 * (1.0 + tau_2) * (n_fib + delta_n)
    )

    parallel_from_fibre = alpha_y * (tau_y + 2.0 * eta)
    alpha_y_recomposed = alpha_y * (1.0 + tau_y) + alpha_sum * delta_n
    artifact_residuals = [
        float(chart[key]) - float(value_law["repair_chart"][key])
        for key in ("tau2_tree_exact", "delta_n_tree_exact", "tauY_fiber")
    ]
    artifact_residuals.extend(
        float(couplings[key]) - float(value_law["repaired_couplings"][key])
        for key in (
            "delta_alpha2",
            "delta_alphaY_parallel",
            "delta_alphaY_perp",
            "alpha2_prime",
            "alphaY_prime",
        )
    )
    artifact_residuals.extend(
        float(outputs[key]) - float(value_law["coherent_emitted_quintet"][key])
        for key in ("MW_pole", "MZ_pole", "alpha_em_eff_inv", "sin2w_eff")
    )

    premises = [
        {
            "id": "QT1",
            "name": "two_channel_quotient_exhaustion",
            "certificate_evidence_status": "assumed_not_enumerated",
            "source_emitted_by_current_repo": False,
            "missing": "finite quotient canonicalizer and proof that no third or cross channel survives",
        },
        {
            "id": "QT2",
            "name": "primitive_mixed_cumulant_normalization",
            "certificate_evidence_status": "assumed_not_enumerated",
            "source_emitted_by_current_repo": False,
            "missing": "finite path measure and equiprobable four-slot proof",
        },
        {
            "id": "QT3",
            "name": "oriented_root_path_incidence_and_Z6_subtraction",
            "certificate_evidence_status": "assumed_not_enumerated",
            "source_emitted_by_current_repo": False,
            "missing": (
                "explicit charged/neutral path lists, exact 1/3 color measure, and central "
                "projector trace 1/6; without QT3 the displayed value-law coefficients remain free"
            ),
        },
        {
            "id": "QT4",
            "name": "sign_and_least_norm_fibre_closure",
            "certificate_evidence_status": "assumed_not_enumerated",
            "source_emitted_by_current_repo": False,
            "missing": "finite-carrier Gram and residual-pairing certificate",
        },
        {
            "id": "QT5",
            "name": "path_completeness_and_rigidity",
            "certificate_evidence_status": "assumed_not_enumerated",
            "source_emitted_by_current_repo": False,
            "missing": (
                "deformation enumeration or positive target-free MAR gap; without QT5 additional "
                "output-changing paths or deformations remain unexcluded even if QT3 is assumed"
            ),
        },
    ]

    return {
        "artifact": "oph_d10_ew_quotient_transport_receipt",
        "generated_utc": _now_utc(),
        "status": "exact_conditional_implication_QT1_QT5_source_entailment_open",
        "certificate_status": "certificate_assumed_not_enumerated",
        "certificate_verified": False,
        "conditional_theorem_algebra_verified": True,
        "unconditional_d10_source_theorem_closed": False,
        "source_entailment_status": "QT1_QT5_are_assumptions_not_current_source_emissions",
        "promotion_allowed": False,
        "certificate_contract": {
            "schema_path": "code/particles/calibration/d10_ew_quotient_path_certificate.schema.json",
            "expected_input_path": (
                "code/particles/runs/calibration/d10_ew_quotient_path_certificate.json"
            ),
            "schema_present": CERTIFICATE_SCHEMA.is_file(),
            "certificate_input_present": EXPECTED_CERTIFICATE.is_file(),
            "certificate_input_status": (
                "present_but_not_verified"
                if EXPECTED_CERTIFICATE.is_file()
                else "missing_no_path_enumeration_supplied"
            ),
            "schema_conformance_is_sufficient_for_source_evidence": False,
            "realized_carrier_enumeration_verifier_present": False,
            "synthetic_path_fixture_eligible_as_source_evidence": False,
            "simulation_or_exact_enumeration_required": True,
            "required_sections": [
                "realized_carrier_and_quotient_canonicalizer_hashes",
                "explicit_charged_and_neutral_path_lists",
                "exact_rational_color_measure_and_incidence_counts",
                "Z6_projector_and_normalized_trace",
                "fibre_Gram_and_residual_pairing",
                "complete_deformation_enumeration_and_positive_MAR_gap",
                "same_branch_target_free_source_DAG",
                "exact_verifier_receipt_bound_to_the_certificate_hash",
            ],
            "promotion_allowed": False,
        },
        "source_artifacts": {
            "source_pair": "code/particles/runs/calibration/d10_ew_source_transport_pair.json",
            "canonical_value_law": (
                "code/particles/runs/calibration/d10_ew_target_free_repair_value_law.json"
            ),
            "canonical_evaluator": (
                "code/particles/calibration/derive_d10_ew_target_free_repair_value_law.py::"
                "evaluate_candidate_from_source_basis"
            ),
        },
        "premises": premises,
        "domain_checks": {
            "alpha2_positive": alpha_2 > 0.0,
            "alphaY_positive": alpha_y > 0.0,
            "eta_positive": eta > 0.0,
            "v_positive": v_value > 0.0,
            "rho_in_open_unit_interval": 0.0 < float(basis["beta_EW"]) < 1.0,
            "repaired_couplings_positive": (
                float(couplings["alpha2_prime"]) > 0.0
                and float(couplings["alphaY_prime"]) > 0.0
            ),
        },
        "numerical_specialization": {
            "source_tuple_scope": "archived_legacy_D10_calibration_tuple_not_strict_P_cand_branch",
            "rho": basis["beta_EW"],
            "alpha_U": basis["alpha_u_seed"],
            "lambda_EW": basis["lambda_EW"],
            "tau2": tau_2,
            "delta_n": delta_n,
            "tauY": tau_y,
            "alpha2_hat": couplings["alpha2_prime"],
            "alphaY_hat": couplings["alphaY_prime"],
            "W_GeV": outputs["MW_pole"],
            "Z_GeV": outputs["MZ_pole"],
        },
        "value_law_consistency_checks": {
            "canonical_artifact_max_abs_residual": _max_abs(artifact_residuals),
            "parallel_fibre_identity_residual": float(couplings["delta_alphaY_parallel"])
            - parallel_from_fibre,
            "alphaY_recomposition_residual": float(couplings["alphaY_prime"])
            - alpha_y_recomposed,
            "mass_chart_jacobian_det": jacobian_det,
            "mass_chart_nondegenerate": jacobian_det > 0.0,
        },
        "relation_to_color_balanced_candidate": {
            "same_law": False,
            "models_are_alternatives_not_cumulative_gates": True,
            "statement": (
                "The sqrt(N_c)/2 and raw N_c color-balanced rule is a distinct leading-quadratic "
                "alternative. The complete archived law is not the simple c=1/(4*rho) competitor, "
                "and the present data comparison does not exclude the complete archived law."
            ),
        },
        "same_branch_source_packet_closed": False,
        "physical_mass_semantics": "canonically_normalized_D10_mass_chart_not_complex_pole",
        "physical_pole_attachment_closed": False,
        "physical_pole_attachment": {
            "closed": False,
            "required_convention": "s_B=(M_B-i*Gamma_B/2)^2",
            "analytic_continuation_and_riemann_sheet_fixed": False,
            "self_energy_and_pole_enclosure_present": False,
        },
        "notes": [
            "The correspondence verifier validates the conditional algebra and local chart; it proves none of QT1--QT5.",
            "QT3 is the coefficient-emission premise, whereas QT5 separately excludes additional admissible deformations.",
            "The numerical specialization uses the archived legacy D10 tuple, not a demonstrated same-branch evaluation at P_cand.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the conditional D10 quotient-transport receipt.")
    parser.add_argument("--source-pair", default=str(SOURCE_PAIR_JSON))
    parser.add_argument("--value-law", default=str(VALUE_LAW_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    artifact = build_artifact(_load_json(Path(args.source_pair)), _load_json(Path(args.value_law)))
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
