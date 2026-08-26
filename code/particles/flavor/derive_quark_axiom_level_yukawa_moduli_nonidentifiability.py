#!/usr/bin/env python3
"""Emit the current-signature quark-Yukawa non-identifiability argument.

This artifact records an explicit algebraic counterfamily of physically
inequivalent one-Higgs, three-generation Yukawa packages.  Its application to
OPH is conditional on the current registered A1--A3 source signature and
declared structural packet containing no typed constraint, dynamics, or
output map that separates or excludes that family.  The emitter serializes
the argument and its registry-audit assertions; it is not a machine proof of
the completeness of the registered source signature.

No quark reference value, fitted spread, selected-family target, or numerical
flavor template is loaded.  The result uses no additional axiom.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = (
    ROOT
    / "particles"
    / "runs"
    / "flavor"
    / "quark_axiom_level_yukawa_moduli_nonidentifiability.json"
)


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact() -> dict[str, Any]:
    return {
        "artifact": "oph_quark_axiom_level_yukawa_moduli_nonidentifiability",
        "generated_utc": _timestamp(),
        "proof_status": "conditional_current_registered_signature_nonidentifiability_argument",
        "claim_tier": "conditional_current_registered_source_signature_obstruction",
        "scope": "current_registered_A1_A3_signature_plus_declared_structural_SM_branch_and_fixed_P",
        "scope_kind": "current_registered_signature_not_all_future_OPH_completions",
        "machine_proof": False,
        "receipt_role": "serialization_and_regression_guard_for_a_conditional_paper_argument",
        "registered_signature_completeness_verified_by_emitter": False,
        "additional_axioms_used": False,
        "conditional_scope_premises_used": True,
        "theorem_grade_obstruction": True,
        "source_only_numeric_quark_spectrum_emitted": False,
        "public_numeric_quark_rows_allowed": False,
        "theorem_statement": (
            "Conditional on the current registered A1--A3 source signature, a fixed pixel closure P, "
            "and a declared structural Standard-Model package containing no orbit-separating output "
            "map or constraint, the registered antecedents do not define the quark Yukawa eigenvalues. "
            "For every generic admissible one-Higgs "
            "three-generation package there is a continuous positive-rescaling family with the same screen, "
            "agreement, reference and aggregation, gauge, anomaly, hypercharge, Higgs, generation, CKM, "
            "CP-capability, and weak-UV data but different Yukawa singular values and quark masses. Hence a unique "
            "quark mass spectrum does not factor through that current registered source signature."
        ),
        "counterfamily": {
            "baseline": (
                "Y_q = U_qL diag(exp(mu_q + sigma_q*v_q1), exp(mu_q + sigma_q*v_q2), "
                "exp(mu_q + sigma_q*v_q3)) U_qR^dagger, q in {u,d}, sum_i v_qi = 0"
            ),
            "family": (
                "Y_q(lambda_q) = U_qL diag(exp(mu_q + lambda_q*sigma_q*v_q1), "
                "exp(mu_q + lambda_q*sigma_q*v_q2), exp(mu_q + lambda_q*sigma_q*v_q3)) "
                "U_qR^dagger"
            ),
            "parameter_space": "(lambda_u,lambda_d) in (R_{>0})^2",
            "free_action": (
                "(lambda_u,lambda_d).(sigma_u,sigma_d) = "
                "(lambda_u*sigma_u,lambda_d*sigma_d)"
            ),
            "simple_spectrum_preserved": True,
            "left_and_right_frames_preserved": True,
            "CKM_matrix_preserved": True,
            "CP_capability_preserved": True,
            "normalized_determinants_preserved": (
                "det(exp(lambda_q*sigma_q*v_q)) = "
                "exp(lambda_q*sigma_q*sum_i(v_qi)) = 1"
            ),
            "quark_spreads_changed": True,
            "quark_mass_singular_values_changed": True,
        },
        "axiom_invariance_audit": {
            "classification": "conditional_registry_audit_assertions_not_machine_verified_by_this_emitter",
            "all_registered_axiom_data_preserved": True,
            "Axiom_1_screen_net": "unchanged by Yukawa-eigenvalue rescaling",
            "Axiom_2_observer_agreement": "the declared accepted-data meaning diagrams are unchanged by Yukawa-eigenvalue rescaling",
            "Axiom_3_information_projection": {
                "gauge_invariant_local_Yukawa_densities_allowed": True,
                "constraint_values_or_multipliers_numerically_fixed_by_axiom": False,
                "map_from_P_to_Yukawa_multipliers_supplied": False,
                "Yukawa_output_map_supplied": False,
                "optimizer_to_Yukawa_output_map_supplied": False,
                "model_space_selection_allowed": False,
            },
            "declared_structural_packet": {
                "Yukawa_completability_preserved": True,
                "intrinsic_CKM_CP_capability_preserved": True,
                "weak_sector_UV_counting_clause_preserved": True,
                "Yukawa_eigenvalues_are_declared_structural_coordinates": False,
                "different_Yukawa_invariants_are_physically_equivalent": False,
                "counterfamily_members_share_all_registered_source_data": True,
                "counterfamily_members_remain_physically_distinct": True,
            },
            "fixed_P": "unchanged across the counterfamily",
        },
        "why_A3_does_not_select_Yukawa_moduli": {
            "A3_selects_states_inside_one_fixed_feasible_space": True,
            "A3_contains_no_typed_optimizer_to_Yukawa_map": True,
            "positive_rescaling_has_no_smallest_positive_element": True,
            "infimum": "lambda_u=lambda_d=0",
            "zero_limit_effect": "massless or degenerate quarks rather than the observed spectrum",
            "adding_a_norm_entropy_cost_description_length_or_RG_functional": (
                "would add a new output-selection premise unless derived within a valid completion"
            ),
        },
        "proof_steps": [
            "Fix any generic admissible one-Higgs three-generation Yukawa package with simple spectra.",
            "Apply independent positive rescalings to the two centered log-spectrum profiles while fixing their frames.",
            "Gauge representations, anomaly cancellation, hypercharges, Higgs content, CKM data, and CP capability are unchanged.",
            "Conditional on the current registry audit, the registered A1--A3 data and declared structural checks are unchanged because no typed orbit-separating constraint or output map is present.",
            "The packages are physically inequivalent because their Yukawa singular values differ and physical equivalence preserves Yukawa invariants.",
            "Therefore the current registered signature, declared structural packet, and fixed P admit multiple quark spectra, contradicting uniqueness of a mass map on that registered signature.",
        ],
        "corollaries": {
            "unique_source_map_P_to_sigma_u_sigma_d_exists": False,
            "unique_source_map_P_to_six_quark_masses_exists": False,
            "selected_frame_representative_independence_breaks_counterfamily": False,
            "refinement_naturality_breaks_counterfamily": False,
            "A3_breaks_counterfamily_under_registered_contract": False,
            "artifact_level_two_modulus_obstruction_is_axiomatically_expected": True,
        },
        "reference_data_policy": {
            "direct_input_artifacts": [],
            "quark_reference_values_consumed": False,
            "PDG_or_API_rows_consumed": False,
            "current_family_targets_consumed": False,
            "selected_class_target_witnesses_consumed": False,
            "fitted_spreads_consumed": False,
            "numerical_flavor_template_consumed": False,
            "no_target_leak_by_construction": True,
        },
        "public_policy": {
            "numeric_prediction_status": "not_defined_by_current_axioms",
            "allowed_public_result": "conditional_current_registered_signature_nonidentifiability_argument",
            "target_reconstruction_tables": "audit_only_not_predictions",
        },
        "notes": [
            "The positive-rescaling algebra is exact; its OPH application is conditional on the current registered source signature and is not a theorem about every future completion.",
            "This emitter serializes assertions and regression fields; it does not machine-verify registry completeness or the paper argument.",
            "No additional axiom or numerical normalization principle is introduced.",
            "A completion retaining this family must separate or exclude its orbit. A derived constraint, dynamics or action, source functional, stronger carrier or representation, or independently justified premise can do so; the current registry contains none with that role.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build the conditional current-signature quark-Yukawa non-identifiability artifact."
    )
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    payload = build_artifact()
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
