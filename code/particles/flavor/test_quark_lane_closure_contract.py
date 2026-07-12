#!/usr/bin/env python3
"""Validate the quark lane's source-only obstruction contract."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SIGMA_OBSTRUCTION_SCRIPT = (
    ROOT / "particles" / "flavor" / "derive_quark_sigma_source_nonidentifiability_obstruction.py"
)
AXIOM_LEVEL_OBSTRUCTION_SCRIPT = (
    ROOT / "particles" / "flavor" / "derive_quark_axiom_level_yukawa_moduli_nonidentifiability.py"
)
SCHEME_OBSTRUCTION_SCRIPT = (
    ROOT / "particles" / "flavor" / "derive_quark_running_mass_scheme_convention_obstruction.py"
)
PUBLIC_EXACT_YUKAWA_THEOREM_SCRIPT = (
    ROOT / "particles" / "flavor" / "derive_quark_public_exact_yukawa_end_to_end_theorem.py"
)
PUBLIC_EXACT_YUKAWA_PROMOTION_FRONTIER_SCRIPT = (
    ROOT / "particles" / "flavor" / "derive_quark_public_exact_yukawa_promotion_frontier.py"
)
SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_lane_closure_contract.py"
OUTPUT = ROOT / "particles" / "runs" / "flavor" / "quark_lane_closure_contract.json"


def test_quark_lane_contract_records_two_independent_obstructions_and_audit_only_sidecars() -> None:
    for script in (
        SIGMA_OBSTRUCTION_SCRIPT,
        AXIOM_LEVEL_OBSTRUCTION_SCRIPT,
        SCHEME_OBSTRUCTION_SCRIPT,
        PUBLIC_EXACT_YUKAWA_THEOREM_SCRIPT,
        PUBLIC_EXACT_YUKAWA_PROMOTION_FRONTIER_SCRIPT,
        SCRIPT,
    ):
        subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT)

    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_quark_lane_closure_contract"
    assert payload["base_theorem_emitted_package_artifact"] == "oph_quark_maximal_theorem_emitted_package"
    assert payload["proof_status"] == "closed_sharper_obstructions_numeric_quark_predictions_withheld"
    assert payload["public_promotion_allowed"] is False
    assert payload["numeric_quark_prediction_rows_allowed"] is False

    source_obstruction = payload["source_spread_nonidentifiability_obstruction"]
    assert source_obstruction["artifact"] == "oph_quark_sigma_source_nonidentifiability_obstruction"
    assert source_obstruction["proof_status"] == "closed_exact_current_corpus_obstruction"
    assert source_obstruction["theorem_grade_obstruction"] is True
    assert source_obstruction["compatible_source_spread_fiber"] == "(R_{>0})^2"
    assert source_obstruction["fiber_dimension"] == 2
    assert source_obstruction["independent_coordinates"] == ["sigma_u", "sigma_d"]
    assert source_obstruction["source_only_sigma_emitted"] is False
    assert source_obstruction["numeric_quark_rows_allowed"] is False
    assert source_obstruction["github_issues"] == [377, 379, 380]
    assert source_obstruction["dependency_audit"]["no_target_leak"] is True
    assert source_obstruction["dependency_audit"]["allowed_forbidden_disjoint"] is True
    assert source_obstruction["minimal_future_extension"]["required_independent_scalar_count"] == 2
    assert "sigma_u_target" not in source_obstruction
    assert "sigma_d_target" not in source_obstruction

    axiom_obstruction = payload["axiom_level_yukawa_moduli_nonidentifiability"]
    assert axiom_obstruction["artifact"] == (
        "oph_quark_axiom_level_yukawa_moduli_nonidentifiability"
    )
    assert axiom_obstruction["proof_status"] == "closed_axiom_level_nondefinability_theorem"
    assert axiom_obstruction["additional_axioms_used"] is False
    assert axiom_obstruction["counterfamily"]["parameter_space"] == (
        "(lambda_u,lambda_d) in (R_{>0})^2"
    )
    assert axiom_obstruction["MAR_audit"]["counterfamily_members_have_equal_MAR_score"] is True
    assert axiom_obstruction["reference_data_policy"]["no_target_leak_by_construction"] is True
    assert axiom_obstruction["public_numeric_quark_rows_allowed"] is False

    scheme_obstruction = payload["running_scheme_and_physical_yukawa_obstruction"]
    assert scheme_obstruction["artifact"] == "oph_quark_running_mass_scheme_convention_obstruction"
    assert scheme_obstruction["proof_status"] == (
        "closed_structural_finite_renormalization_nonidentifiability_obstruction"
    )
    assert scheme_obstruction["github_issues"] == [381, 382]
    assert scheme_obstruction["reference_data_policy"]["no_target_leak_by_construction"] is True
    assert scheme_obstruction["row_partition"]["light_running_coordinates"] == ["u", "d", "s"]
    assert scheme_obstruction["row_partition"]["heavy_running_coordinates"] == ["c", "b"]
    assert scheme_obstruction["row_partition"]["separate_extraction_coordinates"] == ["t"]
    assert scheme_obstruction["row_partition"]["single_running_quark_sextet_claim_allowed"] is False
    matrix_audit = scheme_obstruction["stored_matrix_dimensional_audit"]
    assert matrix_audit["current_classification"] == "mixed_scheme_GeV_mass_texture_matrices"
    assert matrix_audit["stored_entry_dimension"] == "GeV"
    assert matrix_audit["certified_physical_yukawa_matrices"] is False
    assert matrix_audit["dimensionless_normalization_supplied"] is False

    exact_target = payload["exact_pdg_derivation_target"]
    assert exact_target["target_name"] == (
        "exact_mixed_convention_quark_target_packet_on_declared_current_family_transport_frame"
    )
    assert exact_target["status"] == "closed_target_anchored_algebraic_audit_only"
    assert exact_target["artifact"] == "oph_quark_current_family_end_to_end_exact_pdg_derivation_chain"
    assert exact_target["target_anchored"] is True
    assert exact_target["source_only_prediction"] is False
    assert exact_target["single_running_quark_sextet_claim_allowed"] is False

    exact_yukawa_target = payload["exact_yukawa_derivation_target"]
    assert exact_yukawa_target["target_name"] == "mixed_scheme_dimensionful_quark_mass_texture_audit"
    assert exact_yukawa_target["status"] == "closed_algebraic_mass_texture_audit_not_physical_yukawa"
    assert exact_yukawa_target["forward_certified"] is True
    assert exact_yukawa_target["matrix_kind"] == "mixed_scheme_GeV_mass_texture_matrices"
    assert exact_yukawa_target["certified_physical_yukawa_matrices"] is False

    public_target = payload["public_exact_yukawa_derivation_target"]
    assert public_target["status"] == "blocked_mixed_scheme_dimensionful_mass_textures_not_physical_yukawas"
    assert "QUARK_COMMON_SCALE_DIMENSIONLESS_YUKAWA_CERTIFICATE" in public_target[
        "minimal_exact_blocker_set"
    ]
    assert public_target["numeric_values_role"] == "target_anchored_audit_only"
    assert public_target["physical_yukawa_claim_allowed"] is False

    public_frontier = payload["public_final_theorem_frontier"]
    assert public_frontier["proof_status"] == "closed_exact_current_corpus_obstruction"
    assert public_frontier["closure_kind"] == "sharper_current_corpus_nonidentifiability_obstruction"
    assert public_frontier["public_promotion_allowed"] is False
    assert public_frontier["resolved_by_theorem_artifact"] == (
        "oph_quark_sigma_source_nonidentifiability_obstruction"
    )
    assert public_frontier["compatible_source_spread_fiber"] == "(R_{>0})^2"

    yukawa_frontier = payload["public_exact_yukawa_promotion_frontier"]
    assert yukawa_frontier["target_name"] == "selected_class_mixed_scheme_mass_texture_audit"
    assert yukawa_frontier["proof_status"] == (
        "blocked_mixed_scheme_dimensionful_mass_textures_not_physical_yukawas"
    )
    assert yukawa_frontier["public_promotion_allowed"] is False
    assert yukawa_frontier["source_spread_obstruction_artifact"] == (
        "oph_quark_sigma_source_nonidentifiability_obstruction"
    )
    assert yukawa_frontier["scheme_obstruction_artifact"] == (
        "oph_quark_running_mass_scheme_convention_obstruction"
    )
    assert yukawa_frontier["certified_physical_yukawa_matrices"] is False

    compression = payload["candidate_one_theorem_physical_compression"]
    assert compression["status"] == "closed_two_modulus_nonidentifiability_obstruction"
    assert compression["current_corpus_source_spread_fiber"] == "(R_{>0})^2"
    assert compression["remaining_source_dimension"] == 2
    assert compression["remaining_nonalgebraic_theorem"] == (
        "QUARK_SOURCE_SPREAD_PAIR_ACTION_BREAKING_THEOREM"
    )

    missing_ids = [item["id"] for item in payload["exact_missing_theorems"]]
    assert missing_ids == [
        "quark_source_spread_pair_action_breaking_theorem",
        "quark_rg_covariant_scheme_readout_or_invariant",
        "quark_common_scale_dimensionless_yukawa_certificate",
    ]

    frontier = payload["public_current_family_yukawa_frontier"]
    assert frontier["target_1_status"] == "closed_profile_ray_only_absolute_spreads_nonidentifiable"
    sigma_branch = frontier["theorem_grade_sigma_branch"]
    assert sigma_branch["artifact"] == "oph_quark_sigma_source_nonidentifiability_obstruction"
    assert sigma_branch["spread_emitter_status"] == "source_only_spread_pair_nonidentifiable"
    assert sigma_branch["compatible_fiber"] == "(R_{>0})^2"
    assert sigma_branch["source_only_sigma_emitted"] is False
    assert frontier["transport_reduction"]["conditional_remaining_scalar_after_external_sigma_fixing"] == (
        "Delta_ud_overlap"
    )
    assert "QUARK_COMMON_SCALE_DIMENSIONLESS_YUKAWA_CERTIFICATE" in frontier[
        "minimal_exact_blocker_set"
    ]
    assert frontier["scheme_and_dimensional_yukawa_boundary"]["matrix_kind"] == (
        "mixed_scheme_GeV_mass_texture_matrices"
    )
