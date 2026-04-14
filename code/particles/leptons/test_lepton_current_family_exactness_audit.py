#!/usr/bin/env python3
"""Validate the current-family charged-lepton exactness audit artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
READOUT_SCRIPT = ROOT / "particles" / "leptons" / "derive_lepton_log_spectrum_readout.py"
FORWARD_SCRIPT = ROOT / "particles" / "leptons" / "build_forward_charged_leptons.py"
OBSTRUCTION_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_sector_local_current_support_obstruction_certificate.py"
SUPPORT_EXTENSION_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_sector_local_minimal_source_support_extension_emitter.py"
COMPLETION_LAW_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_sector_local_support_extension_completion_law.py"
PAIR_READBACK_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_sector_local_support_extension_source_scalar_pair_readback.py"
ENDPOINT_RATIO_BREAKER_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_sector_local_support_extension_endpoint_ratio_breaker.py"
ABSOLUTE_SCALE_GAP_IDENTITY_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_absolute_scale_transport_gap_identity.py"
ABSOLUTE_SCALE_UNDERDETERMINATION_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_absolute_scale_underdetermination_theorem.py"
TRACE_LIFT_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_uncentered_trace_lift_scaffold.py"
DETERMINANT_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_determinant_line_section_extension.py"
ANCHOR_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_absolute_anchor_section.py"
COCYCLE_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_uncentered_trace_lift_cocycle_reduction.py"
EQUALIZER_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_physical_identity_mode_equalizer.py"
DESCENT_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_mu_physical_descent_reduction.py"
NO_GO_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_centered_operator_mu_phys_no_go.py"
ROUTE_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_post_promotion_absolute_closure_route.py"
FRONTIER_FACTORIZATION_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_absolute_frontier_factorization.py"
END_TO_END_IMPOSSIBILITY_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_end_to_end_impossibility_theorem.py"
AUDIT_SCRIPT = ROOT / "particles" / "leptons" / "derive_lepton_current_family_exactness_audit.py"
OUTPUT = ROOT / "particles" / "runs" / "leptons" / "lepton_current_family_exactness_audit.json"


def test_lepton_exactness_audit_identifies_common_shift_as_insufficient() -> None:
    subprocess.run([sys.executable, str(READOUT_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(FORWARD_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(OBSTRUCTION_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SUPPORT_EXTENSION_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(COMPLETION_LAW_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(PAIR_READBACK_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(ENDPOINT_RATIO_BREAKER_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(ABSOLUTE_SCALE_GAP_IDENTITY_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(ABSOLUTE_SCALE_UNDERDETERMINATION_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(TRACE_LIFT_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(DETERMINANT_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(ANCHOR_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(COCYCLE_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(EQUALIZER_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(DESCENT_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(NO_GO_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(ROUTE_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(FRONTIER_FACTORIZATION_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(END_TO_END_IMPOSSIBILITY_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(AUDIT_SCRIPT)], check=True, cwd=ROOT)

    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_lepton_current_family_exactness_audit"
    assert payload["centered_hierarchy_audit"]["ratio_invariant_under_common_shift"] is True
    assert payload["centered_hierarchy_audit"]["residual_norm"] > 1.0
    assert payload["centered_hierarchy_audit"]["target_sigma_total_log_per_side"] > payload["centered_hierarchy_audit"]["current_sigma_total_log_per_side"]
    assert payload["centered_hierarchy_audit"]["target_gap_pair"]["rho_gap_ratio"] > payload["centered_hierarchy_audit"]["current_gap_pair"]["rho_gap_ratio"]
    assert payload["active_readout_contract"]["hierarchy_mode"] == "two_scalar_ordered_gap_pair_candidate"
    assert abs(payload["active_readout_contract"]["eta_e_split_log_per_side"]) < 1.0e-12
    assert payload["two_scalar_on_current_sigma_audit"]["residual_norm_after_best_eta_on_current_sigma"] > 1.0
    rel_errors = payload["common_shift_audit"]["relative_errors_after_best_shift"]
    assert max(abs(value) for value in rel_errors) > 0.5
    assert payload["same_support_obstruction_audit"]["same_support_exhausted"] is True
    assert payload["same_support_obstruction_audit"]["sigma_support_gap"] > 0.0
    assert payload["smallest_constructive_missing_object"] == "eta_source_support_extension_log_per_side"
    assert payload["current_support_obstruction_certificate"]["proof_status"] == "current_support_obstruction_certificate_closed"
    assert payload["support_extension_emitter"]["artifact"] == "oph_charged_sector_local_minimal_source_support_extension_emitter"
    assert payload["support_extension_completion_law"]["artifact"] == "oph_charged_sector_local_support_extension_completion_law"
    assert payload["support_extension_endpoint_ratio_breaker"]["artifact"] == "oph_charged_sector_local_support_extension_endpoint_ratio_breaker"
    assert payload["support_extension_source_scalar_pair_readback"]["artifact"] == "oph_charged_sector_local_support_extension_source_scalar_pair_readback"
    assert payload["absolute_scale_gap_identity"]["artifact"] == "oph_charged_absolute_scale_transport_gap_identity"
    assert payload["absolute_scale_underdetermination_theorem"]["artifact"] == "oph_charged_absolute_scale_underdetermination_theorem"
    assert abs(payload["absolute_scale_gap_identity"]["identity_residual"]) < 1.0e-12
    assert payload["absolute_scale_closure_status"]["present_chain_under_determines_g_e"] is True
    assert payload["absolute_scale_closure_status"]["charged_absolute_equalizer_status"] == "NO_GO_COMMON_SHIFT"
    assert payload["absolute_scale_closure_status"]["honest_missing_transport_scalar"] == "A_ch"
    assert payload["absolute_scale_closure_status"]["honest_post_promotion_single_slot"] == "refinement_stable_uncentered_trace_lift"
    assert payload["absolute_scale_closure_status"]["honest_post_promotion_internal_carrier"] == (
        "scalar_affine_cocycle_primitive"
    )
    assert payload["absolute_scale_closure_status"]["honest_post_promotion_exact_descended_scalar"] == (
        "charged_physical_affine_scalar_mu"
    )
    assert payload["absolute_scale_closure_status"]["promotion_only_centered_operator_no_go"] == (
        "charged_centered_operator_cannot_emit_mu_phys"
    )
    assert payload["absolute_scale_closure_status"]["absolute_frontier_factorization_artifact"] == "oph_charged_absolute_frontier_factorization"
    assert payload["absolute_scale_closure_status"]["hard_reject"]["g_e"] == 0.6822819838027987
    assert payload["trace_lift_cocycle_reduction"]["artifact"] == "oph_charged_uncentered_trace_lift_cocycle_reduction"
    assert payload["trace_lift_cocycle_reduction"]["single_slot_preserved"] == "refinement_stable_uncentered_trace_lift"
    assert payload["trace_lift_cocycle_reduction"]["irreducible_new_degree_of_freedom"] == (
        "scalar affine cocycle primitive mu"
    )
    assert payload["trace_lift_physical_descent"]["artifact"] == "oph_charged_mu_physical_descent_reduction"
    assert payload["trace_lift_physical_descent"]["exact_smaller_missing_object"] == (
        "charged_physical_affine_scalar_mu"
    )
    assert payload["centered_operator_mu_no_go"]["artifact"] == "oph_charged_centered_operator_mu_phys_no_go"
    assert payload["centered_operator_mu_no_go"]["theorem_id"] == "charged_centered_operator_cannot_emit_mu_phys"
    assert payload["absolute_frontier_factorization"]["artifact"] == "oph_charged_absolute_frontier_factorization"
    assert payload["absolute_frontier_factorization"]["current_surface_missing_object"] == "charged_absolute_anchor_A_ch"
    assert payload["absolute_frontier_factorization"]["post_promotion_single_slot"] == "refinement_stable_uncentered_trace_lift"
    assert payload["absolute_frontier_factorization"]["post_promotion_internal_carrier"] == (
        "scalar_affine_cocycle_primitive"
    )
    assert payload["absolute_frontier_factorization"]["post_promotion_exact_descended_scalar"] == (
        "charged_physical_affine_scalar_mu"
    )
    assert payload["absolute_frontier_factorization"]["promotion_only_no_go"] == (
        "charged_centered_operator_cannot_emit_mu_phys"
    )
    assert payload["end_to_end_closure_decision"]["artifact"] == "oph_charged_end_to_end_impossibility_theorem"
    assert payload["end_to_end_closure_decision"]["verdict"] == "no_current_corpus_end_to_end_closure"
    assert payload["end_to_end_closure_decision"]["closure_now"] is False
    assert payload["end_to_end_closure_decision"]["post_promotion_route_artifact"] == "oph_charged_post_promotion_absolute_closure_route"
    assert payload["end_to_end_closure_decision"]["exact_irreducible_chain"][1]["id"] == "refinement_stable_uncentered_trace_lift"
    assert payload["end_to_end_closure_decision"]["induced_after_irreducible_chain"]["charged_absolute_anchor"]["id"] == "charged_absolute_anchor_A_ch"
    assert payload["charged_sector_response_operator_candidate"]["name"] == "C_hat_e^{cand}"
    assert payload["charged_sector_response_operator_candidate"]["smallest_missing_clause"] == "compression_descendant_commutator_vanishes_or_is_uniformly_quadratic_small_after_central_split"
    assert payload["charged_sector_response_operator_candidate"]["latent_in_flavor_chain"] is True
    assert payload["red_team_branch_verdict"]["smallest_missing_theorem_object"] == "oph_generation_bundle_branch_generator_splitting"
    assert payload["red_team_branch_verdict"]["post_promotion_single_slot"] == "refinement_stable_uncentered_trace_lift"
    waiting = payload["exact_waiting_set"]
    assert waiting["mandatory_package_a"]["id"] == "charged_sector_response_pushforward_to_C_hat_e"
    assert waiting["mandatory_package_a"]["status"] == "resolved_as_exact_current_corpus_obstruction_certificate"
    assert waiting["mandatory_package_a"]["closed_issue"] == "papers.compact.e.30a-resolve-branch-generator-splitting-or-no-go"
    assert waiting["mandatory_package_a"]["exact_obstruction_certificate_artifact"] == (
        "oph_generation_bundle_branch_generator_splitting_obstruction_certificate"
    )
    assert waiting["mandatory_package_a"]["blocked_candidate_object"] == "C_hat_e^{cand}"
    assert waiting["mandatory_package_b"]["id"] == "refinement_stable_uncentered_trace_lift"
    assert waiting["mandatory_package_b"]["replaces_invalid_route"] == "charged_common_refinement_transport_equalizer"
    assert waiting["mandatory_package_c"]["id"] == "charged_absolute_anchor_A_ch"
    assert waiting["mandatory_package_c"]["status"] == "derived_once_package_b_exists"
    assert waiting["optional_package_d"]["id"] == "charged_holonomy_bridge_for_legacy_delta_2_over_9"
