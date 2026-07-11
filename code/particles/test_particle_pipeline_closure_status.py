#!/usr/bin/env python3
"""Smoke tests for the simplified particle pipeline closure status."""

from __future__ import annotations

import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "scripts"))

from build_particle_pipeline_closure_status import build_status  # noqa: E402


def test_particle_pipeline_closure_status_scope_locks_hadrons_and_workers() -> None:
    status = build_status()

    assert status["scope"]["source_only_hadrons_in_current_local_scope"] is False
    assert status["scope"]["empirical_hadron_closure_surface"] is True
    assert status["scope"]["chrome_workers_needed"] is False
    assert "GLORB/Echosahedron" in status["scope"]["hadron_scope_reason"]
    assert status["finalization_gates"]["obstruction_only_worker_result_allowed"] is True
    assert status["finalization_gates"]["empirical_hadron_closure_policy_documented"] is True
    assert status["finalization_gates"]["empirical_hadron_spectral_dataset_integrated"] is False
    assert status["finalization_gates"]["hierarchy_local_global_resonance_closed"] is True
    assert status["finalization_gates"]["higgs_naturality_defect_closed"] is True
    assert status["finalization_gates"]["pixel_screen_resonance_summary_closed"] is True
    assert status["finalization_gates"]["symmetry_only_particle_promotion_blocked"] is True
    assert status["artifacts"]["measured_endpoint_calibration"]["status"] == (
        "oph_plus_empirical_hadron_closure_endpoint"
    )
    assert status["artifacts"]["hierarchy_pixel_screen_resonance"]["exists"] is True
    assert status["artifacts"]["hierarchy_pixel_screen_resonance"]["status"] == (
        "closed_receipt_summary_from_existing_hierarchy_certificates"
    )
    assert status["artifacts"]["empirical_ee_hadrons_source_registry"]["exists"] is True
    assert status["artifacts"]["empirical_ee_hadronic_spectral_measure_schema"]["exists"] is True
    gates = {gate["issue"]: gate for gate in status["issue_gates"]}
    assert gates[536]["state"] == "closed_claim_scope_repaired_quantum_particle_gate_fail_closed"
    assert gates[536]["closable_now"] is True
    assert gates[536]["zero_gev_particle_rows_emitted"] is False
    assert gates[536]["quantum_particle_promotion_allowed"] is False
    assert gates[153]["state"] == "closed_out_of_scope_computationally_blocked"
    assert gates[153]["closable_now"] is True
    assert gates[153]["requires_oph_hardware_backend"] is True
    assert gates[153]["empirical_hadron_closure_allowed"] is True
    assert gates[153]["closed_as_out_of_scope"] is True
    assert gates[153]["chrome_workers"] == "do_not_use_for_backend_execution"
    assert gates[157]["state"] == "closed_out_of_scope_computationally_blocked"
    assert gates[157]["closable_now"] is True
    assert gates[157]["requires_oph_hardware_backend"] is True
    assert gates[157]["empirical_hadron_closure_allowed"] is True
    assert gates[157]["closed_as_out_of_scope"] is True
    assert gates[157]["chrome_workers"] == "do_not_use_for_backend_execution"
    assert gates[223]["state"] == "closed_blocker_isolated_endpoint_package"
    assert gates[332]["state"] == "closed_exact_selected_branch"
    assert gates[332]["closable_now"] is True
    assert gates[337]["state"] == "closed_projection_bridge_with_exact_residual"
    assert gates[337]["closable_now"] is True
    assert gates[335]["state"] == "closed_full_local_global_hierarchy_resonance"
    assert gates[335]["closable_now"] is True
    assert gates[223]["closable_now"] is True
    assert gates[223]["successor_issue"] == 235
    assert gates[223]["chrome_workers"] == "not_needed_for_closed_package"
    assert gates[235]["state"] == "closed_blocker_isolated_source_residual_no_go"
    assert gates[235]["closable_now"] is True
    assert gates[235]["closed_as_first_missing_lemma_isolated"] is True
    assert gates[235]["promotion_allowed"] is False
    assert gates[235]["source_spectral_reduction"] == "source_spectral_reduction_theorem_emitted_measure_payload_absent"
    assert gates[235]["minimal_new_payload"] == "oph_qcd_ward_projected_hadronic_spectral_measure"
    assert gates[235]["chrome_workers"] == "not_needed_until_source_spectral_measure_payload_exists"
    companion = {branch["label"]: branch for branch in status["companion_status_branches"]}
    assert companion["Strong CP"]["state"] == "open_theta_qcd_bar_theta_vanishing_gap"
    assert gates[224]["state"] == "closed_canonical_guarded_trunk_adoption"
    assert gates[32]["state"] == "closed_declared_convention_contract"
    assert gates[32]["closable_now"] is True
    assert gates[207]["state"] == "closed_current_corpus_codomain_no_go"
    assert gates[207]["closable_now"] is True
    assert gates[201]["state"] == "closed_current_corpus_charged_end_to_end_no_go"
    assert gates[201]["closable_now"] is True
    assert gates[199]["state"] == "closed_current_corpus_global_classification_no_go"
    assert gates[199]["closable_now"] is True
    assert gates[234]["state"] == "closed_provenance_ledger_and_declared_sensitivity_taxonomy"
    assert gates[155]["state"] == "open_theta_qcd_bar_theta_vanishing_gap"
    assert gates[155]["closable_now"] is False
    assert gates[117]["state"] == "rejected_candidate_source_basis_and_kernel_open"
    assert gates[117]["closable_now"] is False
    assert gates[198]["closable_now"] is True
    assert status["latest_nonhadron_predictions"]["higgs"]["value"] == 125.1995304097179
    assert status["latest_nonhadron_predictions"]["higgs"]["unit"] == "GeV"
    assert "photon" not in status["latest_nonhadron_predictions"]
    carrier_modes = {row["carrier_id"]: row for row in status["classical_carrier_modes"]}
    assert set(carrier_modes) == {"photon", "gluon", "graviton"}
    assert carrier_modes["graviton"]["quantum_particle_gate"]["passed"] is False
    assert "top_quark" not in status["latest_nonhadron_predictions"]
    assert "electron_neutrino" not in status["latest_nonhadron_predictions"]
    withheld = {row["particle_id"]: row for row in status["withheld_non_prediction_rows"]}
    assert withheld["top_quark"]["reason"] == "target_anchored_witness_kept_in_exact_fit_audit_not_public_prediction"
    assert withheld["electron_neutrino"]["reason"] == (
        "target_informed_candidate_rejected_by_correlated_profile"
    )
