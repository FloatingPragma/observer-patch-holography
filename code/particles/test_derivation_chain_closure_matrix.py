#!/usr/bin/env python3
"""Smoke tests for the derivation-chain closure matrix."""

from __future__ import annotations

import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "scripts"))

from build_derivation_chain_closure_matrix import build_payload  # noqa: E402


def test_derivation_chain_closure_matrix_keeps_stage_gates_explicit() -> None:
    payload = build_payload()

    assert payload["artifact"] == "oph_particle_derivation_chain_closure_matrix"
    assert payload["status"] == "executable_nonhadron_chain_matrix_emitted"
    assert payload["closure_summary"]["all_derivation_chains_claimed_closed"] is False
    assert payload["closure_summary"]["source_backend_absent_chains"] == ["hadrons"]
    assert payload["closure_summary"]["empirical_closure_policy_chains"] == ["hadrons"]
    assert payload["worker_policy"]["chrome_pro_workers_needed_now"] is False
    assert payload["provenance_status"] == "closed_provenance_ledger_and_declared_sensitivity_taxonomy"
    rows = {row["chain"]: row for row in payload["rows"]}
    assert set(rows) == {
        "p_closure_root",
        "conditional_classical_carrier_modes",
        "electroweak_massive_bosons",
        "hierarchy_naturality_bridge",
        "higgs_top_declared_surface",
        "charged_leptons",
        "selected_class_quarks",
        "neutrino_absolute_attachment",
        "hadrons",
    }
    assert rows["p_closure_root"]["promotable"] is False
    assert rows["p_closure_root"]["open_gates"] == []
    assert rows["p_closure_root"]["closed_issue_refs"] == [224]
    assert rows["p_closure_root"]["stage_gate"] == (
        "populated source spectral measure payload + full interval certificate"
    )
    carrier_row = rows["conditional_classical_carrier_modes"]
    assert carrier_row["status"] == "closed_scoped_classical_modes_quantum_particle_gate_open"
    assert carrier_row["promotable"] is False
    assert carrier_row["open_gates"]
    assert carrier_row["outputs"]["photon"]["hard_quadratic_mass_parameter_squared"] == 0
    assert carrier_row["outputs"]["photon"]["quantum_particle_gate"] == "not_passed"
    assert rows["electroweak_massive_bosons"]["outputs"] == {}
    assert rows["electroweak_massive_bosons"]["promotable"] is False
    assert rows["electroweak_massive_bosons"]["status"] == (
        "no_public_prediction_row_target_free_d10_repair_not_promoted"
    )
    assert rows["hierarchy_naturality_bridge"]["status"] == (
        "closed_selected_branch_local_global_hierarchy_naturality"
    )
    assert rows["hierarchy_naturality_bridge"]["promotable"] is True
    assert rows["hierarchy_naturality_bridge"]["open_gates"] == []
    assert rows["hierarchy_naturality_bridge"]["closed_issue_refs"] == [332, 335, 337]
    assert rows["hierarchy_naturality_bridge"]["outputs"]["epsilon_H"] == "0"
    assert rows["hierarchy_naturality_bridge"]["full_theorem_grade_resonance_promoted"] is True
    assert rows["charged_leptons"]["status"] == "closed_current_corpus_charged_end_to_end_no_go"
    assert rows["charged_leptons"]["promotable"] is False
    assert rows["charged_leptons"]["open_gates"] == []
    assert rows["charged_leptons"]["closed_issue_refs"] == [201]
    assert rows["higgs_top_declared_surface"]["status"] == "conditional_declared_surface_higgs_top_candidate"
    assert rows["higgs_top_declared_surface"]["promotable"] is False
    assert rows["higgs_top_declared_surface"]["open_gates"] == [
        "closed_promotable_EWTargetFreeRepairValueLaw_D10"
    ]
    assert rows["higgs_top_declared_surface"]["closed_issue_refs"] == [207]
    assert rows["selected_class_quarks"]["status"] == "selected_class_target_anchored_exact_witness_not_strict_source"
    assert rows["selected_class_quarks"]["promotable"] is False
    assert rows["selected_class_quarks"]["open_gates"] == [
        "QUARK_SIGMA_SOURCE_SELECTOR",
        "NO_TARGET_LEAK_DAG_QUARK_SIGMA_SOURCE",
    ]
    assert rows["selected_class_quarks"]["closed_issue_refs"] == [199, 207, 212]
    assert rows["neutrino_absolute_attachment"]["status"] == (
        "rejected_target_informed_weighted_cycle_candidate"
    )
    assert rows["neutrino_absolute_attachment"]["promotable"] is False
    assert rows["neutrino_absolute_attachment"]["open_gates"] == [
        "source_emitted_family_transport_kernel",
        "source_derived_weight_exponent_and_cycle_matrix_law",
        "source_derived_basis_permutation_and_holonomy_orientation",
        "pre_reference_hash_lock",
    ]
    assert rows["hadrons"]["status"] == "source_backend_absent_empirical_closure_policy_emitted"
    assert rows["hadrons"]["open_gates"] == []
    assert rows["hadrons"]["closed_issue_refs"] == [153, 157]
    assert "p_closure_root" in payload["closure_summary"]["remaining_nonclosed_chains"]
    assert "higgs_top_declared_surface" in payload["closure_summary"]["remaining_nonclosed_chains"]
    assert "selected_class_quarks" in payload["closure_summary"]["remaining_nonclosed_chains"]
    assert "neutrino_absolute_attachment" in payload["closure_summary"]["remaining_nonclosed_chains"]
    assert "charged_leptons" not in payload["closure_summary"]["remaining_nonclosed_chains"]
    assert "hadrons" not in payload["closure_summary"]["remaining_nonclosed_chains"]
    assert payload["closure_summary"]["source_backend_absent_chains"] == ["hadrons"]
    assert payload["particle_five_gates"]["199"]["state"] == "closed_current_corpus_global_classification_no_go"
    assert payload["particle_five_gates"]["201"]["state"] == "closed_current_corpus_charged_end_to_end_no_go"
    assert payload["particle_five_gates"]["153"]["state"] == "closed_out_of_scope_computationally_blocked"
    assert payload["particle_five_gates"]["223"]["state"] == "closed_blocker_isolated_endpoint_package"
    assert payload["particle_five_gates"]["235"]["state"] == "closed_blocker_isolated_source_residual_no_go"
    assert payload["particle_five_gates"]["32"]["state"] == "closed_declared_convention_contract"
    assert payload["particle_five_gates"]["224"]["state"] == "closed_canonical_guarded_trunk_adoption"
    assert payload["particle_five_gates"]["225"]["state"] == "closed_material_sync_no_live_publish"
    assert payload["particle_five_gates"]["234"]["state"] == "closed_provenance_ledger_and_declared_sensitivity_taxonomy"
    assert payload["particle_five_gates"]["332"]["state"] == "closed_exact_selected_branch"
    assert payload["particle_five_gates"]["335"]["state"] == "closed_full_local_global_hierarchy_resonance"
    assert payload["particle_five_gates"]["337"]["state"] == "closed_projection_bridge_with_exact_residual"
