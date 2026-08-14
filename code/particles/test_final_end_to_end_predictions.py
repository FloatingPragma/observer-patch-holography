#!/usr/bin/env python3
"""Smoke tests for the final current end-to-end prediction bundle."""

from __future__ import annotations

import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "scripts"))

from build_final_end_to_end_predictions import build_payload  # noqa: E402


def test_final_end_to_end_predictions_preserve_scientific_boundaries_and_values() -> None:
    payload = build_payload()

    assert payload["artifact"] == "oph_final_current_end_to_end_particle_predictions"
    assert payload["p_closure"]["may_feed_live_particle_predictions"] is False
    assert payload["hadron_policy"]["source_only_hadron_predictions_emitted"] is False
    assert payload["hadron_policy"]["empirical_hadron_closure_allowed_for_display"] is True
    assert "github_issues" not in payload["hadron_policy"]
    fine = payload["fine_structure"]
    assert fine["source_side_no_hadron_near_endpoint"]["row_class"] == (
        "source_side_no_hadron_near_endpoint"
    )
    assert fine["source_side_no_hadron_near_endpoint"]["alpha_inv"].startswith("137.035959513608567790")
    assert fine["source_side_no_hadron_near_endpoint"]["missing_hadronic_correction_alpha_inv"].startswith(
        "0.000039663391432209"
    )
    assert fine["source_side_no_hadron_near_endpoint"]["status"] == (
        "certified_source_root_plus_public_width_diagnostic"
    )
    assert fine["root_only_audit"]["alpha_inv"].startswith("136.994835177412937295")
    assert payload["fine_structure"]["source_only_oph"]["status"] == "not_emitted"
    assert payload["fine_structure"]["oph_plus_empirical_hadron_closure"]["row_class"] == (
        "oph_plus_empirical_hadron_closure"
    )
    assert payload["fine_structure"]["oph_plus_empirical_hadron_closure"]["alpha_inv"] == "137.035999177"
    assert payload["fine_structure"]["empirical_payload_policy"]["source_only_theorem_status"] == "not_promoted"
    assert (
        payload["fine_structure"]["empirical_payload_policy"][
            "raw_cross_section_data_used_in_this_calibration"
        ]
        is False
    )
    assert "finalization_gates" not in payload
    assert "particle_five_issue_gates" not in payload
    assert "companion_open_branches" not in payload
    assert "pipeline_status" not in payload["source_surfaces"]
    hierarchy = payload["hierarchy_and_naturality"]
    assert hierarchy["status"]["resonance_status"] == "exact_conditional_local_global_hierarchy_resonance"
    assert hierarchy["status"]["full_theorem_grade_resonance_promoted"] is False
    assert hierarchy["status"]["work_in_progress_receipts"] == [
        "source production of the unit cost, inverse pairing, and D-optimal selector",
        "HIERARCHY-SCREEN-READOUT: identify log(E_cell/v)=Gamma_screen and attach it to the alpha_U/B_EW branch",
    ]
    assert "independently of the screen sieve" in hierarchy["claim_boundary"][
        "count_independence"
    ]
    assert "HIERARCHY-SCREEN-READOUT" in hierarchy["claim_boundary"]["missing_bridge"]
    assert "Pi_EW=4P" in hierarchy["claim_boundary"]["missing_bridge"]
    assert hierarchy["factor_origins"]["higgs_naturality_defect"] == "0"
    assert hierarchy["local_global_bridge"]["bridge_residual"].strip("0.") == ""
    assert "discharge CP-1, CP-2, and CP-3" in hierarchy["claim_boundary"][
        "capacity_conditions"
    ]
    pixel_screen = hierarchy["pixel_screen_resonance"]
    assert pixel_screen["accepted"] is True
    assert pixel_screen["tile_identity"]["cell_count_formula"] == (
        "K_cell = N_CRC^EW / (P_star/4) = 4*N_CRC^EW/P_star"
    )
    assert pixel_screen["tile_identity"]["relative_reconstruction_error"] == "0"
    assert pixel_screen["dimensionless_de_sitter_coordinate"]["Lambda_a_cell_formula"] == (
        "Lambda_CRC*a_cell = 3*pi*P_star/N_CRC^EW"
    )
    assert pixel_screen["dimensionless_de_sitter_coordinate"]["relative_cell_coordinate_error"] == "0"
    assert pixel_screen["shared_12_24_port_alignment"]["identification_premise"] == (
        "HIERARCHY-SCREEN-READOUT"
    )
    alignment = pixel_screen["shared_12_24_port_alignment"]
    assert alignment["repair_register_slots"] == 24
    assert alignment["screen_ports"] == 12
    assert alignment["identification_status"] == "work_in_progress"
    predictions = {entry["particle_id"]: entry for entry in payload["predictions"]}
    assert "photon" not in predictions
    assert "gluon" not in predictions
    assert "graviton" not in predictions
    carrier_modes = {row["carrier_id"]: row for row in payload["classical_carrier_modes"]}
    assert set(carrier_modes) == {"photon", "gluon", "graviton"}
    assert carrier_modes["photon"]["hard_quadratic_mass_parameter_squared"] == 0
    assert carrier_modes["photon"]["particle_promotion_allowed"] is False
    assert "w_boson" not in predictions
    assert "z_boson" not in predictions
    assert "higgs" not in predictions
    conditional = {
        entry["particle_id"]: entry
        for entry in payload["conditional_nonpromotable_rows"]
    }
    assert conditional["higgs"]["value"] == 125.1995304097179
    assert conditional["higgs"]["promotable"] is False
    assert "electron" not in predictions
    assert "top_quark" not in predictions
    assert "electron_neutrino" not in predictions
    withheld = {row["particle_id"]: row for row in payload["withheld_non_prediction_rows"]}
    assert withheld["electron"]["reason"] == "target_anchored_witness_kept_in_exact_fit_audit_not_public_prediction"
    assert withheld["electron"]["public_theorem_value"] is None
    assert withheld["electron"]["formula_if_anchor_exists"] == "m_e(P)=exp(A_ch(P)-4.495209645475038)"
    assert "charged_determinant_trace_lift_attachment" in withheld["electron"]["missing_for_promotion"]
    charged_boundary = payload["charged_lepton_anchor_boundary"]
    assert charged_boundary["status"] == "missing_theorem"
    assert charged_boundary["required_identity"] == "3*A_ch(P)=sum_psi M_ch[psi]*log(q_psi(P))"
    assert charged_boundary["current_closed_chain"]["A_ch_to_charged_masses"] is True
    assert charged_boundary["current_closed_chain"]["P_to_A_ch"] is False
    axiom_boundary = payload["quark_axiom_level_nondefinability"]
    assert axiom_boundary["artifact"] == "oph_quark_axiom_level_yukawa_moduli_nonidentifiability"
    assert axiom_boundary["proof_status"] == "closed_axiom_level_nondefinability_theorem"
    assert axiom_boundary["additional_axioms_used"] is False
    assert axiom_boundary["counterfamily"]["parameter_space"] == (
        "(lambda_u,lambda_d) in (R_{>0})^2"
    )
    assert axiom_boundary["three_axiom_audit"]["all_registered_axiom_data_preserved"] is True
    assert axiom_boundary["three_axiom_audit"]["Axiom_3_information_projection"][
        "Yukawa_output_map_supplied"
    ] is False
    assert axiom_boundary["public_numeric_quark_rows_allowed"] is False
    assert withheld["top_quark"]["reason"] == "target_anchored_witness_kept_in_exact_fit_audit_not_public_prediction"
    assert withheld["electron_neutrino"]["reason"] == (
        "target_informed_candidate_rejected_by_correlated_profile"
    )
    assert payload["direct_top_auxiliary_comparison"]["value_policy"] == (
        "compare_only_codomain_values_withheld_from_final_prediction_output"
    )
    assert payload["direct_top_auxiliary_comparison"]["primary_top_codomain"] == "Q007TP4"
    assert payload["direct_top_auxiliary_comparison"]["bridge_status"] == (
        "hard_no_go_current_corpus_compare_only_direct_top_codomain"
    )
