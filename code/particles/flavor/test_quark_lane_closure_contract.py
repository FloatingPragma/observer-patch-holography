#!/usr/bin/env python3
"""Validate the exact three-step quark closure contract artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SELECTOR_VALUE_SCRIPT = ROOT / "particles" / "flavor" / "derive_light_quark_overlap_defect_value_law.py"
T1_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_d12_t1_value_law.py"
PHYSICAL_BRANCH_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_physical_branch_repair_theorem.py"
SELECTED_SHEET_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_current_family_selected_sheet_closure.py"
EXACT_READOUT_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_current_family_exact_readout.py"
BACKREAD_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_d12_internal_backread_continuation_closure.py"
SPREAD_MAP_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_spread_map.py"
OVERLAP_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_d12_overlap_transport_law.py"
AFFINE_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_current_family_affine_anchor_theorem.py"
SIGMA_TARGET_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_current_family_exact_sigma_target.py"
CURRENT_PDG_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_current_family_exact_pdg_theorem.py"
ABSOLUTE_COLLAPSE_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_absolute_readout_algebraic_collapse.py"
CURRENT_FAMILY_TRANSPORT_LIFT_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_current_family_transport_frame_sector_attached_lift.py"
CURRENT_FAMILY_PHYSICAL_SIGMA_THEOREM_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_current_family_transport_frame_physical_sigma_lift_theorem.py"
CURRENT_FAMILY_STRENGTHENED_PHYSICAL_SIGMA_THEOREM_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_current_family_transport_frame_strengthened_physical_sigma_lift_theorem.py"
CURRENT_FAMILY_ABSOLUTE_THEOREM_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_current_family_transport_frame_absolute_sector_readout_theorem.py"
CURRENT_FAMILY_TRANSPORT_LIGHT_RATIO_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_current_family_transport_frame_light_ratio_theorem.py"
CURRENT_FAMILY_TRANSPORT_D12_VALUE_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_current_family_transport_frame_d12_value_package.py"
CURRENT_FAMILY_TRANSPORT_COMPLETION_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_current_family_transport_frame_exact_pdg_completion.py"
CURRENT_FAMILY_TRANSPORT_FORWARD_YUKAWAS_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_current_family_transport_frame_exact_forward_yukawas.py"
CURRENT_FAMILY_TRANSPORT_YUKAWA_THEOREM_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_current_family_transport_frame_exact_yukawa_theorem.py"
EXACT_YUKAWA_END_TO_END_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_exact_yukawa_end_to_end_theorem.py"
SIGMA_SOURCE_REQUIRED_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_sigma_source_datum_no_target_leak_required.py"
PUBLIC_SIGMA_THEOREM_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_public_physical_sigma_datum_descent.py"
PUBLIC_EXACT_YUKAWA_THEOREM_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_public_exact_yukawa_end_to_end_theorem.py"
PUBLIC_EXACT_YUKAWA_PROMOTION_FRONTIER_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_public_exact_yukawa_promotion_frontier.py"
PUBLIC_STRENGTHENED_FRONTIER_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_public_strengthened_physical_sigma_lift_frontier.py"
CURRENT_FAMILY_END_TO_END_CHAIN_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_current_family_end_to_end_exact_pdg_derivation_chain.py"
SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_lane_closure_contract.py"
OUTPUT = ROOT / "particles" / "runs" / "flavor" / "quark_lane_closure_contract.json"


def test_quark_lane_closure_contract_records_closed_exact_pdg_target_and_open_public_promotion() -> None:
    for script in (
        SELECTOR_VALUE_SCRIPT,
        T1_SCRIPT,
        PHYSICAL_BRANCH_SCRIPT,
        EXACT_READOUT_SCRIPT,
        SELECTED_SHEET_SCRIPT,
        BACKREAD_SCRIPT,
        SPREAD_MAP_SCRIPT,
        OVERLAP_SCRIPT,
        AFFINE_SCRIPT,
        SIGMA_TARGET_SCRIPT,
        CURRENT_PDG_SCRIPT,
        ABSOLUTE_COLLAPSE_SCRIPT,
        CURRENT_FAMILY_TRANSPORT_LIFT_SCRIPT,
        CURRENT_FAMILY_PHYSICAL_SIGMA_THEOREM_SCRIPT,
        CURRENT_FAMILY_STRENGTHENED_PHYSICAL_SIGMA_THEOREM_SCRIPT,
        CURRENT_FAMILY_ABSOLUTE_THEOREM_SCRIPT,
        CURRENT_FAMILY_TRANSPORT_LIGHT_RATIO_SCRIPT,
        CURRENT_FAMILY_TRANSPORT_D12_VALUE_SCRIPT,
        CURRENT_FAMILY_TRANSPORT_COMPLETION_SCRIPT,
        CURRENT_FAMILY_TRANSPORT_FORWARD_YUKAWAS_SCRIPT,
        CURRENT_FAMILY_TRANSPORT_YUKAWA_THEOREM_SCRIPT,
        EXACT_YUKAWA_END_TO_END_SCRIPT,
        SIGMA_SOURCE_REQUIRED_SCRIPT,
        PUBLIC_SIGMA_THEOREM_SCRIPT,
        PUBLIC_STRENGTHENED_FRONTIER_SCRIPT,
        PUBLIC_EXACT_YUKAWA_THEOREM_SCRIPT,
        PUBLIC_EXACT_YUKAWA_PROMOTION_FRONTIER_SCRIPT,
        CURRENT_FAMILY_END_TO_END_CHAIN_SCRIPT,
        SCRIPT,
    ):
        subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT)

    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_quark_lane_closure_contract"
    assert payload["base_theorem_emitted_package_artifact"] == "oph_quark_maximal_theorem_emitted_package"
    assert payload["proof_status"] == "selected_class_exact_witness_blocked_by_target_derived_sigma_datum"
    assert payload["public_promotion_allowed"] is False
    assert payload["non_circularity_status"]["missing_source_object"] == (
        "quark_sigma_source_datum_no_target_leak_required"
    )
    exact_target = payload["exact_pdg_derivation_target"]
    assert exact_target["target_name"] == "exact_running_quark_sextet_on_declared_current_family_transport_frame"
    assert exact_target["status"] == "closed"
    assert exact_target["minimal_exact_blocker_set"] == []
    assert exact_target["artifact"] == "oph_quark_current_family_end_to_end_exact_pdg_derivation_chain"
    assert exact_target["wrapper_theorem"] == "oph_quark_exact_pdg_end_to_end_theorem"
    assert exact_target["exact_running_values_gev"]["t"] == 172.35235532883115
    assert exact_target["not_the_same_as"] == "target_free_public_physical_sheet_promotion"
    exact_yukawa_target = payload["exact_yukawa_derivation_target"]
    assert exact_yukawa_target["target_name"] == "exact_forward_quark_yukawas_on_declared_current_family_transport_frame"
    assert exact_yukawa_target["status"] == "closed"
    assert exact_yukawa_target["artifact"] == "oph_quark_current_family_transport_frame_exact_yukawa_theorem"
    assert exact_yukawa_target["wrapper_theorem"] == "oph_quark_exact_yukawa_end_to_end_theorem"
    assert exact_yukawa_target["forward_certified"] is True
    public_exact_yukawa_target = payload["public_exact_yukawa_derivation_target"]
    assert public_exact_yukawa_target["status"] == "blocked_target_derived_sigma_source_missing"
    assert public_exact_yukawa_target["artifact"] == "oph_quark_public_exact_yukawa_end_to_end_theorem"
    assert public_exact_yukawa_target["forward_yukawa_artifact"]["artifact"] == "oph_quark_current_family_transport_frame_exact_forward_yukawas"
    assert payload["selected_local_sheet_status"]["sigma_id"] == "sigma_ref"
    assert payload["exact_sidecar_mass_surface"]["scope"] == "current_family_only"
    assert payload["exact_sidecar_mass_surface"]["current_family_affine_anchor_theorem"] == "oph_quark_current_family_affine_anchor_theorem"
    assert payload["exact_sidecar_mass_surface"]["current_family_exact_pdg_theorem"] == "oph_quark_current_family_exact_pdg_theorem"
    assert payload["exact_sidecar_mass_surface"]["exact_outputs_gev"]["u"] == 0.0021599999999999996
    assert payload["exact_sidecar_mass_surface"]["exact_outputs_gev"]["d"] == 0.004700000000000002
    target_surface = payload["current_family_physical_target_surface"]
    assert target_surface["affine_anchor_theorem"]["artifact"] == "oph_quark_current_family_affine_anchor_theorem"
    assert target_surface["exact_sigma_target"]["artifact"] == "oph_quark_current_family_exact_sigma_target"
    assert target_surface["exact_pdg_reconstruction"]["artifact"] == "oph_quark_current_family_exact_pdg_theorem"
    assert target_surface["absolute_readout_algebraic_collapse"]["artifact"] == "oph_quark_absolute_readout_algebraic_collapse"
    assert target_surface["transport_frame_sector_attached_lift"]["artifact"] == "oph_quark_current_family_transport_frame_sector_attached_lift"
    assert target_surface["transport_frame_sector_attached_lift"]["sigma_id"] == "sigma_phys_transport_frame_current_family"
    assert target_surface["restricted_physical_sigma_lift_theorem"]["artifact"] == "oph_quark_current_family_transport_frame_physical_sigma_lift_theorem"
    assert target_surface["restricted_physical_sigma_lift_theorem"]["corresponds_to_global_contract"]["id"] == "quark_physical_sigma_ud_lift"
    assert target_surface["restricted_strengthened_physical_sigma_lift_theorem"]["artifact"] == "oph_quark_current_family_transport_frame_strengthened_physical_sigma_lift_theorem"
    assert target_surface["restricted_strengthened_physical_sigma_lift_theorem"]["compressed_global_contract"]["id"] == "strengthened_quark_physical_sigma_ud_lift"
    assert target_surface["restricted_absolute_sector_readout_theorem"]["artifact"] == "oph_quark_current_family_transport_frame_absolute_sector_readout_theorem"
    assert target_surface["restricted_absolute_sector_readout_theorem"]["corresponds_to_global_contract"]["id"] == "quark_absolute_sector_readout_theorem"
    assert target_surface["restricted_light_ratio_theorem"]["artifact"] == "oph_quark_current_family_transport_frame_light_ratio_theorem"
    assert abs(float(target_surface["restricted_light_ratio_theorem"]["ell_ud"]) - 0.7774542870199399) < 1.0e-12
    assert target_surface["restricted_d12_value_package"]["artifact"] == "oph_quark_current_family_transport_frame_d12_value_package"
    assert abs(float(target_surface["restricted_d12_value_package"]["closed_d12_scalars"]["t1"]) - 0.6478785725166165) < 1.0e-12
    assert target_surface["transport_frame_exact_pdg_completion"]["artifact"] == "oph_quark_current_family_transport_frame_exact_pdg_completion"
    assert target_surface["transport_frame_exact_pdg_completion"]["exact_running_values_gev"]["u"] == 0.0021600000000000005
    assert target_surface["exact_forward_yukawas"]["artifact"] == "oph_quark_current_family_transport_frame_exact_forward_yukawas"
    assert target_surface["exact_forward_yukawas"]["forward_certified"] is True
    assert target_surface["exact_yukawa_theorem"]["artifact"] == "oph_quark_current_family_transport_frame_exact_yukawa_theorem"
    assert target_surface["end_to_end_exact_pdg_derivation_chain"]["artifact"] == "oph_quark_current_family_end_to_end_exact_pdg_derivation_chain"
    assert target_surface["end_to_end_exact_pdg_derivation_chain"]["exact_running_values_gev"]["d"] == 0.004699999999999999
    public_frontier = payload["public_final_theorem_frontier"]
    assert public_frontier["artifact"] == "oph_quark_public_strengthened_physical_sigma_lift_frontier"
    assert public_frontier["proof_status"] == "blocked_target_derived_sigma_datum_descent"
    assert public_frontier["public_promotion_allowed"] is False
    assert public_frontier["resolved_by_theorem_artifact"] == "oph_quark_public_physical_sigma_datum_descent"
    assert public_frontier["final_public_theorem_candidate"]["id"] == "selected_bridge_fiber_sigma_descent_not_source_selection"
    assert public_frontier["final_public_theorem_candidate"]["induces_global_contract"]["id"] == (
        "strengthened_quark_physical_sigma_ud_lift"
    )
    assert public_frontier["alternate_upstream_route"]["id"] == "oph_generation_bundle_branch_generator_splitting"
    assert public_frontier["alternate_upstream_route"]["status"] == "upstream_alternative_route_currently_deprioritized"
    public_yukawa_frontier = payload["public_exact_yukawa_promotion_frontier"]
    assert public_yukawa_frontier["artifact"] == "oph_quark_public_exact_yukawa_promotion_frontier"
    assert public_yukawa_frontier["target_name"] == "selected_class_conditional_forward_quark_yukawas"
    assert public_yukawa_frontier["proof_status"] == "blocked_target_derived_sigma_source_missing"
    assert public_yukawa_frontier["public_promotion_allowed"] is False
    assert public_yukawa_frontier["resolved_by_theorem_artifact"] == "oph_quark_public_exact_yukawa_end_to_end_theorem"
    assert public_yukawa_frontier["final_public_theorem_candidate"]["id"] == "selected_bridge_fiber_sigma_descent_not_source_selection"
    assert public_yukawa_frontier["closed_public_endpoint"]["public_exact_outputs"]["forward_yukawa_artifact"]["artifact"] == (
        "oph_quark_current_family_transport_frame_exact_forward_yukawas"
    )
    compression = payload["candidate_one_theorem_physical_compression"]
    assert compression["status"] == "blocked_target_derived_sigma_source_missing"
    assert compression["artifact"] == "oph_quark_public_physical_sigma_datum_descent"
    assert compression["supporting_algebraic_collapse_artifact"] == "oph_quark_absolute_readout_algebraic_collapse"
    assert compression["remaining_nonalgebraic_theorem"] == (
        "QUARK_SIGMA_SOURCE_SELECTOR"
    )
    assert compression["remaining_exact_gap"] == "target_derived_sigma_datum_used_for_selected_class_exact_witness"
    assert payload["continuation_only_mass_sidecar"]["closed_mass_side_package"]["t1"] == 0.6715870378831591
    assert [item["id"] for item in payload["internalized_theorems"]] == [
        "light_quark_overlap_defect_value_law",
        "quark_d12_t1_value_law",
    ]
    assert payload["exact_missing_theorems"] == []
    frontier = payload["public_current_family_yukawa_frontier"]
    assert frontier["sharper_target_1_primitive"]["artifact"] == "oph_light_quark_overlap_defect_value_law"
    assert frontier["sharper_target_1_primitive"]["equivalent_ray_coordinate_presentation"] == "quark_d12_t1_value_law"
    assert frontier["minimal_exact_blocker_set"] == []
    assert frontier["target_1_status"] == "closed"
    assert frontier["theorem_grade_sigma_branch"]["artifact"] == "oph_family_excitation_spread_map"
    assert frontier["theorem_grade_sigma_branch"]["spread_emitter_status"] == "closed"
    assert frontier["transport_reduction"]["remaining_scalar_on_fixed_sigma_branch"] == "Delta_ud_overlap"
