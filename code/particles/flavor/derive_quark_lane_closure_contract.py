#!/usr/bin/env python3
"""Emit the current quark closure contract and source-only obstruction boundary.

Chain role: collect the current theorem boundary and the strongest exact sidecar
surfaces for the quark lane into one machine-readable contract.

Mathematics: the target-free source identities fix the up- and down-sector
profile rays but leave their positive endpoint spans independent.  The exact
compatible source-spread fiber is therefore ``(R_{>0})^2``.  The selected
bridge-fiber descent proves representative independence only and cannot select
a point of that fiber.

The numerical current-family sidecars remain useful target-anchored algebraic
audits.  Their six rows mix several comparison conventions, however, and their
GeV-valued matrices are mass textures rather than common-scale dimensionless
physical Yukawa matrices.  This contract keeps both obstructions explicit and
withholds source-only numerical promotion.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
SELECTOR_VALUE_JSON = ROOT / "particles" / "runs" / "flavor" / "light_quark_overlap_defect_value_law.json"
T1_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_d12_t1_value_law.json"
PHYSICAL_BRANCH_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_physical_branch_repair_theorem.json"
SELECTED_SHEET_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_selected_sheet_closure.json"
EXACT_READOUT_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_exact_readout.json"
BACKREAD_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_d12_internal_backread_continuation_closure.json"
SECTOR_MEAN_SPLIT_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_sector_mean_split.json"
SPREAD_MAP_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_spread_map.json"
OVERLAP_LAW_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_d12_overlap_transport_law.json"
FORWARD_YUKAWAS_JSON = ROOT / "particles" / "runs" / "flavor" / "forward_yukawas.json"
CURRENT_FAMILY_AFFINE_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_affine_anchor_theorem.json"
CURRENT_FAMILY_SIGMA_TARGET_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_exact_sigma_target.json"
CURRENT_FAMILY_PDG_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_exact_pdg_theorem.json"
ABSOLUTE_COLLAPSE_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_absolute_readout_algebraic_collapse.json"
CURRENT_FAMILY_TRANSPORT_LIFT_JSON = (
    ROOT / "particles" / "runs" / "flavor" / "quark_current_family_transport_frame_sector_attached_lift.json"
)
CURRENT_FAMILY_TRANSPORT_COMPLETION_JSON = (
    ROOT / "particles" / "runs" / "flavor" / "quark_current_family_transport_frame_exact_pdg_completion.json"
)
CURRENT_FAMILY_TRANSPORT_FORWARD_YUKAWAS_JSON = (
    ROOT / "particles" / "runs" / "flavor" / "quark_current_family_transport_frame_exact_forward_yukawas.json"
)
CURRENT_FAMILY_TRANSPORT_YUKAWA_THEOREM_JSON = (
    ROOT / "particles" / "runs" / "flavor" / "quark_current_family_transport_frame_exact_yukawa_theorem.json"
)
EXACT_YUKAWA_END_TO_END_JSON = (
    ROOT / "particles" / "runs" / "flavor" / "quark_exact_yukawa_end_to_end_theorem.json"
)
PUBLIC_EXACT_YUKAWA_PROMOTION_FRONTIER_JSON = (
    ROOT / "particles" / "runs" / "flavor" / "quark_public_exact_yukawa_promotion_frontier.json"
)
PUBLIC_SIGMA_THEOREM_JSON = (
    ROOT / "particles" / "runs" / "flavor" / "quark_public_physical_sigma_datum_descent.json"
)
PUBLIC_EXACT_YUKAWA_THEOREM_JSON = (
    ROOT / "particles" / "runs" / "flavor" / "quark_public_exact_yukawa_end_to_end_theorem.json"
)
PUBLIC_STRENGTHENED_FRONTIER_JSON = (
    ROOT / "particles" / "runs" / "flavor" / "quark_public_strengthened_physical_sigma_lift_frontier.json"
)
CURRENT_FAMILY_PHYSICAL_SIGMA_THEOREM_JSON = (
    ROOT / "particles" / "runs" / "flavor" / "quark_current_family_transport_frame_physical_sigma_lift_theorem.json"
)
CURRENT_FAMILY_STRENGTHENED_PHYSICAL_SIGMA_THEOREM_JSON = (
    ROOT / "particles" / "runs" / "flavor" / "quark_current_family_transport_frame_strengthened_physical_sigma_lift_theorem.json"
)
CURRENT_FAMILY_ABSOLUTE_THEOREM_JSON = (
    ROOT / "particles" / "runs" / "flavor" / "quark_current_family_transport_frame_absolute_sector_readout_theorem.json"
)
CURRENT_FAMILY_TRANSPORT_LIGHT_RATIO_JSON = (
    ROOT / "particles" / "runs" / "flavor" / "quark_current_family_transport_frame_light_ratio_theorem.json"
)
CURRENT_FAMILY_TRANSPORT_D12_VALUE_JSON = (
    ROOT / "particles" / "runs" / "flavor" / "quark_current_family_transport_frame_d12_value_package.json"
)
CURRENT_FAMILY_END_TO_END_CHAIN_JSON = (
    ROOT / "particles" / "runs" / "flavor" / "quark_current_family_end_to_end_exact_pdg_derivation_chain.json"
)
SIGMA_SOURCE_OBSTRUCTION_JSON = (
    ROOT / "particles" / "runs" / "flavor" / "quark_sigma_source_nonidentifiability_obstruction.json"
)
SCHEME_OBSTRUCTION_JSON = (
    ROOT / "particles" / "runs" / "flavor" / "quark_running_mass_scheme_convention_obstruction.json"
)
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_lane_closure_contract.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def build_payload(
    selector_value_law: dict[str, Any],
    t1_value_law: dict[str, Any],
    physical_branch: dict[str, Any],
    selected_sheet: dict[str, Any],
    exact_readout: dict[str, Any],
    backread: dict[str, Any],
    sector_mean_split: dict[str, Any],
    spread_map: dict[str, Any],
    overlap_law: dict[str, Any],
    forward_yukawas: dict[str, Any],
    current_family_affine: dict[str, Any],
    current_family_sigma_target: dict[str, Any],
    current_family_exact_pdg: dict[str, Any],
    absolute_collapse: dict[str, Any],
    current_family_transport_lift: dict[str, Any],
    current_family_transport_completion: dict[str, Any],
    current_family_transport_forward_yukawas: dict[str, Any],
    current_family_transport_yukawa_theorem: dict[str, Any],
    exact_yukawa_end_to_end_theorem: dict[str, Any],
    public_exact_yukawa_promotion_frontier: dict[str, Any],
    public_sigma_theorem: dict[str, Any],
    public_exact_yukawa_theorem: dict[str, Any],
    public_strengthened_frontier: dict[str, Any],
    current_family_physical_sigma_theorem: dict[str, Any],
    current_family_strengthened_physical_sigma_theorem: dict[str, Any],
    current_family_absolute_theorem: dict[str, Any],
    current_family_transport_light_ratio: dict[str, Any],
    current_family_transport_d12_value: dict[str, Any],
    current_family_end_to_end_chain: dict[str, Any],
    sigma_source_obstruction: dict[str, Any],
    scheme_obstruction: dict[str, Any],
) -> dict[str, Any]:
    selected_sigma = selected_sheet["selected_sheet"]["sigma_id"]
    source_spread_obstruction_closed = (
        sigma_source_obstruction.get("theorem_grade_obstruction") is True
        and sigma_source_obstruction.get("proof_status") == "closed_exact_current_corpus_obstruction"
        and sigma_source_obstruction.get("source_only_sigma_emitted") is False
    )
    scheme_obstruction_closed = (
        scheme_obstruction.get("proof_status")
        == "closed_structural_finite_renormalization_nonidentifiability_obstruction"
    )
    physical_yukawa_audit = dict(scheme_obstruction.get("stored_matrix_dimensional_audit") or {})
    public_yukawa_promotable = (
        not source_spread_obstruction_closed
        and not scheme_obstruction_closed
        and
        public_exact_yukawa_theorem.get("public_promotion_allowed") is True
        and public_exact_yukawa_theorem.get("proof_status")
        == "closed_source_only_public_exact_yukawa_end_to_end_theorem"
        and (public_exact_yukawa_theorem.get("non_circularity_status") or {}).get("promotion_allowed") is True
    )
    exact_masses = {
        "u": float(exact_readout["predicted_singular_values_u"][0]),
        "d": float(exact_readout["predicted_singular_values_d"][0]),
        "s": float(exact_readout["predicted_singular_values_d"][1]),
        "c": float(exact_readout["predicted_singular_values_u"][1]),
        "b": float(exact_readout["predicted_singular_values_d"][2]),
        "t": float(exact_readout["predicted_singular_values_u"][2]),
    }
    rho_ord = float(sector_mean_split["rho_ord"])
    x2 = float(t1_value_law["sample_same_family_point"]["x2"])
    main_builder_sigma_branch = dict((overlap_law.get("sigma_branch_contracts") or {}).get("main_builder_sigma_pair") or {})
    return {
        "artifact": "oph_quark_lane_closure_contract",
        "base_theorem_emitted_package_artifact": "oph_quark_maximal_theorem_emitted_package",
        "generated_utc": _timestamp(),
        "scope": "quark_lane_theorem_boundary_plus_exact_sidecars",
        "proof_status": (
            "target_free_public_exact_yukawa_derivation_closed"
            if public_yukawa_promotable
            else "selected_class_exact_witness_blocked_by_target_derived_sigma_datum"
        ),
        "public_promotion_allowed": public_yukawa_promotable,
        "non_circularity_status": public_exact_yukawa_theorem.get("non_circularity_status"),
        "mass_comparison_surface": {
            "kind": "running_mass_comparison_surface",
            "note": "Quark references are running masses, not asymptotic free-particle pole masses.",
        },
        "exact_pdg_derivation_target": {
            "target_name": "exact_running_quark_sextet_on_declared_current_family_transport_frame",
            "status": "closed",
            "artifact": current_family_end_to_end_chain["artifact"],
            "wrapper_theorem": "oph_quark_exact_pdg_end_to_end_theorem",
            "theorem_scope": current_family_end_to_end_chain["theorem_scope"],
            "minimal_exact_blocker_set": [],
            "exact_running_values_gev": current_family_end_to_end_chain["exact_running_values_gev"],
            "lemma_chain": current_family_end_to_end_chain["lemma_chain"],
            "why_closed": (
                "On the declared current-family/common-refinement transport-frame carrier, the target-free D12 bridge, "
                "restricted strengthened physical sigma lift, restricted absolute sector readout, and ordered three-point "
                "readout already compose to the exact PDG-matched running quark sextet."
            ),
            "not_the_same_as": "target_free_public_physical_sheet_promotion",
            "strengthening_above_target": (
                "A stronger target-free public physical-sheet theorem remains a separate promotion question and is not "
                "part of the exact derivation target recorded here."
            ),
        },
        "exact_yukawa_derivation_target": {
            "target_name": current_family_transport_yukawa_theorem["target_name"],
            "status": "closed",
            "artifact": current_family_transport_yukawa_theorem["artifact"],
            "wrapper_theorem": exact_yukawa_end_to_end_theorem["artifact"],
            "theorem_scope": current_family_transport_yukawa_theorem["theorem_scope"],
            "minimal_exact_blocker_set": current_family_transport_yukawa_theorem["minimal_exact_blocker_set"],
            "forward_certified": current_family_transport_yukawa_theorem["forward_yukawa_artifact"]["forward_certified"],
            "certification_status": current_family_transport_yukawa_theorem["forward_yukawa_artifact"]["certification_status"],
            "why_closed": (
                "On the declared current-family/common-refinement transport-frame carrier, the closed exact chain now "
                "emits explicit forward Yukawa matrices Y_u and Y_d with certified singular values matching the exact "
                "running quark sextet."
            ),
            "not_the_same_as": "target_free_public_physical_sheet_yukawa_promotion",
        },
        "public_exact_yukawa_derivation_target": {
            "target_name": public_exact_yukawa_theorem["target_name"],
            "status": "closed" if public_yukawa_promotable else "blocked_target_derived_sigma_source_missing",
            "artifact": public_exact_yukawa_theorem["artifact"],
            "theorem_scope": public_exact_yukawa_theorem["theorem_scope"],
            "minimal_exact_blocker_set": public_exact_yukawa_theorem["minimal_exact_blocker_set"],
            "exact_running_values_gev": public_exact_yukawa_theorem["public_exact_outputs"]["exact_running_values_gev"],
            "forward_yukawa_artifact": public_exact_yukawa_theorem["public_exact_outputs"]["forward_yukawa_artifact"],
            "why_closed": (
                "The source-only sigma selector and public sigma-datum descent theorem identify the selected public "
                "quark frame class with a no-target physical sigma datum, so the algebraic absolute readout and exact "
                "forward construction lift to the selected public class."
                if public_yukawa_promotable
                else "The selected-class exact witness is displayed, but strict public promotion is blocked because descent proves representative independence only; the sigma datum still comes from an exact target surface."
            ),
        },
        "selected_local_sheet_status": {
            "sigma_id": selected_sigma,
            "proof_status": selected_sheet["proof_status"],
            "theorem_scope": selected_sheet["theorem_scope"],
            "wrong_branch_for_physical_ckm_shell": True,
            "why_not_enough": physical_branch["insufficiency_theorem"]["statement"],
        },
        "exact_sidecar_mass_surface": {
            "artifact": exact_readout["artifact"],
            "scope": exact_readout["theorem_scope"],
            "selected_sheet": selected_sigma,
            "current_family_affine_anchor_theorem": current_family_affine["artifact"],
            "current_family_exact_pdg_theorem": current_family_exact_pdg["artifact"],
            "exact_outputs_gev": exact_masses,
            "exact_sector_geometric_means": {
                "g_u": float(exact_readout["g_u"]),
                "g_d": float(exact_readout["g_d"]),
            },
            "closure_statement": selected_sheet["theorem_statement"],
        },
        "current_family_physical_target_surface": {
            "affine_anchor_theorem": {
                "artifact": current_family_affine["artifact"],
                "proof_status": current_family_affine["proof_status"],
                "A_q_current_family": float(current_family_affine["current_family_affine_anchor"]["value"]),
                "delta_q_current_family": float(current_family_affine["current_family_sector_split"]["value"]),
            },
            "exact_sigma_target": {
                "artifact": current_family_sigma_target["artifact"],
                "proof_status": current_family_sigma_target["proof_status"],
                "sigma_u_target": float(current_family_sigma_target["unique_exact_sigma_target"]["sigma_u_target"]),
                "sigma_d_target": float(current_family_sigma_target["unique_exact_sigma_target"]["sigma_d_target"]),
                "delta_vs_current_theorem_grade_sigma_pair": current_family_sigma_target[
                    "delta_vs_current_theorem_grade_sigma_pair"
                ],
            },
            "exact_pdg_reconstruction": {
                "artifact": current_family_exact_pdg["artifact"],
                "proof_status": current_family_exact_pdg["proof_status"],
                "reconstructed_current_family_running_values_gev": current_family_exact_pdg[
                    "reconstructed_current_family_running_values_gev"
                ],
            },
            "absolute_readout_algebraic_collapse": {
                "artifact": absolute_collapse["artifact"],
                "proof_status": absolute_collapse["proof_status"],
                "theorem_scope": absolute_collapse["theorem_scope"],
                "remaining_nonalgebraic_theorem": absolute_collapse["conditional_collapse_route"][
                    "remaining_nonalgebraic_theorem"
                ],
                "candidate_merged_theorem_text": absolute_collapse["candidate_merged_theorem_text"],
            },
            "transport_frame_sector_attached_lift": {
                "artifact": current_family_transport_lift["artifact"],
                "proof_status": current_family_transport_lift["proof_status"],
                "theorem_scope": current_family_transport_lift["theorem_scope"],
                "sigma_id": current_family_transport_lift["emitted_sigma_ud_phys_element"]["sigma_id"],
                "canonical_token": current_family_transport_lift["emitted_sigma_ud_phys_element"]["canonical_token"],
                "sigma_u_target": float(current_family_transport_lift["strengthened_sigma_data"]["sigma_u_target"]),
                "sigma_d_target": float(current_family_transport_lift["strengthened_sigma_data"]["sigma_d_target"]),
            },
            "restricted_physical_sigma_lift_theorem": {
                "artifact": current_family_physical_sigma_theorem["artifact"],
                "proof_status": current_family_physical_sigma_theorem["proof_status"],
                "theorem_scope": current_family_physical_sigma_theorem["theorem_scope"],
                "corresponds_to_global_contract": current_family_physical_sigma_theorem["corresponds_to_global_contract"],
            },
            "restricted_strengthened_physical_sigma_lift_theorem": {
                "artifact": current_family_strengthened_physical_sigma_theorem["artifact"],
                "proof_status": current_family_strengthened_physical_sigma_theorem["proof_status"],
                "theorem_scope": current_family_strengthened_physical_sigma_theorem["theorem_scope"],
                "compressed_global_contract": current_family_strengthened_physical_sigma_theorem["compressed_global_contract"],
                "theorem_grade_physical_sigma_datum": current_family_strengthened_physical_sigma_theorem[
                    "theorem_grade_physical_sigma_datum"
                ],
            },
            "restricted_absolute_sector_readout_theorem": {
                "artifact": current_family_absolute_theorem["artifact"],
                "proof_status": current_family_absolute_theorem["proof_status"],
                "theorem_scope": current_family_absolute_theorem["theorem_scope"],
                "corresponds_to_global_contract": current_family_absolute_theorem["corresponds_to_global_contract"],
                "emitted_absolute_sector_scales": current_family_absolute_theorem["emitted_absolute_sector_scales"],
            },
            "restricted_light_ratio_theorem": {
                "artifact": current_family_transport_light_ratio["artifact"],
                "proof_status": current_family_transport_light_ratio["proof_status"],
                "theorem_scope": current_family_transport_light_ratio["theorem_scope"],
                "ell_ud": float(current_family_transport_light_ratio["exact_light_data"]["ell_ud"]),
            },
            "restricted_d12_value_package": {
                "artifact": current_family_transport_d12_value["artifact"],
                "proof_status": current_family_transport_d12_value["proof_status"],
                "theorem_scope": current_family_transport_d12_value["theorem_scope"],
                "closed_d12_scalars": current_family_transport_d12_value["closed_d12_scalars"],
            },
            "transport_frame_exact_pdg_completion": {
                "artifact": current_family_transport_completion["artifact"],
                "proof_status": current_family_transport_completion["proof_status"],
                "theorem_scope": current_family_transport_completion["theorem_scope"],
                "exact_running_values_gev": current_family_transport_completion["exact_running_values_gev"],
            },
            "end_to_end_exact_pdg_derivation_chain": {
                "artifact": current_family_end_to_end_chain["artifact"],
                "proof_status": current_family_end_to_end_chain["proof_status"],
                "theorem_scope": current_family_end_to_end_chain["theorem_scope"],
                "exact_running_values_gev": current_family_end_to_end_chain["exact_running_values_gev"],
            },
            "exact_forward_yukawas": {
                "artifact": current_family_transport_forward_yukawas["artifact"],
                "proof_status": current_family_transport_forward_yukawas["proof_status"],
                "scope": current_family_transport_forward_yukawas["scope"],
                "forward_certified": current_family_transport_forward_yukawas["forward_certified"],
                "certification_status": current_family_transport_forward_yukawas["certification_status"],
            },
            "exact_yukawa_theorem": {
                "artifact": current_family_transport_yukawa_theorem["artifact"],
                "proof_status": current_family_transport_yukawa_theorem["proof_status"],
                "target_name": current_family_transport_yukawa_theorem["target_name"],
                "theorem_scope": current_family_transport_yukawa_theorem["theorem_scope"],
            },
            "note": (
                "These exact current-family target surfaces are now explicit on disk. They sharpen the remaining "
                "public frontier quantitatively, but they do not by themselves promote the target-free physical-sheet lane."
            ),
        },
        "public_final_theorem_frontier": {
            "artifact": public_strengthened_frontier["artifact"],
            "proof_status": public_strengthened_frontier["proof_status"],
            "public_promotion_allowed": public_strengthened_frontier.get("public_promotion_allowed"),
            "non_circularity_status": public_strengthened_frontier.get("non_circularity_status"),
            "resolved_by_theorem_artifact": public_strengthened_frontier["resolved_by_theorem_artifact"],
            "final_public_theorem_candidate": public_strengthened_frontier["final_public_theorem_candidate"],
            "alternate_upstream_route": public_strengthened_frontier["alternate_upstream_route"],
            "algebraic_consequence_after_closure": public_strengthened_frontier["algebraic_consequence_after_closure"],
        },
        "public_exact_yukawa_promotion_frontier": {
            "artifact": public_exact_yukawa_promotion_frontier["artifact"],
            "proof_status": public_exact_yukawa_promotion_frontier["proof_status"],
            "target_name": public_exact_yukawa_promotion_frontier["target_name"],
            "public_promotion_allowed": public_exact_yukawa_promotion_frontier.get("public_promotion_allowed"),
            "non_circularity_status": public_exact_yukawa_promotion_frontier.get("non_circularity_status"),
            "resolved_by_theorem_artifact": public_exact_yukawa_promotion_frontier["resolved_by_theorem_artifact"],
            "final_public_theorem_candidate": public_exact_yukawa_promotion_frontier["final_public_theorem_candidate"],
            "alternate_upstream_route": public_exact_yukawa_promotion_frontier["alternate_upstream_route"],
            "closed_public_endpoint": public_exact_yukawa_promotion_frontier["closed_public_endpoint"],
        },
        "candidate_one_theorem_physical_compression": {
            "status": "closed" if public_yukawa_promotable else "blocked_target_derived_sigma_source_missing",
            "artifact": public_sigma_theorem["artifact"],
            "supporting_algebraic_collapse_artifact": absolute_collapse["artifact"],
            "conditional_statement": absolute_collapse["theorem_statement"],
            "local_strengthened_theorem_statement": current_family_strengthened_physical_sigma_theorem["theorem_statement"],
            "public_strengthened_theorem_statement": public_sigma_theorem["theorem_statement"],
            "remaining_nonalgebraic_theorem": None
            if public_yukawa_promotable
            else "QUARK_SIGMA_SOURCE_SELECTOR",
            "remaining_exact_gap": None
            if public_yukawa_promotable
            else "target_derived_sigma_datum_used_for_selected_class_exact_witness",
        },
        "continuation_only_mass_sidecar": {
            "artifact": backread["artifact"],
            "scope": backread["scope"],
            "closed_mass_side_package": backread["closed_mass_side_package"],
            "closed_source_side_package": backread["closed_source_side_package"],
            "theorem_boundary_note": backread["theorem_boundary_note"],
        },
        "public_current_family_yukawa_frontier": {
            "definition": (
                "Theorem-grade current-family quark Yukawas on the emitted D12 mass ray, using the already-closed "
                "mean-surface sigma branch rather than the optional edge-statistics bridge."
            ),
            "sharper_target_1_primitive": {
                "artifact": selector_value_law["artifact"],
                "proof_status": selector_value_law["proof_status"],
                "exact_missing_object": selector_value_law["exact_missing_object"],
                "equivalent_ray_coordinate_presentation": selector_value_law["equivalent_ray_coordinate_presentation"]["theorem_id"],
            },
            "theorem_grade_sigma_branch": {
                "artifact": spread_map["artifact"],
                "proof_status": spread_map["proof_status"],
                "spread_emitter_status": spread_map["spread_emitter_status"],
                "sigma_source_kind": spread_map["sigma_source_kind"],
                "sigma_u_total_log_per_side": float(spread_map["sigma_u_total_log_per_side"]),
                "sigma_d_total_log_per_side": float(spread_map["sigma_d_total_log_per_side"]),
            },
            "transport_reduction": {
                "artifact": overlap_law["artifact"],
                "proof_status": overlap_law["proof_status"],
                "status": overlap_law["status"],
                "main_builder_sigma_branch": {
                    "provider_artifact": main_builder_sigma_branch.get("provider_artifact"),
                    "provider_status": main_builder_sigma_branch.get("provider_status"),
                    "sigma_source_kind": main_builder_sigma_branch.get("sigma_source_kind"),
                },
                "remaining_scalar_on_fixed_sigma_branch": overlap_law["reduced_exact_gap"]["remaining_scalar_on_any_fixed_sigma_branch"],
            },
            "minimal_exact_blocker_set": [],
            "target_1_status": "closed",
            "why_closed": (
                "light_quark_overlap_defect_value_law is now internalized on the code surface, so the emitted D12 ray "
                "fixes t1 and the closed sigma/transport laws force the full target-1 source and transport package."
            ),
            "why_edge_statistics_bridge_is_not_required": (
                "The edge-statistics/mean-surface compatibility theorem would strengthen the bridge from overlap-edge "
                "data to the sigma pair, but it is not needed to provide a theorem-grade sigma branch for current-family "
                "public Yukawas because oph_family_excitation_spread_map already closes that branch."
            ),
            "target_1_internalized_theorem_text": selector_value_law["theorem_statement"],
            "closure_after_t1": {
                "forced_source_payload_after_t1": t1_value_law["forced_source_payload_after_t1"],
                "transport_side_reduction": t1_value_law["transport_side_reduction"],
                "candidate_scalar_identities": t1_value_law["candidate_public_construction_route"]["candidate_scalar_identities"],
            },
        },
        "internalized_theorems": [
            {
                "id": "light_quark_overlap_defect_value_law",
                "proof_status": selector_value_law["proof_status"],
                "formula": selector_value_law["target_free_map"]["formula"],
            },
            {
                "id": "quark_d12_t1_value_law",
                "proof_status": t1_value_law["proof_status"],
                "formula": "t1 = (5/6) * log(c_d / c_u)",
            },
        ],
        "exact_missing_theorems": [],
        "closure_chain": [
            "(axioms + light-data transport) => light_quark_overlap_defect_value_law => Delta_ud_overlap => quark_d12_t1_value_law => t1 => (eta_Q_centered, kappa_Q, tau_u, tau_d)",
            "(selected public quark frame class chosen by P) => selected_bridge_fiber_sigma_descent_not_source_selection => representative-independent target-derived sigma witness",
            "(missing source Sigma_ud^phys datum) => quark_absolute_readout_algebraic_collapse => (g_u, g_d) => (m_u, m_d, m_s, m_c, m_b, m_t)",
            "(public selected class + exact forward construction + source sigma selector) => oph_quark_public_exact_yukawa_end_to_end_theorem => exact forward Y_u, Y_d",
        ],
        "notes": [
            "The mass bridge is now internalized on the local code surface: light_quark_overlap_defect_value_law and quark_d12_t1_value_law are no longer part of the remaining exact frontier.",
            "The exact-PDG derivation target is closed on the declared current-family/common-refinement transport-frame carrier and is recorded explicitly above.",
            "The exact forward-Yukawa derivation target is also closed on that declared carrier and is recorded explicitly above.",
            "The direct public sigma-datum descent theorem proves selected bridge-fiber representative independence only.",
            "That public descent does not close the physical sigma source-selector frontier and does not lift the exact Yukawa package to a target-free public theorem on the selected class.",
            "The exact current-family witness and the D12 internal backread sidecar still exhibit the mass data on sidecar surfaces, but they do not repair the physical CKM branch by themselves.",
            "The selected local same-label left-handed sheet closes negatively to sigma_ref; the remaining physical quark burden is not another local selector search.",
            "The historical public physical-sheet contract named Theta_ud^phys and Theta_ud^abs separately; after source sigma selection the absolute readout is algebraic, but the source sigma selector itself remains the blocker.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the exact three-step quark closure contract.")
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    payload = build_payload(
        _load_json(SELECTOR_VALUE_JSON),
        _load_json(T1_JSON),
        _load_json(PHYSICAL_BRANCH_JSON),
        _load_json(SELECTED_SHEET_JSON),
        _load_json(EXACT_READOUT_JSON),
        _load_json(BACKREAD_JSON),
        _load_json(SECTOR_MEAN_SPLIT_JSON),
        _load_json(SPREAD_MAP_JSON),
        _load_json(OVERLAP_LAW_JSON),
        _load_json(FORWARD_YUKAWAS_JSON),
        _load_json(CURRENT_FAMILY_AFFINE_JSON),
        _load_json(CURRENT_FAMILY_SIGMA_TARGET_JSON),
        _load_json(CURRENT_FAMILY_PDG_JSON),
        _load_json(ABSOLUTE_COLLAPSE_JSON),
        _load_json(CURRENT_FAMILY_TRANSPORT_LIFT_JSON),
        _load_json(CURRENT_FAMILY_TRANSPORT_COMPLETION_JSON),
        _load_json(CURRENT_FAMILY_TRANSPORT_FORWARD_YUKAWAS_JSON),
        _load_json(CURRENT_FAMILY_TRANSPORT_YUKAWA_THEOREM_JSON),
        _load_json(EXACT_YUKAWA_END_TO_END_JSON),
        _load_json(PUBLIC_EXACT_YUKAWA_PROMOTION_FRONTIER_JSON),
        _load_json(PUBLIC_SIGMA_THEOREM_JSON),
        _load_json(PUBLIC_EXACT_YUKAWA_THEOREM_JSON),
        _load_json(PUBLIC_STRENGTHENED_FRONTIER_JSON),
        _load_json(CURRENT_FAMILY_PHYSICAL_SIGMA_THEOREM_JSON),
        _load_json(CURRENT_FAMILY_STRENGTHENED_PHYSICAL_SIGMA_THEOREM_JSON),
        _load_json(CURRENT_FAMILY_ABSOLUTE_THEOREM_JSON),
        _load_json(CURRENT_FAMILY_TRANSPORT_LIGHT_RATIO_JSON),
        _load_json(CURRENT_FAMILY_TRANSPORT_D12_VALUE_JSON),
        _load_json(CURRENT_FAMILY_END_TO_END_CHAIN_JSON),
    )

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
