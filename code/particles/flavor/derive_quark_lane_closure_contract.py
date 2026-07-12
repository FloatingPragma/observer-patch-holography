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
AXIOM_LEVEL_YUKAWA_OBSTRUCTION_JSON = (
    ROOT
    / "particles"
    / "runs"
    / "flavor"
    / "quark_axiom_level_yukawa_moduli_nonidentifiability.json"
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
    axiom_level_yukawa_obstruction: dict[str, Any],
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
        and public_exact_yukawa_theorem.get("public_promotion_allowed") is True
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
    return {
        "artifact": "oph_quark_lane_closure_contract",
        "base_theorem_emitted_package_artifact": "oph_quark_maximal_theorem_emitted_package",
        "generated_utc": _timestamp(),
        "scope": "quark_lane_source_obstructions_plus_target_anchored_audit_sidecars",
        "proof_status": (
            "target_free_public_exact_yukawa_derivation_closed"
            if public_yukawa_promotable
            else "closed_sharper_obstructions_numeric_quark_predictions_withheld"
        ),
        "public_promotion_allowed": public_yukawa_promotable,
        "numeric_quark_prediction_rows_allowed": public_yukawa_promotable,
        "non_circularity_status": {
            **dict(public_exact_yukawa_theorem.get("non_circularity_status") or {}),
            "source_spread_fiber": sigma_source_obstruction["exact_ray_classification"]["fiber"],
            "source_spread_fiber_dimension": sigma_source_obstruction["exact_ray_classification"][
                "fiber_dimension"
            ],
            "source_spread_obstruction_closed": source_spread_obstruction_closed,
            "scheme_coordinate_obstruction_closed": scheme_obstruction_closed,
            "physical_yukawa_normalization_closed": physical_yukawa_audit.get(
                "certified_physical_yukawa_matrices"
            )
            is True,
        },
        "mass_comparison_surface": {
            "kind": "mixed_declared_coordinate_target_audit",
            "row_partition": scheme_obstruction["row_partition"],
            "single_running_quark_sextet_claim_allowed": False,
            "note": (
                "The light rows, heavy self-scale rows, and separate top extraction coordinate do not form one "
                "source-emitted running-mass sextet."
            ),
        },
        "exact_pdg_derivation_target": {
            "target_name": current_family_end_to_end_chain["target_name"],
            "status": "closed_target_anchored_algebraic_audit_only",
            "artifact": current_family_end_to_end_chain["artifact"],
            "wrapper_theorem": "oph_quark_exact_pdg_end_to_end_theorem",
            "theorem_scope": current_family_end_to_end_chain["theorem_scope"],
            "minimal_exact_blocker_set": [],
            "target_anchored": True,
            "source_only_prediction": False,
            "single_running_quark_sextet_claim_allowed": False,
            "comparison_coordinate_partition": current_family_end_to_end_chain.get(
                "comparison_coordinate_partition", scheme_obstruction["row_partition"]
            ),
            "exact_running_values_gev": current_family_end_to_end_chain["exact_running_values_gev"],
            "lemma_chain": current_family_end_to_end_chain["lemma_chain"],
            "audit_scope_statement": (
                "The declared current-family carrier reconstructs its chosen target packet exactly after the "
                "target-attached sigma datum is supplied. This certifies the conditional algebra, not source emission."
            ),
            "not_the_same_as": "target_free_public_physical_sheet_promotion",
            "strengthening_above_target": (
                "A source map must first break the independent positive spread action and emit an RG-covariant "
                "trajectory or invariant before any declared comparison chart is applied."
            ),
        },
        "exact_yukawa_derivation_target": {
            "target_name": "mixed_scheme_dimensionful_quark_mass_texture_audit",
            "status": "closed_algebraic_mass_texture_audit_not_physical_yukawa",
            "artifact": current_family_transport_yukawa_theorem["artifact"],
            "wrapper_theorem": exact_yukawa_end_to_end_theorem["artifact"],
            "theorem_scope": current_family_transport_yukawa_theorem["theorem_scope"],
            "minimal_exact_blocker_set": current_family_transport_yukawa_theorem["minimal_exact_blocker_set"],
            "forward_certified": current_family_transport_yukawa_theorem["forward_yukawa_artifact"]["forward_certified"],
            "certification_status": current_family_transport_yukawa_theorem["forward_yukawa_artifact"]["certification_status"],
            "matrix_kind": physical_yukawa_audit["current_classification"],
            "stored_entry_dimension": physical_yukawa_audit["stored_entry_dimension"],
            "certified_physical_yukawa_matrices": False,
            "audit_scope_statement": (
                "The stored matrices reproduce the target packet as GeV-valued singular spectra. They are mass "
                "textures, not common-scale dimensionless physical Yukawa matrices."
            ),
            "not_the_same_as": "target_free_public_physical_sheet_yukawa_promotion",
        },
        "public_exact_yukawa_derivation_target": {
            "target_name": public_exact_yukawa_theorem["target_name"],
            "status": "closed" if public_yukawa_promotable else public_exact_yukawa_theorem["proof_status"],
            "artifact": public_exact_yukawa_theorem["artifact"],
            "theorem_scope": public_exact_yukawa_theorem["theorem_scope"],
            "minimal_exact_blocker_set": public_exact_yukawa_theorem["minimal_exact_blocker_set"],
            "exact_running_values_gev": public_exact_yukawa_theorem["public_exact_outputs"]["exact_running_values_gev"],
            "forward_yukawa_artifact": public_exact_yukawa_theorem["public_exact_outputs"]["forward_yukawa_artifact"],
            "numeric_values_role": "target_anchored_audit_only",
            "physical_yukawa_claim_allowed": False,
            "status_statement": (
                "A future source theorem breaks the two-modulus action and emits a representative-independent spread "
                "pair on the selected class. After common-scale RG transport and Higgs normalization, the physical "
                "dimensionless Yukawa construction then closes."
                if public_yukawa_promotable
                else (
                    "The source spread fiber retains two independent positive moduli. In addition, the selected-class "
                    "numerical matrices are mixed-scheme GeV mass textures and lack a common-scale dimensionless "
                    "Yukawa conversion."
                )
            ),
        },
        "selected_local_sheet_status": {
            "sigma_id": selected_sigma,
            "proof_status": selected_sheet["proof_status"],
            "theorem_scope": selected_sheet["theorem_scope"],
            "wrong_branch_for_physical_ckm_shell": True,
            "why_not_enough": physical_branch["insufficiency_theorem"]["statement"],
        },
        "source_spread_nonidentifiability_obstruction": {
            "artifact": sigma_source_obstruction["artifact"],
            "proof_status": sigma_source_obstruction["proof_status"],
            "claim_tier": sigma_source_obstruction["claim_tier"],
            "theorem_grade_obstruction": sigma_source_obstruction["theorem_grade_obstruction"],
            "theorem_statement": sigma_source_obstruction["theorem_statement"],
            "compatible_source_spread_fiber": sigma_source_obstruction["exact_ray_classification"]["fiber"],
            "fiber_dimension": sigma_source_obstruction["exact_ray_classification"]["fiber_dimension"],
            "independent_coordinates": sigma_source_obstruction["exact_ray_classification"][
                "independent_coordinates"
            ],
            "free_action": sigma_source_obstruction["free_action_certificate"]["action"],
            "source_only_sigma_emitted": sigma_source_obstruction["source_only_sigma_emitted"],
            "numeric_quark_rows_allowed": sigma_source_obstruction["numeric_quark_rows_allowed"],
            "github_issues": sigma_source_obstruction["github_issues"],
            "dependency_audit": sigma_source_obstruction["dependency_audit"],
            "minimal_future_extension": sigma_source_obstruction["minimal_future_extension"],
        },
        "axiom_level_yukawa_moduli_nonidentifiability": {
            "artifact": axiom_level_yukawa_obstruction["artifact"],
            "proof_status": axiom_level_yukawa_obstruction["proof_status"],
            "theorem_statement": axiom_level_yukawa_obstruction["theorem_statement"],
            "additional_axioms_used": axiom_level_yukawa_obstruction["additional_axioms_used"],
            "counterfamily": axiom_level_yukawa_obstruction["counterfamily"],
            "MAR_audit": axiom_level_yukawa_obstruction["axiom_invariance_audit"]["Axiom_5_MAR"],
            "reference_data_policy": axiom_level_yukawa_obstruction["reference_data_policy"],
            "public_numeric_quark_rows_allowed": axiom_level_yukawa_obstruction[
                "public_numeric_quark_rows_allowed"
            ],
        },
        "running_scheme_and_physical_yukawa_obstruction": {
            "artifact": scheme_obstruction["artifact"],
            "proof_status": scheme_obstruction["proof_status"],
            "closure_kind": scheme_obstruction["closure_kind"],
            "theorem_statement": scheme_obstruction["theorem_statement"],
            "github_issues": scheme_obstruction["github_issues"],
            "reference_data_policy": scheme_obstruction["reference_data_policy"],
            "row_partition": scheme_obstruction["row_partition"],
            "formal_quotient_obstruction": scheme_obstruction["formal_quotient_obstruction"],
            "rg_covariant_source_output": scheme_obstruction["rg_covariant_source_output"],
            "stored_matrix_dimensional_audit": physical_yukawa_audit,
            "public_numeric_quark_rows_allowed": scheme_obstruction["closure_effect"][
                "public_numeric_quark_rows_allowed"
            ],
        },
        "exact_sidecar_mass_surface": {
            "artifact": exact_readout["artifact"],
            "scope": exact_readout["theorem_scope"],
            "selected_sheet": selected_sigma,
            "role": "target_anchored_current_family_audit",
            "source_only_prediction": False,
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
                "These values and matrices are retained as explicit target-anchored audit surfaces. They neither "
                "select a point of the source spread fiber nor define a single common-scale physical Yukawa packet."
            ),
        },
        "public_final_theorem_frontier": {
            "artifact": public_strengthened_frontier["artifact"],
            "proof_status": sigma_source_obstruction["proof_status"],
            "closure_kind": "sharper_current_corpus_nonidentifiability_obstruction",
            "public_promotion_allowed": False,
            "non_circularity_status": public_strengthened_frontier.get("non_circularity_status"),
            "resolved_by_theorem_artifact": sigma_source_obstruction["artifact"],
            "legacy_selected_class_descent_artifact": public_strengthened_frontier[
                "resolved_by_theorem_artifact"
            ],
            "compatible_source_spread_fiber": sigma_source_obstruction["exact_ray_classification"]["fiber"],
            "fiber_dimension": sigma_source_obstruction["exact_ray_classification"]["fiber_dimension"],
            "alternate_upstream_route": public_strengthened_frontier["alternate_upstream_route"],
            "reopen_requirement": sigma_source_obstruction["minimal_future_extension"],
            "status_statement": (
                "Selected-class descent removes representative dependence but does not break the independent "
                "positive rescaling action on the up- and down-sector spreads."
            ),
        },
        "public_exact_yukawa_promotion_frontier": {
            "artifact": public_exact_yukawa_promotion_frontier["artifact"],
            "proof_status": public_exact_yukawa_promotion_frontier["proof_status"],
            "target_name": public_exact_yukawa_promotion_frontier["target_name"],
            "public_promotion_allowed": False,
            "non_circularity_status": public_exact_yukawa_promotion_frontier.get("non_circularity_status"),
            "resolved_by_theorem_artifact": public_exact_yukawa_promotion_frontier["resolved_by_theorem_artifact"],
            "closed_public_endpoint": public_exact_yukawa_promotion_frontier["closed_public_endpoint"],
            "source_spread_obstruction_artifact": sigma_source_obstruction["artifact"],
            "scheme_obstruction_artifact": scheme_obstruction["artifact"],
            "matrix_kind": physical_yukawa_audit["current_classification"],
            "certified_physical_yukawa_matrices": False,
        },
        "candidate_one_theorem_physical_compression": {
            "status": (
                "closed" if public_yukawa_promotable else "closed_two_modulus_nonidentifiability_obstruction"
            ),
            "artifact": sigma_source_obstruction["artifact"],
            "legacy_selected_class_descent_artifact": public_sigma_theorem["artifact"],
            "supporting_algebraic_collapse_artifact": absolute_collapse["artifact"],
            "conditional_statement": absolute_collapse["theorem_statement"],
            "current_corpus_source_spread_fiber": sigma_source_obstruction["exact_ray_classification"]["fiber"],
            "remaining_source_dimension": sigma_source_obstruction["exact_ray_classification"]["fiber_dimension"],
            "remaining_nonalgebraic_theorem": (
                None if public_yukawa_promotable else "QUARK_SOURCE_SPREAD_PAIR_ACTION_BREAKING_THEOREM"
            ),
            "remaining_exact_gap": (
                None
                if public_yukawa_promotable
                else "two_independent_positive_spread_moduli_plus_scheme_and_dimensionless_yukawa_chart"
            ),
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
                "Target-free source identities close the ordered D12 profile rays only. Absolute up/down spreads "
                "remain a two-modulus fiber, and the target-side GeV matrices remain mixed-scheme mass textures."
            ),
            "sharper_target_1_primitive": {
                "artifact": selector_value_law["artifact"],
                "proof_status": selector_value_law["proof_status"],
                "exact_missing_object": selector_value_law["exact_missing_object"],
                "equivalent_ray_coordinate_presentation": selector_value_law["equivalent_ray_coordinate_presentation"]["theorem_id"],
            },
            "theorem_grade_sigma_branch": {
                "artifact": sigma_source_obstruction["artifact"],
                "proof_status": sigma_source_obstruction["proof_status"],
                "spread_emitter_status": "source_only_spread_pair_nonidentifiable",
                "sigma_source_kind": "independent_positive_moduli_not_selected_by_current_source_corpus",
                "compatible_fiber": sigma_source_obstruction["exact_ray_classification"]["fiber"],
                "fiber_dimension": sigma_source_obstruction["exact_ray_classification"]["fiber_dimension"],
                "source_only_sigma_emitted": False,
                "target_anchored_audit_provider": spread_map["artifact"],
            },
            "transport_reduction": {
                "artifact": overlap_law["artifact"],
                "proof_status": overlap_law["proof_status"],
                "status": overlap_law["status"],
                "conditional_remaining_scalar_after_external_sigma_fixing": overlap_law["reduced_exact_gap"][
                    "remaining_scalar_on_any_fixed_sigma_branch"
                ],
                "why_not_a_source_closure": (
                    "Fixing a sigma branch is an extra premise. The target-free corpus does not select the two "
                    "positive sector spans to which the conditional reduction is applied."
                ),
            },
            "minimal_exact_blocker_set": [
                "QUARK_SOURCE_SPREAD_PAIR_ACTION_BREAKING_THEOREM",
                "QUARK_RG_COVARIANT_TRAJECTORY_OR_INVARIANT",
                "QUARK_OPERATIONAL_SCHEME_AND_SCALE_SECTION",
                "QUARK_THRESHOLD_AND_TOP_CONVERSION",
                "QUARK_COMMON_SCALE_DIMENSIONLESS_YUKAWA_CERTIFICATE",
            ],
            "target_1_status": "closed_profile_ray_only_absolute_spreads_nonidentifiable",
            "why_not_closed": (
                "The internalized overlap-defect law fixes a ray coordinate, not its two absolute sector spans. "
                "Independent positive rescaling preserves all current source identities and changes the mass readout."
            ),
            "why_edge_statistics_candidate_does_not_close_source": (
                "Even after granting the edge packet, sigma_u=S_13+c_u*delta21 and "
                "sigma_d=S_23+c_d*delta21 retain independent free coefficients c_u and c_d."
            ),
            "scheme_and_dimensional_yukawa_boundary": {
                "artifact": scheme_obstruction["artifact"],
                "row_partition": scheme_obstruction["row_partition"],
                "matrix_kind": physical_yukawa_audit["current_classification"],
                "physical_dimensionless_relation": physical_yukawa_audit["physical_dimensionless_relation"],
                "promotion_requirements": physical_yukawa_audit["promotion_requirements"],
            },
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
        "exact_missing_theorems": [
            {
                "id": "quark_source_spread_pair_action_breaking_theorem",
                "current_obstruction": sigma_source_obstruction["artifact"],
                "required_independent_scalar_count": sigma_source_obstruction["minimal_future_extension"][
                    "required_independent_scalar_count"
                ],
            },
            {
                "id": "quark_rg_covariant_scheme_readout_or_invariant",
                "current_obstruction": scheme_obstruction["artifact"],
                "missing_section": scheme_obstruction["formal_quotient_obstruction"]["missing_section"],
            },
            {
                "id": "quark_common_scale_dimensionless_yukawa_certificate",
                "current_matrix_kind": physical_yukawa_audit["current_classification"],
                "promotion_requirements": physical_yukawa_audit["promotion_requirements"],
            },
        ],
        "closure_chain": [
            "(axioms + light-data transport) => light_quark_overlap_defect_value_law => Delta_ud_overlap => quark_d12_t1_value_law => t1 => (eta_Q_centered, kappa_Q, tau_u, tau_d)",
            "(target-free source identities) => (v_u,v_d) with E_u=sigma_u*v_u and E_d=sigma_d*v_d => compatible fiber (R_{>0})^2",
            "(selected public quark frame class chosen by P) => representative-independent descent, but no selection of (sigma_u,sigma_d)",
            "(externally supplied target sigma pair) => conditional affine/readout algebra => mixed-convention six-row target packet",
            "(common-scale running coordinates + threshold/top conversion + v(mu)) => dimensionless physical Yukawa matrices; these inputs are not emitted by the current source corpus",
        ],
        "notes": [
            "The overlap-defect and D12 value laws close shape information, not the two absolute spread moduli.",
            "The theorem-grade current-corpus obstruction closes issues 377, 379, and 380 only in the accepted sharper-obstruction mode; it does not emit numerical quark rows.",
            "The target-anchored current-family algebra remains available for audit and regression testing.",
            "The declared six-row target packet mixes light MSbar coordinates at 2 GeV, charm and bottom self-scale coordinates, and a separate top extraction coordinate.",
            "The stored GeV-valued matrices are mass textures. Physical Yukawa promotion needs a common scale, scheme and threshold transport, top conversion, a running Higgs expectation value, and dimensionless normalization.",
            "The scheme obstruction closes issues 381 and 382 only in the accepted sharper-obstruction mode and likewise does not authorize numerical prediction rows.",
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
        _load_json(SIGMA_SOURCE_OBSTRUCTION_JSON),
        _load_json(AXIOM_LEVEL_YUKAWA_OBSTRUCTION_JSON),
        _load_json(SCHEME_OBSTRUCTION_JSON),
    )

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
