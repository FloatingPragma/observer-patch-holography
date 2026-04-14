#!/usr/bin/env python3
"""Shared helpers for charged absolute-side scaffold emitters."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
REPO_ROOT = ROOT.parent

RUNS_DIR = ROOT / "particles" / "runs" / "leptons"
GENERATION_BUNDLE_JSON = ROOT / "particles" / "runs" / "flavor" / "generation_bundle_branch_generator.json"
GENERATION_BUNDLE_OBSTRUCTION_JSON = (
    ROOT / "particles" / "runs" / "flavor" / "generation_bundle_branch_generator_splitting_obstruction.json"
)
UNDERDETERMINATION_JSON = RUNS_DIR / "charged_absolute_scale_underdetermination_theorem.json"
TRACE_LIFT_JSON = RUNS_DIR / "charged_uncentered_trace_lift_scaffold.json"
TRACE_LIFT_COCYCLE_JSON = RUNS_DIR / "charged_uncentered_trace_lift_cocycle_reduction.json"
TRACE_LIFT_PHYSICAL_DESCENT_JSON = RUNS_DIR / "charged_mu_physical_descent_reduction.json"
PHYSICAL_EQUALIZER_JSON = RUNS_DIR / "charged_physical_identity_mode_equalizer.json"
CENTERED_OPERATOR_MU_NO_GO_JSON = RUNS_DIR / "charged_centered_operator_mu_phys_no_go.json"
DETERMINANT_LINE_JSON = RUNS_DIR / "charged_determinant_line_section_extension.json"
ANCHOR_SECTION_JSON = RUNS_DIR / "charged_absolute_anchor_section.json"
POST_PROMOTION_ROUTE_JSON = RUNS_DIR / "charged_post_promotion_absolute_closure_route.json"
ABSOLUTE_FRONTIER_FACTORIZATION_JSON = RUNS_DIR / "charged_absolute_frontier_factorization.json"
P_TO_AFFINE_REDUCTION_JSON = RUNS_DIR / "charged_p_to_affine_anchor_reduction.json"
PHYSICAL_CLASS_AFFINE_SCALAR_REDUCTION_JSON = RUNS_DIR / "charged_physical_class_affine_scalar_reduction.json"
PHYSICAL_CANONICAL_LIFT_JSON = RUNS_DIR / "charged_physical_determinant_line_canonical_uncentered_lift.json"
MASS_READOUT_FROM_AFFINE_JSON = RUNS_DIR / "charged_mass_readout_from_affine_anchor.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def artifact_ref(path: Path) -> str:
    try:
        return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()
    except ValueError:
        return str(path)


def anchor_input_contract() -> dict[str, list[str]]:
    return {
        "must_use": [
            "theorem-grade charged sector-response artifact after C_hat_e promotion",
            "lepton_current_family_exactness_audit.json",
        ],
        "must_not_use": [
            "measured charged masses",
            "compare-only D12 continuation targets",
            "PMNS or CKM target fits",
            "shared-budget seed alone",
            "centered-shape-only functionals",
        ],
    }


def anchor_hard_rejections(underdetermination: dict[str, Any]) -> dict[str, Any]:
    compare = underdetermination.get("compare_only_continuation_target", {})
    hard_reject = underdetermination.get("hard_reject", {})
    return {
        "common_shift_invariant_functionals": (
            "cannot emit A_ch because they satisfy F(logm + c*(1,1,1)) = F(logm)"
        ),
        "gamma_min_restore_pick": {
            "g_e": hard_reject.get("g_e"),
            "Delta_e_abs": hard_reject.get("Delta_e_abs"),
            "status": "invalid_orbit_pick",
        },
        "d12_compare_only_target": {
            "g_e_star": compare.get("g_e_star"),
            "Delta_e_abs_star": compare.get("delta_e_abs_star"),
            "status": "compare_only_not_theorem",
        },
    }


def trace_lift_scalar_cocycle_contract() -> dict[str, Any]:
    return {
        "ambient_identification": (
            "after theorem-grade centered identification of the promoted charged response "
            "across the live refinement family"
        ),
        "pairwise_difference_rule": "C_tilde_e(r') - C_tilde_e(r) = delta(r,r') I",
        "cocycle_laws": [
            "delta(r,r) = 0",
            "delta(r,r') = -delta(r',r)",
            "delta(r,r'') = delta(r,r') + delta(r',r'')",
        ],
        "primitive_required_on_fill": "exists mu(r) with delta(r,r') = mu(r') - mu(r)",
        "canonical_readout_on_fill": {
            "uncentered_trace_lift": "C_tilde_e(r) = C_hat_e(r) + mu(r) I",
            "determinant_line_section": "s_det(r) = 3 * mu(r)",
            "affine_anchor": "A_ch(r) = mu(r)",
        },
        "do_not_use": [
            "measured charged masses",
            "compare-only D12 continuation targets",
            "extra matrix-valued refinement data beyond scalar identity mode",
        ],
    }


def trace_lift_physical_descent_contract() -> dict[str, Any]:
    return {
        "precondition": (
            "theorem-grade centered promotion plus refinement stability on theorem-grade physical Y_e"
        ),
        "forced_refinement_identity_mode": "delta(r,r') = 0 for refinement representatives of the same physical Y_e",
        "descended_scalar": "exists unique mu_phys(Y_e) with mu(r) = mu_phys(Y_e) on every refinement representative r",
        "canonical_readout_on_fill": {
            "uncentered_trace_lift": "C_tilde_e(Y_e) = C_hat_e(Y_e) + mu_phys(Y_e) I",
            "determinant_line_section": "s_det(Y_e) = 3 * mu_phys(Y_e)",
            "affine_anchor": "A_ch(Y_e) = mu_phys(Y_e)",
        },
        "exact_smaller_missing_object": "charged_physical_affine_scalar_mu",
    }


def charged_waiting_set(generation_bundle: dict[str, Any]) -> dict[str, Any]:
    promotion_gate = dict(generation_bundle.get("promotion_gate", {}))
    candidate = dict(generation_bundle.get("charged_sector_response_operator_candidate", {}))
    return {
        "mandatory_package_a": {
            "id": "charged_sector_response_pushforward_to_C_hat_e",
            "closed_issue": "papers.compact.e.30a-resolve-branch-generator-splitting-or-no-go",
            "summary": (
                "The current corpus emits an exact obstruction certificate against promoting "
                "the latent charged operator candidate C_hat_e^{cand} from the displayed "
                "projector and proxy-generator data alone."
            ),
            "immediate_downstream_effect": (
                "No theorem-grade C_hat_e follows from the displayed current-corpus transfer "
                "data alone, so the live open charged issues begin at the post-promotion lift "
                "and affine-anchor chain rather than at the exact obstruction certificate."
            ),
            "status": "resolved_as_exact_current_corpus_obstruction_certificate",
            "exact_obstruction_certificate_artifact": (
                "oph_generation_bundle_branch_generator_splitting_obstruction_certificate"
            ),
            "exact_obstruction_certificate_artifact_ref": artifact_ref(
                GENERATION_BUNDLE_OBSTRUCTION_JSON
            ),
            "blocked_candidate_object": candidate.get("name", "C_hat_e^{cand}"),
            "upstream_missing_theorem": candidate.get(
                "declaration_missing_theorem",
                generation_bundle.get("remaining_missing_theorem"),
            ),
            "smallest_missing_clause": candidate.get(
                "smallest_missing_clause",
                promotion_gate.get("smaller_exact_missing_clause"),
            ),
            "exact_vanishing_proved": promotion_gate.get("exact_vanishing_proved"),
            "uniform_quadratic_smallness_proved": promotion_gate.get("uniform_quadratic_smallness_proved"),
            "current_strength_statement": promotion_gate.get("current_strength_statement"),
            "reopen_condition": (
                "a stronger transfer theorem beyond the displayed current-corpus projector and "
                "commutator data"
            ),
        },
        "mandatory_package_b": {
            "id": "refinement_stable_uncentered_trace_lift",
            "linked_issue": "papers.compact.e.30b-derive-post-promotion-uncentered-trace-lift",
            "summary": (
                "Derive the refinement-stable uncentered trace lift of the promoted charged "
                "response on theorem-grade physical Y_e or an equivalent determinant line."
            ),
            "immediate_downstream_effect": (
                "The determinant-line section and A_ch are then induced canonically by "
                "A_ch = (1/3) log det(Y_e) = (1/3) tr(log Y_e)."
            ),
            "status": "open_future_single_slot_only",
            "replaces_invalid_route": "charged_common_refinement_transport_equalizer",
        },
        "mandatory_package_c": {
            "id": "charged_absolute_anchor_A_ch",
            "linked_issue": "papers.compact.e.30c-derive-d10-p-to-charged-affine-anchor",
            "summary": (
                "Read out the affine charged absolute coordinate A_ch from the induced "
                "determinant-line section, with A_ch(logm + c*(1,1,1)) = A_ch(logm) + c."
            ),
            "immediate_downstream_effect": (
                "Once A_ch exists, the absolute charged scale is emitted by g_e = exp(A_ch), "
                "and Delta_e_abs follows as log(g_ch_shared) - A_ch."
            ),
            "status": "derived_once_package_b_exists",
        },
        "optional_package_d": {
            "id": "charged_holonomy_bridge_for_legacy_delta_2_over_9",
            "summary": (
                "Retain a charged holonomy bridge only if the older delta = 2/9 D12 benchmark "
                "is kept as a theorem-grade bridge instead of a diagnostic continuation."
            ),
            "required_only_if": "legacy_continuation_bridge_kept_as_theorem_grade",
        },
    }
