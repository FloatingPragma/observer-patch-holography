#!/usr/bin/env python3
"""Audit the exactness gap on the current local charged-lepton family."""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path

from charged_absolute_route_common import charged_waiting_set


ROOT = Path(__file__).resolve().parents[2]
REFERENCE_JSON = ROOT / "particles" / "data" / "particle_reference_values.json"
FORWARD_JSON = ROOT / "particles" / "runs" / "leptons" / "forward_charged_leptons.json"
READOUT_JSON = ROOT / "particles" / "runs" / "leptons" / "lepton_log_spectrum_readout.json"
ORDERED_PACKAGE_READBACK_JSON = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_ordered_package_readback.json"
CURRENT_SUPPORT_OBSTRUCTION_JSON = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_current_support_obstruction_certificate.json"
SUPPORT_EXTENSION_EMITTER_JSON = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_minimal_source_support_extension_emitter.json"
SUPPORT_EXTENSION_COMPLETION_LAW_JSON = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_support_extension_completion_law.json"
ETA_SOURCE_READBACK_JSON = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_support_extension_eta_source_readback.json"
ENDPOINT_RATIO_BREAKER_JSON = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_support_extension_endpoint_ratio_breaker.json"
SOURCE_SCALAR_PAIR_READBACK_JSON = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_support_extension_source_scalar_pair_readback.json"
CHARGED_D12_CONTINUATION_JSON = ROOT / "particles" / "runs" / "leptons" / "charged_d12_continuation_followup.json"
ABSOLUTE_SCALE_GAP_IDENTITY_JSON = ROOT / "particles" / "runs" / "leptons" / "charged_absolute_scale_transport_gap_identity.json"
ABSOLUTE_SCALE_UNDERDETERMINATION_JSON = ROOT / "particles" / "runs" / "leptons" / "charged_absolute_scale_underdetermination_theorem.json"
GENERATION_BUNDLE_JSON = ROOT / "particles" / "runs" / "flavor" / "generation_bundle_branch_generator.json"
END_TO_END_IMPOSSIBILITY_JSON = ROOT / "particles" / "runs" / "leptons" / "charged_end_to_end_impossibility_theorem.json"
POST_PROMOTION_ROUTE_JSON = (
    ROOT / "particles" / "runs" / "leptons" / "charged_post_promotion_absolute_closure_route.json"
)
ABSOLUTE_FRONTIER_FACTORIZATION_JSON = (
    ROOT / "particles" / "runs" / "leptons" / "charged_absolute_frontier_factorization.json"
)
TRACE_LIFT_COCYCLE_REDUCTION_JSON = (
    ROOT / "particles" / "runs" / "leptons" / "charged_uncentered_trace_lift_cocycle_reduction.json"
)
TRACE_LIFT_PHYSICAL_DESCENT_JSON = (
    ROOT / "particles" / "runs" / "leptons" / "charged_mu_physical_descent_reduction.json"
)
CENTERED_OPERATOR_MU_NO_GO_JSON = (
    ROOT / "particles" / "runs" / "leptons" / "charged_centered_operator_mu_phys_no_go.json"
)
DEFAULT_OUT = ROOT / "particles" / "runs" / "leptons" / "lepton_current_family_exactness_audit.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _centered_logs(values: list[float]) -> tuple[list[float], float]:
    logs = [math.log(value) for value in values]
    mean_log = sum(logs) / len(logs)
    return [value - mean_log for value in logs], mean_log


def _residual_norm(values: list[float]) -> float:
    return math.sqrt(sum(value * value for value in values))


def main() -> int:
    parser = argparse.ArgumentParser(description="Build a charged-lepton current-family exactness audit artifact.")
    parser.add_argument("--forward", default=str(FORWARD_JSON))
    parser.add_argument("--readout", default=str(READOUT_JSON))
    parser.add_argument("--ordered-package-readback", default=str(ORDERED_PACKAGE_READBACK_JSON))
    parser.add_argument("--current-support-obstruction", default=str(CURRENT_SUPPORT_OBSTRUCTION_JSON))
    parser.add_argument("--support-extension-emitter", default=str(SUPPORT_EXTENSION_EMITTER_JSON))
    parser.add_argument("--support-extension-completion-law", default=str(SUPPORT_EXTENSION_COMPLETION_LAW_JSON))
    parser.add_argument("--eta-source-readback", default=str(ETA_SOURCE_READBACK_JSON))
    parser.add_argument("--endpoint-ratio-breaker", default=str(ENDPOINT_RATIO_BREAKER_JSON))
    parser.add_argument("--source-scalar-pair-readback", default=str(SOURCE_SCALAR_PAIR_READBACK_JSON))
    parser.add_argument("--charged-d12-continuation", default=str(CHARGED_D12_CONTINUATION_JSON))
    parser.add_argument("--absolute-scale-gap-identity", default=str(ABSOLUTE_SCALE_GAP_IDENTITY_JSON))
    parser.add_argument("--absolute-scale-underdetermination", default=str(ABSOLUTE_SCALE_UNDERDETERMINATION_JSON))
    parser.add_argument("--generation-bundle", default=str(GENERATION_BUNDLE_JSON))
    parser.add_argument("--post-promotion-route", default=str(POST_PROMOTION_ROUTE_JSON))
    parser.add_argument("--absolute-frontier-factorization", default=str(ABSOLUTE_FRONTIER_FACTORIZATION_JSON))
    parser.add_argument("--trace-lift-cocycle-reduction", default=str(TRACE_LIFT_COCYCLE_REDUCTION_JSON))
    parser.add_argument("--trace-lift-physical-descent", default=str(TRACE_LIFT_PHYSICAL_DESCENT_JSON))
    parser.add_argument("--centered-operator-mu-no-go", default=str(CENTERED_OPERATOR_MU_NO_GO_JSON))
    parser.add_argument("--end-to-end-impossibility", default=str(END_TO_END_IMPOSSIBILITY_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    references = json.loads(REFERENCE_JSON.read_text(encoding="utf-8"))["entries"]
    forward = json.loads(Path(args.forward).read_text(encoding="utf-8"))
    readout = json.loads(Path(args.readout).read_text(encoding="utf-8"))
    ordered_package_readback_path = Path(args.ordered_package_readback)
    ordered_package_readback = (
        json.loads(ordered_package_readback_path.read_text(encoding="utf-8"))
        if ordered_package_readback_path.exists()
        else None
    )
    obstruction_path = Path(args.current_support_obstruction)
    current_support_obstruction = (
        json.loads(obstruction_path.read_text(encoding="utf-8"))
        if obstruction_path.exists()
        else None
    )
    support_extension_path = Path(args.support_extension_emitter)
    support_extension_emitter = (
        json.loads(support_extension_path.read_text(encoding="utf-8"))
        if support_extension_path.exists()
        else None
    )
    support_extension_completion_law_path = Path(args.support_extension_completion_law)
    support_extension_completion_law = (
        json.loads(support_extension_completion_law_path.read_text(encoding="utf-8"))
        if support_extension_completion_law_path.exists()
        else None
    )
    eta_source_readback_path = Path(args.eta_source_readback)
    eta_source_readback = (
        json.loads(eta_source_readback_path.read_text(encoding="utf-8"))
        if eta_source_readback_path.exists()
        else None
    )
    endpoint_ratio_breaker_path = Path(args.endpoint_ratio_breaker)
    endpoint_ratio_breaker = (
        json.loads(endpoint_ratio_breaker_path.read_text(encoding="utf-8"))
        if endpoint_ratio_breaker_path.exists()
        else None
    )
    source_scalar_pair_readback_path = Path(args.source_scalar_pair_readback)
    source_scalar_pair_readback = (
        json.loads(source_scalar_pair_readback_path.read_text(encoding="utf-8"))
        if source_scalar_pair_readback_path.exists()
        else None
    )
    charged_d12_continuation_path = Path(args.charged_d12_continuation)
    charged_d12_continuation = (
        json.loads(charged_d12_continuation_path.read_text(encoding="utf-8"))
        if charged_d12_continuation_path.exists()
        else None
    )
    absolute_scale_gap_identity_path = Path(args.absolute_scale_gap_identity)
    absolute_scale_gap_identity = (
        json.loads(absolute_scale_gap_identity_path.read_text(encoding="utf-8"))
        if absolute_scale_gap_identity_path.exists()
        else None
    )
    absolute_scale_underdetermination_path = Path(args.absolute_scale_underdetermination)
    absolute_scale_underdetermination = (
        json.loads(absolute_scale_underdetermination_path.read_text(encoding="utf-8"))
        if absolute_scale_underdetermination_path.exists()
        else None
    )
    generation_bundle_path = Path(args.generation_bundle)
    generation_bundle = (
        json.loads(generation_bundle_path.read_text(encoding="utf-8"))
        if generation_bundle_path.exists()
        else None
    )
    post_promotion_route_path = Path(args.post_promotion_route)
    post_promotion_route = (
        json.loads(post_promotion_route_path.read_text(encoding="utf-8"))
        if post_promotion_route_path.exists()
        else None
    )
    absolute_frontier_factorization_path = Path(args.absolute_frontier_factorization)
    absolute_frontier_factorization = (
        json.loads(absolute_frontier_factorization_path.read_text(encoding="utf-8"))
        if absolute_frontier_factorization_path.exists()
        else None
    )
    trace_lift_cocycle_reduction_path = Path(args.trace_lift_cocycle_reduction)
    trace_lift_cocycle_reduction = (
        json.loads(trace_lift_cocycle_reduction_path.read_text(encoding="utf-8"))
        if trace_lift_cocycle_reduction_path.exists()
        else None
    )
    trace_lift_physical_descent_path = Path(args.trace_lift_physical_descent)
    trace_lift_physical_descent = (
        json.loads(trace_lift_physical_descent_path.read_text(encoding="utf-8"))
        if trace_lift_physical_descent_path.exists()
        else None
    )
    centered_operator_mu_no_go_path = Path(args.centered_operator_mu_no_go)
    centered_operator_mu_no_go = (
        json.loads(centered_operator_mu_no_go_path.read_text(encoding="utf-8"))
        if centered_operator_mu_no_go_path.exists()
        else None
    )
    end_to_end_impossibility_path = Path(args.end_to_end_impossibility)
    end_to_end_impossibility = (
        json.loads(end_to_end_impossibility_path.read_text(encoding="utf-8"))
        if end_to_end_impossibility_path.exists()
        else None
    )

    current = [float(value) for value in forward["singular_values_abs"]]
    target = [
        float(references["electron"]["value_gev"]),
        float(references["muon"]["value_gev"]),
        float(references["tau"]["value_gev"]),
    ]
    centered_current, mean_log_current = _centered_logs(current)
    centered_target, mean_log_target = _centered_logs(target)
    centered_residual = [centered_target[idx] - centered_current[idx] for idx in range(3)]
    current_gamma21 = centered_current[1] - centered_current[0]
    current_gamma32 = centered_current[2] - centered_current[1]
    target_gamma21 = centered_target[1] - centered_target[0]
    target_gamma32 = centered_target[2] - centered_target[1]
    current_sigma = centered_current[2] - centered_current[0]
    target_sigma = centered_target[2] - centered_target[0]
    x2 = float(readout["ordered_family_coordinate"][1])
    target_eta = ((1.0 + x2) * target_gamma32) - ((1.0 - x2) * target_gamma21)
    diagnostic_eta = float(readout.get("eta_e_rigid_fallback"))
    best_eta_on_current_sigma = (target_gamma32 - target_gamma21) + (x2 * current_sigma)
    best_gap_fit = {
        "gamma21_log_per_side": (((1.0 + x2) * current_sigma) - best_eta_on_current_sigma) / 2.0,
        "gamma32_log_per_side": (((1.0 - x2) * current_sigma) + best_eta_on_current_sigma) / 2.0,
    }
    best_centered_two_scalar = [
        -((2.0 * best_gap_fit["gamma21_log_per_side"]) + best_gap_fit["gamma32_log_per_side"]) / 3.0,
        (best_gap_fit["gamma21_log_per_side"] - best_gap_fit["gamma32_log_per_side"]) / 3.0,
        (best_gap_fit["gamma21_log_per_side"] + (2.0 * best_gap_fit["gamma32_log_per_side"])) / 3.0,
    ]
    best_two_scalar_residual = [
        centered_target[idx] - best_centered_two_scalar[idx]
        for idx in range(3)
    ]

    common_shift = mean_log_target - mean_log_current
    common_scale = math.exp(common_shift)
    best_common_shift_fit = [common_scale * value for value in current]
    rel_errors = [(best_common_shift_fit[idx] - target[idx]) / target[idx] for idx in range(3)]
    support_extension_candidate = None
    if support_extension_emitter is not None and support_extension_emitter.get("eta_source_support_extension_log_per_side_candidate") is not None:
        candidate_shape = [
            float(value)
            for value in support_extension_emitter.get("singular_values_abs_ext_candidate", [])
        ]
        support_extension_candidate = {
            "eta_source_support_extension_log_per_side_candidate": float(
                support_extension_emitter["eta_source_support_extension_log_per_side_candidate"]
            ),
            "kappa_ext_candidate": float(support_extension_emitter["kappa_ext_candidate"]),
            "shape_singular_values_ext_candidate": [
                float(value) for value in support_extension_emitter.get("shape_singular_values_ext_candidate", [])
            ],
            "singular_values_abs_ext_candidate": candidate_shape,
            "relative_errors_against_reference": (
                [
                    (candidate_shape[idx] - target[idx]) / target[idx]
                    for idx in range(3)
                ]
                if len(candidate_shape) == 3
                else None
            ),
            "candidate_next_single_residual_object": support_extension_emitter.get(
                "candidate_next_single_residual_object"
            ),
            "candidate_status": support_extension_emitter.get("candidate_support_extension_status"),
        }

    artifact = {
        "artifact": "oph_lepton_current_family_exactness_audit",
        "generated_utc": _timestamp(),
        "scope": "current_family_only",
        "current_candidate": {
            "singular_values_abs": current,
            "closure_state": forward.get("closure_state"),
            "shape_shift_missing": bool(forward.get("shape_shift_missing", False)),
        },
        "reference_targets": {
            "singular_values_abs": target,
        },
        "centered_hierarchy_audit": {
            "current_centered_log": centered_current,
            "target_centered_log": centered_target,
            "residual": centered_residual,
            "residual_norm": _residual_norm(centered_residual),
            "ratio_invariant_under_common_shift": True,
            "current_sigma_total_log_per_side": current_sigma,
            "target_sigma_total_log_per_side": target_sigma,
            "current_gap_pair": {
                "gamma21_log_per_side": current_gamma21,
                "gamma32_log_per_side": current_gamma32,
                "rho_gap_ratio": current_gamma21 / current_gamma32,
                "eta_e_split_log_per_side_diagnostic": diagnostic_eta,
            },
            "target_gap_pair": {
                "gamma21_log_per_side": target_gamma21,
                "gamma32_log_per_side": target_gamma32,
                "rho_gap_ratio": target_gamma21 / target_gamma32,
                "eta_e_split_log_per_side": target_eta,
            },
        },
        "two_scalar_on_current_sigma_audit": {
            "sigma_e_total_log_per_side_current": current_sigma,
            "eta_e_rigid_fallback": diagnostic_eta,
            "eta_e_best_fit_on_current_sigma": best_eta_on_current_sigma,
            "best_gap_pair_on_current_sigma": best_gap_fit,
            "best_centered_log_on_current_sigma": best_centered_two_scalar,
            "residual_after_best_eta_on_current_sigma": best_two_scalar_residual,
            "residual_norm_after_best_eta_on_current_sigma": _residual_norm(best_two_scalar_residual),
        },
        "common_shift_audit": {
            "shape_log_shift_formula": "common additive shift on E_e_log_centered",
            "best_common_shift": common_shift,
            "best_common_scale": common_scale,
            "best_common_shift_fit": best_common_shift_fit,
            "relative_errors_after_best_shift": rel_errors,
        },
        "active_readout_contract": {
            "hierarchy_mode": readout.get("hierarchy_mode"),
            "sigma_e_total_log_per_side": float(readout["sigma_e_total_log_per_side"]),
            "eta_e_split_log_per_side": readout.get("eta_e_split_log_per_side"),
            "eta_e_rigid_fallback": readout.get("eta_e_rigid_fallback"),
            "rho_ord": float(readout["rho_ord"]),
            "shape_log_shift_e": readout.get("shape_log_shift_e"),
            "smallest_constructive_missing_object": (
                support_extension_completion_law.get("smallest_constructive_missing_object")
                if support_extension_completion_law is not None
                else support_extension_emitter.get("smallest_constructive_missing_object")
                if support_extension_emitter is not None
                else current_support_obstruction.get("smallest_constructive_missing_object")
                if current_support_obstruction is not None
                else readout.get("smallest_constructive_missing_object")
            ),
        },
        "current_package_readback": None if ordered_package_readback is None else {
            "artifact": ordered_package_readback.get("artifact"),
            "proof_status": ordered_package_readback.get("proof_status"),
            "source_side_quadratic_midpoint_defect_log_per_side_emitted": ordered_package_readback.get("source_side_quadratic_midpoint_defect_log_per_side_emitted"),
            "current_package_linear_subray_only": ordered_package_readback.get("current_package_linear_subray_only"),
            "same_support_exhausted": ordered_package_readback.get("same_support_exhausted"),
            "a_e_log_coeff": ordered_package_readback.get("a_e_log_coeff"),
            "b_e_log_coeff": ordered_package_readback.get("b_e_log_coeff"),
        },
        "same_support_obstruction_audit": {
            "same_support_exhausted": True,
            "sigma_support_gap": target_sigma - current_sigma,
            "ordered_gap_ratio_current": current_gamma21 / current_gamma32,
            "ordered_gap_ratio_reference": target_gamma21 / target_gamma32,
            "best_same_sigma_two_scalar_residual_norm": _residual_norm(best_two_scalar_residual),
        },
        "current_support_obstruction_certificate": None if current_support_obstruction is None else {
            "artifact": current_support_obstruction.get("artifact"),
            "proof_status": current_support_obstruction.get("proof_status"),
            "same_support_exhausted": current_support_obstruction.get("same_support_exhausted"),
            "current_support_linear_subray_only": current_support_obstruction.get("current_support_linear_subray_only"),
            "same_support_transverse_coeff_closed": current_support_obstruction.get("same_support_transverse_coeff_closed"),
            "smallest_constructive_missing_object": current_support_obstruction.get("smallest_constructive_missing_object"),
        },
        "support_extension_emitter": None if support_extension_emitter is None else {
            "artifact": support_extension_emitter.get("artifact"),
            "proof_status": support_extension_emitter.get("proof_status"),
            "extension_basis_kind": support_extension_emitter.get("extension_basis_kind"),
            "new_beyond_support_scalar_slots_required": support_extension_emitter.get("new_beyond_support_scalar_slots_required"),
            "smallest_constructive_missing_object": support_extension_emitter.get("smallest_constructive_missing_object"),
            "fixed_current_span_certificate": support_extension_emitter.get("fixed_current_span_certificate"),
            "eta_source_support_extension_log_per_side_candidate": support_extension_emitter.get("eta_source_support_extension_log_per_side_candidate"),
            "candidate_support_extension_status": support_extension_emitter.get("candidate_support_extension_status"),
            "candidate_next_single_residual_object": support_extension_emitter.get("candidate_next_single_residual_object"),
        },
        "support_extension_completion_law": None if support_extension_completion_law is None else {
            "artifact": support_extension_completion_law.get("artifact"),
            "proof_status": support_extension_completion_law.get("proof_status"),
            "smallest_constructive_missing_object": support_extension_completion_law.get("smallest_constructive_missing_object"),
            "next_single_residual_object": support_extension_completion_law.get("next_single_residual_object"),
            "next_residual_after_debug_eta_promotion": support_extension_completion_law.get("next_residual_after_debug_eta_promotion"),
            "sigma_source_support_extension_total_log_per_side": support_extension_completion_law.get("sigma_source_support_extension_total_log_per_side"),
            "eta_source_support_extension_log_per_side": support_extension_completion_law.get("eta_source_support_extension_log_per_side"),
        },
        "support_extension_eta_source_readback": None if eta_source_readback is None else {
            "artifact": eta_source_readback.get("artifact"),
            "status": eta_source_readback.get("status"),
            "proof_status": eta_source_readback.get("proof_status"),
            "eta_readback_invariant_name": eta_source_readback.get("eta_readback_invariant_name"),
            "eta_equivalence_formula": eta_source_readback.get("eta_equivalence_formula"),
            "next_single_residual_object_after_eta": eta_source_readback.get("next_single_residual_object_after_eta"),
            "downstream_sigma_artifact": eta_source_readback.get("downstream_sigma_artifact"),
        },
        "support_extension_endpoint_ratio_breaker": None if endpoint_ratio_breaker is None else {
            "artifact": endpoint_ratio_breaker.get("artifact"),
            "status": endpoint_ratio_breaker.get("status"),
            "proof_status": endpoint_ratio_breaker.get("proof_status"),
            "precondition_residual_object": endpoint_ratio_breaker.get("precondition_residual_object"),
            "smallest_constructive_missing_object_within_primitive": endpoint_ratio_breaker.get("smallest_constructive_missing_object_within_primitive"),
            "endpoint_ratio_breaker_invariant_name": endpoint_ratio_breaker.get("endpoint_ratio_breaker_invariant_name"),
        },
        "support_extension_source_scalar_pair_readback": None if source_scalar_pair_readback is None else {
            "artifact": source_scalar_pair_readback.get("artifact"),
            "status": source_scalar_pair_readback.get("status"),
            "proof_status": source_scalar_pair_readback.get("proof_status"),
            "scalar_order": source_scalar_pair_readback.get("scalar_order"),
            "next_single_residual_object": source_scalar_pair_readback.get("next_single_residual_object"),
            "next_single_residual_object_after_eta": source_scalar_pair_readback.get("next_single_residual_object_after_eta"),
            "eta_readback_invariant_name": source_scalar_pair_readback.get("eta_readback_invariant_name"),
            "sigma_readback_invariant_name": source_scalar_pair_readback.get("sigma_readback_invariant_name"),
        },
        "absolute_scale_gap_identity": None if absolute_scale_gap_identity is None else {
            "artifact": absolute_scale_gap_identity.get("artifact"),
            "proof_status": absolute_scale_gap_identity.get("proof_status"),
            "identity_formula": absolute_scale_gap_identity.get("identity_formula"),
            "identity_residual": absolute_scale_gap_identity.get("identity_residual"),
            "typed_restore_formulas": absolute_scale_gap_identity.get("typed_restore_formulas"),
            "typed_restore_values": absolute_scale_gap_identity.get("typed_restore_values"),
        },
        "absolute_scale_underdetermination_theorem": None if absolute_scale_underdetermination is None else {
            "artifact": absolute_scale_underdetermination.get("artifact"),
            "proof_status": absolute_scale_underdetermination.get("proof_status"),
            "theorem_statement": absolute_scale_underdetermination.get("theorem_statement"),
            "same_carrier_mass_formulas": absolute_scale_underdetermination.get("same_carrier_mass_formulas"),
            "centered_sum_rule": absolute_scale_underdetermination.get("centered_sum_rule"),
            "determinant_rules": absolute_scale_underdetermination.get("determinant_rules"),
            "no_go_theorem": absolute_scale_underdetermination.get("no_go_theorem"),
            "shared_budget_seed": absolute_scale_underdetermination.get("shared_budget_seed"),
            "compare_only_continuation_target": absolute_scale_underdetermination.get("compare_only_continuation_target"),
            "next_exact_missing_object": absolute_scale_underdetermination.get("next_exact_missing_object"),
            "minimal_new_theorem": absolute_scale_underdetermination.get("minimal_new_theorem"),
        },
        "absolute_scale_closure_status": None if absolute_scale_underdetermination is None else {
            "present_chain_under_determines_g_e": True,
            "current_theorem_output": "E_e_log_centered mod common shift",
            "charged_absolute_equalizer_status": absolute_scale_underdetermination.get("charged_absolute_equalizer"),
            "compare_only_g_e_star": absolute_scale_underdetermination.get("compare_only_continuation_target", {}).get("g_e_star"),
            "compare_only_delta_e_abs_star": absolute_scale_underdetermination.get("compare_only_continuation_target", {}).get("delta_e_abs_star"),
            "hard_reject": absolute_scale_underdetermination.get("hard_reject"),
            "honest_missing_transport_scalar": absolute_scale_underdetermination.get("minimal_new_theorem", {}).get("required_new_scalar"),
            "honest_post_promotion_single_slot": (
                post_promotion_route.get("post_promotion_single_slot", {}).get("id")
                if post_promotion_route is not None
                else None
            ),
            "honest_post_promotion_internal_carrier": (
                post_promotion_route.get("post_promotion_single_slot", {}).get("internal_carrier")
                if post_promotion_route is not None
                else None
            ),
            "honest_post_promotion_exact_descended_scalar": (
                post_promotion_route.get("post_promotion_single_slot", {}).get("exact_descended_scalar", {}).get("id")
                if post_promotion_route is not None
                else None
            ),
            "promotion_only_centered_operator_no_go": (
                post_promotion_route.get("promotion_only_no_go", {}).get("theorem_id")
                if post_promotion_route is not None
                else None
            ),
            "absolute_frontier_factorization_artifact": (
                absolute_frontier_factorization.get("artifact")
                if absolute_frontier_factorization is not None
                else None
            ),
        },
        "trace_lift_cocycle_reduction": None if trace_lift_cocycle_reduction is None else {
            "artifact": trace_lift_cocycle_reduction.get("artifact"),
            "status": trace_lift_cocycle_reduction.get("status"),
            "single_slot_preserved": trace_lift_cocycle_reduction.get("single_slot_preserved"),
            "irreducible_new_degree_of_freedom": trace_lift_cocycle_reduction.get(
                "matrix_vs_scalar_content", {}
            ).get("irreducible_new_degree_of_freedom"),
            "pairwise_difference_rule": trace_lift_cocycle_reduction.get("scalar_cocycle_contract", {}).get(
                "pairwise_difference_rule"
            ),
            "primitive_required_on_fill": trace_lift_cocycle_reduction.get("scalar_cocycle_contract", {}).get(
                "primitive_required_on_fill"
            ),
        },
        "trace_lift_physical_descent": None if trace_lift_physical_descent is None else {
            "artifact": trace_lift_physical_descent.get("artifact"),
            "status": trace_lift_physical_descent.get("status"),
            "exact_smaller_missing_object": trace_lift_physical_descent.get("exact_smaller_missing_object"),
            "forced_refinement_identity_mode": trace_lift_physical_descent.get("forced_vanishing", {}).get(
                "on_same_physical_Y_e"
            ),
        },
        "centered_operator_mu_no_go": None if centered_operator_mu_no_go is None else {
            "artifact": centered_operator_mu_no_go.get("artifact"),
            "status": centered_operator_mu_no_go.get("status"),
            "theorem_id": centered_operator_mu_no_go.get("no_go_theorem", {}).get("id"),
            "forbidden_target": centered_operator_mu_no_go.get("target_scalar", {}).get("id"),
            "trace_zero_by_construction": centered_operator_mu_no_go.get("input_surface", {}).get(
                "trace_zero_by_construction"
            ),
        },
        "absolute_frontier_factorization": None if absolute_frontier_factorization is None else {
            "artifact": absolute_frontier_factorization.get("artifact"),
            "status": absolute_frontier_factorization.get("status"),
            "current_surface_missing_object": absolute_frontier_factorization.get("current_surface_layer", {}).get("exact_missing_object"),
            "post_promotion_single_slot": absolute_frontier_factorization.get("post_promotion_layer", {}).get("irreducible_single_slot", {}).get("id"),
            "post_promotion_internal_carrier": absolute_frontier_factorization.get("post_promotion_layer", {}).get(
                "irreducible_single_slot", {}
            ).get("internal_carrier"),
            "post_promotion_exact_descended_scalar": absolute_frontier_factorization.get(
                "frontier_ledger", {}
            ).get("post_promotion_exact_descended_scalar"),
            "promotion_only_no_go": absolute_frontier_factorization.get("post_promotion_layer", {}).get(
                "promotion_only_no_go", {}
            ).get("theorem_id"),
            "reduction_theorem_id": absolute_frontier_factorization.get("frontier_ledger", {}).get("reduction_theorem_id"),
            "theorem_statement": absolute_frontier_factorization.get("theorem_statement"),
        },
        "end_to_end_closure_decision": None if end_to_end_impossibility is None else {
            "artifact": end_to_end_impossibility.get("artifact"),
            "verdict": end_to_end_impossibility.get("verdict"),
            "closure_now": end_to_end_impossibility.get("closure_now"),
            "exact_irreducible_chain": end_to_end_impossibility.get("exact_irreducible_chain"),
            "induced_after_irreducible_chain": end_to_end_impossibility.get("induced_after_irreducible_chain"),
            "post_promotion_route_artifact": (
                post_promotion_route.get("artifact")
                if post_promotion_route is not None
                else None
            ),
            "theorem_forbid_emit_now": end_to_end_impossibility.get("theorem_forbid_emit_now"),
        },
        "charged_sector_response_operator_candidate": None if generation_bundle is None else {
            "name": generation_bundle.get("charged_sector_response_operator_candidate", {}).get("name", "C_hat_e^{cand}"),
            "artifact": generation_bundle.get("artifact"),
            "status": generation_bundle.get("charged_sector_response_operator_candidate", {}).get("declaration_status", "candidate_only"),
            "declaration_missing_theorem": generation_bundle.get("charged_sector_response_operator_candidate", {}).get("declaration_missing_theorem", generation_bundle.get("remaining_missing_theorem")),
            "smallest_missing_clause": generation_bundle.get("charged_sector_response_operator_candidate", {}).get(
                "smallest_missing_clause",
                generation_bundle.get("promotion_gate", {}).get("smaller_exact_missing_clause"),
            ),
            "exact_vanishing_proved": generation_bundle.get("promotion_gate", {}).get("exact_vanishing_proved"),
            "uniform_quadratic_smallness_proved": generation_bundle.get("promotion_gate", {}).get("uniform_quadratic_smallness_proved"),
            "current_strength_statement": generation_bundle.get("promotion_gate", {}).get("current_strength_statement"),
            "matrix": generation_bundle.get("charged_sector_response_operator_candidate", {}).get("matrix"),
            "ordered_spectrum": generation_bundle.get("charged_sector_response_operator_candidate", {}).get("ordered_spectrum"),
            "same_label_overlap_amplitudes": generation_bundle.get("projective_readout_certificate", {}).get("same_label_overlap_amplitudes"),
            "sigma_formula": generation_bundle.get("charged_sector_response_operator_candidate", {}).get("sigma_formula"),
            "eta_formula_on_current_family": "x2 * sigma(C_hat_e^{cand}) - 3 * lambda_mid(C_hat_e^{cand})",
            "ordered_family_coordinate_x2": x2,
            "sigma_current": current_sigma,
            "eta_current": float(readout["eta_e_split_log_per_side"]),
            "latent_in_flavor_chain": True,
            "declaration_status": "candidate_only",
            "declared_operator_name": None,
        },
        "exact_waiting_set": charged_waiting_set(generation_bundle or {}),
        "red_team_branch_verdict": {
            "status": "current_branch_cannot_be_closed_as_stated",
            "smallest_wrong_frontier": [
                "eta_source_support_extension_log_per_side",
                "sigma_source_support_extension_total_log_per_side",
            ],
            "smallest_missing_theorem_object": "oph_generation_bundle_branch_generator_splitting",
            "blocked_candidate_object": "C_hat_e^{cand}",
            "smallest_missing_clause": (
                generation_bundle.get("promotion_gate", {}).get("smaller_exact_missing_clause")
                if generation_bundle is not None
                else None
            ),
            "next_exact_object_after_that": "refinement_stable_uncentered_trace_lift",
            "post_promotion_single_slot": (
                post_promotion_route.get("post_promotion_single_slot", {}).get("id")
                if post_promotion_route is not None
                else None
            ),
            "next_exact_object_after_that_if_closed": "charged_absolute_anchor_A_ch",
            "do_not_promote": [
                "eta_source_support_extension_log_per_side",
                "sigma_source_support_extension_total_log_per_side",
                "compare_only_g_e_star",
                "compare_only_delta_e_abs_star",
                "Delta_e_abs = 0.30236566025890826",
                "g_e = 0.6822819838027987",
            ],
        },
        "charged_d12_continuation_followup": None if charged_d12_continuation is None else {
            "artifact": charged_d12_continuation.get("artifact"),
            "status": charged_d12_continuation.get("status"),
            "theorem_tier": charged_d12_continuation.get("theorem_tier"),
            "d12_continuation_assumptions": charged_d12_continuation.get("d12_continuation_assumptions"),
            "eta_source_support_extension_log_per_side": (
                charged_d12_continuation.get("d12_continuation_pair", {}).get("eta_source_support_extension_log_per_side")
            ),
            "sigma_source_support_extension_total_log_per_side": (
                charged_d12_continuation.get("d12_continuation_pair", {}).get("sigma_source_support_extension_total_log_per_side")
            ),
            "centered_log_residual_norm": (
                charged_d12_continuation.get("compare_only_shape_check_against_reference_masses", {}).get("centered_log_residual_norm")
            ),
            "g_e_compare_only_needed_for_exact_absolute_masses": (
                charged_d12_continuation.get("compare_only_shape_check_against_reference_masses", {}).get("g_e_compare_only_needed_for_exact_absolute_masses")
            ),
        },
        "support_extension_candidate_audit": support_extension_candidate,
        "sigma_support_gap": target_sigma - current_sigma,
        "ordered_gap_ratio_current": current_gamma21 / current_gamma32,
        "ordered_gap_ratio_reference": target_gamma21 / target_gamma32,
        "candidate_residual_object_if_rigid_eta_promoted": (
            support_extension_emitter.get("candidate_next_single_residual_object")
            if support_extension_emitter is not None
            else None
        ),
        "smallest_exact_obstruction": "the current ordered package is explicitly read back and its midpoint defect closes to zero on the present support, so the same-support family is exhausted and still fails to reproduce the charged hierarchy",
        "smallest_constructive_missing_object": (
            support_extension_completion_law.get("smallest_constructive_missing_object")
            if support_extension_completion_law is not None
            else support_extension_emitter.get("smallest_constructive_missing_object")
            if support_extension_emitter is not None
            else current_support_obstruction.get("smallest_constructive_missing_object")
            if current_support_obstruction is not None
            else "oph_charged_sector_local_current_support_obstruction_certificate"
        ),
        "notes": [
            "The current charged-lepton builder now exposes the two-scalar ordered-gap family in formula form plus one common shift.",
            "A common shift rescales all three masses together and cannot fix a large hierarchy mismatch.",
            "The current emitted sigma_e is far too small to reproduce the electron-tau spread even if eta_e is optimized on top of it.",
            "The ordered support and affine-quadratic parameterization are already fixed, and on the current support the midpoint-defect emitter closes to zero.",
            "The current-support obstruction certificate is now on disk, and the next charged mover is the minimal support-extension emitter on the canonical quadratic ordered direction.",
            "The full two-scalar support-extension completion law is now explicit on disk; the live same-carrier primitive is the eta source-readback, followed by the sigma endpoint-ratio breaker.",
            "The stronger same-carrier source-scalar pair readback is also now explicit on disk, collecting those eta and sigma invariants into one ordered primitive beneath the full completion shell.",
            "At theorem level, eta and sigma are no longer the deepest honest waiting set. The live builder still exposes eta then sigma as the first same-carrier residuals, but the paper-facing exact burden is first to promote the latent candidate C_hat_e^{cand} by closing the branch-generator splitting theorem, then to restore the lost affine mode through a refinement-stable uncentered trace lift of the charged response, from which the determinant-line section and the affine absolute coordinate A_ch are induced.",
            "That post-promotion lift slot is now reduced more sharply too: after centered promotion the only remaining ambiguity is a scalar affine cocycle primitive mu on the refinement family, not an extra matrix-valued theorem beyond the uncentered trace lift.",
            "And because the lift is already required to be refinement-stable on theorem-grade physical Y_e, that primitive descends further to one physical affine scalar mu_phys(Y_e).",
            "A sharper impossibility theorem is now on disk too: even a future theorem-grade centered C_hat_e cannot emit mu_phys(Y_e) by itself, because centered operator data stays common-shift invariant.",
            (
                "The charged sector-response operator remains undeclared: only the latent candidate C_hat_e^{cand} is on disk, and its promotion is blocked by the upstream theorem oph_generation_bundle_branch_generator_splitting together with the smaller clause compression_descendant_commutator_vanishes_or_is_uniformly_quadratic_small_after_central_split."
                if generation_bundle is not None
                else "No latent charged sector-response candidate is attached to this audit yet."
            ),
            (
                "On the live corpus, the commutator-transfer bridge proves neither exact vanishing nor uniform quadratic smallness after the central split; only the desired conditional bridge is recorded."
                if generation_bundle is not None
                else "No commutator-transfer strength statement is attached to this audit yet."
            ),
            (
                "The present charged theorem determines only the centered charged log class modulo a common additive shift, so no theorem-grade g_e, Delta_e_abs, or charged absolute equalizer exists on the live theorem lane."
                if absolute_scale_underdetermination is not None
                else "No explicit charged absolute-scale underdetermination theorem is attached to this audit yet."
            ),
            "The eta-only extension acts at fixed current span, so it preserves the endpoint ratio tau/e and can only move the middle state against that fixed endpoint pair.",
            "A rigid eta candidate can be written from the current ordered-gap ratio alone, but it lands far from the charged targets; the remaining charged burden then shifts to the still-open total source span.",
            (
                "A D12 continuation bridge is now explicit too: under the extra assumptions A1-A3 it emits eta = -6.729586682888832 and sigma = 8.154061112725994, with centered-log residual norm about 2.13e-05 against charged references, but the required absolute scale remains compare-only and therefore nonpromotable."
                if charged_d12_continuation is not None
                else "No D12 charged continuation bridge is attached to this audit."
            ),
            (
                "The current-family absolute-scale restore candidate is also cleaner now: the common gap subtracted from log(g_e_raw) matches the emitted overlap-edge theorem gap gamma on the current family, so mu_e_absolute_log_candidate = log(g_e_raw) - gamma_gap is explicit. But that restore shell is not theorem-grade, because it only chooses a representative on the common-shift orbit."
                if absolute_scale_gap_identity is not None
                else "No current-family charged absolute-scale gap identity is attached to this audit."
            ),
            (
                "The D12 continuation bridge remains compare-only on the absolute side: it would need g_e_star = 0.04577885783568762, equivalently Delta_e_abs_star = 3.003986333402356, but that value is not OPH-emitted."
                if absolute_scale_underdetermination is not None
                else "No compare-only charged absolute-scale target is attached to this audit."
            ),
        ],
    }

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
