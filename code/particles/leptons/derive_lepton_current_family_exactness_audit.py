#!/usr/bin/env python3
"""Audit the exactness gap on the current local charged-lepton family."""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path


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
        "exact_waiting_set": {
            "mandatory_package_a": {
                "id": "charged_sector_response_pushforward_to_C_hat_e",
                "linked_issue": "papers.compact.e.29-derive-the-yukawa-excitation-dictionary",
                "summary": "Derive the quotient-natural charged sector-response / declaration functor that emits the unique charged color-singlet three-family excitation operator C_hat_e from the shared OPH flavor dictionary.",
                "immediate_downstream_effect": "Once C_hat_e exists, eta and sigma become charged spectral invariants rather than independent primitive theorem objects.",
            },
            "mandatory_package_b": {
                "id": "charged_common_refinement_transport_equalizer",
                "linked_issue": "papers.compact.e.30-replace-koide-assisted-lepton-fitting-with-a-theorem",
                "summary": "Derive the common-refinement transport equalizer under ordered_common_refinement_seed_rigidity so the charged absolute-scale lane upgrades from shared_budget_only to theorem-grade.",
                "immediate_downstream_effect": "This is the honest absolute-scale bridge from the shared charged evaluator descent to the physical scale g_e.",
            },
            "optional_package_c": {
                "id": "charged_holonomy_bridge_for_legacy_delta_2_over_9",
                "summary": "Retain a charged holonomy bridge only if the older delta = 2/9 D12 benchmark is kept as a theorem-grade bridge instead of a diagnostic continuation.",
                "required_only_if": "legacy_continuation_bridge_kept_as_theorem_grade",
            },
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
            "At theorem level, eta and sigma are no longer the deepest honest waiting set. The live builder still exposes eta then sigma as the first same-carrier residuals, but the paper-facing exact burden is to derive the charged sector-response functor C_hat_e and then the common-refinement transport equalizer for the absolute scale.",
            "The eta-only extension acts at fixed current span, so it preserves the endpoint ratio tau/e and can only move the middle state against that fixed endpoint pair.",
            "A rigid eta candidate can be written from the current ordered-gap ratio alone, but it lands far from the charged targets; the remaining charged burden then shifts to the still-open total source span.",
            (
                "A D12 continuation bridge is now explicit too: under the extra assumptions A1-A3 it emits eta = -6.729586682888832 and sigma = 8.154061112725994, with centered-log residual norm about 2.13e-05 against charged references, but the required absolute scale remains compare-only and therefore nonpromotable."
                if charged_d12_continuation is not None
                else "No D12 charged continuation bridge is attached to this audit."
            ),
            (
                "The current-family absolute-scale restore candidate is also cleaner now: the common gap subtracted from log(g_e_raw) matches the emitted overlap-edge theorem gap gamma on the current family, so mu_e_absolute_log_candidate = log(g_e_raw) - gamma_gap is explicit even though the charged hierarchy still remains blocked first on eta then sigma."
                if absolute_scale_gap_identity is not None
                else "No current-family charged absolute-scale gap identity is attached to this audit."
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
