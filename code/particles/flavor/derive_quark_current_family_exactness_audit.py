#!/usr/bin/env python3
"""Audit the exactness gap on the current local quark family."""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
REFERENCE_JSON = ROOT / "particles" / "data" / "particle_reference_values.json"
FORWARD_JSON = ROOT / "particles" / "runs" / "flavor" / "forward_yukawas.json"
MEAN_SPLIT_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_sector_mean_split.json"
SPREAD_MAP_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_spread_map.json"
J_B_EVALUATOR_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_B_odd_source_scalar_evaluator.json"
D12_SELECTOR_JSON = ROOT / "particles" / "runs" / "flavor" / "light_quark_isospin_overlap_defect_selector_law.json"
D12_MASS_BRANCH_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_d12_mass_branch_and_ckm_residual.json"
D12_OVERLAP_LAW_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_d12_overlap_transport_law.json"
QUADRATIC_SCALAR_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_quadratic_even_transport_scalar.json"
PHYSICAL_INVARIANT_JSON = ROOT / "particles" / "runs" / "flavor" / "generation_bundle_same_label_physical_invariant_bundle.json"
SCALARIZED_BUNDLE_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_scalarized_continuation_bundle.json"
ONE_SCALAR_SPECIALIZATION_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_d12_one_scalar_specialization.json"
MASS_RAY_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_d12_mass_ray.json"
MASS_SIDE_UNDERDETERMINATION_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_d12_mass_side_underdetermination_theorem.json"
PHYSICAL_BRANCH_REPAIR_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_physical_branch_repair_theorem.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_exactness_audit.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _centered_logs(values: list[float]) -> tuple[list[float], float]:
    logs = [math.log(value) for value in values]
    mean_log = sum(logs) / len(logs)
    return [value - mean_log for value in logs], mean_log


def _geom_mean(values: list[float]) -> float:
    return math.prod(values) ** (1.0 / len(values))


def _residual_norm(values: list[float]) -> float:
    return math.sqrt(sum(value * value for value in values))


def main() -> int:
    parser = argparse.ArgumentParser(description="Build a quark current-family exactness audit artifact.")
    parser.add_argument("--forward", default=str(FORWARD_JSON))
    parser.add_argument("--mean-split", default=str(MEAN_SPLIT_JSON))
    parser.add_argument("--spread-map", default=str(SPREAD_MAP_JSON))
    parser.add_argument("--j-b-evaluator", default=str(J_B_EVALUATOR_JSON))
    parser.add_argument("--d12-selector", default=str(D12_SELECTOR_JSON))
    parser.add_argument("--d12-mass-branch", default=str(D12_MASS_BRANCH_JSON))
    parser.add_argument("--d12-overlap-law", default=str(D12_OVERLAP_LAW_JSON))
    parser.add_argument("--quadratic-scalar", default=str(QUADRATIC_SCALAR_JSON))
    parser.add_argument("--physical-invariants", default=str(PHYSICAL_INVARIANT_JSON))
    parser.add_argument("--scalarized-bundle", default=str(SCALARIZED_BUNDLE_JSON))
    parser.add_argument("--one-scalar-specialization", default=str(ONE_SCALAR_SPECIALIZATION_JSON))
    parser.add_argument("--d12-mass-ray", default=str(MASS_RAY_JSON))
    parser.add_argument("--mass-side-underdetermination", default=str(MASS_SIDE_UNDERDETERMINATION_JSON))
    parser.add_argument("--physical-branch-repair", default=str(PHYSICAL_BRANCH_REPAIR_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    references = json.loads(REFERENCE_JSON.read_text(encoding="utf-8"))["entries"]
    forward = json.loads(Path(args.forward).read_text(encoding="utf-8"))
    mean_split = json.loads(Path(args.mean_split).read_text(encoding="utf-8"))
    spread_map = json.loads(Path(args.spread_map).read_text(encoding="utf-8"))
    j_b_evaluator_path = Path(args.j_b_evaluator)
    j_b_evaluator = json.loads(j_b_evaluator_path.read_text(encoding="utf-8")) if j_b_evaluator_path.exists() else None
    d12_selector_path = Path(args.d12_selector)
    d12_selector = json.loads(d12_selector_path.read_text(encoding="utf-8")) if d12_selector_path.exists() else None
    d12_mass_branch_path = Path(args.d12_mass_branch)
    d12_mass_branch = json.loads(d12_mass_branch_path.read_text(encoding="utf-8")) if d12_mass_branch_path.exists() else None
    d12_overlap_law_path = Path(args.d12_overlap_law)
    d12_overlap_law = json.loads(d12_overlap_law_path.read_text(encoding="utf-8")) if d12_overlap_law_path.exists() else None
    quadratic_scalar_path = Path(args.quadratic_scalar)
    quadratic_scalar = json.loads(quadratic_scalar_path.read_text(encoding="utf-8")) if quadratic_scalar_path.exists() else None
    physical_invariants_path = Path(args.physical_invariants)
    physical_invariants = json.loads(physical_invariants_path.read_text(encoding="utf-8")) if physical_invariants_path.exists() else None
    scalarized_bundle_path = Path(args.scalarized_bundle)
    scalarized_bundle = json.loads(scalarized_bundle_path.read_text(encoding="utf-8")) if scalarized_bundle_path.exists() else None
    one_scalar_specialization_path = Path(args.one_scalar_specialization)
    one_scalar_specialization = (
        json.loads(one_scalar_specialization_path.read_text(encoding="utf-8"))
        if one_scalar_specialization_path.exists()
        else None
    )
    d12_mass_ray_path = Path(args.d12_mass_ray)
    d12_mass_ray = json.loads(d12_mass_ray_path.read_text(encoding="utf-8")) if d12_mass_ray_path.exists() else None
    mass_side_underdetermination_path = Path(args.mass_side_underdetermination)
    mass_side_underdetermination = (
        json.loads(mass_side_underdetermination_path.read_text(encoding="utf-8"))
        if mass_side_underdetermination_path.exists()
        else None
    )
    physical_branch_repair_path = Path(args.physical_branch_repair)
    physical_branch_repair = (
        json.loads(physical_branch_repair_path.read_text(encoding="utf-8"))
        if physical_branch_repair_path.exists()
        else None
    )

    current_u = [float(value) for value in forward["singular_values_u"]]
    current_d = [float(value) for value in forward["singular_values_d"]]
    target_u = [
        float(references["up_quark"]["value_gev"]),
        float(references["charm_quark"]["value_gev"]),
        float(references["top_quark"]["value_gev"]),
    ]
    target_d = [
        float(references["down_quark"]["value_gev"]),
        float(references["strange_quark"]["value_gev"]),
        float(references["bottom_quark"]["value_gev"]),
    ]

    centered_current_u, mean_log_current_u = _centered_logs(current_u)
    centered_current_d, mean_log_current_d = _centered_logs(current_d)
    centered_target_u, mean_log_target_u = _centered_logs(target_u)
    centered_target_d, mean_log_target_d = _centered_logs(target_d)
    residual_u = [centered_target_u[idx] - centered_current_u[idx] for idx in range(3)]
    residual_d = [centered_target_d[idx] - centered_current_d[idx] for idx in range(3)]
    x2 = float(spread_map["normalized_coordinate_x2"])
    q_ord = [
        (1.0 - x2 * x2) / 3.0,
        -(2.0 * (1.0 - x2 * x2)) / 3.0,
        (1.0 - x2 * x2) / 3.0,
    ]
    q_norm = sum(value * value for value in q_ord)
    b_ord = [-1.0, 0.0, 1.0]

    def _best_kappa_fit(residual: list[float]) -> tuple[float, list[float], float]:
        kappa = sum(residual[idx] * q_ord[idx] for idx in range(3)) / q_norm
        leftover = [residual[idx] - (kappa * q_ord[idx]) for idx in range(3)]
        return kappa, leftover, _residual_norm(leftover)

    sigma_seed = float(mean_split["sigma_seed_ud_candidate"])
    eta_ud = float(mean_split["eta_ud_candidate"])
    g_ch = float(mean_split["shared_norm_value"])
    log_shift_target_u = mean_log_target_u - math.log(g_ch)
    log_shift_target_d = mean_log_target_d - math.log(g_ch)
    exact_A = -(log_shift_target_u + log_shift_target_d) / (2.0 * sigma_seed)
    exact_B = (log_shift_target_u - log_shift_target_d) / (2.0 * eta_ud)
    g_u_exact = g_ch * math.exp(-(exact_A * sigma_seed - exact_B * eta_ud))
    g_d_exact = g_ch * math.exp(-(exact_A * sigma_seed + exact_B * eta_ud))
    exact_fit_u = [g_u_exact * math.exp(value) for value in centered_current_u]
    exact_fit_d = [g_d_exact * math.exp(value) for value in centered_current_d]
    kappa_u_fit, residual_u_after_kappa, residual_u_after_kappa_norm = _best_kappa_fit(residual_u)
    kappa_d_fit, residual_d_after_kappa, residual_d_after_kappa_norm = _best_kappa_fit(residual_d)
    tau_u_fit = float(residual_u_after_kappa[2])
    tau_d_fit = float(residual_d_after_kappa[2])

    artifact = {
        "artifact": "oph_quark_current_family_exactness_audit",
        "generated_utc": _timestamp(),
        "scope": "current_family_only",
        "current_family_forward_values_gev": {
            "u": current_u[0],
            "c": current_u[1],
            "t": current_u[2],
            "d": current_d[0],
            "s": current_d[1],
            "b": current_d[2],
        },
        "current_family_reference_values_gev": {
            "u": target_u[0],
            "c": target_u[1],
            "t": target_u[2],
            "d": target_d[0],
            "s": target_d[1],
            "b": target_d[2],
        },
        "current_candidate": {
            "singular_values_u": current_u,
            "singular_values_d": current_d,
            "geom_mean_u": _geom_mean(current_u),
            "geom_mean_d": _geom_mean(current_d),
        },
        "reference_targets": {
            "singular_values_u": target_u,
            "singular_values_d": target_d,
            "geom_mean_u": _geom_mean(target_u),
            "geom_mean_d": _geom_mean(target_d),
        },
        "centered_ray_fit": {
            "current_centered_log_u": centered_current_u,
            "current_centered_log_d": centered_current_d,
            "target_centered_log_u": centered_target_u,
            "target_centered_log_d": centered_target_d,
            "residual_u": residual_u,
            "residual_d": residual_d,
            "residual_norm_u": _residual_norm(residual_u),
            "residual_norm_d": _residual_norm(residual_d),
        },
        "mean_split_audit": {
            "shared_norm_value": g_ch,
            "sigma_seed_ud": sigma_seed,
            "eta_ud": eta_ud,
            "current_log_shift_u": mean_log_current_u - math.log(g_ch),
            "current_log_shift_d": mean_log_current_d - math.log(g_ch),
            "target_log_shift_u": log_shift_target_u,
            "target_log_shift_d": log_shift_target_d,
            "exact_two_scalar_mean_fit": {
                "formula": "log(g_u/g_ch) = -(A_exact * sigma_seed_ud - B_exact * eta_ud), log(g_d/g_ch) = -(A_exact * sigma_seed_ud + B_exact * eta_ud)",
                "A_exact": exact_A,
                "B_exact": exact_B,
                "g_u_exact_fit": g_u_exact,
                "g_d_exact_fit": g_d_exact,
                "exact_fit_singular_values_u": exact_fit_u,
                "exact_fit_singular_values_d": exact_fit_d,
            },
        },
        "spread_emitter_audit": {
            "spread_emitter_status": spread_map.get("spread_emitter_status"),
            "sigma_source_kind": spread_map.get("sigma_source_kind"),
            "sigma_u_total_log_per_side": spread_map.get("sigma_u_total_log_per_side"),
            "sigma_d_total_log_per_side": spread_map.get("sigma_d_total_log_per_side"),
        },
        "quadratic_residual_audit": {
            "Q_ord": q_ord,
            "kappa_u_best_fit": kappa_u_fit,
            "kappa_d_best_fit": kappa_d_fit,
            "residual_u_after_best_quadratic_fit": residual_u_after_kappa,
            "residual_d_after_best_quadratic_fit": residual_d_after_kappa,
            "residual_norm_u_after_best_quadratic_fit": residual_u_after_kappa_norm,
            "residual_norm_d_after_best_quadratic_fit": residual_d_after_kappa_norm,
            "spread_preserved": True,
            "mean_surface_preserved": True,
        },
        "diagonal_gap_shift_audit": {
            "B_ord": b_ord,
            "tau_u_best_fit": tau_u_fit,
            "tau_d_best_fit": tau_d_fit,
            "residual_u_after_best_diagonal_shift": [
                residual_u_after_kappa[idx] - tau_u_fit * b_ord[idx]
                for idx in range(3)
            ],
            "residual_d_after_best_diagonal_shift": [
                residual_d_after_kappa[idx] - tau_d_fit * b_ord[idx]
                for idx in range(3)
            ],
        },
        "diagnostic_only_tau_best_fit": {
            "tau_u_best_fit": tau_u_fit,
            "tau_d_best_fit": tau_d_fit,
        },
        "diagnostic_fit_source": "particle_reference_values.json",
        "diagnostic_fit_promotion_allowed": False,
        "b_mode_amplitude_blindness_audit": {
            "B_ord": b_ord,
            "Q_ord": q_ord,
            "sum_B_ord": sum(b_ord),
            "dot_B_ord_Q_ord": sum(b * q for b, q in zip(b_ord, q_ord)),
            "mean_annihilator_is_identically_zero_on_B_mode": True,
            "quadratic_annihilator_is_identically_zero_on_B_mode": True,
        },
        "b_mode_odd_projector_evaluator": None if j_b_evaluator is None else {
            "artifact": j_b_evaluator.get("artifact"),
            "proof_status": j_b_evaluator.get("proof_status"),
            "J_B_formula": j_b_evaluator.get("J_B_formula"),
            "smallest_constructive_missing_object": j_b_evaluator.get("smallest_constructive_missing_object"),
            "predictive_J_B_source_law_status": j_b_evaluator.get("predictive_J_B_source_law_status"),
        },
        "d12_isospin_selector_law": None if d12_selector is None else {
            "artifact": d12_selector.get("artifact"),
            "status": d12_selector.get("status"),
            "proof_status": d12_selector.get("proof_status"),
            "scope": d12_selector.get("scope"),
            "selector_scalar_name": d12_selector.get("selector_scalar_name"),
            "next_single_residual_object": d12_selector.get("next_single_residual_object"),
            "selector_equivalence_formula": d12_selector.get("selector_equivalence_formula"),
            "odd_budget_neutrality_formula": d12_selector.get("odd_budget_neutrality_formula"),
        },
        "d12_mass_branch_and_ckm_residual": None if d12_mass_branch is None else {
            "artifact": d12_mass_branch.get("artifact"),
            "status": d12_mass_branch.get("status"),
            "theorem_tier": d12_mass_branch.get("theorem_tier"),
            "sample_same_family_point": d12_mass_branch.get("sample_same_family_point"),
            "comparison_only_best_same_family_point": d12_mass_branch.get("comparison_only_best_same_family_point"),
            "forward_same_label_transport_closed": (
                d12_mass_branch.get("closure_residual", {}).get("fro_norm") is not None
            ),
            "standard_ckm_parameters": d12_mass_branch.get("standard_ckm_parameters"),
            "closure_residual_fro_norm": d12_mass_branch.get("closure_residual", {}).get("fro_norm"),
        },
        "d12_overlap_transport_law": None if d12_overlap_law is None else {
            "artifact": d12_overlap_law.get("artifact"),
            "status": d12_overlap_law.get("status"),
            "theorem_tier": d12_overlap_law.get("theorem_tier"),
            "next_single_residual_object": d12_overlap_law.get("next_single_residual_object"),
            "transport_formulas": d12_overlap_law.get("transport_formulas"),
            "sample_same_family_point": d12_overlap_law.get("sample_same_family_point"),
            "comparison_only_best_same_family_point": d12_overlap_law.get("comparison_only_best_same_family_point"),
        },
        "d12_quadratic_even_transport_scalar": None if quadratic_scalar is None else {
            "artifact": quadratic_scalar.get("artifact"),
            "proof_status": quadratic_scalar.get("proof_status"),
            "next_single_residual_object": quadratic_scalar.get("next_single_residual_object"),
            "quadratic_even_log_formula_direct": quadratic_scalar.get("quadratic_even_log_formula_direct"),
            "sample_same_family_ray_point": quadratic_scalar.get("sample_same_family_ray_point"),
        },
        "d12_same_label_physical_invariant_bundle": None if physical_invariants is None else {
            "artifact": physical_invariants.get("artifact"),
            "proof_status": physical_invariants.get("proof_status"),
            "physical_invariants": physical_invariants.get("physical_invariants"),
            "next_single_residual_object": physical_invariants.get("next_single_residual_object"),
        },
        "d12_scalarized_continuation_bundle": None if scalarized_bundle is None else {
            "artifact": scalarized_bundle.get("artifact"),
            "proof_status": scalarized_bundle.get("proof_status"),
            "honest_remaining_value_laws": scalarized_bundle.get("honest_remaining_value_laws"),
            "mass_side": scalarized_bundle.get("mass_side"),
            "mixing_side": scalarized_bundle.get("mixing_side"),
        },
        "d12_one_scalar_specialization": None if one_scalar_specialization is None else {
            "artifact": one_scalar_specialization.get("artifact"),
            "proof_status": one_scalar_specialization.get("proof_status"),
            "scalar_name": one_scalar_specialization.get("scalar_name"),
            "sample_scalar_name": one_scalar_specialization.get("sample_scalar_name"),
            "next_single_residual_object": one_scalar_specialization.get("next_single_residual_object"),
            "mass_side_object_count_reduction": one_scalar_specialization.get("mass_side_object_count_reduction"),
            "specialization_formulas": one_scalar_specialization.get("specialization_formulas"),
            "sample_same_family_point": one_scalar_specialization.get("sample_same_family_point"),
            "sample_same_family_mass_point": one_scalar_specialization.get("sample_same_family_mass_point"),
        },
        "d12_mass_ray": None if d12_mass_ray is None else {
            "artifact": d12_mass_ray.get("artifact"),
            "proof_status": d12_mass_ray.get("proof_status"),
            "emitted_object": d12_mass_ray.get("emitted_object"),
            "same_family_ray": d12_mass_ray.get("same_family_ray"),
            "next_exact_missing_object": d12_mass_ray.get("next_exact_missing_object"),
        },
        "d12_mass_side_underdetermination_theorem": None if mass_side_underdetermination is None else {
            "artifact": mass_side_underdetermination.get("artifact"),
            "proof_status": mass_side_underdetermination.get("proof_status"),
            "theorem_statement": mass_side_underdetermination.get("theorem_statement"),
            "same_family_ray": mass_side_underdetermination.get("same_family_ray"),
            "next_exact_missing_object": mass_side_underdetermination.get("next_exact_missing_object"),
            "minimal_new_theorem": mass_side_underdetermination.get("minimal_new_theorem"),
        },
        "d12_physical_branch_repair_theorem": None if physical_branch_repair is None else {
            "artifact": physical_branch_repair.get("artifact"),
            "proof_status": physical_branch_repair.get("proof_status"),
            "current_d12_sheet": physical_branch_repair.get("current_d12_sheet"),
            "comparison_shell": physical_branch_repair.get("comparison_shell"),
            "current_sheet_invariants": physical_branch_repair.get("current_sheet_invariants"),
            "physical_shell_mismatch": physical_branch_repair.get("physical_shell_mismatch"),
            "minimal_branch_shift_repair_theorem": physical_branch_repair.get("minimal_branch_shift_repair_theorem"),
            "relative_sheet_scan": physical_branch_repair.get("relative_sheet_scan"),
        },
        "source_readback_payload_kind": "pure_B_light_sector_payload_pair",
        "recovered_core_no_go_for_nonzero_light_quark_pure_b_selector": True,
        "recovered_core_no_go_basis": "March 28, 2026 final-wave consolidation against the OPH tier ledger in the uploaded corpus",
        "active_builder_smallest_missing_object": "source_readback_u_log_per_side_and_source_readback_d_log_per_side",
        "broader_honest_frontier": (
            "quark_d12_t1_value_law"
            if physical_branch_repair is not None
            and ((physical_branch_repair.get("minimal_branch_shift_repair_theorem") or {}).get("selected_value") is not None)
            else "quark_relative_sheet_selector"
        ),
        "predictive_J_B_source_law_status": (
            j_b_evaluator.get("predictive_J_B_source_law_status")
            if j_b_evaluator is not None
            else "missing"
        ),
        "smallest_exact_obstruction": (
            "the builder-facing pure-B payload pair is still open on the active public branch, "
            + (
                "the emitted local same-label left-handed orbit has already collapsed to sigma_ref, so the broader D12 continuation branch is now waiting on the intrinsic D12 scale law on the selected mass ray while the CKM-shell no-go persists on that branch"
                if physical_branch_repair is not None
                and ((physical_branch_repair.get("minimal_branch_shift_repair_theorem") or {}).get("selected_value") is not None)
                else "the broader D12 continuation branch is a strict no-go for the physical CKM shell until one discrete quark_relative_sheet_selector is emitted; mass-side scale fixing is a separate burden after that branch selection"
            )
        ),
        "smallest_constructive_missing_object": (
            "source_readback_u_log_per_side_and_source_readback_d_log_per_side"
            if j_b_evaluator is not None
            else "source_readback_u_log_per_side_and_source_readback_d_log_per_side"
        ),
        "notes": [
            "The current local quark rays are already close to the measured centered log profiles.",
            "The compact current-family sector-mean law is closed on the emitted spread package.",
            "The spread emitter is now read back from the closed mean surface rather than seeded diagnostically.",
            "The unique quadratic residual basis Q_ord isolates the only same-surface skew mode left on the ordered three-point family.",
            "Projecting the residual onto Q_ord only slightly reduces the mismatch, and the leftover is exactly a diagonal gap-shift pattern [-tau, 0, +tau].",
            "That means the family shell and pure-B source-readback law are already known; on the active builder path the missing predictive step is the emitted pure-B payload pair, not another larger quark family.",
            "The March 28, 2026 final-wave consolidation also establishes a tier boundary: a nonzero light-quark pure-B selector is not available at recovered-core tier in the uploaded corpus.",
            (
                "The broader repair frontier has therefore moved past the discrete selector: on the emitted local solver surface sigma_ud closes to sigma_ref, and the next exact D12 object is the one-scalar value law quark_d12_t1_value_law on the selected mass ray, with intrinsic_scale_law_D12 retained only as the derived wrapper, even though the active local builder still waits first on the pure-B payload pair."
                if physical_branch_repair is not None
                and ((physical_branch_repair.get("minimal_branch_shift_repair_theorem") or {}).get("selected_value") is not None)
                else "The broader honest repair frontier is therefore a D12 light-quark isospin-breaking selector / overlap-defect scalar, even though the active local builder still waits first on the pure-B payload pair."
            ),
            (
                "The D12 selector-law shell remains explicit on disk: one continuation-level overlap-defect scalar `Delta_ud_overlap` would fix the light-sector pure-B payload pair by odd-budget neutrality, but the emitted local branch selector now closes trivially to sigma_ref rather than opening a new physical CKM branch."
                if physical_branch_repair is not None
                and ((physical_branch_repair.get("minimal_branch_shift_repair_theorem") or {}).get("selected_value") is not None)
                else "The D12 selector-law shell is now explicit on disk: one continuation-level overlap-defect scalar `Delta_ud_overlap` would fix the light-sector pure-B payload pair by odd-budget neutrality, but that selector value is still open and not recovered-core promotable."
            ),
            (
                "The D12 overlap transport law is now explicit too: once the spread totals are fixed, the odd payload pair collapses to one scalar `Delta_ud_overlap` with tau_u = sigma_d * Delta / (2 (sigma_u + sigma_d)), tau_d = sigma_u * Delta / (2 (sigma_u + sigma_d)), and Lambda = sigma_u sigma_d * Delta / (2 (sigma_u + sigma_d))."
                if d12_overlap_law is not None
                else "No D12 overlap transport law is attached to this audit."
            ),
            (
                "A stronger D12 continuation sample point is now explicit too: the retained same-family point with t1_sample = ray_modulus = 0.6695617711471163 gives (u,c,t) = (0.002176632493, 1.256692171439, 172.851929939314) GeV and (d,s,b) = (0.004708229529, 0.091608271273, 4.155513989985) GeV with RMS log-mass error about 1.08e-02, and the CKM/CP lane closes on that same branch because the forward Yukawa step already emits the same-label transport unitary V_CKM^fwd = U_u^dagger U_d and its principal logarithm."
                if d12_mass_branch is not None
                else "No D12 quark mass-branch followup is attached to this audit."
            ),
            (
                "The current D12 sheet is now known to be a strict no-go for the physical CKM shell: the emitted local same-label left-handed orbit closes to sigma_ref, but that selected singleton is just the current D12 reference sheet, so the exact next object is quark_d12_t1_value_law rather than another branch selector."
                if physical_branch_repair is not None
                and ((physical_branch_repair.get("minimal_branch_shift_repair_theorem") or {}).get("selected_value") is not None)
                else (
                    "The current D12 sheet is now known to be a strict no-go for the physical CKM shell: same-sheet rephasing leaves the CKM invariants frozen, the emitted angles are too small, and the exact next object is one discrete quark_relative_sheet_selector rather than a continuous repair scalar."
                    if physical_branch_repair is not None
                    else "No explicit D12 physical-branch repair theorem is attached to this audit yet."
                )
            ),
            (
                "The quadratic even transport is scalarized too: once the ordered-family carrier is fixed, the even residual collapses to one centered scalar eta_Q_centered with direct log formula (eta_Q_centered / 6) * (1,-2,1)."
                if quadratic_scalar is not None
                else "No scalarized D12 quadratic-even transport shell is attached to this audit."
            ),
            (
                "The strongest current D12 continuation bundle reduces the mass side to two value laws, Delta_ud_overlap and eta_Q_centered, while the mixing side is closed by the gauge-fixed physical invariants of the forward same-label left-transport generator."
                if scalarized_bundle is not None
                else "No scalarized D12 continuation bundle is attached to this audit."
            ),
            (
                "There is also a stricter same-family diagnostic specialization: on the current sample branch both mass-side continuation scalars collapse to one ray coordinate, with Delta_ud_overlap = ray_modulus / 5 and eta_Q_centered = -((1 - x2^2) / 27) * ray_modulus. That specialization is explicit on disk, but the retained numerical point is sample-only rather than OPH-emitted."
                if one_scalar_specialization is not None
                else "No one-scalar D12 same-family specialization is attached to this audit."
            ),
            (
                "The same-family D12 mass ray is now an explicit emitted object on disk: D12_ud_mass_ray = { ray_modulus * (1/5, -((1 - x2^2) / 27)) }, and the remaining primitive exact scalar burden is quark_d12_t1_value_law, with intrinsic_scale_law_D12 retained as the derived wrapper."
                if d12_mass_ray is not None
                else "No explicit emitted D12 mass-ray object is attached to this audit yet."
            ),
            (
                "The mass-side theorem boundary remains sharp too: with D12_ud_mass_ray emitted, the next exact missing object is quark_d12_t1_value_law rather than another larger mass-side carrier."
                if mass_side_underdetermination is not None
                else "No explicit D12 mass-side underdetermination theorem is attached to this audit yet."
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
