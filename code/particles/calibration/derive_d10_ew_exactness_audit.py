#!/usr/bin/env python3
"""Audit the exactness status of the D10 electroweak family and repair surface."""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
REFERENCE_JSON = ROOT / "particles" / "data" / "particle_reference_values.json"
FAMILY_JSON = ROOT / "particles" / "runs" / "calibration" / "d10_ew_observable_family.json"
SOURCE_PAIR_JSON = ROOT / "particles" / "runs" / "calibration" / "d10_ew_source_transport_pair.json"
READOUT_JSON = ROOT / "particles" / "runs" / "calibration" / "d10_ew_source_transport_readout.json"
EXACT_CLOSURE_JSON = ROOT / "particles" / "runs" / "calibration" / "d10_ew_exact_closure_beyond_current_carrier.json"
EXACT_WZ_COORDINATE_JSON = ROOT / "particles" / "runs" / "calibration" / "d10_ew_exact_wz_coordinate_beyond_single_tree_identity.json"
AFFINE_GERM_JSON = ROOT / "particles" / "runs" / "calibration" / "d10_ew_fixed_eta_post_transport_affine_germ.json"
FIBERWISE_TREE_LAW_JSON = ROOT / "particles" / "runs" / "calibration" / "d10_ew_fiberwise_population_tree_law_beneath_single_tree_identity.json"
TAU2_OBSTRUCTION_JSON = ROOT / "particles" / "runs" / "calibration" / "d10_ew_tau2_current_carrier_obstruction.json"
EXACT_MASS_PAIR_CHART_JSON = ROOT / "particles" / "runs" / "calibration" / "d10_ew_exact_mass_pair_chart_current_carrier.json"
REPAIR_BRANCH_JSON = ROOT / "particles" / "runs" / "calibration" / "d10_ew_repair_branch_beyond_current_carrier.json"
TARGET_POINT_DIAGNOSTIC_JSON = ROOT / "particles" / "runs" / "calibration" / "d10_ew_repair_target_point_diagnostic.json"
W_ANCHOR_FACTORIZATION_JSON = ROOT / "particles" / "runs" / "calibration" / "d10_ew_w_anchor_neutral_shear_factorization_official_pdg_2025_update.json"
W_ANCHOR_BOX_JSON = ROOT / "particles" / "runs" / "calibration" / "d10_ew_w_anchor_neutral_shear_box_dominance.json"
TARGET_FREEZE_SPLIT_JSON = ROOT / "particles" / "runs" / "calibration" / "d10_ew_target_freeze_and_subobject_split.json"
MINIMAL_CONDITIONAL_JSON = ROOT / "particles" / "runs" / "calibration" / "d10_ew_minimal_conditional_theorem.json"
TARGET_EMITTER_JSON = ROOT / "particles" / "runs" / "calibration" / "d10_ew_target_emitter_candidate.json"
TARGET_FREE_REPAIR_JSON = ROOT / "particles" / "runs" / "calibration" / "d10_ew_target_free_repair_value_law.json"
FORWARD_TRANSMUTATION_JSON = ROOT / "particles" / "runs" / "calibration" / "d10_ew_forward_transmutation_certificate.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_exactness_audit.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _coherent_from_alpha_sin(alpha_em_inv: float, sin2w: float, v_value: float) -> dict[str, float]:
    alpha_sum = 1.0 / (alpha_em_inv * sin2w * (1.0 - sin2w))
    alpha_y = sin2w * alpha_sum
    alpha2 = (1.0 - sin2w) * alpha_sum
    return {
        "alphaY_prime": alpha_y,
        "alpha2_prime": alpha2,
        "MW_pole": v_value * math.sqrt(math.pi * alpha2),
        "MZ_pole": v_value * math.sqrt(math.pi * alpha_sum),
        "alpha_em_eff_inv": alpha_em_inv,
        "sin2w_eff": sin2w,
        "v_report": v_value,
    }


def _coherent_from_reference_wz(mw_ref: float, mz_ref: float, v_value: float) -> dict[str, float]:
    alpha2 = (mw_ref / v_value) ** 2 / math.pi
    alpha_sum = (mz_ref / v_value) ** 2 / math.pi
    alpha_y = alpha_sum - alpha2
    return {
        "alphaY_prime": alpha_y,
        "alpha2_prime": alpha2,
        "MW_pole": mw_ref,
        "MZ_pole": mz_ref,
        "alpha_em_eff_inv": (alpha_y + alpha2) / (alpha_y * alpha2),
        "sin2w_eff": alpha_y / (alpha_y + alpha2),
        "v_report": v_value,
    }


def _solve_bisect(
    fn,
    target: float,
    *,
    lower: float,
    upper: float,
    iterations: int = 120,
) -> float | None:
    f_lower = fn(lower) - target
    f_upper = fn(upper) - target
    if f_lower == 0.0:
        return lower
    if f_upper == 0.0:
        return upper
    if f_lower * f_upper > 0.0:
        return None
    lo = lower
    hi = upper
    for _ in range(iterations):
        mid = 0.5 * (lo + hi)
        f_mid = fn(mid) - target
        if f_mid == 0.0:
            return mid
        if f_lower * f_mid <= 0.0:
            hi = mid
            f_upper = f_mid
        else:
            lo = mid
            f_lower = f_mid
    return 0.5 * (lo + hi)


def main() -> int:
    parser = argparse.ArgumentParser(description="Build a D10 electroweak exactness audit artifact.")
    parser.add_argument("--family", default=str(FAMILY_JSON))
    parser.add_argument("--source-pair", default=str(SOURCE_PAIR_JSON))
    parser.add_argument("--readout", default=str(READOUT_JSON))
    parser.add_argument("--exact-closure", default=str(EXACT_CLOSURE_JSON))
    parser.add_argument("--exact-wz-coordinate", default=str(EXACT_WZ_COORDINATE_JSON))
    parser.add_argument("--affine-germ", default=str(AFFINE_GERM_JSON))
    parser.add_argument("--fiberwise-tree-law", default=str(FIBERWISE_TREE_LAW_JSON))
    parser.add_argument("--tau2-obstruction", default=str(TAU2_OBSTRUCTION_JSON))
    parser.add_argument("--exact-mass-pair-chart", default=str(EXACT_MASS_PAIR_CHART_JSON))
    parser.add_argument("--repair-branch", default=str(REPAIR_BRANCH_JSON))
    parser.add_argument("--target-point-diagnostic", default=str(TARGET_POINT_DIAGNOSTIC_JSON))
    parser.add_argument("--w-anchor-factorization", default=str(W_ANCHOR_FACTORIZATION_JSON))
    parser.add_argument("--w-anchor-box", default=str(W_ANCHOR_BOX_JSON))
    parser.add_argument("--target-freeze-split", default=str(TARGET_FREEZE_SPLIT_JSON))
    parser.add_argument("--minimal-conditional-promotion", default=str(MINIMAL_CONDITIONAL_JSON))
    parser.add_argument("--target-emitter", default=str(TARGET_EMITTER_JSON))
    parser.add_argument("--target-free-repair", default=str(TARGET_FREE_REPAIR_JSON))
    parser.add_argument("--forward-transmutation", default=str(FORWARD_TRANSMUTATION_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    references = json.loads(REFERENCE_JSON.read_text(encoding="utf-8"))["entries"]
    family = json.loads(Path(args.family).read_text(encoding="utf-8"))
    source_pair = json.loads(Path(args.source_pair).read_text(encoding="utf-8"))
    readout = json.loads(Path(args.readout).read_text(encoding="utf-8"))
    exact_closure_path = Path(args.exact_closure)
    exact_closure = json.loads(exact_closure_path.read_text(encoding="utf-8")) if exact_closure_path.exists() else None
    exact_wz_coordinate_path = Path(args.exact_wz_coordinate)
    exact_wz_coordinate = (
        json.loads(exact_wz_coordinate_path.read_text(encoding="utf-8"))
        if exact_wz_coordinate_path.exists()
        else None
    )
    affine_germ_path = Path(args.affine_germ)
    affine_germ = json.loads(affine_germ_path.read_text(encoding="utf-8")) if affine_germ_path.exists() else None
    fiberwise_tree_law_path = Path(args.fiberwise_tree_law)
    fiberwise_tree_law = (
        json.loads(fiberwise_tree_law_path.read_text(encoding="utf-8"))
        if fiberwise_tree_law_path.exists()
        else None
    )
    tau2_obstruction_path = Path(args.tau2_obstruction)
    tau2_obstruction = (
        json.loads(tau2_obstruction_path.read_text(encoding="utf-8"))
        if tau2_obstruction_path.exists()
        else None
    )
    exact_mass_pair_chart_path = Path(args.exact_mass_pair_chart)
    exact_mass_pair_chart = (
        json.loads(exact_mass_pair_chart_path.read_text(encoding="utf-8"))
        if exact_mass_pair_chart_path.exists()
        else None
    )
    repair_branch_path = Path(args.repair_branch)
    repair_branch = json.loads(repair_branch_path.read_text(encoding="utf-8")) if repair_branch_path.exists() else None
    target_point_diagnostic_path = Path(args.target_point_diagnostic)
    target_point_diagnostic = (
        json.loads(target_point_diagnostic_path.read_text(encoding="utf-8"))
        if target_point_diagnostic_path.exists()
        else None
    )
    w_anchor_factorization_path = Path(args.w_anchor_factorization)
    w_anchor_factorization = (
        json.loads(w_anchor_factorization_path.read_text(encoding="utf-8"))
        if w_anchor_factorization_path.exists()
        else None
    )
    w_anchor_box_path = Path(args.w_anchor_box)
    w_anchor_box = (
        json.loads(w_anchor_box_path.read_text(encoding="utf-8"))
        if w_anchor_box_path.exists()
        else None
    )
    target_freeze_split_path = Path(args.target_freeze_split)
    target_freeze_split = (
        json.loads(target_freeze_split_path.read_text(encoding="utf-8"))
        if target_freeze_split_path.exists()
        else None
    )
    minimal_conditional_path = Path(args.minimal_conditional_promotion)
    minimal_conditional = (
        json.loads(minimal_conditional_path.read_text(encoding="utf-8"))
        if minimal_conditional_path.exists()
        else None
    )
    target_emitter_path = Path(args.target_emitter)
    target_emitter = (
        json.loads(target_emitter_path.read_text(encoding="utf-8"))
        if target_emitter_path.exists()
        else None
    )
    target_free_repair_path = Path(args.target_free_repair)
    target_free_repair = (
        json.loads(target_free_repair_path.read_text(encoding="utf-8"))
        if target_free_repair_path.exists()
        else None
    )
    forward_transmutation_path = Path(args.forward_transmutation)
    forward_transmutation = (
        json.loads(forward_transmutation_path.read_text(encoding="utf-8"))
        if forward_transmutation_path.exists()
        else None
    )

    reported = dict(family["reported_outputs"])
    public_quintet = dict(readout.get("public_emitted_quintet") or {})
    current_compact_quintet = dict(readout.get("current_compact_emitted_quintet") or {})
    v_value = float(reported["v"])
    mw_ref = float(references["w_boson"]["value_gev"])
    mz_ref = float(references["z_boson"]["value_gev"])
    alpha_y0 = float(source_pair["source_pair"]["alphaY_mz"])
    alpha2_0 = float(source_pair["source_pair"]["alpha2_mz"])
    eta_fixed = float(source_pair["compact_hypercharge_only_mass_slice"]["eta_EW"])
    reference_wz_slice = _coherent_from_reference_wz(mw_ref, mz_ref, v_value)
    exact_mass_ratio_sin2 = 1.0 - (mw_ref / mz_ref) ** 2
    running_alpha = float(reported["alpha_em_inv_mz"])
    running_sin2 = float(reported["sin2w_mz"])
    alpha_sin_quintet = _coherent_from_alpha_sin(running_alpha, running_sin2, v_value)
    sigma_domain_lower = -1.0 + abs(eta_fixed) + 1.0e-12
    sigma_probe_upper = 10.0

    def _alpha_inv_on_fixed_eta(sigma_ew: float) -> float:
        alpha_y = alpha_y0 * (1.0 + sigma_ew - eta_fixed)
        alpha2 = alpha2_0 * (1.0 + sigma_ew + eta_fixed)
        return (alpha_y + alpha2) / (alpha_y * alpha2)

    def _sin2_on_fixed_eta(sigma_ew: float) -> float:
        alpha_y = alpha_y0 * (1.0 + sigma_ew - eta_fixed)
        alpha2 = alpha2_0 * (1.0 + sigma_ew + eta_fixed)
        return alpha_y / (alpha_y + alpha2)

    sigma_from_mw = (mw_ref / v_value) ** 2 / (math.pi * alpha2_0) - 1.0 - eta_fixed
    sigma_from_mz = (
        (mz_ref / v_value) ** 2 / (math.pi * (alpha_y0 + alpha2_0))
        - 1.0
        - eta_fixed * ((alpha2_0 - alpha_y0) / (alpha_y0 + alpha2_0))
    )
    sigma_from_running_alpha = _solve_bisect(
        _alpha_inv_on_fixed_eta,
        running_alpha,
        lower=sigma_domain_lower,
        upper=sigma_probe_upper,
    )
    sin2_low = _sin2_on_fixed_eta(sigma_domain_lower)
    sin2_high = _sin2_on_fixed_eta(sigma_probe_upper)
    sigma_from_running_sin2 = _solve_bisect(
        _sin2_on_fixed_eta,
        running_sin2,
        lower=sigma_domain_lower,
        upper=sigma_probe_upper,
    )
    sin2_min = min(sin2_low, sin2_high)
    sin2_max = max(sin2_low, sin2_high)
    target_free_closed = target_free_repair is not None and target_free_repair.get("status") == "closed"
    target_free_validation = (
        dict(target_free_repair.get("compare_only_validation_against_frozen_surface") or {})
        if target_free_repair is not None
        else {}
    )
    frozen_repair_quintet = (
        dict(w_anchor_factorization.get("coherent_repaired_quintet") or {})
        if w_anchor_factorization is not None
        else {}
    )

    artifact = {
        "artifact": "oph_d10_ew_exactness_audit",
        "generated_utc": _timestamp(),
        "scope": "current_family_only",
        "reference_inputs": {
            "MW_pole_gev": mw_ref,
            "MZ_pole_gev": mz_ref,
            "alpha_em_eff_inv_running_surface": running_alpha,
            "sin2w_eff_running_surface": running_sin2,
            "v_report_gev": v_value,
        },
        "current_running_family_quintet": {
            "MW_pole": float(reported["m_w_run"]),
            "MZ_pole": float(reported["m_z_run"]),
            "alpha_em_eff_inv": running_alpha,
            "sin2w_eff": running_sin2,
            "v_report": v_value,
        },
        "reference_wz_audit_slice": {
            "status": "audit_only_inverse_slice",
            **reference_wz_slice,
        },
        "coherent_alpha_sin_candidate": {
            "status": "derived_from_running_alpha_and_running_sin2",
            **alpha_sin_quintet,
        },
        "tree_level_identity_audit": {
            "required_identity": "sin2w_eff = 1 - (MW_pole / MZ_pole)^2 on one coherent family",
            "mass_ratio_sin2_from_reference_WZ": exact_mass_ratio_sin2,
            "running_surface_sin2": running_sin2,
            "identity_residual_on_mixed_reference_surface": running_sin2 - exact_mass_ratio_sin2,
        },
        "reference_wz_audit_slice_residual_vs_running_alpha_sin": {
            "delta_alpha_em_eff_inv": reference_wz_slice["alpha_em_eff_inv"] - running_alpha,
            "delta_sin2w_eff": reference_wz_slice["sin2w_eff"] - running_sin2,
        },
        "alpha_sin_candidate_residual_vs_reference_masses": {
            "delta_MW_pole": alpha_sin_quintet["MW_pole"] - mw_ref,
            "delta_MZ_pole": alpha_sin_quintet["MZ_pole"] - mz_ref,
        },
        "fixed_eta_single_sigma_audit": {
            "family": {
                "eta_EW_fixed": eta_fixed,
                "tau_Y_formula": "sigma_EW - eta_EW_fixed",
                "tau_2_formula": "sigma_EW + eta_EW_fixed",
                "domain_lower_bound_sigma_EW": sigma_domain_lower,
            },
            "sigma_required_by_reference_mass_pair": {
                "sigma_from_reference_MW": sigma_from_mw,
                "sigma_from_reference_MZ": sigma_from_mz,
                "mass_pair_sigma_split": sigma_from_mw - sigma_from_mz,
                "mass_pair_sigma_midpoint": 0.5 * (sigma_from_mw + sigma_from_mz),
            },
            "sigma_required_by_running_readouts": {
                "sigma_from_running_alpha_em_eff_inv": sigma_from_running_alpha,
                "sigma_from_running_sin2w_eff": sigma_from_running_sin2,
                "running_sin2w_eff_admissible_on_fixed_eta_family": sigma_from_running_sin2 is not None,
                "sin2w_eff_range_on_fixed_eta_family": {
                    "sigma_low": sigma_domain_lower,
                    "sigma_high": sigma_probe_upper,
                    "sin2w_eff_low": sin2_low,
                    "sin2w_eff_high": sin2_high,
                    "sin2w_eff_min": sin2_min,
                    "sin2w_eff_max": sin2_max,
                },
            },
            "verdict": {
                "mass_pair_nearly_coherent": abs(sigma_from_mw - sigma_from_mz) < 5.0e-4,
                "running_alpha_conflicts_with_mass_pair": (
                    sigma_from_running_alpha is not None
                    and abs(sigma_from_running_alpha - 0.5 * (sigma_from_mw + sigma_from_mz)) > 1.0e-2
                ),
                "running_sin2_forces_family_escape": sigma_from_running_sin2 is None,
            },
        },
        "current_carrier_closure_summary": {
            "exact_current_carrier_chart_closed": (
                exact_mass_pair_chart.get("proof_status")
                == "exact_mass_pair_chart_closed_and_current_selector_pullback_has_unique_zero_at_current_point"
                if exact_mass_pair_chart is not None
                else False
            ),
            "current_carrier_exact_mass_pair": {
                "MW_pole": float(current_compact_quintet.get("MW_pole", reported["m_w_run"])),
                "MZ_pole": float(current_compact_quintet.get("MZ_pole", reported["m_z_run"])),
            },
            "current_carrier_obstruction_certified": (
                tau2_obstruction.get("proof_status")
                == "no_single_tau2_on_closed_current_carrier_can_hit_exact_W_and_exact_Z"
                if tau2_obstruction is not None
                else False
            ),
            "current_carrier_builder_local_frontier": (
                readout.get("current_carrier_builder_local_frontier")
                or (
                    exact_mass_pair_chart.get("next_single_residual_object")
                    if exact_mass_pair_chart and exact_mass_pair_chart.get("next_single_residual_object") is not None
                    else "EWExactMassPairSelector_D10"
                )
            ),
            "broader_honest_repair_frontier": readout.get("broader_honest_repair_frontier"),
            "final_wave_consolidation_verdict": (
                "The present selected/current carrier closes its own exact W/Z chart at its local pair, while the public target-free source-only theorem emits the coherent electroweak quintet from the D10 source basis alone."
                if target_free_closed
                else
                "The present selected/current carrier closes its own exact W/Z chart at its local pair, while the public frozen-target repair surface emits the authoritative exact W/Z pair from one coherent repaired coupling pair."
                if readout.get("w_anchor_neutral_shear_factorization_status") == "closed_freeze_once_coherent_repair_law"
                else
                "The present selected/current carrier closes its own exact W/Z chart, but exact PDG W/Z "
                "does not honestly follow on that carrier from the current emitted selector package."
            ),
        },
        "d10_repair_branch_beyond_current_carrier": None if repair_branch is None else {
            "artifact": repair_branch.get("artifact"),
            "status": repair_branch.get("status"),
            "proof_status": repair_branch.get("proof_status"),
            "object_id": repair_branch.get("object_id"),
            "replaces_builder_local_frontier": repair_branch.get("replaces_builder_local_frontier"),
            "stronger_residual_object": repair_branch.get("stronger_residual_object"),
            "required_closure_kind": repair_branch.get("required_closure_kind"),
            "operative_primitive": repair_branch.get("operative_primitive"),
        },
        "repair_target_point_diagnostic": None if target_point_diagnostic is None else {
            "artifact": target_point_diagnostic.get("artifact"),
            "status": target_point_diagnostic.get("status"),
            "spec_id": target_point_diagnostic.get("spec_id"),
            "MW_pole_target_gev": target_point_diagnostic.get("MW_pole_target_gev"),
            "MZ_pole_target_gev": target_point_diagnostic.get("MZ_pole_target_gev"),
            "tau2_tree_exact_target": target_point_diagnostic.get("tau2_tree_exact_target"),
            "delta_n_tree_exact_target": target_point_diagnostic.get("delta_n_tree_exact_target"),
            "delta_alpha2_tree_target": target_point_diagnostic.get("delta_alpha2_tree_target"),
            "delta_alphaY_tree_target": target_point_diagnostic.get("delta_alphaY_tree_target"),
        },
        "w_anchor_neutral_shear_factorization": None if w_anchor_factorization is None else {
            "artifact": w_anchor_factorization.get("artifact"),
            "status": w_anchor_factorization.get("status"),
            "exact_missing_law": w_anchor_factorization.get("exact_missing_law"),
            "new_smaller_primitive": w_anchor_factorization.get("new_smaller_primitive"),
            "central_target_point": w_anchor_factorization.get("central_target_point"),
            "coherent_repaired_quintet": w_anchor_factorization.get("coherent_repaired_quintet"),
            "conclusion": w_anchor_factorization.get("conclusion"),
        },
        "w_anchor_neutral_shear_box_dominance": None if w_anchor_box is None else {
            "artifact": w_anchor_box.get("artifact"),
            "status": w_anchor_box.get("status"),
            "bounds": w_anchor_box.get("bounds"),
            "verdict": w_anchor_box.get("verdict"),
        },
        "target_freeze_and_subobject_split": None if target_freeze_split is None else {
            "artifact": target_freeze_split.get("artifact"),
            "status": target_freeze_split.get("status"),
            "target_freeze_required": target_freeze_split.get("target_freeze_required"),
            "subobject_split": target_freeze_split.get("subobject_split"),
        },
        "target_free_source_only_underdetermination": None if minimal_conditional is None else {
            "artifact": minimal_conditional.get("artifact"),
            "status": (
                "superseded_by_target_free_repair_theorem"
                if target_free_closed
                else minimal_conditional.get("status")
            ),
            "superseded_by": minimal_conditional.get("superseded_by"),
            "unconditional_theorem": minimal_conditional.get("unconditional_theorem"),
            "conditional_principle": minimal_conditional.get("conditional_principle"),
            "conditional_theorem": minimal_conditional.get("conditional_theorem"),
            "n_c_3_specialization": minimal_conditional.get("n_c_3_specialization"),
        },
        "target_free_source_only_candidate": None if target_emitter is None else {
            "artifact": target_emitter.get("artifact"),
            "status": "promoted" if target_free_closed else target_emitter.get("status"),
            "promoted_to": target_emitter.get("promoted_to"),
            "object_id": target_emitter.get("object_id"),
            "emitter_scalar": target_emitter.get("emitter_scalar"),
            "target_emitter_law": target_emitter.get("target_emitter_law"),
            "coherent_emitted_quintet": target_emitter.get("coherent_emitted_quintet"),
            "comparison_to_frozen_local_reference_surface": target_emitter.get("comparison_to_frozen_local_reference_surface"),
        },
        "target_free_source_only_repair_theorem": None if target_free_repair is None else {
            "artifact": target_free_repair.get("artifact"),
            "status": target_free_repair.get("status"),
            "object_id": target_free_repair.get("object_id"),
            "theorem": target_free_repair.get("theorem"),
            "repair_chart": target_free_repair.get("repair_chart"),
            "repaired_couplings": target_free_repair.get("repaired_couplings"),
            "coherent_emitted_quintet": target_free_repair.get("coherent_emitted_quintet"),
            "compare_only_validation_against_frozen_surface": target_free_repair.get("compare_only_validation_against_frozen_surface"),
        },
        "forward_transmutation_certificate": None if forward_transmutation is None else {
            "artifact": forward_transmutation.get("artifact"),
            "status": forward_transmutation.get("status"),
            "object_id": forward_transmutation.get("object_id"),
            "notation_split": forward_transmutation.get("notation_split"),
            "forward_core_solution": forward_transmutation.get("forward_core_solution"),
            "source_only_reconstruction": forward_transmutation.get("source_only_reconstruction"),
            "forward_checks": forward_transmutation.get("forward_checks"),
        },
        "exact_closure_beyond_current_carrier": None if exact_closure is None else {
            "artifact": exact_closure.get("artifact"),
            "status": exact_closure.get("status"),
            "proof_status": exact_closure.get("proof_status"),
            "exactness_surface_kind": exact_closure.get("exactness_surface_kind"),
            "completion_kind": exact_closure.get("completion_kind"),
            "stronger_residual_object": exact_closure.get("stronger_residual_object"),
            "exact_outputs": exact_closure.get("exact_outputs"),
        },
        "exact_wz_coordinate_beyond_single_tree_identity": None if exact_wz_coordinate is None else {
            "artifact": exact_wz_coordinate.get("artifact"),
            "status": exact_wz_coordinate.get("status"),
            "proof_status": exact_wz_coordinate.get("proof_status"),
            "depends_on_object": exact_wz_coordinate.get("depends_on_object"),
            "coordinate_symbol": exact_wz_coordinate.get("coordinate_symbol"),
            "next_residual_object_if_open": exact_wz_coordinate.get("next_residual_object_if_open"),
        },
        "fiberwise_population_tree_law_beneath_single_tree_identity": None if fiberwise_tree_law is None else {
            "artifact": fiberwise_tree_law.get("artifact"),
            "status": fiberwise_tree_law.get("status"),
            "proof_status": fiberwise_tree_law.get("proof_status"),
            "strictly_smaller_than": fiberwise_tree_law.get("strictly_smaller_than"),
            "next_single_residual_object": fiberwise_tree_law.get("next_single_residual_object"),
        },
        "fixed_eta_post_transport_affine_germ": None if affine_germ is None else {
            "artifact": affine_germ.get("artifact"),
            "status": affine_germ.get("status"),
            "proof_status": affine_germ.get("proof_status"),
            "strictly_smaller_than": affine_germ.get("strictly_smaller_than"),
            "diagnostic_only": affine_germ.get("diagnostic_only"),
            "next_residual_object": affine_germ.get("next_residual_object"),
        },
        "tau2_current_carrier_obstruction": None if tau2_obstruction is None else {
            "artifact": tau2_obstruction.get("artifact"),
            "status": tau2_obstruction.get("status"),
            "proof_status": tau2_obstruction.get("proof_status"),
            "strictly_smaller_than": tau2_obstruction.get("strictly_smaller_than"),
            "next_single_residual_object": tau2_obstruction.get("next_single_residual_object"),
        },
        "exact_mass_pair_chart_current_carrier": None if exact_mass_pair_chart is None else {
            "artifact": exact_mass_pair_chart.get("artifact"),
            "status": exact_mass_pair_chart.get("status"),
            "proof_status": exact_mass_pair_chart.get("proof_status"),
            "strictly_smaller_than": exact_mass_pair_chart.get("strictly_smaller_than"),
            "next_single_residual_object": exact_mass_pair_chart.get("next_single_residual_object"),
        },
        "freeze_once_coherent_repair_summary": None if w_anchor_factorization is None else {
            "status": w_anchor_factorization.get("status"),
            "repair_law_id": w_anchor_factorization.get("exact_missing_law", {}).get("object_id"),
            "historical_compare_only": target_free_closed,
            "frozen_surface_mass_pair": {
                "MW_pole": float(frozen_repair_quintet.get("MW_pole", public_quintet.get("MW_pole", reported["m_w_run"]))),
                "MZ_pole": float(frozen_repair_quintet.get("MZ_pole", public_quintet.get("MZ_pole", reported["m_z_run"]))),
            },
            "target_free_predictive_emission_closed": w_anchor_factorization.get("target_free_predictive_emission_closed"),
            "stricter_still_open_object": w_anchor_factorization.get("conclusion", {}).get("stricter_still_open_object"),
            "compare_only_validation_against_target_free_surface": (
                {
                    "delta_MW_gev": target_free_validation.get("delta_MW_gev"),
                    "delta_MZ_gev": target_free_validation.get("delta_MZ_gev"),
                    "delta_alpha_em_eff_inv": target_free_validation.get("delta_alpha_em_eff_inv"),
                    "delta_sin2w_eff": target_free_validation.get("delta_sin2w_eff"),
                }
                if target_free_closed
                else None
            ),
        },
        "smallest_exact_obstruction": (
            None
            if target_free_closed
            else
            "On one frozen authoritative D10 target pair, the coherent repair law is closed and emits exact W/Z on one family; the only stricter remaining D10 step is the target-free repair value law from P alone."
            if w_anchor_factorization and w_anchor_factorization.get("status") == "closed_freeze_once_coherent_repair_law"
            else
            "the exact two-coordinate current-carrier chart is closed, but the current selector still has a unique zero at the current point, so exact W/Z closure now needs one nonzero selector on that chart"
            if exact_mass_pair_chart and exact_mass_pair_chart.get("status") == "closed_smaller_primitive"
            else
            "current closed one-variable carrier moves W and Z with the same sign, so exact mass-pair closure needs one additional neutral-leg scalar beyond tau2_tree_exact"
            if tau2_obstruction and tau2_obstruction.get("status") == "closed_smaller_primitive"
            else
            "the selected D10 carrier point admits a closed split exactness law, and the fiberwise population tree law reduces exact W/Z closure to one remaining charged-leg scalar"
            if exact_wz_coordinate and exact_wz_coordinate.get("next_residual_object_if_open") == "tau2_tree_exact"
            else "the selected D10 carrier point admits a closed split exactness law, but the stronger single post-transport tree identity remains open"
            if exact_closure and exact_closure.get("status") == "closed"
            else "the reopened two-scalar D10 carrier is now canonically selected, but the selected current-carrier point still does not close the full electroweak target surface exactly"
        ),
        "active_builder_smallest_missing_object": (
            None
            if target_free_closed
            else
            w_anchor_factorization.get("conclusion", {}).get("stricter_still_open_object")
            if w_anchor_factorization and w_anchor_factorization.get("status") == "closed_freeze_once_coherent_repair_law"
            else
            exact_mass_pair_chart.get("next_single_residual_object")
            if exact_mass_pair_chart and exact_mass_pair_chart.get("next_single_residual_object") is not None
            else
            tau2_obstruction.get("next_single_residual_object")
            if tau2_obstruction and tau2_obstruction.get("next_single_residual_object") is not None
            else
            exact_wz_coordinate.get("next_residual_object_if_open")
            if exact_wz_coordinate and exact_wz_coordinate.get("next_residual_object_if_open") is not None
            else
            exact_closure.get("stronger_residual_object")
            if exact_closure and exact_closure.get("status") == "closed"
            else "EWTransportExactClosureBeyondCurrentCarrier_D10"
        ),
        "broader_honest_repair_frontier": (
            None
            if target_free_closed
            else
            w_anchor_factorization.get("conclusion", {}).get("stricter_still_open_object")
            if w_anchor_factorization and w_anchor_factorization.get("status") == "closed_freeze_once_coherent_repair_law"
            else
            repair_branch.get("object_id")
            if repair_branch and repair_branch.get("object_id") is not None
            else "D10RepairBranchBeyondCurrentCarrier"
        ),
        "exact_pdg_wz_frontier": (
            target_free_repair.get("object_id")
            if target_free_closed
            else
            w_anchor_factorization.get("conclusion", {}).get("stricter_still_open_object")
            if w_anchor_factorization and w_anchor_factorization.get("status") == "closed_freeze_once_coherent_repair_law"
            else
            repair_branch.get("object_id")
            if repair_branch and repair_branch.get("object_id") is not None
            else "D10RepairBranchBeyondCurrentCarrier"
        ),
        "smallest_constructive_missing_object": (
            None
            if target_free_closed
            else
            w_anchor_factorization.get("conclusion", {}).get("stricter_still_open_object")
            if w_anchor_factorization and w_anchor_factorization.get("status") == "closed_freeze_once_coherent_repair_law"
            else
            repair_branch.get("object_id")
            if repair_branch and repair_branch.get("object_id") is not None
            else "D10RepairBranchBeyondCurrentCarrier"
        ),
        "exactness_precondition": {
            "frozen_d10_target_spec_required": True,
            "reason": "exact W/Z closure is not well-posed until one authoritative calibration target pair is frozen at declared precision",
            "current_target_point_diagnostic_kind": (
                target_point_diagnostic.get("spec_id")
                if target_point_diagnostic is not None
                else None
            ),
        },
        "notes": [
            "This audit computes a purely inverse reference-W/Z slice for diagnosis only; that slice is not a predictive artifact and is not emitted by the live D10 source-pair builder.",
            "Any coherent family matching those reference W/Z values forces sin2(theta_W) = 1 - (MW/MZ)^2, so the current mixed reference surface cannot be hit exactly by one simple tree-level family.",
            "On the fixed-eta current family, the reference W and Z masses already point to nearly the same sigma_EW, but that sigma_EW is incompatible with the running alpha_em^-1 target and the running sin^2(theta_W) target is not even admissible on that one-parameter slice.",
            (
                "The carrier-level selector, split exact-closure law, fiberwise tree law, exact current-carrier mass chart, freeze-once coherent repair law, and target-free source-only repair theorem are all closed. The W/Z mass lane comes from the target-free theorem, the frozen authoritative repair surface is compare-only validation, and the fine-structure constant is read separately on the Ward-projected U(1)_Q transport family."
                if target_free_closed
                else
                "The carrier-level selector, split exact-closure law, fiberwise tree law, exact current-carrier mass chart, and freeze-once coherent repair law are all closed. Exact public W/Z now come from the frozen authoritative repair surface, while the stricter still-open D10 step is the target-free repair value law from P alone."
                if w_anchor_factorization and w_anchor_factorization.get("status") == "closed_freeze_once_coherent_repair_law"
                else
                "The carrier-level selector, split exact-closure law, fiberwise tree law, and exact current-carrier mass chart are all closed on the existing carrier; the remaining builder-local D10 mass-side object is the selector `EWExactMassPairSelector_D10` on that chart, but the broader exact-PDG frontier is the repair branch beyond the current carrier."
                if exact_mass_pair_chart and exact_mass_pair_chart.get("status") == "closed_smaller_primitive"
                else
                "The carrier-level selector, split exact-closure law, and fiberwise tree law are closed on the existing carrier, and the live obstruction is narrower: a single tau2 move cannot hit exact W and exact Z simultaneously on the current carrier, so the next mass-side scalar is delta_n_tree_exact."
                if tau2_obstruction and tau2_obstruction.get("status") == "closed_smaller_primitive"
                else
                "The carrier-level selector and the split exact closure law are both closed on the existing reopened two-scalar carrier; the fiberwise J_pop_EW minimizer removes the placeholder tree law, leaving only tau2_tree_exact as the remaining D10 scalar."
                if exact_wz_coordinate and exact_wz_coordinate.get("next_residual_object_if_open") == "tau2_tree_exact"
                else "The carrier-level selector and the split exact closure law are both closed on the existing reopened two-scalar carrier; the stronger remaining question is whether one unsplit post-transport tree identity can replace that split readout."
                if exact_closure and exact_closure.get("status") == "closed"
                else "The carrier-level selector is now closed on the existing reopened two-scalar carrier, so the exact local question is whether exact electroweak closure requires a new invariant beyond the current selected carrier point."
            ),
            (
                "The broader D10 geometry is the repair branch beyond the present current carrier. That branch remains explicit, the W/Z mass lane uses the target-free source-only repair law, and the electromagnetic readout sits on the Ward-projected U(1)_Q transport family."
                if target_free_closed
                else "The broader D10 geometry is still the repair branch beyond the present current carrier, but on one frozen authoritative target pair the coherent repair law is already explicit and exact."
            ),
            (
                "The forward transmutation certificate now records the non-circular P -> alpha_U -> t map explicitly and separates the calibration beta_ratio_EW from the transmutation counting factor beta_transmutation_EW = N_c + 1."
                if forward_transmutation is not None
                else "No explicit forward transmutation certificate is attached to this audit."
            ),
            "The current-carrier exact pair remains a distinct object from the public repaired pair; they should not be conflated in audits or status surfaces.",
            (
                "No stricter D10 electroweak mass-side theorem object remains open on the active Phase II calibration surface."
                if target_free_closed
                else "The only stricter remaining D10 theorem step is a target-free law that emits the same nonzero repair directly from P alone with no frozen authoritative W/Z input."
            ),
            (
                "The source-only underdetermination theorem, the minimal conditional route through ColorBalancedQuadraticRepairDescent_D10, and the former candidate EWTargetEmitter_D10 all sit beneath the target-free mass-side theorem; they do not define the public electromagnetic readout."
                if target_free_closed and (minimal_conditional is not None or target_emitter is not None)
                else
                "That target-free step is now split sharply on disk: the current source-only D10 corpus underdetermines the repair coefficients, the smallest honest conditional theorem route goes through ColorBalancedQuadraticRepairDescent_D10, and the strongest current source-only candidate is EWTargetEmitter_D10."
                if (minimal_conditional is not None or target_emitter is not None)
                else "No sharper target-free split is attached to this audit."
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
