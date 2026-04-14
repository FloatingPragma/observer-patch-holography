#!/usr/bin/env python3
"""Map the selected D10 carrier point and promoted repair laws into the electroweak readout.

Chain role: turn the chosen D10 two-scalar transport point into the emitted
`W/Z` mass pair and the associated electroweak quintet diagnostics.

Mathematics: factorized shared-scalar readout on top of the baseline D10
quintet, with explicit separation between the selected carrier point, the
historical freeze-once coherent repair surface, and the promoted target-free
source-only repair theorem.

OPH-derived inputs: the D10 observable family plus the selected two-scalar
carrier from the source transport pair.

Output: the public electroweak quintet plus the residual data that separates
the closed current-carrier chart from the historical frozen-target validation
surface.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_FAMILY = ROOT / "particles" / "runs" / "calibration" / "d10_ew_observable_family.json"
DEFAULT_SOURCE_PAIR = ROOT / "particles" / "runs" / "calibration" / "d10_ew_source_transport_pair.json"
DEFAULT_EXACT_CLOSURE = ROOT / "particles" / "runs" / "calibration" / "d10_ew_exact_closure_beyond_current_carrier.json"
DEFAULT_EXACT_WZ_COORDINATE = ROOT / "particles" / "runs" / "calibration" / "d10_ew_exact_wz_coordinate_beyond_single_tree_identity.json"
DEFAULT_TAU2_OBSTRUCTION = ROOT / "particles" / "runs" / "calibration" / "d10_ew_tau2_current_carrier_obstruction.json"
DEFAULT_EXACT_MASS_PAIR_CHART = ROOT / "particles" / "runs" / "calibration" / "d10_ew_exact_mass_pair_chart_current_carrier.json"
DEFAULT_W_ANCHOR_FACTORIZATION = ROOT / "particles" / "runs" / "calibration" / "d10_ew_w_anchor_neutral_shear_factorization_official_pdg_2025_update.json"
DEFAULT_MINIMAL_CONDITIONAL = ROOT / "particles" / "runs" / "calibration" / "d10_ew_minimal_conditional_theorem.json"
DEFAULT_TARGET_EMITTER = ROOT / "particles" / "runs" / "calibration" / "d10_ew_target_emitter_candidate.json"
DEFAULT_TARGET_FREE_REPAIR = ROOT / "particles" / "runs" / "calibration" / "d10_ew_target_free_repair_value_law.json"
DEFAULT_FORWARD_TRANSMUTATION = ROOT / "particles" / "runs" / "calibration" / "d10_ew_forward_transmutation_certificate.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_source_transport_readout.json"


BUILDER_LOCAL_FRONTIER = "EWExactMassPairSelector_D10"
GLOBAL_REPAIR_FRONTIER = "D10RepairBranchBeyondCurrentCarrier"
TARGET_FREE_REPAIR_FRONTIER = "EWTargetFreeRepairValueLaw_D10"
WARD_PROJECTED_TRANSPORT_THEOREM = "WardProjectedU1QTransportLaw_D10"
THOMSON_ALPHA_INV = 137.035999177


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(
    family: dict,
    source_pair: dict,
    exact_closure: dict | None = None,
    exact_wz_coordinate: dict | None = None,
    tau2_obstruction: dict | None = None,
    exact_mass_pair_chart: dict | None = None,
    w_anchor_factorization: dict | None = None,
    minimal_conditional: dict | None = None,
    target_emitter: dict | None = None,
    target_free_repair: dict | None = None,
    forward_transmutation: dict | None = None,
) -> dict:
    reported = dict(family["reported_outputs"])
    source_slots = dict(source_pair["source_pair"])
    seed_trial = dict(source_pair["first_nonzero_oph_seed_trial"])
    transport_quintet = dict(seed_trial["coherent_output_quintet"])
    compact_mass_slice = dict(source_pair["compact_hypercharge_only_mass_slice"])
    compact_mass_quintet = dict(compact_mass_slice["coherent_output_quintet"])
    alpha_y0 = float(source_slots["alphaY_mz"])
    alpha2_0 = float(source_slots["alpha2_mz"])
    sin2_theta_w0 = alpha_y0 / (alpha_y0 + alpha2_0)
    cos2_theta_w0 = (alpha2_0 - alpha_y0) / (alpha_y0 + alpha2_0)
    base_quintet = {
        "MW_pole": float(reported["m_w_run"]),
        "MZ_pole": float(reported["m_z_run"]),
        "alpha_em_eff_inv": float(reported["alpha_em_inv_mz"]),
        "sin2w_eff": float(reported["sin2w_mz"]),
        "v_report": float(reported["v"]),
    }
    selected_point = dict(source_pair.get("predictive_population_point", {}))
    selected_population_closed = bool(source_pair.get("predictive_population_closed", False))
    if selected_population_closed and compact_mass_quintet:
        public_quintet = {
            "MW_pole": float(compact_mass_quintet["MW_pole"]),
            "MZ_pole": float(compact_mass_quintet["MZ_pole"]),
            "alpha_em_eff_inv": float(compact_mass_quintet["alpha_em_eff_inv"]),
            "sin2w_eff": float(compact_mass_quintet["sin2w_eff"]),
            "v_report": float(compact_mass_quintet["v_report"]),
        }
    else:
        public_quintet = {
            "MW_pole": float(transport_quintet["MW_pole"]),
            "MZ_pole": float(transport_quintet["MZ_pole"]),
            "alpha_em_eff_inv": float(transport_quintet["alpha_em_eff_inv"]),
            "sin2w_eff": float(transport_quintet["sin2w_eff"]),
            "v_report": float(transport_quintet["v_report"]),
        }
    shared_scalar_values = {
        "delta_alpha": public_quintet["alpha_em_eff_inv"] / base_quintet["alpha_em_eff_inv"] - 1.0,
        "delta_kappa": public_quintet["sin2w_eff"] / base_quintet["sin2w_eff"] - 1.0,
        "delta_rho": public_quintet["MZ_pole"] / base_quintet["MZ_pole"] - 1.0,
        "delta_rW": (
            public_quintet["MW_pole"] / base_quintet["MW_pole"]
            - public_quintet["MZ_pole"] / base_quintet["MZ_pole"]
        ),
    }
    compact_shared_scalar_values = {
        "delta_alpha": compact_mass_quintet["alpha_em_eff_inv"] / base_quintet["alpha_em_eff_inv"] - 1.0,
        "delta_kappa": compact_mass_quintet["sin2w_eff"] / base_quintet["sin2w_eff"] - 1.0,
        "delta_rho": compact_mass_quintet["MZ_pole"] / base_quintet["MZ_pole"] - 1.0,
        "delta_rW": (
            compact_mass_quintet["MW_pole"] / base_quintet["MW_pole"]
            - compact_mass_quintet["MZ_pole"] / base_quintet["MZ_pole"]
        ),
    }
    compact_alpha_y = alpha_y0 * (1.0 + float(compact_mass_slice["tau_Y"]))
    compact_alpha2 = alpha2_0 * (1.0 + float(compact_mass_slice["tau_2"]))
    exact_closure = exact_closure or {}
    exact_closure_closed = exact_closure.get("status") == "closed"
    exact_outputs = dict(exact_closure.get("exact_outputs", {})) if exact_closure_closed else None
    exact_wz_coordinate = exact_wz_coordinate or {}
    tau2_obstruction = tau2_obstruction or {}
    exact_mass_pair_chart = exact_mass_pair_chart or {}
    w_anchor_factorization = w_anchor_factorization or {}
    minimal_conditional = minimal_conditional or {}
    target_emitter = target_emitter or {}
    target_free_repair = target_free_repair or {}
    forward_transmutation = forward_transmutation or {}
    freeze_once_repair_closed = (
        w_anchor_factorization.get("status") == "closed_freeze_once_coherent_repair_law"
    )
    target_free_repair_closed = target_free_repair.get("status") == "closed"
    factorization_point = dict(w_anchor_factorization.get("central_target_point") or {})
    factorized_repair_package = dict(w_anchor_factorization.get("repair_law", {}).get("coherent_repair_package", {}))
    factorized_repair_chart = dict(w_anchor_factorization.get("repair_law", {}).get("chart_coordinates", {}))
    factorized_repair_quintet = dict(w_anchor_factorization.get("coherent_repaired_quintet") or {})
    target_free_chart = dict(target_free_repair.get("repair_chart") or {})
    target_free_couplings = dict(target_free_repair.get("repaired_couplings") or {})
    target_free_quintet = dict(target_free_repair.get("coherent_emitted_quintet") or {})

    if exact_mass_pair_chart.get("status") == "closed_smaller_primitive":
        current_carrier_builder_local_frontier = (
            exact_mass_pair_chart.get("next_single_residual_object") or BUILDER_LOCAL_FRONTIER
        )
    elif exact_closure_closed and tau2_obstruction.get("next_single_residual_object") is not None:
        current_carrier_builder_local_frontier = tau2_obstruction.get("next_single_residual_object")
    elif exact_closure_closed and exact_wz_coordinate.get("next_residual_object_if_open") is not None:
        current_carrier_builder_local_frontier = exact_wz_coordinate.get("next_residual_object_if_open")
    elif exact_closure_closed:
        current_carrier_builder_local_frontier = exact_closure.get("stronger_residual_object")
    elif selected_population_closed:
        current_carrier_builder_local_frontier = "EWTransportExactClosureBeyondCurrentCarrier_D10"
    else:
        current_carrier_builder_local_frontier = "EWGaugeSourceTransportPairPopulationEvaluator_D10"

    if target_free_repair_closed and target_free_quintet:
        public_quintet = {
            "MW_pole": float(target_free_quintet["MW_pole"]),
            "MZ_pole": float(target_free_quintet["MZ_pole"]),
            "alpha_em_eff_inv": float(target_free_quintet["alpha_em_eff_inv"]),
            "sin2w_eff": float(target_free_quintet["sin2w_eff"]),
            "v_report": float(target_free_quintet["v_report"]),
        }
        active_builder_smallest_missing_object = None
        broader_honest_repair_frontier = None
        exact_pdg_wz_frontier = TARGET_FREE_REPAIR_FRONTIER
    elif freeze_once_repair_closed and factorized_repair_quintet:
        public_quintet = {
            "MW_pole": float(factorized_repair_quintet["MW_pole"]),
            "MZ_pole": float(factorized_repair_quintet["MZ_pole"]),
            "alpha_em_eff_inv": float(factorized_repair_quintet["alpha_em_eff_inv"]),
            "sin2w_eff": float(factorized_repair_quintet["sin2w_eff"]),
            "v_report": float(factorized_repair_quintet["v_report"]),
        }
        active_builder_smallest_missing_object = TARGET_FREE_REPAIR_FRONTIER
        broader_honest_repair_frontier = TARGET_FREE_REPAIR_FRONTIER
        exact_pdg_wz_frontier = TARGET_FREE_REPAIR_FRONTIER
    else:
        active_builder_smallest_missing_object = current_carrier_builder_local_frontier
        broader_honest_repair_frontier = (
            GLOBAL_REPAIR_FRONTIER if exact_mass_pair_chart.get("status") == "closed_smaller_primitive" else None
        )
        exact_pdg_wz_frontier = (
            GLOBAL_REPAIR_FRONTIER
            if exact_mass_pair_chart.get("status") == "closed_smaller_primitive"
            else active_builder_smallest_missing_object
        )
    public_surface_scope = (
        ["MW_pole", "MZ_pole", "alpha_em_eff_inv", "sin2w_eff", "v_report"]
        if target_free_repair_closed
        else ["MW_pole", "MZ_pole"]
    )
    public_mass_lane_quintet = dict(public_quintet)
    public_readout_quintet = {
        "MW_pole": float(public_quintet["MW_pole"]),
        "MZ_pole": float(public_quintet["MZ_pole"]),
        "alpha_em_eff_inv": None,
        "sin2w_eff": None,
        "v_report": float(public_quintet["v_report"]),
    }
    return {
        "artifact": "oph_d10_ew_source_transport_readout",
        "generated_utc": _timestamp(),
        "theorem_candidate": "EWTransportReadoutCoherence_D10",
        "proof_status": (
            "source_locked_wz_mass_lane_closed__ward_projected_em_transport_family"
            if target_free_repair_closed
            else
            "freeze_once_coherent_repair_law_closed"
            if freeze_once_repair_closed
            else
            "selected_current_carrier_split_exact_closure_closed"
            if exact_closure_closed
            else "population_selected_current_carrier"
            if selected_population_closed
            else "coherent_current_family_candidate_only"
        ),
        "smallest_predictive_missing_object": active_builder_smallest_missing_object,
        "active_builder_smallest_missing_object": active_builder_smallest_missing_object,
        "current_carrier_builder_local_frontier": current_carrier_builder_local_frontier,
        "broader_honest_repair_frontier": broader_honest_repair_frontier,
        "exact_pdg_wz_frontier": exact_pdg_wz_frontier,
        "predictive_promotion_allowed": target_free_repair_closed,
        "public_surface_candidate_allowed": True,
        "public_surface_candidate_scope": (
            ["MW_pole", "MZ_pole", "v_report"] if target_free_repair_closed else public_surface_scope
        ),
        "public_surface_policy": (
            "source_locked_wz_mass_lane_plus_ward_projected_em_transport_family"
            if target_free_repair_closed
            else
            "freeze_once_authoritative_target_coherent_repair_surface"
            if freeze_once_repair_closed
            else "best_available_reference_free_mass_pair_candidate"
        ),
        "population_evaluator_candidate": "oph_d10_ew_population_evaluator",
        "predictive_trace_evaluator_candidate": "oph_d10_ew_transport_trace_evaluator",
        "exact_closure_beyond_current_carrier_artifact": exact_closure.get("artifact"),
        "exact_closure_beyond_current_carrier_status": exact_closure.get("status"),
        "exact_closure_beyond_current_carrier_object_id": exact_closure.get("object_id"),
        "exact_wz_coordinate_beyond_single_tree_identity_artifact": exact_wz_coordinate.get("artifact"),
        "exact_wz_coordinate_beyond_single_tree_identity_status": exact_wz_coordinate.get("status"),
        "exact_wz_coordinate_beyond_single_tree_identity_object_id": exact_wz_coordinate.get("object_id"),
        "exact_wz_coordinate_beyond_single_tree_identity_depends_on": exact_wz_coordinate.get("depends_on_object"),
        "tau2_current_carrier_obstruction_artifact": tau2_obstruction.get("artifact"),
        "tau2_current_carrier_obstruction_status": tau2_obstruction.get("status"),
        "tau2_current_carrier_obstruction_object_id": tau2_obstruction.get("object_id"),
        "exact_mass_pair_chart_current_carrier_artifact": exact_mass_pair_chart.get("artifact"),
        "exact_mass_pair_chart_current_carrier_status": exact_mass_pair_chart.get("status"),
        "exact_mass_pair_chart_current_carrier_object_id": exact_mass_pair_chart.get("object_id"),
        "w_anchor_neutral_shear_factorization_artifact": w_anchor_factorization.get("artifact"),
        "w_anchor_neutral_shear_factorization_status": w_anchor_factorization.get("status"),
        "freeze_once_repair_law_id": w_anchor_factorization.get("exact_missing_law", {}).get("object_id"),
        "freeze_once_repair_target_spec": dict(w_anchor_factorization.get("target_spec", {})),
        "freeze_once_coherent_repair_package": factorized_repair_package or None,
        "freeze_once_chart_coordinates": factorized_repair_chart or None,
        "freeze_once_coherent_repaired_quintet": factorized_repair_quintet or None,
        "freeze_once_target_free_residual_object": w_anchor_factorization.get("conclusion", {}).get("stricter_still_open_object"),
        "target_free_repair_value_law_artifact": target_free_repair.get("artifact"),
        "target_free_repair_value_law_status": target_free_repair.get("status"),
        "target_free_repair_value_law_object_id": target_free_repair.get("object_id"),
        "target_free_repair_chart": target_free_chart or None,
        "target_free_repaired_couplings": target_free_couplings or None,
        "target_free_coherent_emitted_quintet": target_free_quintet or None,
        "forward_transmutation_certificate_artifact": forward_transmutation.get("artifact"),
        "forward_transmutation_certificate_status": forward_transmutation.get("status"),
        "forward_transmutation_certificate_object_id": forward_transmutation.get("object_id"),
        "forward_transmutation_certificate": (
            {
                "notation_split": forward_transmutation.get("notation_split"),
                "forward_core_solution": forward_transmutation.get("forward_core_solution"),
                "source_only_reconstruction": forward_transmutation.get("source_only_reconstruction"),
                "forward_checks": forward_transmutation.get("forward_checks"),
            }
            if forward_transmutation
            else None
        ),
        "target_free_repair_status_split": {
            "status": "closed" if target_free_repair_closed else "open",
            "theorem": target_free_repair.get("object_id") if target_free_repair_closed else None,
            "unconditional_source_only_status": minimal_conditional.get("unconditional_theorem", {}).get("name"),
            "minimal_conditional_principle": minimal_conditional.get("conditional_principle", {}).get("name"),
            "minimal_conditional_theorem": minimal_conditional.get("conditional_theorem", {}).get("name"),
            "strongest_source_only_candidate": target_emitter.get("object_id"),
        },
        "minimal_conditional_promotion_artifact": minimal_conditional.get("artifact"),
        "minimal_conditional_promotion_status": minimal_conditional.get("status"),
        "minimal_conditional_promotion": (
            {
                "unconditional_theorem": minimal_conditional.get("unconditional_theorem"),
                "conditional_principle": minimal_conditional.get("conditional_principle"),
                "conditional_theorem": minimal_conditional.get("conditional_theorem"),
                "n_c_3_specialization": minimal_conditional.get("n_c_3_specialization"),
            }
            if minimal_conditional
            else None
        ),
        "target_emitter_candidate_artifact": target_emitter.get("artifact"),
        "target_emitter_candidate_status": target_emitter.get("status"),
        "target_emitter_candidate": (
            {
                "object_id": target_emitter.get("object_id"),
                "emitter_scalar": target_emitter.get("emitter_scalar"),
                "target_emitter_law": target_emitter.get("target_emitter_law"),
                "coherent_emitted_quintet": target_emitter.get("coherent_emitted_quintet"),
                "comparison_to_frozen_local_reference_surface": target_emitter.get("comparison_to_frozen_local_reference_surface"),
            }
            if target_emitter
            else None
        ),
        "source_pair_reduction": "reopened_two_scalar_source_family",
        "coherent_readout_law_id": "EWTransportReadoutCoherence_D10",
        "best_mass_pair_law_id": (
            TARGET_FREE_REPAIR_FRONTIER
            if target_free_repair_closed
            else
            "freeze_once_coherent_d10_electroweak_repair_law"
            if freeze_once_repair_closed
            else compact_mass_slice["name"]
        ),
        "predictive_mass_promotion_allowed": target_free_repair_closed,
        "source_transport_pair_artifact": source_pair["artifact"],
        "family_source_id": "d10_running_tree",
        "scheme_id": "freeze_once",
        "origin_kernel_id": "EWTransportKernel_D10",
        "shared_scalar_package_id": "Sigma_EW_D10",
        "base_running_quintet": dict(base_quintet),
        "fixed_eta_slice": {
            "alphaY_0": alpha_y0,
            "alpha2_0": alpha2_0,
            "sin2_thetaW0": sin2_theta_w0,
            "cos2_thetaW0": cos2_theta_w0,
            "eta_EW_formula": "alpha_u * cos(2*theta_W0)",
            "eta_EW": float(compact_mass_slice["eta_EW"]),
            "sigma_EW_symbol": "free",
            "tau_Y_formula": "sigma_EW - eta_EW",
            "tau_2_formula": "sigma_EW + eta_EW",
        },
        "two_scalar_population_status": source_pair.get("two_scalar_population_status"),
        "population_selector_status": source_pair.get("population_selector_status"),
        "predictive_population_verdict": source_pair.get("predictive_population_verdict"),
        "predictive_population_closed": selected_population_closed,
        "selected_population_point": selected_point or None,
        "population_driven_quintet": (
            {
                "MW_pole": float(compact_mass_quintet["MW_pole"]),
                "MZ_pole": float(compact_mass_quintet["MZ_pole"]),
                "alpha_em_eff_inv": float(compact_mass_quintet["alpha_em_eff_inv"]),
                "sin2w_eff": float(compact_mass_quintet["sin2w_eff"]),
                "v_report": float(compact_mass_quintet["v_report"]),
            }
            if selected_population_closed
            else None
        ),
        "population_forward_map": dict(source_pair.get("population_forward_map", {})),
        "population_basis": dict(source_pair.get("population_basis", {})),
        "population_reconstruction": dict(source_pair.get("population_reconstruction", {})),
        "population_atomic_quartet": dict(source_pair.get("population_atomic_quartet", {})),
        "population_minimality_certificate": dict(source_pair.get("population_minimality_certificate", {})),
        "population_nonuniqueness_certificate": dict(source_pair.get("population_nonuniqueness_certificate", {})),
        "one_scalar_reduction_certificate": dict(source_pair.get("one_scalar_reduction_certificate", {})),
        "forbidden_inverse_witness_formulas": dict(source_pair.get("forbidden_inverse_witness_formulas", {})),
        "quartet_inverse_diagnostic_only": {
            "u_EW_from_mW": "mW^2 / (pi * v_report^2 * alpha2_0)",
            "n_EW_from_mZ": "mZ^2 / (pi * v_report^2 * (alphaY_0 + alpha2_0))",
            "eta_EW_from_u_n": "(u_EW - n_EW) / (1 - beta_EW)",
            "sigma_EW_from_u_n": "(n_EW - 1 - beta_EW * (u_EW - 1)) / (1 - beta_EW)",
        },
        "source_pair": {
            "sigma_EW": float(source_pair["special_slices"]["current_one_seed_slice"]["sigma_EW"]),
            "eta_EW": float(source_pair["special_slices"]["current_one_seed_slice"]["eta_EW"]),
            "tau_Y": float(seed_trial["tau_Y"]),
            "tau_2": float(seed_trial["tau_2"]),
        },
        "transport_entry_values": dict(seed_trial["transport_entry_values"]),
        "reported_readout_assignment": {
            "MW_pole": (
                "coherent_target_free_repair_couplings"
                if target_free_repair_closed
                else "coherent_frozen_target_repair_couplings"
                if freeze_once_repair_closed
                else "shared_scalar_package"
            ),
            "MZ_pole": (
                "coherent_target_free_repair_couplings"
                if target_free_repair_closed
                else "coherent_frozen_target_repair_couplings"
                if freeze_once_repair_closed
                else "shared_scalar_package"
            ),
            "alpha_em_eff_inv": (
                "ward_projected_u1q_transport_family"
                if target_free_repair_closed
                else
                "coherent_frozen_target_repair_couplings"
                if freeze_once_repair_closed
                else
                "source_normalized_hypercharge_readout"
                if exact_closure_closed
                else "shared_scalar_package"
            ),
            "sin2w_eff": (
                "ward_projected_u1q_transport_family"
                if target_free_repair_closed
                else
                "coherent_frozen_target_repair_couplings"
                if freeze_once_repair_closed
                else
                "source_normalized_hypercharge_readout"
                if exact_closure_closed
                else "shared_scalar_package"
            ),
            "v_report": "inherit_running_core",
        },
        "public_readout_split": {
            "wz_mass_lane_surface": (
                TARGET_FREE_REPAIR_FRONTIER
                if target_free_repair_closed
                else "FreezeOnceCoherentD10ElectroweakRepairLaw_D10"
                if freeze_once_repair_closed
                else current_carrier_builder_local_frontier
            ),
            "electromagnetic_source_anchor": "source_locked_running_family_anchor",
            "electromagnetic_transport_surface": WARD_PROJECTED_TRANSPORT_THEOREM,
            "electromagnetic_transport_kernel": "EWTransportKernel_D10",
            "compact_hypercharge_slice_supplies_public_alpha_surface": False,
        },
        "ward_projected_transport_family": {
            "theorem_id": WARD_PROJECTED_TRANSPORT_THEOREM,
            "charge_operator": "Q = T3 + Y",
            "projector": "Ward projection to the unbroken U(1)_Q channel",
            "transport_kernel_id": "EWTransportKernel_D10",
            "transport_readout_clause": "EWTransportReadoutCoherence_D10",
            "scalar_provenance_clause": "EWScalarProvenanceEquality_D10",
            "anchor_scale": "m_Z^2",
            "anchor_alpha_em_eff_inv": float(base_quintet["alpha_em_eff_inv"]),
            "thomson_endpoint_alpha_em_eff_inv": THOMSON_ALPHA_INV,
            "delta_alpha_from_anchor_to_thomson": THOMSON_ALPHA_INV / float(base_quintet["alpha_em_eff_inv"]) - 1.0,
            "transport_ratio_tQ_0_over_tQ_mZ2": float(base_quintet["alpha_em_eff_inv"]) / THOMSON_ALPHA_INV,
            "physical_readout_formula": "alpha_em^-1(q^2;P) = 8*pi^2 / t_Q(q^2;P)",
            "thomson_limit_formula": "alpha_Th^-1(P) = lim_{q^2 -> 0} alpha_em^-1(q^2;P)",
            "shared_provenance": {
                "family_source_id": "d10_running_tree",
                "scheme_id": "freeze_once",
                "origin_kernel_id": "EWTransportKernel_D10",
            },
        },
        "shared_scalar_values_reported": dict(shared_scalar_values),
        "public_emitted_quintet": public_readout_quintet,
        "public_mass_lane_quintet": public_mass_lane_quintet,
        "exact_closure_emitted_quintet": exact_outputs,
        "coherent_quintet_family_formula": {
            "alphaY_star": "alphaY_0 * (1 + sigma_EW - eta_EW)",
            "alpha2_star": "alpha2_0 * (1 + sigma_EW + eta_EW)",
            "MW_pole": "v_report * sqrt(pi * alpha2_star)",
            "MZ_pole": "v_report * sqrt(pi * (alphaY_star + alpha2_star))",
            "alpha_em_eff_inv": "(alphaY_star + alpha2_star) / (alphaY_star * alpha2_star)",
            "sin2w_eff": "alphaY_star / (alphaY_star + alpha2_star)",
        },
        "quartet_atomicity": {
            "all_four_readouts_share_one_population_point": True,
            "independent_post_population_readout_scalar_remaining": False,
            "candidate_status": "all_or_none",
        },
        "current_compact_point": {
            "constraint": "sigma_EW = -eta_EW",
            "sigma_EW_formula": "-eta_EW",
            "sigma_EW": float(compact_mass_slice["sigma_EW"]),
            "eta_EW": float(compact_mass_slice["eta_EW"]),
            "tau_Y": float(compact_mass_slice["tau_Y"]),
            "tau_2": float(compact_mass_slice["tau_2"]),
            "alphaY_star": compact_alpha_y,
            "alpha2_star": compact_alpha2,
        },
        "current_compact_shared_scalar_values": dict(compact_shared_scalar_values),
        "current_compact_emitted_quintet": {
            "MW_pole": float(compact_mass_quintet["MW_pole"]),
            "MZ_pole": float(compact_mass_quintet["MZ_pole"]),
            "alpha_em_eff_inv": float(compact_mass_quintet["alpha_em_eff_inv"]),
            "sin2w_eff": float(compact_mass_quintet["sin2w_eff"]),
            "v_report": float(compact_mass_quintet["v_report"]),
        },
        "tree_identity_residuals": {
            "sin2_from_mass_ratio_minus_reported": (
                1.0
                - (compact_mass_quintet["MW_pole"] / compact_mass_quintet["MZ_pole"]) ** 2
                - compact_mass_quintet["sin2w_eff"]
            ),
            "alpha_from_transported_pair_minus_reported": (
                (compact_alpha_y + compact_alpha2) / (compact_alpha_y * compact_alpha2)
                - compact_mass_quintet["alpha_em_eff_inv"]
            ),
        },
        "mass_pair_predictive_candidate": {
            "status": (
                "target_free_source_only_repair_exact"
                if target_free_repair_closed
                else
                "freeze_once_coherent_repair_exact"
                if freeze_once_repair_closed
                else
                "selected_current_carrier_nonexact_pair"
                if selected_population_closed
                else "candidate_only"
            ),
            "law": (
                {
                    "name": TARGET_FREE_REPAIR_FRONTIER,
                    "equivalent_two_coordinate_chart": "(tau2_tree_exact, delta_n_tree_exact)",
                    "proof_gate": "single_family_single_P_no_mixed_readout",
                    "historical_validation_surface": "FreezeOnceCoherentD10ElectroweakRepairLaw_D10",
                }
                if target_free_repair_closed
                else
                {
                    "name": "freeze_once_coherent_d10_electroweak_repair_law",
                    "equivalent_two_coordinate_chart": "(tau2_w_anchor, delta_n_dagger)",
                    "proof_gate": "single_family_single_P_no_mixed_readout",
                    "stricter_still_open_object": TARGET_FREE_REPAIR_FRONTIER,
                }
                if freeze_once_repair_closed
                else dict(compact_mass_slice["law"])
            ),
            "sigma_EW": float(compact_mass_slice["sigma_EW"]) if not (freeze_once_repair_closed or target_free_repair_closed) else None,
            "eta_EW": float(compact_mass_slice["eta_EW"]) if not (freeze_once_repair_closed or target_free_repair_closed) else None,
            "tau_Y": (
                float(compact_mass_slice["tau_Y"])
                if not (freeze_once_repair_closed or target_free_repair_closed)
                else float(target_free_chart.get("tauY_fiber"))
                if target_free_repair_closed
                else float(factorization_point.get("tauY_fiber_dagger", 0.0))
            ),
            "tau_2": (
                float(compact_mass_slice["tau_2"])
                if not (freeze_once_repair_closed or target_free_repair_closed)
                else float(target_free_chart.get("tau2_tree_exact", 0.0))
                if target_free_repair_closed
                else float(factorized_repair_chart.get("tau2_w_anchor", 0.0))
            ),
            "delta_n_dagger": (
                float(target_free_chart.get("delta_n_tree_exact"))
                if target_free_repair_closed
                else None if not freeze_once_repair_closed else float(factorized_repair_chart.get("delta_n_dagger"))
            ),
            "delta_alpha2_dagger": (
                float(target_free_couplings.get("delta_alpha2"))
                if target_free_repair_closed
                else None if not freeze_once_repair_closed else float(factorized_repair_package.get("delta_alpha2_dagger"))
            ),
            "delta_alphaY_parallel": (
                float(target_free_couplings.get("delta_alphaY_parallel"))
                if target_free_repair_closed
                else None if not freeze_once_repair_closed else float(factorized_repair_package.get("delta_alphaY_parallel"))
            ),
            "delta_alphaY_perp": (
                float(target_free_couplings.get("delta_alphaY_perp"))
                if target_free_repair_closed
                else None if not freeze_once_repair_closed else float(factorized_repair_package.get("delta_alphaY_perp"))
            ),
            "MW_pole": float(public_quintet["MW_pole"]),
            "MZ_pole": float(public_quintet["MZ_pole"]),
        },
        "coherent_emitter_formula": {
            "MW_pole": "W0 * (1 + delta_rho + delta_rW)",
            "MZ_pole": "Z0 * (1 + delta_rho)",
            "alpha_em_eff_inv": "a0 * (1 + delta_alpha)",
            "sin2w_eff": "s0 * (1 + delta_kappa)",
            "v_report": "v0",
        },
        "transport_slice_diagnostics": {
            "transport_quintet": dict(public_mass_lane_quintet),
            "base_running_quintet": dict(base_quintet),
        },
        "proof_gate": {
            "provenance_equality_required": True,
            "single_post_transport_tree_identity_required": not exact_closure_closed,
            "single_remaining_mass_scalar_after_fiber_law": exact_wz_coordinate.get("next_residual_object_if_open") == "tau2_tree_exact",
            "no_run_pole_mix": True,
            "common_provenance_required": True,
        },
        "notes": [
            "One running-family base quintet and one shared scalar package Sigma_EW_D10 organize the D10 lane.",
            "The selected carrier point emits the mass pair directly from the transported D10 couplings.",
            "The reduced two-scalar carrier gives the compact source family on the selected D10 lane.",
            (
                "The exact mass-pair chart on the selected carrier is closed, and the builder-local mass-side residual is the selector `EWExactMassPairSelector_D10` on that chart."
                if exact_mass_pair_chart.get("status") == "closed_smaller_primitive"
                else "The smaller fiberwise population tree law removes the placeholder unsplit tree shell, so the remaining D10 mass-side residual is the single scalar tau2_tree_exact."
                if exact_wz_coordinate.get("next_residual_object_if_open") == "tau2_tree_exact"
                else "The exact W/Z coordinate shell still depends on a stronger unsplit tree identity."
            ),
            (
                "The D10 mass-side surface carries the W/Z pair from the source trunk. The freeze-once coherent pair serves as compare-only validation on the same family and agrees with the source-only mass-side law to machine scale."
                if target_free_repair_closed
                else "The freeze-once coherent repair law is closed on one authoritative frozen target pair, so the W/Z surface comes from one coherent coupling pair on the shared D10 family."
                if freeze_once_repair_closed
                else "The broader honest exact-PDG W/Z frontier is the repair branch `D10RepairBranchBeyondCurrentCarrier` beyond the present current carrier, not just the builder-local selector shell."
            ),
            "The selected electroweak point is the chosen carrier point itself, not a separate transported seed placeholder.",
            "The compact point records the same family on the fixed-eta slice eta_EW = alpha_u * cos(2*theta_W0) with free sigma_EW; the no-new-parameter point is sigma_EW = -eta_EW.",
            (
                "The compact anti-diagonal slice restores alpha_em^-1 and sin^2(theta_W) through a source-normalized hypercharge readout on the selected carrier point, but that slice does not define the public fine-structure theorem."
                if exact_closure_closed
                else "The fixed-eta trace evaluator remains useful diagnostically, but the live predictive blocker is now exact electroweak closure beyond the current selected carrier point rather than another selector on the exhausted compact slice."
            ),
            "The physical electromagnetic readout is anchored at `a0 = alpha_em^-1(m_Z^2) = 128.30576920234813` on the source-locked running family and is read through Ward projection to the unbroken `U(1)_Q` channel, with Thomson endpoint `alpha^-1(0) = 137.035999177`.",
            (
                "The source-only underdetermination theorem, the minimal conditional route through ColorBalancedQuadraticRepairDescent_D10, and the stronger source-only candidate EWTargetEmitter_D10 stay on disk as lower-level mass-side objects beneath the Ward-projected electromagnetic readout."
                if target_free_repair_closed and (minimal_conditional or target_emitter)
                else "The target-free D10 problem splits across the source-only underdetermination theorem, the conditional route through ColorBalancedQuadraticRepairDescent_D10, and the source-only candidate EWTargetEmitter_D10."
                if minimal_conditional or target_emitter
                else "No sharper source-only target-free D10 split is attached to this readout."
            ),
            (
                "The forward transmutation certificate `EWForwardTransmutationCertificate_D10` records the non-circular P -> alpha_U -> t map explicitly and keeps the source-ratio beta_ratio_EW separate from the transmutation counting factor beta_transmutation_EW = N_c + 1."
                if forward_transmutation
                else "No explicit forward transmutation certificate is attached to this readout."
            ),
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the factorized D10 electroweak readout candidate.")
    parser.add_argument("--family", default=str(DEFAULT_FAMILY))
    parser.add_argument("--source-pair", default=str(DEFAULT_SOURCE_PAIR))
    parser.add_argument("--exact-closure", default=str(DEFAULT_EXACT_CLOSURE))
    parser.add_argument("--exact-wz-coordinate", default=str(DEFAULT_EXACT_WZ_COORDINATE))
    parser.add_argument("--tau2-obstruction", default=str(DEFAULT_TAU2_OBSTRUCTION))
    parser.add_argument("--exact-mass-pair-chart", default=str(DEFAULT_EXACT_MASS_PAIR_CHART))
    parser.add_argument("--w-anchor-factorization", default=str(DEFAULT_W_ANCHOR_FACTORIZATION))
    parser.add_argument("--minimal-conditional-promotion", default=str(DEFAULT_MINIMAL_CONDITIONAL))
    parser.add_argument("--target-emitter", default=str(DEFAULT_TARGET_EMITTER))
    parser.add_argument("--target-free-repair", default=str(DEFAULT_TARGET_FREE_REPAIR))
    parser.add_argument("--forward-transmutation", default=str(DEFAULT_FORWARD_TRANSMUTATION))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    family = json.loads(Path(args.family).read_text(encoding="utf-8"))
    source_pair = json.loads(Path(args.source_pair).read_text(encoding="utf-8"))
    exact_closure_path = Path(args.exact_closure)
    exact_closure = json.loads(exact_closure_path.read_text(encoding="utf-8")) if exact_closure_path.exists() else None
    exact_wz_coordinate_path = Path(args.exact_wz_coordinate)
    exact_wz_coordinate = (
        json.loads(exact_wz_coordinate_path.read_text(encoding="utf-8"))
        if exact_wz_coordinate_path.exists()
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
    w_anchor_factorization_path = Path(args.w_anchor_factorization)
    w_anchor_factorization = (
        json.loads(w_anchor_factorization_path.read_text(encoding="utf-8"))
        if w_anchor_factorization_path.exists()
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
    artifact = build_artifact(
        family,
        source_pair,
        exact_closure,
        exact_wz_coordinate,
        tau2_obstruction,
        exact_mass_pair_chart,
        w_anchor_factorization,
        minimal_conditional,
        target_emitter,
        target_free_repair,
        forward_transmutation,
    )

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
