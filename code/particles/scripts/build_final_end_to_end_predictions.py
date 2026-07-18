#!/usr/bin/env python3
"""Build the final end-to-end particle prediction bundle."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from decimal import Decimal, getcontext
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
PARTICLES_ROOT = ROOT / "particles"
P_ROOT = ROOT / "P_derivation"
P_TRUNK = P_ROOT / "runtime" / "p_closure_trunk_current.json"
MEASURED_ENDPOINT = P_ROOT / "runtime" / "measured_endpoint_calibration_current.json"
PIPELINE_STATUS = PARTICLES_ROOT / "runs" / "status" / "particle_pipeline_closure_status.json"
EXACT_NONHADRON = PARTICLES_ROOT / "exact_nonhadron_masses.json"
RESULTS_STATUS = PARTICLES_ROOT / "results_status.json"
DIRECT_TOP = PARTICLES_ROOT / "runs" / "calibration" / "direct_top_bridge_contract.json"
CHARGED_TRACE_LIFT_REQUIRED = (
    PARTICLES_ROOT / "runs" / "leptons" / "charged_determinant_trace_lift_attachment_required.json"
)
QUARK_SIGMA_REQUIRED = PARTICLES_ROOT / "runs" / "flavor" / "quark_sigma_source_datum_no_target_leak_required.json"
QUARK_SIGMA_OBSTRUCTION = (
    PARTICLES_ROOT / "runs" / "flavor" / "quark_sigma_source_nonidentifiability_obstruction.json"
)
QUARK_AXIOM_LEVEL_OBSTRUCTION = (
    PARTICLES_ROOT
    / "runs"
    / "flavor"
    / "quark_axiom_level_yukawa_moduli_nonidentifiability.json"
)
QUARK_SCHEME_OBSTRUCTION = (
    PARTICLES_ROOT / "runs" / "flavor" / "quark_running_mass_scheme_convention_obstruction.json"
)
EMPIRICAL_EE_REGISTRY = PARTICLES_ROOT / "hadron" / "empirical_ee_hadrons_sources.yaml"
EMPIRICAL_EE_SCHEMA = PARTICLES_ROOT / "hadron" / "empirical_ee_hadronic_spectral_measure.schema.json"
HIERARCHY_ROOT = PARTICLES_ROOT / "hierarchy"
HIERARCHY_RESONANCE = HIERARCHY_ROOT / "certificates" / "R_local_global_hierarchy_resonance_closeout_335.json"
HIERARCHY_EW_CAPACITY = HIERARCHY_ROOT / "certificates" / "R_EW_global_capacity_certificate.json"
HIERARCHY_READBACK = HIERARCHY_ROOT / "certificates" / "R_readback_resolution_certificate.json"
HIERARCHY_M_REP = HIERARCHY_ROOT / "certificates" / "R_m_rep_24_certificate.json"
HIERARCHY_SCREEN_SIEVE = HIERARCHY_ROOT / "certificates" / "R_screen_sieve_icosahedral_certificate.json"
HIERARCHY_NATURALITY = HIERARCHY_ROOT / "issue_332_rg_naturality_certificate.json"
HIERARCHY_PIXEL_SCREEN = HIERARCHY_ROOT / "certificates" / "R_pixel_screen_resonance_summary.json"
DEFAULT_JSON_OUT = PARTICLES_ROOT / "runs" / "status" / "final_end_to_end_predictions.json"
DEFAULT_MD_OUT = PARTICLES_ROOT / "PARTICLE_PIPELINE_STATUS.md"
getcontext().prec = 90


PARTICLE_ORDER = [
    "photon",
    "gluon",
    "graviton",
    "higgs",
    "electron",
    "muon",
    "tau",
    "up_quark",
    "down_quark",
    "strange_quark",
    "charm_quark",
    "bottom_quark",
    "top_quark",
    "electron_neutrino",
    "muon_neutrino",
    "tau_neutrino",
]


def _now_utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _load_optional_json(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    return _load_json(path)


def _rel(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def _display_status(status: str) -> str:
    if status is None:
        return "missing"
    return status.replace("current_corpus", "corpus_limited")


def _prediction_entry(entry: dict[str, Any]) -> dict[str, Any]:
    if entry.get("mass_gev") is not None:
        value = float(entry["mass_gev"])
        unit = "GeV"
    else:
        value = float(entry["mass_eV"])
        unit = "eV"
    return {
        "particle_id": entry["particle_id"],
        "label": entry.get("label"),
        "value": value,
        "unit": unit,
        "exact_kind": entry.get("exact_kind"),
        "scope": entry.get("scope"),
        "promotable": entry.get("promotable"),
        "source_artifact": entry.get("source_artifact"),
        "supporting_scope_closure_artifact": entry.get("supporting_scope_closure_artifact"),
    }


def _fine_structure_surface(measured_endpoint: dict[str, Any]) -> dict[str, Any]:
    calibrated = measured_endpoint["calibrated_values"]
    source_candidate = measured_endpoint["current_source_candidate"]
    endpoint_requirement = measured_endpoint["codata_mapped_endpoint_requirement"]
    alpha_u = _load_json(HIERARCHY_NATURALITY)["source_values"]["alpha_U"]
    root_alpha_inv = Decimal(source_candidate["alpha_inv"])
    alpha_u_dec = Decimal(alpha_u)
    near_endpoint_alpha_inv = root_alpha_inv + alpha_u_dec
    measured_alpha_inv = Decimal(calibrated["alpha_inv_0"])
    missing_hadronic_correction = measured_alpha_inv - near_endpoint_alpha_inv
    relative_shortfall = missing_hadronic_correction / measured_alpha_inv
    near_endpoint = {
        "row_class": "source_side_no_hadron_near_endpoint",
        "status": "certified_source_root_plus_public_width_diagnostic",
        "formula": "alpha_root^-1 + alpha_U(P_star)",
        "alpha_root_inv": source_candidate["alpha_inv"],
        "alpha_U": alpha_u,
        "alpha_inv": str(near_endpoint_alpha_inv),
        "alpha": str(Decimal(1) / near_endpoint_alpha_inv),
        "P": source_candidate["P"],
        "missing_hadronic_correction_alpha_inv": str(missing_hadronic_correction),
        "relative_shortfall": str(relative_shortfall),
        "percent_shortfall": str(relative_shortfall * Decimal(100)),
        "one_part_in": str(measured_alpha_inv / missing_hadronic_correction),
        "open_physical_map_requirements": [
            "source Ward-projected hadronic spectral measure",
            "same-scheme finite remainder",
            "one-coordinate fixed-point evaluation",
            "interval certificate and theory uncertainty",
        ],
        "promotable_as_exact_source_theorem": False,
    }
    root_audit = {
        "row_class": "root_only_audit_trunk",
        "status": "certified_source_root_coordinate",
        "alpha_inv": source_candidate["alpha_inv"],
        "alpha": source_candidate["alpha"],
        "P": source_candidate["P"],
        "missing_inverse_alpha_units_to_thomson_endpoint": source_candidate["missing_inverse_alpha_units"],
        "p_gap_implemented_minus_empirical": source_candidate["p_gap_implemented_minus_calibrated"],
        "promotable_as_exact_source_theorem": False,
    }
    return {
        "source_side_no_hadron_near_endpoint": near_endpoint,
        "source_only_oph": {
            "row_class": "source_only_physical_thomson_endpoint",
            "status": "not_emitted",
            "promotable_as_exact_source_theorem": False,
        },
        "root_only_audit": root_audit,
        "oph_plus_empirical_hadron_closure": {
            "row_class": "oph_plus_empirical_hadron_closure",
            "status": measured_endpoint["status"],
            "alpha_inv": calibrated["alpha_inv_0"],
            "alpha_inv_standard_uncertainty": calibrated["alpha_inv_0_standard_uncertainty"],
            "alpha": calibrated["alpha_0"],
            "P": calibrated["P_from_outer_equation"],
            "P_standard_uncertainty": calibrated["P_standard_uncertainty"],
            "required_transport_delta_alpha_inv": endpoint_requirement["required_transport_delta_alpha_inv"],
            "implemented_transport_delta_alpha_inv": endpoint_requirement["implemented_transport_delta_alpha_inv"],
            "missing_source_transport_delta_alpha_inv": endpoint_requirement[
                "missing_source_transport_delta_alpha_inv"
            ],
            "required_screening_factor": endpoint_requirement["required_screening_factor"],
            "residual_second_order_coefficient": endpoint_requirement[
                "residual_second_order_coefficient"
            ],
            "promotable_as_exact_source_theorem": False,
        },
        "empirical_payload_policy": measured_endpoint["empirical_hadron_closure"],
    }


def _hierarchy_surface() -> dict[str, Any]:
    resonance = _load_optional_json(HIERARCHY_RESONANCE)
    ew_capacity = _load_optional_json(HIERARCHY_EW_CAPACITY)
    readback = _load_optional_json(HIERARCHY_READBACK)
    m_rep = _load_optional_json(HIERARCHY_M_REP)
    screen_sieve = _load_optional_json(HIERARCHY_SCREEN_SIEVE)
    naturality = _load_optional_json(HIERARCHY_NATURALITY)
    pixel_screen = _load_optional_json(HIERARCHY_PIXEL_SCREEN)

    exact_capacity = dict((ew_capacity or {}).get("exact_capacity_fixed_point") or {})
    source_values = dict((ew_capacity or {}).get("source_values") or {})
    exact_statement = dict((resonance or {}).get("exact_surviving_statement") or {})
    target_relation = dict((resonance or {}).get("target_relation") or {})
    readback_resolution = dict((resonance or {}).get("finite_readback_resolution_certificate") or {})
    round_count = dict((resonance or {}).get("round_count_certificate") or {})
    screen_sieve_summary = dict((resonance or {}).get("screen_sieve_certificate") or {})

    return {
        "status": {
            "resonance_status": (resonance or {}).get("status"),
            "accepted": bool((resonance or {}).get("accepted", False)),
            "full_theorem_grade_resonance_promoted": bool(
                (resonance or {}).get("full_theorem_grade_resonance_promoted", False)
            ),
            "work_in_progress_receipts": list((resonance or {}).get("work_in_progress_receipts") or []),
            "ew_capacity_status": (ew_capacity or {}).get("status"),
            "readback_status": (readback or {}).get("status") or readback_resolution.get("status"),
            "m_rep_status": (m_rep or {}).get("status") or round_count.get("status"),
            "screen_sieve_status": (screen_sieve or {}).get("status") or screen_sieve_summary.get("status"),
            "higgs_naturality_status": (naturality or {}).get("status") or (
                "closed_exact_selected_branch" if (naturality or {}).get("accepted") else None
            ),
        },
        "source_values": {
            "P_star": source_values.get("P_star"),
            "alpha_U": source_values.get("alpha_U"),
            "alpha_U_interval": source_values.get("alpha_U_interval"),
        },
        "local_global_bridge": {
            "target_relation": target_relation,
            "projection_map": exact_statement.get("projection_map"),
            "bridge_residual_formula": exact_statement.get("bridge_residual"),
            "exact_bridge_target": exact_statement.get("exact_bridge_target"),
            "N_CRC_EW": exact_capacity.get("N_CRC_EW") or exact_statement.get("N_EW_public_endpoint"),
            "bridge_residual": exact_capacity.get("bridge_residual"),
            "fixed_point_residual_x": exact_capacity.get("fixed_point_residual_x"),
            "v_over_E_cell_source": exact_capacity.get("v_over_E_cell_source"),
            "v_identity_error": exact_capacity.get("v_identity_error"),
        },
        "factor_origins": {
            "screen_ports": screen_sieve_summary.get("orbit_size", 12),
            "curvature_charge": screen_sieve_summary.get("total_curvature_charge", 12),
            "m_rep": round_count.get("m_rep", 24),
            "specialized_exponent": round_count.get("specialized_exponent", "-1/48"),
            "higgs_naturality_defect": (naturality or {}).get("epsilon_H"),
        },
        "claim_boundary": {
            "improves": (
                "The selected branch has zero algebraic bridge residual and epsilon_H=0. "
                "Its hierarchy interpretation is conditional on HIERARCHY-SCREEN-READOUT."
            ),
            "count_independence": (
                "m_rep=24 is derived from the doubled observer-visible product adjoint, "
                "independently of the screen sieve. The screen's twelve ports and oriented "
                "24-slot register are an arithmetic alignment, not a physical identification."
            ),
            "missing_bridge": (
                "HIERARCHY-SCREEN-READOUT must identify log(E_cell/v)=Gamma_screen; "
                "the alpha_U branch must then establish Pi_EW=4P, equivalently B_EW=0."
            ),
            "does_not_promote": [
                "public Thomson endpoint without the missing hadronic spectral payload",
                "electroweak massive-boson mass rows",
                "charged-lepton absolute masses",
                "source-only hadron masses",
                "strong CP",
                "SI G without the full no-G clock stack",
            ],
            "capacity_conditions": [
                "supply HIERARCHY-SCREEN-READOUT and the Pi_EW=4P identification",
                "construct the physical readback map F",
                "discharge CP-1, CP-2, and CP-3",
                "propagate the joint cosmological posterior before comparison",
            ],
        },
        "pixel_screen_resonance": {
            "status": (pixel_screen or {}).get("status"),
            "accepted": bool((pixel_screen or {}).get("accepted", False)),
            "source_pair": (pixel_screen or {}).get("source_pair"),
            "tile_identity": (pixel_screen or {}).get("pixel_screen_tile_identity"),
            "shared_12_24_port_alignment": (pixel_screen or {}).get(
                "shared_12_24_port_alignment"
            ),
            "dimensionless_de_sitter_coordinate": (pixel_screen or {}).get(
                "dimensionless_de_sitter_coordinate"
            ),
            "claim_boundary": (pixel_screen or {}).get("claim_boundary"),
            "checks": (pixel_screen or {}).get("checks"),
        },
        "artifacts": {
            "local_global_resonance": _rel(HIERARCHY_RESONANCE),
            "ew_capacity": _rel(HIERARCHY_EW_CAPACITY),
            "finite_readback": _rel(HIERARCHY_READBACK),
            "m_rep": _rel(HIERARCHY_M_REP),
            "screen_sieve": _rel(HIERARCHY_SCREEN_SIEVE),
            "higgs_naturality": _rel(HIERARCHY_NATURALITY),
            "pixel_screen_resonance": _rel(HIERARCHY_PIXEL_SCREEN),
        },
    }


def build_payload() -> dict[str, Any]:
    p_trunk = _load_json(P_TRUNK)
    measured_endpoint = _load_json(MEASURED_ENDPOINT)
    pipeline = _load_json(PIPELINE_STATUS)
    exact = _load_json(EXACT_NONHADRON)
    results = _load_json(RESULTS_STATUS)
    direct_top = _load_json(DIRECT_TOP)
    charged_trace_required = _load_optional_json(CHARGED_TRACE_LIFT_REQUIRED) or {}
    quark_sigma_required = _load_optional_json(QUARK_SIGMA_REQUIRED) or {}
    quark_sigma_obstruction = _load_optional_json(QUARK_SIGMA_OBSTRUCTION) or {}
    quark_axiom_level_obstruction = _load_optional_json(QUARK_AXIOM_LEVEL_OBSTRUCTION) or {}
    quark_scheme_obstruction = _load_optional_json(QUARK_SCHEME_OBSTRUCTION) or {}
    by_id = {entry["particle_id"]: _prediction_entry(entry) for entry in exact["entries"]}
    predictions = [by_id[particle_id] for particle_id in PARTICLE_ORDER if particle_id in by_id]
    particle_five_gates = [
        gate
        for gate in pipeline["issue_gates"]
        if gate["issue"] in {223, 224, 225, 234, 235, 32, 153, 199, 201, 207}
    ]

    return {
        "artifact": "oph_final_current_end_to_end_particle_predictions",
        "generated_utc": _now_utc(),
        "scope": "nonhadron_particle_pipeline_with_empirical_hadron_closure_policy",
        "claim_status": "final_nonhadron_predictions_with_classical_carriers_and_empirical_hadrons_separated",
        "source_surfaces": {
            "p_trunk": "code/P_derivation/runtime/p_closure_trunk_current.json",
            "measured_endpoint_calibration": (
                "code/P_derivation/runtime/measured_endpoint_calibration_current.json"
            ),
            "thomson_endpoint_package": "code/P_derivation/runtime/thomson_endpoint_package_current.json",
            "pipeline_status": "code/particles/runs/status/particle_pipeline_closure_status.json",
            "exact_nonhadron": "code/particles/exact_nonhadron_masses.json",
            "carrier_mode_acceptance": "code/particles/runs/status/carrier_mode_acceptance.json",
            "results_status": "code/particles/results_status.json",
            "direct_top_bridge": "code/particles/runs/calibration/direct_top_bridge_contract.json",
            "charged_determinant_trace_lift_attachment_required": (
                "code/particles/runs/leptons/charged_determinant_trace_lift_attachment_required.json"
            ),
            "quark_sigma_source_datum_no_target_leak_required": (
                "code/particles/runs/flavor/quark_sigma_source_datum_no_target_leak_required.json"
            ),
            "quark_sigma_source_nonidentifiability_obstruction": (
                "code/particles/runs/flavor/quark_sigma_source_nonidentifiability_obstruction.json"
            ),
            "quark_axiom_level_yukawa_moduli_nonidentifiability": (
                "code/particles/runs/flavor/quark_axiom_level_yukawa_moduli_nonidentifiability.json"
            ),
            "quark_running_mass_scheme_convention_obstruction": (
                "code/particles/runs/flavor/quark_running_mass_scheme_convention_obstruction.json"
            ),
            "hadron_policy": "docs/HADRON.md",
            "empirical_ee_hadrons_source_registry": (
                "code/particles/hadron/empirical_ee_hadrons_sources.yaml"
            ),
            "empirical_ee_hadronic_spectral_measure_schema": (
                "code/particles/hadron/empirical_ee_hadronic_spectral_measure.schema.json"
            ),
            "hierarchy_local_global_resonance": (
                "code/particles/hierarchy/certificates/R_local_global_hierarchy_resonance_closeout_335.json"
            ),
            "hierarchy_ew_capacity": (
                "code/particles/hierarchy/certificates/R_EW_global_capacity_certificate.json"
            ),
            "hierarchy_higgs_naturality": (
                "code/particles/hierarchy/issue_332_rg_naturality_certificate.json"
            ),
            "hierarchy_pixel_screen_resonance": (
                "code/particles/hierarchy/certificates/R_pixel_screen_resonance_summary.json"
            ),
        },
        "p_closure": {
            "P": p_trunk["fixed_point_candidate"]["P"],
            "alpha_inv": p_trunk["fixed_point_candidate"]["alpha_inv"],
            "claim_status": p_trunk["claim_status"],
            "may_feed_live_particle_predictions": p_trunk["consumer_policy"]["may_feed_live_particle_predictions"],
        },
        "runtime_inputs": results.get("inputs", {}),
        "output_classes": [
            "source_side_no_hadron_near_endpoint",
            "source_only_oph",
            "root_only_audit",
            "oph_plus_empirical_hadron_closure",
            "compare_only",
            "work_in_progress",
            "conditional_classical_carrier_mode_not_particle_mass_prediction",
        ],
        "withheld_non_prediction_rows": exact.get("withheld_entries", []),
        "classical_carrier_modes": exact.get("classical_carrier_modes", []),
        "charged_lepton_anchor_boundary": {
            "artifact": charged_trace_required.get("artifact"),
            "status": charged_trace_required.get("status"),
            "required_identity": charged_trace_required.get("required_identity"),
            "equivalent_defect": charged_trace_required.get("equivalent_defect"),
            "missing_for_promotion": charged_trace_required.get("missing_for_promotion", []),
            "forbidden_ancestors": charged_trace_required.get("forbidden_ancestors", []),
            "current_closed_chain": charged_trace_required.get("current_closed_chain", {}),
        },
        "quark_sigma_source_boundary": {
            "artifact": quark_sigma_required.get("artifact"),
            "status": quark_sigma_required.get("status"),
            "claim_tier": quark_sigma_required.get("claim_tier"),
            "required_identity": quark_sigma_required.get("required_identity"),
            "source_only_sigma_emitted": quark_sigma_required.get("source_only_sigma_emitted"),
            "downstream_algebra_closed": quark_sigma_required.get("downstream_algebra_closed"),
            "closure_kind": quark_sigma_required.get("closure_kind"),
            "obstruction_artifact": quark_sigma_required.get("obstruction_artifact"),
            "compatible_spread_fiber": quark_sigma_required.get("compatible_spread_fiber"),
            "compatible_spread_fiber_dimension": quark_sigma_required.get(
                "compatible_spread_fiber_dimension"
            ),
            "independent_unselected_coordinates": quark_sigma_required.get(
                "independent_unselected_coordinates", []
            ),
            "missing_for_promotion": quark_sigma_required.get("missing_for_promotion", []),
            "forbidden_ancestors": quark_sigma_required.get("forbidden_ancestors", []),
            "reopen_requirements": quark_sigma_required.get("reopen_requirements", {}),
            "dependency_audit": quark_sigma_required.get("dependency_audit", {}),
            "issue_acceptance": {
                "377": quark_sigma_obstruction.get("issue_377_acceptance_met_as_obstruction"),
                "379": quark_sigma_obstruction.get("issue_379_acceptance_met_as_obstruction"),
                "380": quark_sigma_obstruction.get("issue_380_acceptance_met_as_obstruction"),
            },
        },
        "quark_scheme_and_yukawa_boundary": {
            "artifact": quark_scheme_obstruction.get("artifact"),
            "proof_status": quark_scheme_obstruction.get("proof_status"),
            "closure_kind": quark_scheme_obstruction.get("closure_kind"),
            "row_partition": quark_scheme_obstruction.get("row_partition", {}),
            "declared_comparison_conventions": quark_scheme_obstruction.get(
                "declared_comparison_conventions", {}
            ),
            "stored_matrix_dimensional_audit": quark_scheme_obstruction.get(
                "stored_matrix_dimensional_audit", {}
            ),
            "issue_acceptance": quark_scheme_obstruction.get("issue_acceptance", {}),
            "closure_effect": quark_scheme_obstruction.get("closure_effect", {}),
            "physical_rg_dynamics_nondefinability": quark_scheme_obstruction.get(
                "physical_rg_dynamics_nondefinability", {}
            ),
            "external_chart_boundary": quark_scheme_obstruction.get(
                "external_chart_boundary", {}
            ),
        },
        "quark_axiom_level_nondefinability": {
            "artifact": quark_axiom_level_obstruction.get("artifact"),
            "proof_status": quark_axiom_level_obstruction.get("proof_status"),
            "additional_axioms_used": quark_axiom_level_obstruction.get(
                "additional_axioms_used"
            ),
            "counterfamily": quark_axiom_level_obstruction.get("counterfamily", {}),
            "MAR_audit": (
                quark_axiom_level_obstruction.get("axiom_invariance_audit", {}).get(
                    "Axiom_5_MAR", {}
                )
            ),
            "public_numeric_quark_rows_allowed": quark_axiom_level_obstruction.get(
                "public_numeric_quark_rows_allowed"
            ),
        },
        "fine_structure": _fine_structure_surface(measured_endpoint),
        "hierarchy_and_naturality": _hierarchy_surface(),
        "finalization_gates": pipeline["finalization_gates"],
        "particle_five_issue_gates": particle_five_gates,
        "companion_open_branches": list(pipeline.get("companion_status_branches", [])),
        "predictions": predictions,
        "hadron_policy": {
            "source_only_hadron_predictions_emitted": False,
            "empirical_hadron_closure_allowed_for_display": True,
            "policy_artifact": "docs/HADRON.md",
            "source_registry": EMPIRICAL_EE_REGISTRY.relative_to(ROOT).as_posix(),
            "empirical_payload_schema": EMPIRICAL_EE_SCHEMA.relative_to(ROOT).as_posix(),
            "reason": (
                "Source-only hadron outputs require a working OPH hadron backend. Empirical "
                "hadron closure values stay in a separate output class; the e+e- spectral payload "
                "has a source registry and schema."
            ),
            "github_issues": [153, 157],
        },
        "direct_top_auxiliary_comparison": {
            "primary_top_codomain": direct_top["current_target_audit_coordinate"]["pdg_summary_id"],
            "auxiliary_direct_top_codomain": direct_top["auxiliary_direct_top_coordinate"]["pdg_summary_id"],
            "value_policy": "compare_only_codomain_values_withheld_from_final_prediction_output",
            "audit_artifact": "code/particles/runs/calibration/direct_top_bridge_contract.json",
            "bridge_status": direct_top["status"],
        },
    }


def render_markdown(payload: dict[str, Any]) -> str:
    lines = [
        "# Particle Pipeline Status",
        "",
        f"Generated: `{payload['generated_utc']}`",
        "",
        f"Scope: `{payload['scope']}`",
        f"Claim label: `{payload['claim_status']}`",
        "",
        "## P Closure",
        "",
        f"- Candidate `P`: `{payload['p_closure']['P']}`",
        f"- Candidate `alpha^-1`: `{payload['p_closure']['alpha_inv']}`",
        f"- Claim label: `{payload['p_closure']['claim_status']}`",
        f"- May feed promoted particle predictions: `{payload['p_closure']['may_feed_live_particle_predictions']}`",
        "",
        "## Particle-Five Receipts",
        "",
        "| Receipt label | Closable | Local artifact | Worker policy |",
        "| --- | --- | --- | --- |",
    ]
    for gate in payload["particle_five_issue_gates"]:
        lines.append(
            f"| `{_display_status(gate['state'])}` | `{gate['closable_now']}` | "
            f"`{gate['local_next_artifact']}` | {gate['chrome_workers']} |"
        )
    companion_open_branches = payload.get("companion_open_branches") or []
    if companion_open_branches:
        lines.extend(
            [
                "",
                "## Companion Claim Boundaries",
                "",
                "| Topic | Claim label | Boundary | Gate |",
                "| --- | --- | --- | --- |",
            ]
        )
        for branch in companion_open_branches:
            lines.append(
                f"| {branch['label']} | `{_display_status(branch['state'])}` | {branch['summary']} | {branch['next_action']} |"
            )
    lines.extend(
        [
            "",
            "## Predictions",
            "",
            "| Particle | Prediction | Claim label | Scope | Promotable |",
            "| --- | ---: | --- | --- | --- |",
        ]
    )
    for entry in payload["predictions"]:
        lines.append(
            f"| `{entry['particle_id']}` | `{entry['value']} {entry['unit']}` | "
            f"`{entry['exact_kind']}` | `{entry['scope']}` | `{entry['promotable']}` |"
        )
    carrier_modes = payload.get("classical_carrier_modes") or []
    if carrier_modes:
        lines.extend(
            [
                "",
                "## Separated Classical Carrier Modes",
                "",
                "These rows are zero hard parameters in declared quadratic actions, not `0 GeV` quantum-particle predictions.",
                "",
                "| Carrier | Hard parameter squared | Classical gate | Quantum gate | Particle promotion |",
                "| --- | ---: | --- | --- | --- |",
            ]
        )
        for row in carrier_modes:
            lines.append(
                f"| `{row['carrier_id']}` | `{row['hard_quadratic_mass_parameter_squared']}` | "
                f"`{row['classical_carrier_gate']['status']}` | "
                f"`{row['quantum_particle_gate']['status']}` | "
                f"`{row['particle_promotion_allowed']}` |"
            )
    withheld = payload.get("withheld_non_prediction_rows") or []
    if withheld:
        lines.extend(
            [
                "",
                "## Withheld Non-Prediction Rows",
                "",
                "These rows are retained in audit surfaces but are not numeric predictions.",
                "",
                "| Particle | Claim label | Reason | Missing gate |",
                "| --- | --- | --- | --- |",
            ]
        )
        for row in withheld:
            missing_gate = ", ".join(row.get("missing_for_promotion") or [])
            if not missing_gate:
                missing_gate = "n/a"
            claim_label = row.get("claim_tier", row["exact_kind"])
            lines.append(
                f"| `{row['particle_id']}` | `{claim_label}` | {row['reason']} | {missing_gate} |"
            )
    charged_boundary = payload.get("charged_lepton_anchor_boundary") or {}
    if charged_boundary.get("artifact"):
        lines.extend(
            [
                "",
                "## Charged-Lepton Anchor Boundary",
                "",
                f"- Artifact: `{charged_boundary['artifact']}`",
                f"- Status: `{charged_boundary['status']}`",
                f"- Required identity: `{charged_boundary['required_identity']}`",
                f"- Equivalent defect: `{charged_boundary['equivalent_defect']}`",
                f"- Closed downstream: `{charged_boundary.get('current_closed_chain', {}).get('A_ch_to_charged_masses')}`",
                f"- Closed upstream from `P`: `{charged_boundary.get('current_closed_chain', {}).get('P_to_A_ch')}`",
            ]
        )
    quark_boundary = payload.get("quark_sigma_source_boundary") or {}
    if quark_boundary.get("artifact"):
        lines.extend(
            [
                "",
                "## Quark Sigma Source Boundary",
                "",
                f"- Artifact: `{quark_boundary['artifact']}`",
                f"- Status: `{quark_boundary['status']}`",
                f"- Claim tier: `{quark_boundary['claim_tier']}`",
                f"- Required identity: `{quark_boundary['required_identity']}`",
                f"- Source-only sigma emitted: `{quark_boundary['source_only_sigma_emitted']}`",
                f"- Closed downstream algebra: `{quark_boundary['downstream_algebra_closed']}`",
                f"- Closure kind: `{quark_boundary.get('closure_kind')}`",
                f"- Obstruction: `{quark_boundary.get('obstruction_artifact')}`",
                f"- Compatible spread fiber: `{quark_boundary.get('compatible_spread_fiber')}` "
                f"(dimension `{quark_boundary.get('compatible_spread_fiber_dimension')}`)",
                f"- Independent unselected coordinates: "
                f"`{quark_boundary.get('independent_unselected_coordinates', [])}`",
                f"- Missing gates: `{quark_boundary.get('missing_for_promotion', [])}`",
            ]
        )
    axiom_boundary = payload.get("quark_axiom_level_nondefinability") or {}
    if axiom_boundary.get("artifact"):
        mar_audit = axiom_boundary.get("MAR_audit") or {}
        counterfamily = axiom_boundary.get("counterfamily") or {}
        lines.extend(
            [
                "",
                "## Quark Axiom-Level Non-Definability",
                "",
                f"- Artifact: `{axiom_boundary['artifact']}`",
                f"- Status: `{axiom_boundary.get('proof_status')}`",
                f"- Additional axioms used: `{axiom_boundary.get('additional_axioms_used')}`",
                f"- Counterfamily: `{counterfamily.get('parameter_space')}`",
                f"- MAR complexity vector: `{mar_audit.get('complexity_vector')}`",
                f"- Equal MAR score across counterfamily: "
                f"`{mar_audit.get('counterfamily_members_have_equal_MAR_score')}`",
                f"- Public numeric quark rows allowed: "
                f"`{axiom_boundary.get('public_numeric_quark_rows_allowed')}`",
            ]
        )
    scheme_boundary = payload.get("quark_scheme_and_yukawa_boundary") or {}
    if scheme_boundary.get("artifact"):
        matrix_audit = scheme_boundary.get("stored_matrix_dimensional_audit") or {}
        lines.extend(
            [
                "",
                "## Quark Scheme and Yukawa Boundary",
                "",
                f"- Artifact: `{scheme_boundary['artifact']}`",
                f"- Status: `{scheme_boundary.get('proof_status')}`",
                f"- Row partition: `{scheme_boundary.get('row_partition', {})}`",
                f"- Matrix classification: `{matrix_audit.get('current_classification')}`",
                f"- Certified physical Yukawa matrices: "
                f"`{matrix_audit.get('certified_physical_yukawa_matrices')}`",
                f"- Public numeric rows allowed: "
                f"`{(scheme_boundary.get('closure_effect') or {}).get('public_numeric_quark_rows_allowed')}`",
            ]
        )
    direct = payload["direct_top_auxiliary_comparison"]
    hierarchy = payload["hierarchy_and_naturality"]
    hierarchy_status = hierarchy["status"]
    hierarchy_bridge = hierarchy["local_global_bridge"]
    hierarchy_factors = hierarchy["factor_origins"]
    pixel_screen = hierarchy["pixel_screen_resonance"]
    tile_identity = pixel_screen.get("tile_identity") or {}
    ds_coordinate = pixel_screen.get("dimensionless_de_sitter_coordinate") or {}
    lines.extend(
        [
            "",
            "## Fine Structure",
            "",
            "| Output class | alpha^-1 coordinate | P | Difference to endpoint | Claim label |",
            "| --- | ---: | ---: | ---: | --- |",
        ]
    )
    fine = payload["fine_structure"]
    source_no_hadron = fine["source_side_no_hadron_near_endpoint"]
    root_audit = fine["root_only_audit"]
    empirical = fine["oph_plus_empirical_hadron_closure"]
    lines.append(
        f"| `source_side_no_hadron_near_endpoint` | `{source_no_hadron['alpha_inv']}` | "
        f"`{source_no_hadron['P']}` | "
        f"`{source_no_hadron['missing_hadronic_correction_alpha_inv']}` | "
        f"`{source_no_hadron['status']}` |"
    )
    lines.append(
        f"| `oph_plus_empirical_hadron_closure` | `{empirical['alpha_inv']}` | "
        f"`{empirical['P']}` | `0` | `{empirical['status']}` |"
    )
    lines.append(
        f"| `root_only_audit` | `{root_audit['alpha_inv']}` | `{root_audit['P']}` | "
        f"`{root_audit['missing_inverse_alpha_units_to_thomson_endpoint']}` | "
        f"`{root_audit['status']}` |"
    )
    lines.extend(
        [
            "",
            f"- Diagnostic formula: `{source_no_hadron['formula']}`",
            f"- Relative coordinate difference to the endpoint: "
            f"`{source_no_hadron['relative_shortfall']}` "
            f"(`{source_no_hadron['percent_shortfall']}` percent)",
            f"- Open physical-map requirements: "
            f"`{source_no_hadron['open_physical_map_requirements']}`",
            "- Source-only physical Thomson endpoint emitted: `False`",
            f"- Empirical payload policy: `{fine['empirical_payload_policy']['dispersion_payload_status']}`",
            "",
            "## Hierarchy And Naturality",
            "",
            f"- Resonance label: `{_display_status(hierarchy_status['resonance_status'])}`",
            f"- Full theorem-grade resonance promoted: `{hierarchy_status['full_theorem_grade_resonance_promoted']}`",
            f"- Work-in-progress receipts: `{hierarchy_status['work_in_progress_receipts']}`",
            f"- Conditional EW bridge capacity (modulo F and CP-1 to CP-3): "
            f"`{hierarchy_bridge['N_CRC_EW']}`",
            f"- Bridge residual: `{hierarchy_bridge['bridge_residual']}`",
            f"- Source `v/E_cell`: `{hierarchy_bridge['v_over_E_cell_source']}`",
            f"- Factor origins: `ports={hierarchy_factors['screen_ports']}`, "
            f"`m_rep={hierarchy_factors['m_rep']}`, "
            f"`exponent={hierarchy_factors['specialized_exponent']}`",
            f"- Count independence: {hierarchy['claim_boundary']['count_independence']}",
            f"- Missing hierarchy bridge: {hierarchy['claim_boundary']['missing_bridge']}",
            f"- Higgs naturality defect: `epsilon_H={hierarchy_factors['higgs_naturality_defect']}`",
            f"- Boundary: {hierarchy['claim_boundary']['improves']}",
            f"- Physical capacity conditions: "
            f"{'; '.join(hierarchy['claim_boundary']['capacity_conditions'])}",
            f"- Not promoted by this bridge: {', '.join(hierarchy['claim_boundary']['does_not_promote'])}",
            "",
            "## Pixel-Screen Resonance",
            "",
            f"- Receipt label: `{_display_status(pixel_screen['status'])}`",
            f"- Accepted: `{pixel_screen['accepted']}`",
            f"- Cell count: `{tile_identity.get('cell_count')}`",
            f"- Cell entropy: `{tile_identity.get('cell_entropy')}`",
            f"- Capacity reconstruction error: `{tile_identity.get('relative_reconstruction_error')}`",
            f"- Dimensionless `Lambda*l_star^2`: `{ds_coordinate.get('Lambda_lstar2')}`",
            f"- Dimensionless `Lambda*a_cell`: `{ds_coordinate.get('Lambda_a_cell')}`",
            f"- Boundary: {pixel_screen.get('claim_boundary', {}).get('closed_here')}",
            "",
            "## Direct-Top Auxiliary Comparison",
            "",
            f"- Primary codomain: `{direct['primary_top_codomain']}`",
            f"- Auxiliary codomain: `{direct['auxiliary_direct_top_codomain']}`",
            f"- Value policy: `{direct['value_policy']}`",
            f"- Audit artifact: `{direct['audit_artifact']}`",
            f"- Bridge label: `{_display_status(direct['bridge_status'])}`",
            "",
            "## Hadrons",
            "",
            f"- Source-only hadron predictions emitted: `{payload['hadron_policy']['source_only_hadron_predictions_emitted']}`",
            f"- Empirical hadron closure allowed for display: `{payload['hadron_policy']['empirical_hadron_closure_allowed_for_display']}`",
            f"- Policy artifact: `{payload['hadron_policy']['policy_artifact']}`",
            f"- Reason: {payload['hadron_policy']['reason']}",
        ]
    )
    return "\n".join(lines)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build final current end-to-end particle predictions.")
    parser.add_argument("--json-out", default=str(DEFAULT_JSON_OUT))
    parser.add_argument("--markdown-out", default=str(DEFAULT_MD_OUT))
    parser.add_argument("--print-json", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    payload = build_payload()
    json_text = json.dumps(payload, indent=2, sort_keys=True) + "\n"

    json_out = Path(args.json_out)
    json_out.parent.mkdir(parents=True, exist_ok=True)
    json_out.write_text(json_text, encoding="utf-8")

    markdown_out = Path(args.markdown_out)
    markdown_out.write_text(render_markdown(payload) + "\n", encoding="utf-8")

    if args.print_json:
        print(json_text, end="")
    else:
        print(f"saved: {json_out}")
        print(f"saved: {markdown_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
