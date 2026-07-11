#!/usr/bin/env python3
"""Build the single closure-status manifest for the particle pipeline.

This issue gate keeps source-only rows, empirical hadron closure rows,
compare-only rows, and work-in-progress rows mechanically separate.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
PARTICLES_ROOT = ROOT / "particles"
P_DERIVATION_ROOT = ROOT / "P_derivation"
DEFAULT_JSON_OUT = PARTICLES_ROOT / "runs" / "status" / "particle_pipeline_closure_status.json"
DEFAULT_MD_OUT = PARTICLES_ROOT / "PARTICLE_PIPELINE_CLOSURE_STATUS.md"

P_TRUNK = P_DERIVATION_ROOT / "runtime" / "p_closure_trunk_current.json"
MEASURED_ENDPOINT = P_DERIVATION_ROOT / "runtime" / "measured_endpoint_calibration_current.json"
THOMSON_CONTRACT = P_DERIVATION_ROOT / "runtime" / "thomson_endpoint_contract_current.json"
THOMSON_PACKAGE = P_DERIVATION_ROOT / "runtime" / "thomson_endpoint_package_current.json"
SCREENING_NO_GO = P_DERIVATION_ROOT / "runtime" / "screening_invariant_no_go_current.json"
INTERVAL_CERTIFICATE = P_DERIVATION_ROOT / "runtime" / "fine_structure_interval_certificate_current.json"
R_Q_CONTRACT = P_DERIVATION_ROOT / "runtime" / "r_q_residual_contract_current.json"
RG_CONTRACT = P_DERIVATION_ROOT / "runtime" / "rg_matching_threshold_contract_current.json"
DIRECT_TOP_CONTRACT = PARTICLES_ROOT / "runs" / "calibration" / "direct_top_bridge_contract.json"
GAP_LEDGER = PARTICLES_ROOT / "runs" / "status" / "particle_derivation_gap_ledger.json"
BLIND_PROVENANCE = PARTICLES_ROOT / "runs" / "status" / "blind_prediction_provenance.json"
QUARK_CONTRACT = PARTICLES_ROOT / "runs" / "flavor" / "quark_lane_closure_contract.json"
QUARK_GLOBAL_OBSTRUCTION = PARTICLES_ROOT / "runs" / "flavor" / "quark_class_uniform_public_frame_descent_obstruction.json"
CHARGED_NONCLOSURE = PARTICLES_ROOT / "runs" / "leptons" / "charged_end_to_end_impossibility_theorem.json"
NEUTRINO_CONTRACT = PARTICLES_ROOT / "runs" / "neutrino" / "neutrino_lane_closure_contract.json"
HADRON_SPECTRAL_CONTRACT = PARTICLES_ROOT / "runs" / "hadron" / "ward_projected_spectral_measure_contract.json"
EMPIRICAL_EE_REGISTRY = PARTICLES_ROOT / "hadron" / "empirical_ee_hadrons_sources.yaml"
EMPIRICAL_EE_SCHEMA = PARTICLES_ROOT / "hadron" / "empirical_ee_hadronic_spectral_measure.schema.json"
EXACT_NONHADRON = PARTICLES_ROOT / "exact_nonhadron_masses.json"
CARRIER_ACCEPTANCE = PARTICLES_ROOT / "runs" / "status" / "carrier_mode_acceptance.json"
RESULTS_STATUS = PARTICLES_ROOT / "results_status.json"
HIERARCHY_ROOT = PARTICLES_ROOT / "hierarchy"
HIERARCHY_NATURALITY = HIERARCHY_ROOT / "issue_332_rg_naturality_certificate.json"
HIERARCHY_EW_PROJECTION = HIERARCHY_ROOT / "certificates" / "R_EW_tick_projection_certificate.json"
HIERARCHY_EW_CAPACITY = HIERARCHY_ROOT / "certificates" / "R_EW_global_capacity_certificate.json"
HIERARCHY_READBACK = HIERARCHY_ROOT / "certificates" / "R_readback_resolution_certificate.json"
HIERARCHY_M_REP = HIERARCHY_ROOT / "certificates" / "R_m_rep_24_certificate.json"
HIERARCHY_SCREEN_SIEVE = HIERARCHY_ROOT / "certificates" / "R_screen_sieve_icosahedral_certificate.json"
HIERARCHY_JOINT_FIXED_POINT = HIERARCHY_ROOT / "certificates" / "R_PN_joint_fixed_point_certificate_report.json"
HIERARCHY_PIXEL_SCREEN = HIERARCHY_ROOT / "certificates" / "R_pixel_screen_resonance_summary.json"
HIERARCHY_LOCAL_GLOBAL_RESONANCE = (
    HIERARCHY_ROOT / "certificates" / "R_local_global_hierarchy_resonance_closeout_335.json"
)


def _now_utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def _load_json(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def _artifact_status(path: Path, payload: dict[str, Any] | None) -> dict[str, Any]:
    status = None
    if payload is not None:
        status = payload.get("status") or payload.get("proof_status") or payload.get("claim_status")
    return {
        "path": _rel(path),
        "exists": path.exists(),
        "artifact": payload.get("artifact") if payload else None,
        "status": status,
        "promotion_allowed": payload.get("promotion_allowed") if payload else None,
    }


def _display_status(status: str) -> str:
    return status.replace("current_corpus", "corpus_limited")


def _companion_status_branches(gap_ledger: dict[str, Any] | None) -> list[dict[str, Any]]:
    if not gap_ledger:
        return []
    branches: list[dict[str, Any]] = []
    for row in gap_ledger.get("rows", []):
        if row.get("id") != "qcd.strong-cp-angle":
            continue
        branches.append(
            {
                "issue": row.get("github_issue"),
                "label": "Strong CP",
                "state": row.get("status"),
                "summary": row.get("current_boundary"),
                "next_action": row.get("next_action"),
            }
        )
    return branches


def _latest_nonhadron_predictions(exact_payload: dict[str, Any] | None) -> dict[str, dict[str, Any]]:
    if not exact_payload:
        return {}
    predictions: dict[str, dict[str, Any]] = {}
    for entry in exact_payload.get("entries", []):
        if entry.get("mass_gev") is not None:
            predictions[entry["particle_id"]] = {
                "value": float(entry["mass_gev"]),
                "unit": "GeV",
                "exact_kind": entry.get("exact_kind"),
                "scope": entry.get("scope"),
                "promotable": entry.get("promotable"),
            }
        elif entry.get("mass_eV") is not None:
            predictions[entry["particle_id"]] = {
                "value": float(entry["mass_eV"]),
                "unit": "eV",
                "exact_kind": entry.get("exact_kind"),
                "scope": entry.get("scope"),
                "promotable": entry.get("promotable"),
            }
    return predictions


def _carrier_summaries(payload: dict[str, Any] | None) -> list[dict[str, Any]]:
    if not payload:
        return []
    return [
        {
            "carrier_id": row["carrier_id"],
            "label": row["label"],
            "claim_kind": row["claim_kind"],
            "branch": row["branch"],
            "hard_quadratic_mass_parameter_squared": row[
                "hard_quadratic_mass_parameter_squared"
            ],
            "classical_carrier_gate": row["classical_carrier_gate"],
            "quantum_particle_gate": row["quantum_particle_gate"],
            "particle_promotion_allowed": row["particle_promotion_allowed"],
            "source_artifact": "code/particles/runs/status/carrier_mode_acceptance.json",
        }
        for row in payload.get("carriers", [])
    ]


def build_status() -> dict[str, Any]:
    p_trunk = _load_json(P_TRUNK)
    measured_endpoint = _load_json(MEASURED_ENDPOINT)
    thomson = _load_json(THOMSON_CONTRACT)
    thomson_package = _load_json(THOMSON_PACKAGE)
    screening_no_go = _load_json(SCREENING_NO_GO)
    interval_certificate = _load_json(INTERVAL_CERTIFICATE)
    r_q_contract = _load_json(R_Q_CONTRACT)
    rg = _load_json(RG_CONTRACT)
    direct_top = _load_json(DIRECT_TOP_CONTRACT)
    gap_ledger = _load_json(GAP_LEDGER)
    blind_provenance = _load_json(BLIND_PROVENANCE)
    quark = _load_json(QUARK_CONTRACT)
    quark_global = _load_json(QUARK_GLOBAL_OBSTRUCTION)
    charged_nonclosure = _load_json(CHARGED_NONCLOSURE)
    neutrino = _load_json(NEUTRINO_CONTRACT)
    hadron_spectral = _load_json(HADRON_SPECTRAL_CONTRACT)
    exact = _load_json(EXACT_NONHADRON)
    carrier_acceptance = _load_json(CARRIER_ACCEPTANCE)
    results_status = _load_json(RESULTS_STATUS)
    hierarchy_naturality = _load_json(HIERARCHY_NATURALITY)
    hierarchy_ew_projection = _load_json(HIERARCHY_EW_PROJECTION)
    hierarchy_ew_capacity = _load_json(HIERARCHY_EW_CAPACITY)
    hierarchy_readback = _load_json(HIERARCHY_READBACK)
    hierarchy_m_rep = _load_json(HIERARCHY_M_REP)
    hierarchy_screen_sieve = _load_json(HIERARCHY_SCREEN_SIEVE)
    hierarchy_joint_fixed_point = _load_json(HIERARCHY_JOINT_FIXED_POINT)
    hierarchy_pixel_screen = _load_json(HIERARCHY_PIXEL_SCREEN)
    hierarchy_local_global = _load_json(HIERARCHY_LOCAL_GLOBAL_RESONANCE)

    return {
        "artifact": "oph_particle_pipeline_closure_status",
        "generated_utc": _now_utc(),
        "purpose": "Single closure gate for source-only rows and empirical hadron closure rows.",
        "scope": {
            "current_pipeline_scope": "source_only_rows_plus_empirical_hadron_closure_policy",
            "source_only_hadrons_in_current_local_scope": False,
            "empirical_hadron_closure_surface": True,
            "hadron_scope_reason": (
                "Source-only hadron rows require a working OPH hardware backend such as "
                "GLORB/Echosahedron. Empirical hadron closure stays in a separate output class; "
                "the e+e- spectral payload has a source registry and schema."
            ),
            "chrome_workers_needed": False,
        },
        "current_surface": {
            "builder": "code/particles/compute_current_output_table.py",
            "default_runtime_root": "temp/particles_runtime",
            "source_repo": "reverse-engineering-reality/code",
            "simplification": (
                "The prediction pipeline keeps source-only OPH, OPH plus empirical hadron closure, "
                "compare-only, and work-in-progress rows mechanically distinct."
            ),
        },
        "artifacts": {
            "p_trunk": _artifact_status(P_TRUNK, p_trunk),
            "measured_endpoint_calibration": _artifact_status(MEASURED_ENDPOINT, measured_endpoint),
            "thomson_endpoint_contract": _artifact_status(THOMSON_CONTRACT, thomson),
            "thomson_endpoint_package": _artifact_status(THOMSON_PACKAGE, thomson_package),
            "screening_invariant_no_go": _artifact_status(SCREENING_NO_GO, screening_no_go),
            "fine_structure_interval_certificate": _artifact_status(INTERVAL_CERTIFICATE, interval_certificate),
            "r_q_residual_contract": _artifact_status(R_Q_CONTRACT, r_q_contract),
            "rg_matching_threshold_contract": _artifact_status(RG_CONTRACT, rg),
            "direct_top_bridge_contract": _artifact_status(DIRECT_TOP_CONTRACT, direct_top),
            "gap_ledger": _artifact_status(GAP_LEDGER, gap_ledger),
            "blind_prediction_provenance": _artifact_status(BLIND_PROVENANCE, blind_provenance),
            "carrier_mode_acceptance": _artifact_status(CARRIER_ACCEPTANCE, carrier_acceptance),
            "quark_lane_closure_contract": _artifact_status(QUARK_CONTRACT, quark),
            "quark_global_classification_obstruction": _artifact_status(QUARK_GLOBAL_OBSTRUCTION, quark_global),
            "charged_end_to_end_impossibility_theorem": _artifact_status(CHARGED_NONCLOSURE, charged_nonclosure),
            "neutrino_lane_closure_contract": _artifact_status(NEUTRINO_CONTRACT, neutrino),
            "hadron_spectral_measure_contract": _artifact_status(HADRON_SPECTRAL_CONTRACT, hadron_spectral),
            "empirical_ee_hadrons_source_registry": {
                "path": _rel(EMPIRICAL_EE_REGISTRY),
                "exists": EMPIRICAL_EE_REGISTRY.exists(),
                "status": "source_registry_present",
                "row_class": "oph_plus_empirical_hadron_closure",
            },
            "empirical_ee_hadronic_spectral_measure_schema": {
                "path": _rel(EMPIRICAL_EE_SCHEMA),
                "exists": EMPIRICAL_EE_SCHEMA.exists(),
                "status": "schema_present",
                "row_class": "oph_plus_empirical_hadron_closure",
            },
            "hierarchy_rg_higgs_naturality": _artifact_status(HIERARCHY_NATURALITY, hierarchy_naturality),
            "hierarchy_ew_projection": _artifact_status(HIERARCHY_EW_PROJECTION, hierarchy_ew_projection),
            "hierarchy_ew_capacity": _artifact_status(HIERARCHY_EW_CAPACITY, hierarchy_ew_capacity),
            "hierarchy_readback_resolution": _artifact_status(HIERARCHY_READBACK, hierarchy_readback),
            "hierarchy_m_rep_24": _artifact_status(HIERARCHY_M_REP, hierarchy_m_rep),
            "hierarchy_screen_sieve": _artifact_status(HIERARCHY_SCREEN_SIEVE, hierarchy_screen_sieve),
            "hierarchy_joint_fixed_point": _artifact_status(
                HIERARCHY_JOINT_FIXED_POINT, hierarchy_joint_fixed_point
            ),
            "hierarchy_pixel_screen_resonance": _artifact_status(
                HIERARCHY_PIXEL_SCREEN, hierarchy_pixel_screen
            ),
            "hierarchy_local_global_resonance": _artifact_status(
                HIERARCHY_LOCAL_GLOBAL_RESONANCE, hierarchy_local_global
            ),
        },
        "issue_gates": [
            {
                "issue": 536,
                "title": "Classical carrier versus quantum-particle gate",
                "state": "closed_claim_scope_repaired_quantum_particle_gate_fail_closed",
                "closable_now": True,
                "local_next_artifact": _rel(CARRIER_ACCEPTANCE),
                "classical_carrier_modes_recorded": True,
                "zero_gev_particle_rows_emitted": False,
                "quantum_particle_promotion_allowed": False,
                "chrome_workers": "not_needed_for_analytic_claim_gate",
            },
            {
                "issue": 332,
                "title": "RG/Higgs naturality closure",
                "state": "closed_exact_selected_branch",
                "closable_now": True,
                "local_next_artifact": _rel(HIERARCHY_NATURALITY),
                "epsilon_H": (hierarchy_naturality or {}).get("epsilon_H"),
                "chrome_workers": "not_needed_for_closed_certificate",
            },
            {
                "issue": 337,
                "title": "Electroweak projection bridge",
                "state": "closed_projection_bridge_with_exact_residual",
                "closable_now": True,
                "local_next_artifact": _rel(HIERARCHY_EW_PROJECTION),
                "chrome_workers": "not_needed_for_closed_certificate",
            },
            {
                "issue": 335,
                "title": "Local/global hierarchy resonance",
                "state": "closed_full_local_global_hierarchy_resonance",
                "closable_now": True,
                "local_next_artifact": _rel(HIERARCHY_LOCAL_GLOBAL_RESONANCE),
                "full_theorem_grade_resonance_promoted": bool(
                    (hierarchy_local_global or {}).get("full_theorem_grade_resonance_promoted", False)
                ),
                "chrome_workers": "not_needed_for_closed_certificate",
            },
            {
                "issue": 223,
                "title": "Ward-projected Thomson endpoint package",
                "state": "closed_blocker_isolated_endpoint_package",
                "closable_now": True,
                "local_next_artifact": _rel(THOMSON_PACKAGE),
                "contract_artifact": _rel(THOMSON_CONTRACT),
                "closed_as_blocker_isolation": True,
                "successor_issue": 235,
                "promotion_allowed": False,
                "chrome_workers": "not_needed_for_closed_package",
            },
            {
                "issue": 235,
                "title": "Source spectral measure payload and interval certificate",
                "state": "closed_blocker_isolated_source_residual_no_go",
                "closable_now": True,
                "local_next_artifact": _rel(THOMSON_CONTRACT),
                "package_artifact": _rel(THOMSON_PACKAGE),
                "no_go_artifact": _rel(SCREENING_NO_GO),
                "interval_certificate_artifact": _rel(INTERVAL_CERTIFICATE),
                "residual_contract_artifact": _rel(R_Q_CONTRACT),
                "closed_as_first_missing_lemma_isolated": True,
                "promotion_allowed": False,
                "first_missing_lemma": "source-emitted same-scheme Ward-projected R_Q(P)",
                "source_spectral_reduction": "source_spectral_reduction_theorem_emitted_measure_payload_absent",
                "minimal_new_payload": "oph_qcd_ward_projected_hadronic_spectral_measure",
                "hadron_dependency_hardware_gated": True,
                "chrome_workers": "not_needed_until_source_spectral_measure_payload_exists",
            },
            {
                "issue": 224,
                "title": "Adopt certified derived P closure root",
                "state": "closed_canonical_guarded_trunk_adoption",
                "closable_now": True,
                "local_next_artifact": _rel(P_TRUNK),
                "closed_as_guarded_candidate_adoption": True,
                "promotion_allowed": bool((p_trunk or {}).get("consumer_policy", {}).get("may_feed_live_particle_predictions", False)),
                "stage_gate": "populated source spectral measure payload plus full interval certificate before live particle promotion",
                "chrome_workers": "not_needed_for_guarded_codepath_closure",
            },
            {
                "issue": 225,
                "title": "Synchronize derived P closure values across material surfaces",
                "state": "closed_material_sync_no_live_publish",
                "closable_now": True,
                "local_next_artifact": "paper/deriving_the_particle_zoo_from_observer_consistency.tex",
                "closed_as_material_update": True,
                "publish_performed": False,
                "chrome_workers": "not_needed_for_material_sync",
            },
            {
                "issue": 32,
                "title": "RG matching and threshold structure",
                "state": "closed_declared_convention_contract",
                "closable_now": True,
                "local_next_artifact": _rel(RG_CONTRACT),
                "closed_as_declared_convention_contract": True,
                "promotion_allowed": False,
                "chrome_workers": "not_needed_for_closed_contract",
            },
            {
                "issue": 153,
                "title": "Hadron backend and systematics",
                "state": "closed_out_of_scope_computationally_blocked",
                "closable_now": True,
                "local_next_artifact": _rel(HADRON_SPECTRAL_CONTRACT),
                "requires_oph_hardware_backend": True,
                "empirical_hadron_closure_allowed": True,
                "closed_as_out_of_scope": True,
                "close_reason": (
                    "Source-only hadron prediction is outside the local environment. Empirical "
                    "hadron closure has a separate output class with an e+e- source registry and schema."
                ),
                "chrome_workers": "do_not_use_for_backend_execution",
            },
            {
                "issue": 157,
                "title": "Nonperturbative QCD hadron branch",
                "state": "closed_out_of_scope_computationally_blocked",
                "closable_now": True,
                "local_next_artifact": _rel(HADRON_SPECTRAL_CONTRACT),
                "requires_oph_hardware_backend": True,
                "empirical_hadron_closure_allowed": True,
                "closed_as_out_of_scope": True,
                "close_reason": (
                    "The compact/paper source-only hadron branch is outside the computational "
                    "scope. The empirical hadron closure surface is the data-driven display path "
                    "with an e+e- source registry and schema."
                ),
                "chrome_workers": "do_not_use_for_backend_execution",
            },
            {
                "issue": 201,
                "title": "Charged determinant trace-lift attachment",
                "state": "closed_current_corpus_charged_end_to_end_no_go",
                "closable_now": True,
                "local_next_artifact": _rel(CHARGED_NONCLOSURE),
                "public_charged_masses_emitted": False,
                "chrome_workers": "not_needed_until_new_uncentered_trace_lift_source_exists",
            },
            {
                "issue": 207,
                "title": "Direct-top codomain bridge",
                "state": "closed_current_corpus_codomain_no_go",
                "closable_now": True,
                "local_next_artifact": _rel(DIRECT_TOP_CONTRACT),
                "chrome_workers": "not_needed_until_new_response_kernel_source_exists",
            },
            {
                "issue": 234,
                "title": "Blind-prediction provenance and convention-sensitivity audits",
                "state": "closed_provenance_ledger_and_declared_sensitivity_taxonomy",
                "closable_now": True,
                "local_next_artifact": _rel(BLIND_PROVENANCE),
                "closed_as_material_audit": True,
                "chrome_workers": "not_needed_for_closed_provenance_taxonomy",
            },
            {
                "issue": 117,
                "title": "Neutrino splittings and mixing",
                "state": (
                    "source_closed_prediction"
                    if neutrino.get("public_promotion_allowed") is True
                    else "rejected_candidate_source_basis_and_kernel_open"
                ),
                "closable_now": neutrino.get("public_promotion_allowed") is True,
                "local_next_artifact": _rel(NEUTRINO_CONTRACT),
                "chrome_workers": "not_needed_for_finite_audit_repairs",
            },
            {
                "issue": 198,
                "title": "Selected public quark-frame sigma descent",
                "state": "closed_selected_class_scope_visible",
                "closable_now": True,
                "local_next_artifact": _rel(QUARK_CONTRACT),
                "chrome_workers": "not_needed",
            },
            {
                "issue": 199,
                "title": "Class-uniform public quark-frame descent",
                "state": "closed_current_corpus_global_classification_no_go",
                "closable_now": True,
                "local_next_artifact": _rel(QUARK_GLOBAL_OBSTRUCTION),
                "selected_class_theorem_preserved": True,
                "chrome_workers": "not_needed_until_new_global_public_frame_classifier_source_exists",
            },
            {
                "issue": 155,
                "title": "Strong-CP branch status",
                "state": "open_theta_qcd_bar_theta_vanishing_gap",
                "closable_now": False,
                "local_next_artifact": _rel(GAP_LEDGER),
                "public_status_only": True,
                "chrome_workers": "not_needed_until_a_concrete_strong_cp_packet_exists",
            },
        ],
        "companion_status_branches": _companion_status_branches(gap_ledger),
        "finalization_gates": {
            "nonhadron_prediction_surface_buildable": True,
            "source_only_hadrons_suppressed_by_default": bool(
                (results_status or {}).get("inputs", {}).get("hadron_profile", "suppressed") == "suppressed"
            ),
            "empirical_hadron_closure_policy_documented": True,
            "empirical_hadron_spectral_dataset_integrated": False,
            "p_trunk_candidate_only": not bool(
                (p_trunk or {}).get("consumer_policy", {}).get("may_feed_live_particle_predictions", False)
            ),
            "obstruction_only_worker_result_allowed": True,
            "paper_material_sync_complete_without_live_publish": True,
            "source_spectral_stage_gate": "populated source spectral measure payload plus interval certificate",
            "hierarchy_local_global_resonance_closed": bool(
                (hierarchy_local_global or {}).get("full_theorem_grade_resonance_promoted", False)
            ),
            "higgs_naturality_defect_closed": (hierarchy_naturality or {}).get("epsilon_H") == "0",
            "pixel_screen_resonance_summary_closed": bool(
                (hierarchy_pixel_screen or {}).get("accepted", False)
            ),
            "symmetry_only_particle_promotion_blocked": bool(
                carrier_acceptance
                and carrier_acceptance.get("abstract_symmetry_group_alone_passes_quantum_gate") is False
                and all(
                    row.get("particle_promotion_allowed") is False
                    for row in carrier_acceptance.get("carriers", [])
                )
            ),
        },
        "latest_nonhadron_predictions": _latest_nonhadron_predictions(exact),
        "withheld_non_prediction_rows": (exact or {}).get("withheld_entries", []),
        "classical_carrier_modes": _carrier_summaries(carrier_acceptance),
    }


def render_markdown(status: dict[str, Any]) -> str:
    lines = [
        "# Particle Pipeline Closure Status",
        "",
        f"Generated: `{status['generated_utc']}`",
        "",
        status["purpose"],
        "",
        "## Scope",
        "",
        f"- Scope: `{status['scope']['current_pipeline_scope']}`",
        f"- Source-only hadrons in local scope: `{status['scope']['source_only_hadrons_in_current_local_scope']}`",
        f"- Empirical hadron closure surface: `{status['scope']['empirical_hadron_closure_surface']}`",
        f"- Chrome workers needed: `{status['scope']['chrome_workers_needed']}`",
        f"- Hadron scope reason: {status['scope']['hadron_scope_reason']}",
        "",
        "## Receipt Gates",
        "",
        "| Receipt label | Closable | Local artifact | Worker policy |",
        "| --- | --- | --- | --- |",
    ]
    for gate in status["issue_gates"]:
        lines.append(
            f"| `{_display_status(gate['state'])}` | `{gate['closable_now']}` | "
            f"`{gate['local_next_artifact']}` | {gate['chrome_workers']} |"
        )

    companion_status_branches = status.get("companion_status_branches") or []
    if companion_status_branches:
        lines.extend(
            [
                "",
                "## Companion Claim Boundaries",
                "",
                "| Topic | Claim label | Boundary | Gate |",
                "| --- | --- | --- | --- |",
            ]
        )
        for branch in companion_status_branches:
            lines.append(
                f"| {branch['label']} | `{_display_status(branch['state'])}` | {branch['summary']} | {branch['next_action']} |"
            )

    lines.extend(
        [
            "",
            "## Latest Non-Hadron Predictions",
            "",
            "| Particle ID | Mass |",
            "| --- | ---: |",
        ]
    )
    for particle_id, prediction in sorted(status["latest_nonhadron_predictions"].items()):
        lines.append(f"| `{particle_id}` | `{prediction['value']} {prediction['unit']}` |")

    carrier_modes = status.get("classical_carrier_modes") or []
    if carrier_modes:
        lines.extend(
            [
                "",
                "## Conditional Classical Carrier Modes",
                "",
                "| Carrier | Hard parameter squared | Classical gate | Quantum gate |",
                "| --- | ---: | --- | --- |",
            ]
        )
        for row in carrier_modes:
            lines.append(
                f"| `{row['carrier_id']}` | `{row['hard_quadratic_mass_parameter_squared']}` | "
                f"`{row['classical_carrier_gate']['status']}` | "
                f"`{row['quantum_particle_gate']['status']}` |"
            )

    withheld = status.get("withheld_non_prediction_rows") or []
    if withheld:
        lines.extend(
            [
                "",
                "## Withheld Non-Prediction Rows",
                "",
                "| Particle ID | Claim label | Reason |",
                "| --- | --- | --- |",
            ]
        )
        for row in withheld:
            lines.append(f"| `{row['particle_id']}` | `{row['exact_kind']}` | {row['reason']} |")

    lines.extend(
        [
            "",
            "## Finalization Gates",
            "",
        ]
    )
    for key, value in status["finalization_gates"].items():
        lines.append(f"- `{key}`: `{value}`")
    return "\n".join(lines)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build the particle pipeline closure status manifest.")
    parser.add_argument("--json-out", default=str(DEFAULT_JSON_OUT))
    parser.add_argument("--markdown-out", default=str(DEFAULT_MD_OUT))
    parser.add_argument("--print-json", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    status = build_status()
    json_text = json.dumps(status, indent=2, sort_keys=True) + "\n"

    json_out = Path(args.json_out)
    json_out.parent.mkdir(parents=True, exist_ok=True)
    json_out.write_text(json_text, encoding="utf-8")

    markdown_out = Path(args.markdown_out)
    markdown_out.write_text(render_markdown(status) + "\n", encoding="utf-8")

    if args.print_json:
        print(json_text, end="")
    else:
        print(f"saved: {json_out}")
        print(f"saved: {markdown_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
