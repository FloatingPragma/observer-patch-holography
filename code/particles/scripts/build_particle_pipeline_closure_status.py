#!/usr/bin/env python3
"""Build the single closure-status manifest for the particle pipeline.

This is the simplified issue gate for the particle pipeline: non-hadron
particle predictions are finalized through the disposable runtime surface,
while hadron production is closed out-of-scope pending OPH backend hardware.
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
EXACT_NONHADRON = PARTICLES_ROOT / "exact_nonhadron_masses.json"
RESULTS_STATUS = PARTICLES_ROOT / "results_status.json"


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


def build_status() -> dict[str, Any]:
    p_trunk = _load_json(P_TRUNK)
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
    results_status = _load_json(RESULTS_STATUS)

    return {
        "artifact": "oph_particle_pipeline_closure_status",
        "generated_utc": _now_utc(),
        "purpose": "Single simplified closure gate for the non-hadron particle pipeline.",
        "scope": {
            "current_pipeline_scope": "nonhadron_particles_plus_candidate_P_root_metadata",
            "hadrons_in_current_local_scope": False,
            "hadron_scope_reason": (
                "Production hadrons require a working OPH hardware backend such as GLORB/Echosahedron. "
                "Issues #153 and #157 are closed as out-of-scope/computationally blocked for the "
                "pipeline; local surrogate code and Chrome workers are non-promoting."
            ),
            "chrome_workers_needed_now": False,
        },
        "current_surface": {
            "builder": "code/particles/compute_current_output_table.py",
            "default_runtime_root": "temp/particles_runtime",
            "source_repo": "reverse-engineering-reality/code",
            "simplification": (
                "The prediction pipeline is one disposable runtime surface plus this closure "
                "manifest, with hadron promotion excluded by scope."
            ),
        },
        "artifacts": {
            "p_trunk": _artifact_status(P_TRUNK, p_trunk),
            "thomson_endpoint_contract": _artifact_status(THOMSON_CONTRACT, thomson),
            "thomson_endpoint_package": _artifact_status(THOMSON_PACKAGE, thomson_package),
            "screening_invariant_no_go": _artifact_status(SCREENING_NO_GO, screening_no_go),
            "fine_structure_interval_certificate": _artifact_status(INTERVAL_CERTIFICATE, interval_certificate),
            "r_q_residual_contract": _artifact_status(R_Q_CONTRACT, r_q_contract),
            "rg_matching_threshold_contract": _artifact_status(RG_CONTRACT, rg),
            "direct_top_bridge_contract": _artifact_status(DIRECT_TOP_CONTRACT, direct_top),
            "gap_ledger": _artifact_status(GAP_LEDGER, gap_ledger),
            "blind_prediction_provenance": _artifact_status(BLIND_PROVENANCE, blind_provenance),
            "quark_lane_closure_contract": _artifact_status(QUARK_CONTRACT, quark),
            "quark_global_classification_obstruction": _artifact_status(QUARK_GLOBAL_OBSTRUCTION, quark_global),
            "charged_end_to_end_impossibility_theorem": _artifact_status(CHARGED_NONCLOSURE, charged_nonclosure),
            "neutrino_lane_closure_contract": _artifact_status(NEUTRINO_CONTRACT, neutrino),
            "hadron_spectral_measure_contract": _artifact_status(HADRON_SPECTRAL_CONTRACT, hadron_spectral),
        },
        "issue_gates": [
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
                "closed_as_out_of_scope": True,
                "close_reason": (
                    "The local environment lacks a working OPH hadron backend; reopen only when "
                    "GLORB/Echosahedron-class backend output and production systematics exist."
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
                "closed_as_out_of_scope": True,
                "close_reason": (
                    "The compact/paper hadron branch is outside the computational scope pending "
                    "a working OPH hadron backend that emits the Ward-projected spectral measure."
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
                "state": "closed_keep_visible_comparison_tension",
                "closable_now": True,
                "local_next_artifact": _rel(NEUTRINO_CONTRACT),
                "chrome_workers": "not_needed",
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
            "hadrons_suppressed_by_default": bool(
                (results_status or {}).get("inputs", {}).get("hadron_profile", "suppressed") == "suppressed"
            ),
            "p_trunk_candidate_only": not bool(
                (p_trunk or {}).get("consumer_policy", {}).get("may_feed_live_particle_predictions", False)
            ),
            "obstruction_only_worker_result_allowed": True,
            "paper_material_sync_complete_without_live_publish": True,
            "source_spectral_stage_gate": "populated source spectral measure payload plus interval certificate",
        },
        "latest_nonhadron_predictions": _latest_nonhadron_predictions(exact),
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
        f"- Hadrons in local scope: `{status['scope']['hadrons_in_current_local_scope']}`",
        f"- Chrome workers needed: `{status['scope']['chrome_workers_needed_now']}`",
        f"- Hadron scope reason: {status['scope']['hadron_scope_reason']}",
        "",
        "## Issue Gates",
        "",
        "| Issue | State | Closable | Local next artifact | Chrome policy |",
        "| --- | --- | --- | --- | --- |",
    ]
    for gate in status["issue_gates"]:
        lines.append(
            f"| #{gate['issue']} | `{_display_status(gate['state'])}` | `{gate['closable_now']}` | "
            f"`{gate['local_next_artifact']}` | {gate['chrome_workers']} |"
        )

    companion_status_branches = status.get("companion_status_branches") or []
    if companion_status_branches:
        lines.extend(
            [
                "",
                "## Companion Status Branches",
                "",
                "| Topic | State | Current boundary | Next action |",
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
