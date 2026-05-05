#!/usr/bin/env python3
"""Build the single closure-status manifest for the particle pipeline.

This is the simplified issue gate for the current pipeline: non-hadron particle
predictions are finalized through the disposable runtime surface, while hadron
production is closed out-of-scope until OPH backend hardware exists.
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
RG_CONTRACT = P_DERIVATION_ROOT / "runtime" / "rg_matching_threshold_contract_current.json"
DIRECT_TOP_CONTRACT = PARTICLES_ROOT / "runs" / "calibration" / "direct_top_bridge_contract.json"
GAP_LEDGER = PARTICLES_ROOT / "runs" / "status" / "particle_derivation_gap_ledger.json"
BLIND_PROVENANCE = PARTICLES_ROOT / "runs" / "status" / "blind_prediction_provenance.json"
QUARK_CONTRACT = PARTICLES_ROOT / "runs" / "flavor" / "quark_lane_closure_contract.json"
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
    rg = _load_json(RG_CONTRACT)
    direct_top = _load_json(DIRECT_TOP_CONTRACT)
    gap_ledger = _load_json(GAP_LEDGER)
    blind_provenance = _load_json(BLIND_PROVENANCE)
    quark = _load_json(QUARK_CONTRACT)
    neutrino = _load_json(NEUTRINO_CONTRACT)
    hadron_spectral = _load_json(HADRON_SPECTRAL_CONTRACT)
    exact = _load_json(EXACT_NONHADRON)
    results_status = _load_json(RESULTS_STATUS)

    return {
        "artifact": "oph_particle_pipeline_closure_status",
        "generated_utc": _now_utc(),
        "purpose": "Single simplified closure gate for the current non-hadron particle pipeline.",
        "scope": {
            "current_pipeline_scope": "nonhadron_particles_plus_candidate_P_root_metadata",
            "hadrons_in_current_local_scope": False,
            "hadron_scope_reason": (
                "Production hadrons require a working OPH hardware backend such as GLORB/Echosahedron. "
                "Issues #153 and #157 are closed as out-of-scope/computationally blocked for the current "
                "pipeline; local surrogate code and Chrome workers are non-promoting."
            ),
            "chrome_workers_needed_now": False,
        },
        "current_surface": {
            "builder": "code/particles/compute_current_output_table.py",
            "default_runtime_root": "temp/particles_runtime",
            "source_repo": "reverse-engineering-reality/code",
            "simplification": (
                "The current prediction pipeline is one disposable runtime surface plus this closure "
                "manifest, with hadron promotion excluded by scope."
            ),
        },
        "artifacts": {
            "p_trunk": _artifact_status(P_TRUNK, p_trunk),
            "thomson_endpoint_contract": _artifact_status(THOMSON_CONTRACT, thomson),
            "rg_matching_threshold_contract": _artifact_status(RG_CONTRACT, rg),
            "direct_top_bridge_contract": _artifact_status(DIRECT_TOP_CONTRACT, direct_top),
            "gap_ledger": _artifact_status(GAP_LEDGER, gap_ledger),
            "blind_prediction_provenance": _artifact_status(BLIND_PROVENANCE, blind_provenance),
            "quark_lane_closure_contract": _artifact_status(QUARK_CONTRACT, quark),
            "neutrino_lane_closure_contract": _artifact_status(NEUTRINO_CONTRACT, neutrino),
            "hadron_spectral_measure_contract": _artifact_status(HADRON_SPECTRAL_CONTRACT, hadron_spectral),
        },
        "issue_gates": [
            {
                "issue": 223,
                "title": "Ward-projected Thomson endpoint",
                "state": "open_constructive_contract",
                "closable_now": False,
                "local_next_artifact": _rel(THOMSON_CONTRACT),
                "hadron_dependency_hardware_gated": True,
                "chrome_workers": "not_useful_until_source_endpoint_packet_exists",
            },
            {
                "issue": 224,
                "title": "Adopt certified derived P closure root",
                "state": "open_waiting_certified_root",
                "closable_now": False,
                "local_next_artifact": _rel(P_TRUNK),
                "promotion_allowed": bool((p_trunk or {}).get("consumer_policy", {}).get("may_feed_live_particle_predictions", False)),
                "chrome_workers": "not_useful_until_endpoint_and_interval_gates_close",
            },
            {
                "issue": 32,
                "title": "RG matching and threshold structure",
                "state": "open_constructive_contract",
                "closable_now": False,
                "local_next_artifact": _rel(RG_CONTRACT),
                "chrome_workers": "only_after_beta_threshold_packet_is_populated",
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
                    "Current environment lacks a working OPH hadron backend; reopen only when "
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
                    "The compact/paper hadron branch is outside the current computational scope "
                    "until a working OPH hadron backend emits the Ward-projected spectral measure."
                ),
                "chrome_workers": "do_not_use_for_backend_execution",
            },
            {
                "issue": 207,
                "title": "Direct-top codomain bridge",
                "state": "open_constructive_conversion_contract",
                "closable_now": False,
                "local_next_artifact": _rel(DIRECT_TOP_CONTRACT),
                "chrome_workers": "only_for_independent_audit_of_a_proposed_response_kernel",
            },
            {
                "issue": 234,
                "title": "Blind-prediction provenance and convention-sensitivity audits",
                "state": "open_provenance_ledger_emitted_sensitivity_open",
                "closable_now": False,
                "local_next_artifact": _rel(BLIND_PROVENANCE),
                "chrome_workers": "only_for_audit_after_convention_sensitivity_sweep_exists",
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
        ],
        "finalization_gates": {
            "nonhadron_prediction_surface_buildable": True,
            "hadrons_suppressed_by_default": bool(
                (results_status or {}).get("inputs", {}).get("hadron_profile", "suppressed") == "suppressed"
            ),
            "p_trunk_candidate_only": not bool(
                (p_trunk or {}).get("consumer_policy", {}).get("may_feed_live_particle_predictions", False)
            ),
            "obstruction_only_worker_result_allowed": False,
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
        f"- Current scope: `{status['scope']['current_pipeline_scope']}`",
        f"- Hadrons in current local scope: `{status['scope']['hadrons_in_current_local_scope']}`",
        f"- Chrome workers needed now: `{status['scope']['chrome_workers_needed_now']}`",
        f"- Hadron scope reason: {status['scope']['hadron_scope_reason']}",
        "",
        "## Issue Gates",
        "",
        "| Issue | State | Closable now | Local next artifact | Chrome policy |",
        "| --- | --- | --- | --- | --- |",
    ]
    for gate in status["issue_gates"]:
        lines.append(
            f"| #{gate['issue']} | `{gate['state']}` | `{gate['closable_now']}` | "
            f"`{gate['local_next_artifact']}` | {gate['chrome_workers']} |"
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
