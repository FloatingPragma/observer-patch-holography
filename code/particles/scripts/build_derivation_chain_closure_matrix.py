#!/usr/bin/env python3
"""Build the particle derivation-chain closure matrix."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
PARTICLES_ROOT = ROOT / "particles"
FINAL_PREDICTIONS = PARTICLES_ROOT / "runs" / "status" / "final_end_to_end_predictions.json"
PIPELINE_STATUS = PARTICLES_ROOT / "runs" / "status" / "particle_pipeline_closure_status.json"
BLIND_PROVENANCE = PARTICLES_ROOT / "runs" / "status" / "blind_prediction_provenance.json"
CHARGED_NONCLOSURE = PARTICLES_ROOT / "runs" / "leptons" / "charged_end_to_end_impossibility_theorem.json"
QUARK_GLOBAL_OBSTRUCTION = (
    PARTICLES_ROOT / "runs" / "flavor" / "quark_class_uniform_public_frame_descent_obstruction.json"
)
DIRECT_TOP_CONTRACT = PARTICLES_ROOT / "runs" / "calibration" / "direct_top_bridge_contract.json"
DEFAULT_JSON_OUT = PARTICLES_ROOT / "runs" / "status" / "derivation_chain_closure_matrix.json"
DEFAULT_MD_OUT = PARTICLES_ROOT / "DERIVATION_CHAIN_CLOSURE_MATRIX.md"


def _now_utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _prediction_map(payload: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {entry["particle_id"]: entry for entry in payload["predictions"]}


def _issue_map(payload: dict[str, Any]) -> dict[int, dict[str, Any]]:
    return {int(gate["issue"]): gate for gate in payload["issue_gates"]}


def _display_status(status: str) -> str:
    return status.replace("current_corpus", "corpus_limited").replace("not_live", "not_certified")


def build_payload() -> dict[str, Any]:
    final_predictions = _load_json(FINAL_PREDICTIONS)
    pipeline = _load_json(PIPELINE_STATUS)
    provenance = _load_json(BLIND_PROVENANCE)
    charged_nonclosure = _load_json(CHARGED_NONCLOSURE)
    quark_global = _load_json(QUARK_GLOBAL_OBSTRUCTION)
    direct_top = _load_json(DIRECT_TOP_CONTRACT)
    predictions = _prediction_map(final_predictions)
    gates = _issue_map(pipeline)

    rows = [
        {
            "chain": "p_closure_root",
            "status": "candidate_not_live_root",
            "claim_level": final_predictions["p_closure"]["claim_status"],
            "outputs": {
                "P": final_predictions["p_closure"]["P"],
                "alpha_inv": final_predictions["p_closure"]["alpha_inv"],
            },
            "promotable": False,
            "open_gates": [32, 235, 224],
            "next_artifact": "code/P_derivation/runtime/rg_matching_threshold_contract_current.json + code/P_derivation/runtime/thomson_endpoint_package_current.json",
        },
        {
            "chain": "structural_massless_bosons",
            "status": "closed_structural_zero",
            "claim_level": "structural",
            "outputs": {key: predictions[key]["value"] for key in ("photon", "gluon", "graviton")},
            "promotable": True,
            "open_gates": [],
            "next_artifact": None,
        },
        {
            "chain": "electroweak_wz",
            "status": "compare_only_reproduction_not_prediction_theorem",
            "claim_level": predictions["w_boson"]["exact_kind"],
            "outputs": {
                "w_boson": predictions["w_boson"]["value"],
                "z_boson": predictions["z_boson"]["value"],
            },
            "promotable": False,
            "open_gates": [32, 235, 224],
            "next_artifact": "code/P_derivation/runtime/rg_matching_threshold_contract_current.json",
        },
        {
            "chain": "higgs_top_declared_surface",
            "status": "closed_on_declared_d10_d11_surface_direct_top_no_go",
            "claim_level": predictions["higgs"]["exact_kind"],
            "outputs": {
                "higgs": predictions["higgs"]["value"],
                "top_companion": predictions["top_quark"]["value"],
            },
            "promotable": True,
            "open_gates": [],
            "closed_issue_refs": [207],
            "next_artifact": "code/particles/runs/calibration/direct_top_bridge_contract.json",
            "codomain_obstruction": direct_top.get("status"),
        },
        {
            "chain": "charged_leptons",
            "status": "closed_current_corpus_charged_end_to_end_no_go",
            "claim_level": predictions["electron"]["exact_kind"],
            "outputs": {
                "electron": predictions["electron"]["value"],
                "muon": predictions["muon"]["value"],
                "tau": predictions["tau"]["value"],
            },
            "promotable": False,
            "open_gates": [],
            "closed_issue_refs": [201],
            "next_artifact": "code/particles/runs/leptons/charged_end_to_end_impossibility_theorem.json",
            "nonclosure_theorem": charged_nonclosure.get("artifact"),
        },
        {
            "chain": "selected_class_quarks",
            "status": "closed_selected_public_class_global_classification_no_go",
            "claim_level": predictions["top_quark"]["exact_kind"],
            "outputs": {
                "up": predictions["up_quark"]["value"],
                "down": predictions["down_quark"]["value"],
                "strange": predictions["strange_quark"]["value"],
                "charm": predictions["charm_quark"]["value"],
                "bottom": predictions["bottom_quark"]["value"],
                "top": predictions["top_quark"]["value"],
            },
            "promotable": True,
            "open_gates": [],
            "closed_issue_refs": [199, 207, 212],
            "next_artifact": "code/particles/runs/flavor/quark_class_uniform_public_frame_descent_obstruction.json",
            "global_classification_obstruction": quark_global.get("proof_status"),
        },
        {
            "chain": "neutrino_absolute_attachment",
            "status": "closed_weighted_cycle_absolute_attachment_with_comparison_tension_visible",
            "claim_level": predictions["electron_neutrino"]["exact_kind"],
            "outputs": {
                "electron_neutrino": predictions["electron_neutrino"]["value"],
                "muon_neutrino": predictions["muon_neutrino"]["value"],
                "tau_neutrino": predictions["tau_neutrino"]["value"],
            },
            "unit": "eV",
            "promotable": True,
            "open_gates": [],
            "next_artifact": None,
        },
        {
            "chain": "hadrons",
            "status": "closed_out_of_scope_computationally_blocked",
            "claim_level": "no_local_prediction_emitted",
            "outputs": {},
            "promotable": False,
            "open_gates": [],
            "closed_issue_refs": [153, 157],
            "closure_reason": (
                "The local environment has no working OPH hadron backend; local surrogate code "
                "and Chrome workers cannot promote hadron predictions."
            ),
            "next_artifact": (
                "none in local scope; reopen only when a GLORB/Echosahedron-class OPH backend "
                "emits production hadron output and systematics"
            ),
        },
    ]
    closed_or_scoped = [
        row["chain"]
        for row in rows
        if row["status"]
        in {
            "closed_structural_zero",
            "closed_on_declared_d10_d11_surface_direct_top_no_go",
            "closed_current_corpus_charged_end_to_end_no_go",
            "closed_selected_public_class_global_classification_no_go",
            "closed_weighted_cycle_absolute_attachment_with_comparison_tension_visible",
        }
    ]
    remaining_nonclosed = [
        row["chain"]
        for row in rows
        if row["chain"] not in closed_or_scoped
        and row["status"] != "closed_out_of_scope_computationally_blocked"
    ]

    return {
        "artifact": "oph_particle_derivation_chain_closure_matrix",
        "generated_utc": _now_utc(),
        "status": "executable_nonhadron_chain_matrix_emitted",
        "closure_summary": {
            "all_derivation_chains_claimed_closed": False,
            "closed_or_scoped_chains": closed_or_scoped,
            "remaining_nonclosed_chains": remaining_nonclosed,
            "closed_out_of_scope_chains": ["hadrons"],
            "hardware_gated_chains": ["hadrons"],
            "policy": (
                "Do not promote candidate, compare-only, witness-only, corpus-limited no-go, or "
                "hardware-gated chains as closed theorem predictions."
            ),
        },
        "source_artifacts": {
            "final_predictions": "code/particles/runs/status/final_end_to_end_predictions.json",
            "pipeline_status": "code/particles/runs/status/particle_pipeline_closure_status.json",
            "blind_provenance": "code/particles/runs/status/blind_prediction_provenance.json",
        },
        "worker_policy": {
            "chrome_pro_workers_needed_now": False,
            "reason": (
                "Hadron issues #153/#157 are closed out-of-scope pending OPH hadron hardware; "
                "the remaining in-scope open P/electroweak chains need a source-emitted residual map and interval certificate before worker audit is meaningful. "
                "The charged absolute-anchor, quark global-classification, and direct-top auxiliary-codomain lanes "
                "have corpus-limited no-go boundaries."
            ),
        },
        "particle_five_gates": {
            str(issue): gates[issue]
            for issue in (32, 153, 199, 201, 207, 223, 224, 235)
            if issue in gates
        },
        "provenance_status": provenance["status"],
        "rows": rows,
    }


def render_markdown(payload: dict[str, Any]) -> str:
    lines = [
        "# Particle Derivation Chain Closure Matrix",
        "",
        f"Generated: `{payload['generated_utc']}`",
        "",
        f"Status: `{payload['status']}`",
        f"All derivation chains claimed closed: `{payload['closure_summary']['all_derivation_chains_claimed_closed']}`",
        f"Remaining nonclosed chains: `{', '.join(payload['closure_summary']['remaining_nonclosed_chains'])}`",
        f"Chrome Pro workers needed: `{payload['worker_policy']['chrome_pro_workers_needed_now']}`",
        f"Reason: {payload['worker_policy']['reason']}",
        "",
        "| Chain | Status | Promotable | Open gates | Outputs | Next artifact |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for row in payload["rows"]:
        outputs = ", ".join(f"`{key}={value}`" for key, value in row["outputs"].items()) or "n/a"
        gates = ", ".join(f"#{issue}" for issue in row["open_gates"]) or "none"
        next_artifact = row["next_artifact"] or "none"
        lines.append(
            f"| `{row['chain']}` | `{_display_status(row['status'])}` | `{row['promotable']}` | {gates} | {outputs} | {next_artifact} |"
        )
    return "\n".join(lines)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build particle derivation-chain closure matrix.")
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
