#!/usr/bin/env python3
"""Build the final end-to-end particle prediction bundle."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
PARTICLES_ROOT = ROOT / "particles"
P_ROOT = ROOT / "P_derivation"
P_TRUNK = P_ROOT / "runtime" / "p_closure_trunk_current.json"
PIPELINE_STATUS = PARTICLES_ROOT / "runs" / "status" / "particle_pipeline_closure_status.json"
EXACT_NONHADRON = PARTICLES_ROOT / "exact_nonhadron_masses.json"
RESULTS_STATUS = PARTICLES_ROOT / "results_status.json"
DIRECT_TOP = PARTICLES_ROOT / "runs" / "calibration" / "direct_top_bridge_contract.json"
DEFAULT_JSON_OUT = PARTICLES_ROOT / "runs" / "status" / "final_end_to_end_predictions.json"
DEFAULT_MD_OUT = PARTICLES_ROOT / "FINAL_END_TO_END_PREDICTIONS.md"


PARTICLE_ORDER = [
    "photon",
    "gluon",
    "graviton",
    "w_boson",
    "z_boson",
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


def _display_status(status: str) -> str:
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


def build_payload() -> dict[str, Any]:
    p_trunk = _load_json(P_TRUNK)
    pipeline = _load_json(PIPELINE_STATUS)
    exact = _load_json(EXACT_NONHADRON)
    results = _load_json(RESULTS_STATUS)
    direct_top = _load_json(DIRECT_TOP)
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
        "scope": "nonhadron_particle_pipeline_with_hadrons_closed_out_of_scope",
        "claim_status": "final_nonhadron_predictions_without_full_hadron_or_certified_P_root_release",
        "source_surfaces": {
            "p_trunk": "code/P_derivation/runtime/p_closure_trunk_current.json",
            "thomson_endpoint_package": "code/P_derivation/runtime/thomson_endpoint_package_current.json",
            "pipeline_status": "code/particles/runs/status/particle_pipeline_closure_status.json",
            "exact_nonhadron": "code/particles/exact_nonhadron_masses.json",
            "results_status": "code/particles/results_status.json",
            "direct_top_bridge": "code/particles/runs/calibration/direct_top_bridge_contract.json",
        },
        "p_closure": {
            "P": p_trunk["fixed_point_candidate"]["P"],
            "alpha_inv": p_trunk["fixed_point_candidate"]["alpha_inv"],
            "claim_status": p_trunk["claim_status"],
            "may_feed_live_particle_predictions": p_trunk["consumer_policy"]["may_feed_live_particle_predictions"],
        },
        "runtime_inputs": results.get("inputs", {}),
        "finalization_gates": pipeline["finalization_gates"],
        "particle_five_issue_gates": particle_five_gates,
        "companion_open_branches": list(pipeline.get("companion_status_branches", [])),
        "predictions": predictions,
        "hadron_policy": {
            "predictions_emitted": False,
            "reason": (
                "Hadrons require a working OPH hadron backend on suitable hardware such as "
                "GLORB/Echosahedron. Issues #153/#157 are closed out-of-scope/computationally "
                "blocked, not solved, and no hadron predictions are emitted."
            ),
            "github_issues": [153, 157],
        },
        "direct_top_auxiliary_comparison": {
            "current_top_coordinate_gev": direct_top["current_theorem_coordinate"]["value_gev"],
            "current_top_codomain": direct_top["current_theorem_coordinate"]["pdg_summary_id"],
            "auxiliary_direct_top_gev": direct_top["auxiliary_direct_top_coordinate"]["value_gev"],
            "auxiliary_direct_top_codomain": direct_top["auxiliary_direct_top_coordinate"]["pdg_summary_id"],
            "direct_minus_current_coordinate_gev": direct_top["comparison_only_readout"]["direct_minus_current_coordinate_gev"],
            "pull_in_combined_sigma": direct_top["comparison_only_readout"]["pull_in_combined_sigma"],
            "bridge_status": direct_top["status"],
        },
    }


def render_markdown(payload: dict[str, Any]) -> str:
    lines = [
        "# Final End-to-End Particle Predictions",
        "",
        f"Generated: `{payload['generated_utc']}`",
        "",
        f"Scope: `{payload['scope']}`",
        f"Claim status: `{payload['claim_status']}`",
        "",
        "## P Closure",
        "",
        f"- Candidate `P`: `{payload['p_closure']['P']}`",
        f"- Candidate `alpha^-1`: `{payload['p_closure']['alpha_inv']}`",
        f"- Claim status: `{payload['p_closure']['claim_status']}`",
        f"- May feed promoted particle predictions: `{payload['p_closure']['may_feed_live_particle_predictions']}`",
        "",
        "## Particle-Five Gates",
        "",
        "| Issue | State | Closable | Local artifact | Worker policy |",
        "| --- | --- | --- | --- | --- |",
    ]
    for gate in payload["particle_five_issue_gates"]:
        lines.append(
            f"| #{gate['issue']} | `{_display_status(gate['state'])}` | `{gate['closable_now']}` | "
            f"`{gate['local_next_artifact']}` | {gate['chrome_workers']} |"
        )
    companion_open_branches = payload.get("companion_open_branches") or []
    if companion_open_branches:
        lines.extend(
            [
                "",
                "## Companion Open Branches",
                "",
                "| Topic | State | Current boundary | Next action |",
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
            "| Particle | Prediction | Status | Scope | Promotable |",
            "| --- | ---: | --- | --- | --- |",
        ]
    )
    for entry in payload["predictions"]:
        lines.append(
            f"| `{entry['particle_id']}` | `{entry['value']} {entry['unit']}` | "
            f"`{entry['exact_kind']}` | `{entry['scope']}` | `{entry['promotable']}` |"
        )
    direct = payload["direct_top_auxiliary_comparison"]
    lines.extend(
        [
            "",
            "## Direct-Top Auxiliary Comparison",
            "",
            f"- Top theorem coordinate: `{direct['current_top_coordinate_gev']} GeV` on `{direct['current_top_codomain']}`",
            f"- Auxiliary direct-top coordinate: `{direct['auxiliary_direct_top_gev']} GeV` on `{direct['auxiliary_direct_top_codomain']}`",
            f"- Difference: `{direct['direct_minus_current_coordinate_gev']} GeV`",
            f"- Pull: `{direct['pull_in_combined_sigma']}` combined sigma",
            f"- Bridge status: `{_display_status(direct['bridge_status'])}`",
            "",
            "## Hadrons",
            "",
            f"- Predictions emitted: `{payload['hadron_policy']['predictions_emitted']}`",
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
