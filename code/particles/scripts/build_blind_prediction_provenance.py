#!/usr/bin/env python3
"""Build the quantitative particle provenance and blind-prediction audit."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
PARTICLES_ROOT = ROOT / "particles"
P_ROOT = ROOT / "P_derivation"
EXACT_NONHADRON = PARTICLES_ROOT / "exact_nonhadron_masses.json"
RESULTS_STATUS = PARTICLES_ROOT / "results_status.json"
PIPELINE_STATUS = PARTICLES_ROOT / "runs" / "status" / "particle_pipeline_closure_status.json"
RG_CONTRACT = P_ROOT / "runtime" / "rg_matching_threshold_contract_current.json"
THOMSON_CONTRACT = P_ROOT / "runtime" / "thomson_endpoint_contract_current.json"
THOMSON_PACKAGE = P_ROOT / "runtime" / "thomson_endpoint_package_current.json"
DEFAULT_JSON_OUT = PARTICLES_ROOT / "runs" / "status" / "blind_prediction_provenance.json"
DEFAULT_MD_OUT = PARTICLES_ROOT / "BLIND_PREDICTION_PROVENANCE.md"


def _now_utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _repo_ref(path: Path) -> str:
    return f"code/{path.relative_to(ROOT)}"


def _classify_entry(entry: dict[str, Any]) -> dict[str, Any]:
    exact_kind = entry.get("exact_kind", "")
    scope = entry.get("scope")
    promotable = bool(entry.get("promotable"))
    particle_id = entry["particle_id"]

    if exact_kind == "structural_zero":
        row_class = "structural_blind_zero"
        target_use = "no_mass_target_used"
        blind_status = "blind_structural"
        convention_sensitivity = "none_for_mass_zero"
    elif "frozen_target" in exact_kind:
        row_class = "compare_only_reproduction"
        target_use = "target_used_as_frozen_reference"
        blind_status = "not_blind"
        convention_sensitivity = "inherits_declared_electroweak_surface"
    elif "target_anchored" in exact_kind:
        row_class = "target_anchored_witness"
        target_use = "target_values_used_to_anchor_family_witness"
        blind_status = "not_blind"
        convention_sensitivity = "not_promotable_until_source_attachment_closes"
    elif "selected_class" in exact_kind:
        row_class = "selected_class_exact_theorem"
        target_use = "reference_codomain_matched_after_selected_class_theorem"
        blind_status = "selected_class_not_global_blind"
        convention_sensitivity = "quark_scheme_and_frame_class_scope_must_remain_visible"
    elif "source_only" in exact_kind:
        row_class = "source_only_declared_surface_theorem"
        target_use = "no_direct_row_target_input_recorded_in_the_forward_split_readout"
        blind_status = "conditionally_blind_on_declared_surface"
        convention_sensitivity = "depends_on_declared_D10_D11_running_matching_threshold_surface"
    elif "weighted_cycle" in exact_kind:
        row_class = "theorem_branch_no_direct_mass_target"
        target_use = "no_absolute_mass_target_input"
        blind_status = "blind_absolute_mass_branch"
        convention_sensitivity = "pmns_comparison_tension_reported_separately"
    else:
        row_class = "unclassified"
        target_use = "requires_manual_audit"
        blind_status = "unknown"
        convention_sensitivity = "requires_manual_audit"

    value = entry.get("mass_gev")
    unit = "GeV"
    if value is None:
        value = entry.get("mass_eV")
        unit = "eV"

    return {
        "particle_id": particle_id,
        "value": value,
        "unit": unit,
        "exact_kind": exact_kind,
        "scope": scope,
        "promotable": promotable,
        "row_class": row_class,
        "target_use": target_use,
        "blind_status": blind_status,
        "source_artifact": entry.get("source_artifact"),
        "supporting_scope_closure_artifact": entry.get("supporting_scope_closure_artifact"),
        "convention_sensitivity": convention_sensitivity,
    }


def build_payload() -> dict[str, Any]:
    exact = _load_json(EXACT_NONHADRON)
    results = _load_json(RESULTS_STATUS)
    pipeline = _load_json(PIPELINE_STATUS)
    rg = _load_json(RG_CONTRACT)
    thomson = _load_json(THOMSON_CONTRACT)
    thomson_package = _load_json(THOMSON_PACKAGE)
    rows = [_classify_entry(entry) for entry in exact.get("entries", [])]

    return {
        "artifact": "oph_blind_prediction_provenance_audit",
        "generated_utc": _now_utc(),
        "github_issue": 234,
        "scope": "public_quantitative_particle_rows",
        "status": "provenance_ledger_emitted_convention_sensitivity_contract_open",
        "promotion_allowed": False,
        "source_surfaces": {
            "exact_nonhadron": _repo_ref(EXACT_NONHADRON),
            "results_status": _repo_ref(RESULTS_STATUS),
            "pipeline_closure_status": _repo_ref(PIPELINE_STATUS),
            "rg_matching_threshold_contract": _repo_ref(RG_CONTRACT),
            "thomson_endpoint_contract": _repo_ref(THOMSON_CONTRACT),
            "thomson_endpoint_package": _repo_ref(THOMSON_PACKAGE),
        },
        "pipeline_inputs": results.get("inputs", {}),
        "finalization_gates": pipeline.get("finalization_gates", {}),
        "row_counts": {
            "total": len(rows),
            "promotable": sum(1 for row in rows if row["promotable"]),
            "not_promotable": sum(1 for row in rows if not row["promotable"]),
            "blind_or_conditionally_blind": sum(
                1 for row in rows if row["blind_status"] in {"blind_structural", "conditionally_blind_on_declared_surface", "blind_absolute_mass_branch"}
            ),
        },
        "rows": rows,
        "convention_sensitivity": {
            "status": "contract_open_not_quantified",
            "rg_contract_status": rg.get("status"),
            "required_objects": [item["id"] for item in rg.get("constructive_objects", [])],
            "endpoint_contract_status": thomson.get("status"),
            "endpoint_package_status": thomson_package.get("claim_status"),
            "next_artifact": "populate beta_provenance_table, threshold_map, scheme_lock, and interval composition certificates",
        },
        "preregistered_blind_workflows": [
            {
                "id": "new_quantity_pre_reference_lock",
                "status": "protocol_emitted_unexercised",
                "rule": (
                    "For any quantitative row outside construction inputs, freeze the source artifacts, "
                    "hash the runtime bundle, record allowed conventions, then fetch or reveal the external reference."
                ),
                "required_evidence": [
                    "source_artifact_hashes",
                    "forbidden_target_inputs",
                    "convention_set",
                    "pre_reference_runtime_output",
                    "post_reference_comparison_only_delta",
                ],
            },
            {
                "id": "convention_sensitivity_sweep",
                "status": "blocked_on_rg_matching_threshold_contract",
                "rule": (
                    "Vary only declared scheme, matching, and threshold choices inside certified intervals; "
                    "report induced intervals for every public quantitative row."
                ),
                "required_evidence": [
                    "scheme_lock",
                    "threshold_map",
                    "matching_interval_composition_certificate",
                    "rowwise_sensitivity_intervals",
                ],
            },
        ],
        "closure_gate": {
            "closable_now": False,
            "reason": (
                "The row provenance ledger and blind workflow protocol exist, but convention sensitivity "
                "is still a contract because #32 has not supplied the scheme/threshold interval certificate."
            ),
            "close_issue_when": [
                "row provenance is regenerated by the live runtime",
                "at least one pre-reference workflow is exercised on a new or withheld output",
                "convention sensitivity intervals are quantified for public rows",
            ],
        },
    }


def render_markdown(payload: dict[str, Any]) -> str:
    lines = [
        "# Blind Prediction Provenance",
        "",
        f"Generated: `{payload['generated_utc']}`",
        "",
        "This ledger records target-use and convention-sensitivity status for the public quantitative particle rows.",
        "",
        "## Closure Gate",
        "",
        f"- Status: `{payload['status']}`",
        f"- Closable: `{payload['closure_gate']['closable_now']}`",
        f"- Reason: {payload['closure_gate']['reason']}",
        "",
        "## Rows",
        "",
        "| Particle | Value | Class | Blind status | Target use | Promotable | Convention sensitivity |",
        "| --- | ---: | --- | --- | --- | --- | --- |",
    ]
    for row in payload["rows"]:
        value = "n/a" if row["value"] is None else f"{row['value']} {row['unit']}"
        lines.append(
            f"| `{row['particle_id']}` | `{value}` | `{row['row_class']}` | `{row['blind_status']}` | "
            f"{row['target_use']} | `{row['promotable']}` | {row['convention_sensitivity']} |"
        )
    lines.extend(
        [
            "",
            "## Preregistered Workflows",
            "",
        ]
    )
    for workflow in payload["preregistered_blind_workflows"]:
        evidence = ", ".join(f"`{item}`" for item in workflow["required_evidence"])
        lines.append(f"- `{workflow['id']}`: `{workflow['status']}`. {workflow['rule']} Required evidence: {evidence}.")
    lines.extend(
        [
            "",
            "## Convention Sensitivity",
            "",
            f"- Status: `{payload['convention_sensitivity']['status']}`",
            f"- RG contract status: `{payload['convention_sensitivity']['rg_contract_status']}`",
            f"- Endpoint contract status: `{payload['convention_sensitivity']['endpoint_contract_status']}`",
            f"- Endpoint package status: `{payload['convention_sensitivity']['endpoint_package_status']}`",
            f"- Next artifact: {payload['convention_sensitivity']['next_artifact']}",
        ]
    )
    return "\n".join(lines)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build the blind-prediction provenance audit.")
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
