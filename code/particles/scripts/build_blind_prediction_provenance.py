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
CARRIER_ACCEPTANCE = PARTICLES_ROOT / "runs" / "status" / "carrier_mode_acceptance.json"
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

    if "frozen_target" in exact_kind:
        row_class = "compare_only_reproduction"
        target_use = "target_used_as_frozen_reference"
        blind_status = "not_blind"
        convention_sensitivity = "inherits_declared_electroweak_surface"
    elif exact_kind == "conditional_declared_surface_higgs_top_candidate":
        row_class = "conditional_declared_surface_candidate"
        target_use = "candidate_upstream_d10_repair_not_source_promoted"
        blind_status = "conditionally_blind_on_declared_surface"
        convention_sensitivity = "depends_on_declared_D10_D11_running_matching_threshold_surface"
    elif exact_kind == "selected_class_target_anchored_exact_witness":
        row_class = "selected_class_target_anchored_witness"
        target_use = "target_derived_sigma_datum_used_for_selected_class_exact_witness"
        blind_status = "selected_class_target_anchored_not_blind"
        convention_sensitivity = "quark_scheme_and_frame_class_scope_must_remain_visible"
    elif "target_anchored" in exact_kind:
        row_class = "target_anchored_witness"
        target_use = "target_values_used_to_anchor_current_family_witness"
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
        row_class = "rejected_target_informed_template_candidate"
        target_use = "target_ranked_selector_development_and_compare_only_absolute_attachment"
        blind_status = "withheld_not_blind_rejected_candidate"
        convention_sensitivity = "all_row_column_and_cycle_orientations_audited_no_nufit61_rescue"
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


def _classify_withheld_entry(entry: dict[str, Any]) -> dict[str, Any]:
    exact_kind = entry.get("exact_kind", "")
    if "target_anchored" in exact_kind:
        target_use = "target_values_or_target_derived_datum_used"
        blind_status = "withheld_not_blind"
    elif "compare_only" in exact_kind:
        target_use = "compare_only_reference_or_absolute_attachment_used"
        blind_status = "withheld_compare_only"
    elif "target_informed" in exact_kind or "rejected" in exact_kind:
        target_use = "target_ranked_selector_development_and_correlated_profile_rejection"
        blind_status = "withheld_not_blind_rejected_candidate"
    else:
        target_use = "withheld_by_public_output_policy"
        blind_status = "withheld"
    return {
        "particle_id": entry["particle_id"],
        "exact_kind": exact_kind,
        "scope": entry.get("scope"),
        "promotable": bool(entry.get("promotable", False)),
        "target_use": target_use,
        "blind_status": blind_status,
        "source_artifact": entry.get("source_artifact"),
        "reason": entry.get("reason", "not_public_prediction_output"),
    }


def _classify_carrier_mode(entry: dict[str, Any]) -> dict[str, Any]:
    return {
        "carrier_id": entry["carrier_id"],
        "row_class": "conditional_classical_carrier_mode_not_particle_mass_prediction",
        "hard_quadratic_mass_parameter_squared": entry["hard_quadratic_mass_parameter_squared"],
        "target_use": "no_particle_mass_target_used",
        "blind_status": "not_a_quantum_particle_prediction",
        "classical_carrier_gate": entry["classical_carrier_gate"]["status"],
        "quantum_particle_gate": entry["quantum_particle_gate"]["status"],
        "particle_promotion_allowed": bool(entry["particle_promotion_allowed"]),
        "branch": entry["branch"],
    }


def build_payload() -> dict[str, Any]:
    exact = _load_json(EXACT_NONHADRON)
    carrier_acceptance = _load_json(CARRIER_ACCEPTANCE)
    results = _load_json(RESULTS_STATUS)
    pipeline = _load_json(PIPELINE_STATUS)
    rg = _load_json(RG_CONTRACT)
    thomson = _load_json(THOMSON_CONTRACT)
    thomson_package = _load_json(THOMSON_PACKAGE)
    rows = [_classify_entry(entry) for entry in exact.get("entries", [])]
    withheld_rows = [_classify_withheld_entry(entry) for entry in exact.get("withheld_entries", [])]
    carrier_mode_rows = [
        _classify_carrier_mode(entry) for entry in carrier_acceptance.get("carriers", [])
    ]

    return {
        "artifact": "oph_blind_prediction_provenance_audit",
        "generated_utc": _now_utc(),
        "github_issue": 234,
        "scope": "public_quantitative_particle_rows",
        "status": "closed_provenance_ledger_and_declared_sensitivity_taxonomy",
        "promotion_allowed": False,
        "source_surfaces": {
            "exact_nonhadron": _repo_ref(EXACT_NONHADRON),
            "carrier_mode_acceptance": _repo_ref(CARRIER_ACCEPTANCE),
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
            "withheld_non_prediction": len(withheld_rows),
            "promotable": sum(1 for row in rows if row["promotable"]),
            "not_promotable": sum(1 for row in rows if not row["promotable"]),
            "separated_classical_carrier_modes": len(carrier_mode_rows),
            "blind_or_conditionally_blind": sum(
                1
                for row in rows
                if row["blind_status"]
                in {"conditionally_blind_on_declared_surface", "blind_absolute_mass_branch"}
            ),
        },
        "rows": rows,
        "withheld_rows": withheld_rows,
        "carrier_mode_rows": carrier_mode_rows,
        "convention_sensitivity": {
            "status": "declared_taxonomy_emitted_numeric_sweep_stage_gated",
            "rg_contract_status": rg.get("status"),
            "required_objects": [item["id"] for item in rg.get("constructive_objects", [])],
            "endpoint_contract_status": thomson.get("status"),
            "endpoint_package_status": thomson_package.get("claim_status"),
            "numeric_sweep_performed": False,
            "next_artifact": "interval composition certificates after the populated source spectral measure payload exists",
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
                "status": "declared_taxonomy_emitted_numeric_sweep_stage_gated",
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
            "closable_now": True,
            "closed_as": "provenance_ledger_and_declared_sensitivity_taxonomy",
            "reason": (
                "The row provenance ledger, blind workflow protocol, and declared convention-sensitivity "
                "taxonomy are emitted. Numeric sensitivity intervals remain tied to the populated source "
                "spectral measure payload and interval certificate."
            ),
            "close_issue_when": [
                "row provenance ledger is emitted",
                "declared sensitivity taxonomy is emitted",
                "numeric sweeps remain stage-gated by the source spectral measure payload",
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
    withheld_rows = payload.get("withheld_rows") or []
    if withheld_rows:
        lines.extend(
            [
                "",
                "## Withheld Non-Prediction Rows",
                "",
                "These rows have audit artifacts but no public prediction value in the output tables.",
                "",
                "| Particle | Claim label | Blind status | Target use | Reason |",
                "| --- | --- | --- | --- | --- |",
            ]
        )
        for row in withheld_rows:
            lines.append(
                f"| `{row['particle_id']}` | `{row['exact_kind']}` | `{row['blind_status']}` | "
                f"{row['target_use']} | {row['reason']} |"
            )
    carrier_rows = payload.get("carrier_mode_rows") or []
    if carrier_rows:
        lines.extend(
            [
                "",
                "## Separated Classical Carrier Modes",
                "",
                "These zero hard quadratic parameters are branch-conditional mode statements, not public quantum-particle mass predictions.",
                "",
                "| Carrier | Hard parameter squared | Classical gate | Quantum gate | Particle promotion |",
                "| --- | ---: | --- | --- | --- |",
            ]
        )
        for row in carrier_rows:
            lines.append(
                f"| `{row['carrier_id']}` | `{row['hard_quadratic_mass_parameter_squared']}` | "
                f"`{row['classical_carrier_gate']}` | `{row['quantum_particle_gate']}` | "
                f"`{row['particle_promotion_allowed']}` |"
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
