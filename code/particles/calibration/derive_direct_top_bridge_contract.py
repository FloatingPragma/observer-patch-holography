#!/usr/bin/env python3
"""Emit the direct-top target-codomain obstruction certificate.

Chain role: keep the target-audit top extraction coordinate explicit while
closing transport to the auxiliary direct-top comparison as a current-corpus no-go.

The current target-audit top coordinate uses the PDG cross-section top entry.
The direct-top PDG entry is a
different extraction codomain, so closing the bridge requires an extraction
response map and uncertainty propagation rather than replacing the theorem row
with the auxiliary central value.  Since the current source corpus emits no such
response map, both coordinates remain comparison data.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
REFERENCE_JSON = ROOT / "particles" / "data" / "particle_reference_values.json"
QUARK_EXACT_READOUT_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_exact_readout.json"
PUBLIC_QUARK_THEOREM_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_public_exact_yukawa_end_to_end_theorem.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "calibration" / "direct_top_bridge_contract.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _symmetric_error(entry: dict[str, Any]) -> float:
    return 0.5 * (float(entry["error_minus_gev"]) + float(entry["error_plus_gev"]))


def _current_top_coordinate(quark_readout: dict[str, Any], public_theorem: dict[str, Any] | None) -> tuple[float, str]:
    if public_theorem is not None:
        public_outputs = public_theorem.get("public_exact_outputs", {})
        exact_values = public_outputs.get("exact_running_values_gev", {})
        if "t" in exact_values:
            return float(exact_values["t"]), "code/particles/runs/flavor/quark_public_exact_yukawa_end_to_end_theorem.json"
    return (
        float(quark_readout["predicted_singular_values_u"][2]),
        "code/particles/runs/flavor/quark_current_family_exact_readout.json",
    )


def build_payload(
    references: dict[str, Any],
    quark_readout: dict[str, Any],
    public_theorem: dict[str, Any] | None,
) -> dict[str, Any]:
    cross_section = references["top_quark"]
    direct = references["top_quark_direct_aux"]
    top_coordinate, source_artifact = _current_top_coordinate(quark_readout, public_theorem)

    cross_section_error = _symmetric_error(cross_section)
    direct_error = _symmetric_error(direct)
    direct_minus_coordinate = float(direct["value_gev"]) - top_coordinate
    combined_sigma = math.sqrt((direct_error * direct_error) + (cross_section_error * cross_section_error))

    return {
        "artifact": "oph_direct_top_bridge_contract",
        "generated_utc": _timestamp(),
        "github_issue": 207,
        "status": "hard_no_go_current_corpus_compare_only_direct_top_codomain",
        "promotion_allowed": False,
        "worker_result_policy": {
            "obstruction_only_result_allowed": True,
            "pro_workers_needed_now": False,
            "reason": "A new worker pass is not useful until a source-side extraction-response object exists.",
        },
        "closure_verdict": {
            "issue_207_acceptance_met_as_obstruction": True,
            "closure_kind": "hard_no_go_current_corpus",
            "closed_target_audit_surface": "Q007TP4_cross_section_top_coordinate",
            "blocked_surface": "Q007TP_auxiliary_direct_top_coordinate",
            "auxiliary_row_policy": "compare_only_not_promotable",
            "why_no_forward_bridge": (
                "The current OPH/D10/D11/quark corpus emits a top coordinate in the PDG cross-section "
                "codomain Q007TP4. It does not emit an extraction-response kernel, mass-scheme conversion, "
                "or uncertainty-propagation certificate from Q007TP4 into the auxiliary direct-measurement "
                "codomain Q007TP. Any such conversion would choose external experimental/MC response data."
            ),
        },
        "current_target_audit_coordinate": {
            "value_gev": top_coordinate,
            "source_artifact": source_artifact,
            "codomain": "PDG cross-section top mass entry",
            "pdg_summary_id": cross_section["source"]["summary_id"],
            "pdg_description": cross_section["description"],
            "pdg_value_gev": float(cross_section["value_gev"]),
            "pdg_symmetric_error_gev": cross_section_error,
        },
        "auxiliary_direct_top_coordinate": {
            "value_gev": float(direct["value_gev"]),
            "codomain": "PDG direct top mass entry",
            "pdg_summary_id": direct["source"]["summary_id"],
            "pdg_description": direct["description"],
            "pdg_symmetric_error_gev": direct_error,
        },
        "comparison_only_readout": {
            "direct_minus_current_coordinate_gev": direct_minus_coordinate,
            "combined_reference_sigma_gev": combined_sigma,
            "pull_in_combined_sigma": direct_minus_coordinate / combined_sigma,
            "within_combined_one_sigma": abs(direct_minus_coordinate) <= combined_sigma,
        },
        "constructive_objects": [
            {
                "id": "top_extraction_codomain_lock",
                "kind": "certificate_interface",
                "target_status": "cross_section_and_direct_top_codomain_metadata_explicit",
                "required_fields": [
                    "cross_section_summary_id",
                    "direct_top_summary_id",
                    "mass_scheme_labels",
                    "observable_extraction_definitions",
                    "row_policy_for_public_tables",
                ],
            },
            {
                "id": "cross_section_to_direct_top_response_kernel",
                "kind": "builder_interface",
                "target_status": "source_or_QFT_emitted_conversion_map",
                "required_fields": [
                    "input_top_coordinate",
                    "source_scheme",
                    "target_direct_top_scheme",
                    "response_kernel_formula",
                    "nonperturbative_or_extraction_systematic_terms",
                    "validity_interval",
                ],
            },
            {
                "id": "direct_top_uncertainty_propagation",
                "kind": "certificate_interface",
                "target_status": "uncertainty_budget_on_converted_coordinate",
                "required_fields": [
                    "coordinate_interval",
                    "scheme_conversion_budget",
                    "experimental_codomain_budget",
                    "correlation_policy",
                    "closure_decision",
                ],
            },
        ],
        "formal_nonidentifiability_witness": {
            "free_response_family": "m_direct(lambda) = m_Q007TP4 + lambda",
            "source_constraints_on_lambda": "none emitted by the current corpus",
            "lambda_matching_auxiliary_central_value_gev": direct_minus_coordinate,
            "why_not_a_theorem": (
                "Choosing lambda from the auxiliary Q007TP central value would use the direct-top row as an "
                "input. Choosing lambda=0 would relabel Q007TP4 as Q007TP. Both choices satisfy the current "
                "target-audit top data equally well because the extraction-response codomain is absent."
            ),
        },
        "closure_gate": {
            "closable_now": True,
            "reason": (
                "The current target-audit top coordinate and the auxiliary direct-top average are compatible "
                "as a comparison, but no source-side extraction-response map has been emitted; therefore "
                "the direct-top theorem bridge is closed as a no-go on the current corpus."
            ),
            "reopen_only_when": [
                "codomain metadata is explicit in the public status surfaces",
                "a response kernel maps the current top coordinate into the direct-top codomain without using the direct-top central value as an input",
                "the propagated interval overlaps or explains the auxiliary direct-top entry",
            ],
        },
        "local_next_steps": [
            "Keep both Q007TP and Q007TP4 as compare-only target coordinates.",
            "Do not launch more workers for this lane unless a source-side or QFT extraction-response map is introduced.",
            "Open a new issue only for a concrete response-kernel artifact, not for further border sharpening.",
        ],
        "forbidden_solver_inputs": [
            "Q007TP_direct_top_central_value_as_calibration_input",
            "post_hoc_shift_to_force_the_direct_top_average",
            "relabeling_Q007TP4_as_Q007TP_without_a_conversion_map",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the direct-top codomain closure certificate.")
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    parser.add_argument("--print-json", action="store_true")
    args = parser.parse_args()

    references = _load_json(REFERENCE_JSON)["entries"]
    public_theorem = _load_json(PUBLIC_QUARK_THEOREM_JSON) if PUBLIC_QUARK_THEOREM_JSON.exists() else None
    payload = build_payload(references, _load_json(QUARK_EXACT_READOUT_JSON), public_theorem)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    out_path.write_text(text, encoding="utf-8")
    if args.print_json:
        print(text, end="")
    else:
        print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
