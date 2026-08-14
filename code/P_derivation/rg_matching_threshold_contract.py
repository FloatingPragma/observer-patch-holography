#!/usr/bin/env python3
"""Emit the scientific contract for RG matching and threshold transport.

The artifact records the exact source-bound precursor, the scientific objects
that are still absent, the corpus-limited non-identifiability result, and the
promotion boundary. Project-control and scheduling metadata belong outside
this repository.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


DEFAULT_OUT = Path(__file__).resolve().parent / "runtime" / "rg_matching_threshold_contract_current.json"
FRONTIER = (
    Path(__file__).resolve().parent
    / "source_rg_frontier"
    / "outputs"
    / "rg_representation_frontier.json"
)


def build_contract() -> dict[str, Any]:
    frontier_raw = FRONTIER.read_bytes()
    frontier = json.loads(frontier_raw)
    if frontier.get("status") != "PARTIAL_EXACT_REPRESENTATION_INDICES__SOURCE_MATCHING_OPEN":
        raise ValueError("the source RG frontier is missing or has an unsafe status")
    frontier_ref = {
        "path": FRONTIER.relative_to(Path(__file__).resolve().parent.parent).as_posix(),
        "bytes": len(frontier_raw),
        "byte_sha256": hashlib.sha256(frontier_raw).hexdigest(),
        "subject_digest": frontier["subject_digest"],
        "status": frontier["status"],
    }
    return {
        "artifact": "oph_rg_matching_threshold_contract",
        "status": "open_source_rg_frontier_partial",
        "promotion_allowed": False,
        "source_frontier": frontier_ref,
        "constructive_objects": [
            {
                "id": "representation_index_frontier",
                "kind": "source_bound_exact_precursor",
                "current_status": "complete_at_finite_representation_scope",
                "source_artifact": frontier_ref["path"],
                "boundary": (
                    "The per-copy quadratic indices are exact. The standard "
                    "four-dimensional one-loop beta functional is imported, "
                    "family/scalar multiplicities are not source-selected, "
                    "and zero gauge indices do not prove full W/Z decoupling."
                ),
            },
            {
                "id": "scheme_lock",
                "kind": "certificate_interface",
                "current_status": "not_emitted",
                "target_status": "one_scheme_used_by_unification_anchor_endpoint_and_mass_readouts",
                "required_fields": [
                    "renormalization_scheme",
                    "normalization_of_U1",
                    "input_scale_definitions",
                    "conversion_maps",
                    "surfaces_using_the_scheme",
                ],
            },
            {
                "id": "threshold_map",
                "kind": "builder_interface",
                "current_status": "not_emitted",
                "target_status": "source_emitted_threshold_placements_and_decoupling_maps",
                "required_fields": [
                    "particle_thresholds",
                    "superpartner_or_effective_thresholds",
                    "matching_scales",
                    "threshold_uncertainty_intervals",
                    "status_per_threshold",
                ],
            },
            {
                "id": "beta_provenance_table",
                "kind": "audit_table",
                "current_status": "partial_parametric_gauge_one_loop_only",
                "target_status": "complete source-derived vector beta provenance at every interval",
                "required_fields": [
                    "gauge_group",
                    "matter_content",
                    "loop_order",
                    "coefficient",
                    "status",
                    "source_artifact_or_reference",
                ],
            },
            {
                "id": "matching_interval_composition_certificate",
                "kind": "certificate_interface",
                "current_status": "not_emitted",
                "target_status": "interval_bound_for_composed_running_map",
                "required_fields": [
                    "input_intervals",
                    "composed_map",
                    "roundoff_budget",
                    "matching_budget",
                    "threshold_budget",
                    "image_interval",
                ],
            },
        ],
        "scientific_boundary": {
            "established": (
                "The source frontier determines exact per-copy representation indices."
            ),
            "not_established": (
                "A source-selected physical action, coupled-field census, threshold "
                "ordering, finite matching transport, or certified vector running map."
            ),
            "corpus_limited_no_go": (
                "The exact representation indices alone do not determine a unique "
                "physical RG and threshold transport packet."
            ),
        },
        "promotion_boundary": {
            "promotion_allowed": False,
            "reason": (
                "The finite source receipts determine exact per-copy representation indices, "
                "but no target-clean OPH source emits the physical action, family and scalar "
                "attachments, complete W/Z-coupled census modulo proved zero-vertex "
                "decoupling, ordered thresholds, finite decoupling and scheme maps, "
                "Jacobians, masks, or certified vector remainders."
            ),
            "required_scientific_objects": [
                "one target-clean OPH source emits the complete W/Z-coupled active and heavy field census modulo proved zero-vertex decoupling",
                "the same source emits every ordered interval and threshold location",
                "finite decoupling and scheme maps, Jacobians, and term masks are executable and independently replayed",
                "the composed vector running map carries certified deterministic remainders",
                "no measured electroweak target or external Standard Model packet enters source construction",
            ],
        },
        "forbidden_promotions": [
            "treating_the_conditional_Ng3_NH1_evaluation_as_source_selected",
            "silently_treating_declared_MSSM_running_as_OPH_derived",
            "using_threshold_choices_as_hidden_fit_parameters",
            "reusing_an_external_validation_packet_as_an_OPH_source",
            "promoting_p_closure_root_without_interval_composition_certificate",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Emit the RG matching and threshold constructive contract.")
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    parser.add_argument("--print-json", action="store_true")
    args = parser.parse_args()

    payload = build_contract()
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text, encoding="utf-8")
    if args.print_json:
        print(text, end="")
    else:
        print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
