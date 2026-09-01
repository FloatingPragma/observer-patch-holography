#!/usr/bin/env python3
"""Verify the compact RER mirror of the OPH-FPE history-family receipt.

This standard-library checker verifies canonical bytes, the projection digest,
the published exact counts, and every promotion/nonclaim gate. It is not the
raw-history replay: the pinned simulator custody record names the producer-free
verifier that reconstructs every cutoff separately from its own embedded raw
observer log.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
RECEIPT = HERE / "source_causal_history_family_publication_projection.json"

EXPECTED_CUTOFFS = [4, 8, 16, 32, 64]
EXPECTED_EVENT_COUNTS = [24, 48, 96, 192, 384]
EXPECTED_EDGE_COUNTS = [38, 78, 158, 318, 638]
EXPECTED_COMPARABLE_COUNTS = [244, 1060, 4420, 18052, 72964]
EXPECTED_HEIGHTS = [12, 24, 48, 96, 192]
EXPECTED_WIDTHS = [2, 2, 2, 2, 2]
EXPECTED_ORDERING_FRACTIONS = [
    0.8840579710144928,
    0.9397163120567376,
    0.9692982456140351,
    0.9845113438045375,
    0.992221496953873,
]


def canonical_bytes(value: object) -> bytes:
    return (
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        )
        + "\n"
    ).encode("ascii")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(message)


def main() -> None:
    raw = RECEIPT.read_bytes()
    report = json.loads(raw.decode("ascii"))
    require(raw == canonical_bytes(report), "receipt is not canonical JSON")

    body = {key: value for key, value in report.items() if key != "projection_sha256"}
    digest = "sha256:" + hashlib.sha256(canonical_bytes(body)).hexdigest()
    require(report.get("projection_sha256") == digest, "projection digest mismatch")
    require(
        report.get("schema")
        == "oph.source-causal-history-family-publication-projection.v1",
        "unexpected schema",
    )
    require(
        report.get("status")
        == "CERTIFIED_INFORMATIONAL_HISTORY_EXTENSION_FAMILY__FIXED_WIDTH_NOT_SPACETIME_REFINEMENT",
        "unexpected status",
    )
    require(report.get("INFORMATIONAL_HISTORY_EXTENSION_FAMILY_RECEIPT") is True,
            "history extension receipt is not attained")
    require(report.get("INFORMATIONAL_INDUCED_PREFIX_REFINEMENT_RECEIPT") is True,
            "induced-prefix custody receipt is not attained")
    require(
        report.get("INFORMATIONAL_INDEPENDENT_CUTOFF_GENERATION_RECEIPT") is True,
        "independent cutoff-generation receipt is not attained",
    )
    require(report.get("all_cutoffs_independently_generated") is True,
            "cutoffs were not independently generated")
    require(report.get("all_induced_order_embeddings_certified") is True,
            "not every induced embedding is certified")
    require(report.get("controls_fail_closed") is True, "controls do not fail closed")

    levels = report.get("levels", [])
    require([row.get("complete_round_cutoff") for row in levels] == EXPECTED_CUTOFFS,
            "round cutoffs changed")
    require([row.get("event_count") for row in levels] == EXPECTED_EVENT_COUNTS,
            "event counts changed")
    require([row.get("direct_edge_count") for row in levels] == EXPECTED_EDGE_COUNTS,
            "direct-edge counts changed")
    require([row.get("comparable_pair_count") for row in levels]
            == EXPECTED_COMPARABLE_COUNTS, "comparable-pair counts changed")
    require([row.get("height") for row in levels] == EXPECTED_HEIGHTS,
            "heights changed")
    require([row.get("width") for row in levels] == EXPECTED_WIDTHS,
            "widths changed")
    require([row.get("ordering_fraction") for row in levels]
            == EXPECTED_ORDERING_FRACTIONS, "ordering fractions changed")
    require(all(row.get("exact_width_certificate") is True for row in levels),
            "an exact-width certificate failed")
    require(all(row.get("observer_chain_cover", {}).get("chain_count") == 2
                for row in levels), "observer chain cover changed")
    evidence_hashes = report.get("cutoff_run_evidence_sha256s", [])
    require(len(evidence_hashes) == len(EXPECTED_CUTOFFS)
            and len(set(evidence_hashes)) == len(EXPECTED_CUTOFFS)
            and all(isinstance(value, str) and value.startswith("sha256:")
                    and len(value) == 71 for value in evidence_hashes),
            "independent cutoff evidence hashes are missing or repeated")
    require(
        all(
            row.get("generated_from_own_cutoff_capture") is True
            and row.get("independent_cutoff_run_evidence_sha256") == evidence_hash
            for row, evidence_hash in zip(levels, evidence_hashes, strict=True)
        ),
        "a level is not bound to its own cutoff capture",
    )

    embeddings = report.get("induced_order_embeddings", [])
    require(
        [(row.get("from_complete_round_cutoff"), row.get("to_complete_round_cutoff"))
         for row in embeddings]
        == list(zip(EXPECTED_CUTOFFS, EXPECTED_CUTOFFS[1:])),
        "embedding ladder changed",
    )
    require(all(row.get("proper_carrier_inclusion") is True
                and row.get("direct_order_is_induced_restriction") is True
                and row.get("transitive_order_is_induced_restriction") is True
                for row in embeddings), "an induced-order certificate failed")

    cone = report.get("prescribed_single_frame_source_port_placement", {})
    require(cone.get("event_count") == 24, "placement carrier changed")
    require(cone.get("comparable_pair_count") == 244, "placement comparables changed")
    require(cone.get("incomparable_pair_count") == 32, "placement incomparables changed")
    require(cone.get("coincident_same_rank_incomparable_pair_count") == 8,
            "same-rank collision count changed")
    require(cone.get("incomparable_zero_spatial_separation_pair_count") == 8,
            "zero-separation incomparable count changed")
    require(cone.get("causal_lower_time_scale_bound") == 1.129775730952839,
            "causal lower bound changed")
    require(cone.get("spacelike_upper_time_scale_bound") == 0.0,
            "spacelike upper bound changed")
    require(cone.get("all_precedence_pairs_future_causal_at_lower_bound") is True,
            "one-way cone certificate failed")
    require(cone.get("injective_four_coordinate_map") is False,
            "prescribed placement unexpectedly became injective")
    require(cone.get("precedence_iff_future_causal") is False,
            "prescribed placement unexpectedly became order-reflecting")
    require(cone.get("global_time_scale_interval_nonempty") is False,
            "prescribed placement unexpectedly gained an admissible scale")
    require(cone.get("FINITE_FAITHFUL_RANK3_CONE_PLACEMENT_RECEIPT") is False,
            "faithful-placement gate changed")
    for field in (
        "inter_carrier_frame_gluing_source_derived",
        "consumed_record_barycentre_rule_source_derived",
        "other_source_selected_placements_excluded",
        "physical_no_go_for_other_source_selected_placements",
    ):
        require(cone.get(field) is False, f"placement boundary changed: {field}")
    require("not a no-go" in cone.get("interpretation", ""),
            "placement interpretation lost its no-go boundary")

    flags = report.get("promotion_and_nonclaim_flags", {})
    require(flags and all(value is False for value in flags.values()),
            "a physical promotion flag is true or missing")
    require(report.get("physical_promotion_allowed") is not True,
            "top-level physical promotion is allowed")
    require(report.get("full_receipt_file_sha256", "").startswith("sha256:"),
            "full receipt file digest is missing")
    require(report.get("full_receipt_report_sha256", "").startswith("sha256:"),
            "full receipt report digest is missing")

    print(json.dumps({
        "status": "VERIFIED_COMPACT_PUBLICATION_PROJECTION",
        "projection_sha256": digest,
        "full_receipt_file_sha256": report["full_receipt_file_sha256"],
        "levels": len(levels),
        "event_counts": EXPECTED_EVENT_COUNTS,
        "widths": EXPECTED_WIDTHS,
        "prescribed_placement_faithful": False,
        "physical_promotion_allowed": False,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
