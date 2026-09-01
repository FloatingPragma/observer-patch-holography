#!/usr/bin/env python3
"""Verify the RER mirror of the exploratory causet-likeness receipt.

This standard-library checker validates canonical bytes, the receipt's
self-digest, the exact published statistics, the replay/control gates, and
every physical-promotion nonclaim.  Reconstructing intervals and the synthetic
Minkowski/de Sitter controls from raw simulator data remains the job of the
producer-free verifier pinned by the paper's simulator revision.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
RECEIPT = HERE / "causet_likeness_receipt.json"


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


def require_statistics(actual: object, expected: dict[str, object], label: str) -> None:
    require(isinstance(actual, dict), f"{label} statistics are missing")
    for field, value in expected.items():
        require(actual.get(field) == value, f"{label} statistic changed: {field}")


def main() -> None:
    raw = RECEIPT.read_bytes()
    report = json.loads(raw.decode("ascii"))
    require(raw == canonical_bytes(report), "receipt is not canonical JSON")

    body = {key: value for key, value in report.items() if key != "report_sha256"}
    digest = "sha256:" + hashlib.sha256(canonical_bytes(body)).hexdigest()
    require(report.get("report_sha256") == digest, "report digest mismatch")
    require(report.get("schema") == "oph.causet-likeness-exploratory.v1",
            "unexpected schema")
    require(
        report.get("status")
        == "NOT_SIMILAR_AT_CURRENT_CUTOFF__INTERVAL_ORDERING_FRACTIONS_OUTSIDE_EXPLORATORY_4D_BAND",
        "unexpected current-cutoff result",
    )
    require(report.get("CAUSET_DIAGNOSTIC_PIPELINE_REPRODUCTION_RECEIPT") is True,
            "diagnostic pipeline was not reproduced")
    require(report.get("OPH_CAUSAL_SET_SIMILARITY_RECEIPT") is False,
            "similarity receipt was unexpectedly promoted")
    require(report.get("CAUSET_MANIFOLDLIKE_RECEIPT") is False,
            "manifoldlikeness receipt was unexpectedly promoted")
    require(report.get("physical_promotion_allowed") is False,
            "physical promotion was unexpectedly allowed")
    require(report.get("controls_fail_closed") is True,
            "controls do not fail closed")

    require_statistics(
        report.get("source_statistics"),
        {
            "event_count": 24,
            "input_edge_count": 38,
            "comparable_pair_count": 244,
            "height": 12,
            "width": 2,
            "maximum_interval_size": 21,
            "adequate_dimension_interval_count": 0,
            "global_ordering_fraction": 0.8840579710144928,
        },
        "source control",
    )
    require(
        report.get("source_control_status")
        == "INCONCLUSIVE__INSUFFICIENT_CERTIFIED_INTERVAL_SIZE",
        "source-control size boundary changed",
    )

    local = report.get("existing_local_domain_diagnostic", {})
    require(isinstance(local, dict), "local-domain diagnostic is missing")
    require(local.get("status") == "EXPLORATORY_INTERVAL_DIAGNOSTIC_READY_NOT_OPH_SIMILARITY",
            "local-domain diagnostic status changed")
    require_statistics(
        local.get("statistics"),
        {
            "event_count": 2304,
            "input_edge_count": 3483,
            "comparable_pair_count": 70661,
            "height": 18,
            "width": 128,
            "maximum_interval_size": 72,
            "adequate_dimension_interval_count": 736,
            "global_ordering_fraction": 0.026633813986587544,
        },
        "local domain",
    )

    comparison = report.get("single_cutoff_matched_interval_comparison", {})
    require(isinstance(comparison, dict), "matched-interval comparison is missing")
    require(comparison.get("epistemic_status") == "POST_HOC_EXPLORATORY_SINGLE_CUTOFF",
            "single-cutoff epistemic status changed")
    require(comparison.get("adequate_interval_count") == 736,
            "matched interval count changed")
    require(comparison.get("inclusive_cardinality_range") == [32, 72],
            "matched cardinality range changed")
    require(comparison.get("exploratory_4d_ordering_fraction_band") == [0.05, 0.2],
            "exploratory four-dimensional band changed")
    require(comparison.get("interval_count_in_exploratory_4d_band") == 0,
            "an interval unexpectedly entered the exploratory band")
    require(comparison.get("interval_fraction_in_exploratory_4d_band") == 0.0,
            "exploratory-band fraction changed")
    require(
        comparison.get("ordering_fraction_quantiles")
        == {
            "maximum": 0.8465909090909091,
            "median": 0.5719512195121951,
            "minimum": 0.33365384615384613,
            "q25": 0.5135135135135135,
            "q75": 0.6369747899159663,
        },
        "matched-interval ordering-fraction quantiles changed",
    )
    require(comparison.get("heldout_negative_controls_distinguished") is True,
            "held-out negative controls were not distinguished")

    invariance = report.get("invariance_controls", {})
    require(isinstance(invariance, dict) and invariance
            and all(value is True for value in invariance.values()),
            "an invariance or replay control failed")
    refinement = report.get("oph_refinement_family_comparison", {})
    require(isinstance(refinement, dict), "refinement comparison is missing")
    require(refinement.get("certified_oph_refinement_map_or_family_available") is False,
            "a refinement family was unexpectedly certified")
    require(refinement.get("similarity_claimed") is False,
            "similarity was unexpectedly claimed")
    require(
        refinement.get("status")
        == "NOT_EVALUATED_NO_CERTIFIED_OPH_REFINEMENT_MAP_OR_FAMILY",
        "refinement-family boundary changed",
    )
    require(report.get("held_out_confirmation_status") == "NOT_RUN_EXPLORATORY_ONLY",
            "held-out confirmation status changed")
    require(
        report.get("flrw_reference_control_status")
        == "IMPLEMENTED_REPLAYED_DE_SITTER_FLAT_PATCH_SPECIAL_FLRW_REFERENCE_ONLY",
        "curved-reference scope changed",
    )

    binding = report.get("source_binding", {})
    require(isinstance(binding, dict), "source binding is missing")
    require(binding.get("source_receipt_path")
            == "data/causal_order/source_derived_causal_order_receipt.json",
            "source receipt path changed")
    for field in (
        "generated_edges_sha256",
        "observer_event_log_sha256",
        "semantic_event_keys_sha256",
        "semantic_events_sha256",
        "semantic_poset_sha256",
        "source_report_sha256",
    ):
        value = binding.get(field)
        require(isinstance(value, str) and value.startswith("sha256:")
                and len(value) == 71, f"source binding is missing: {field}")

    boundary = report.get("claim_boundary", "")
    for phrase in (
        "does not establish physical causal faithfulness",
        "not a physical nonexistence result",
        "not the complete seam-repair event history",
    ):
        require(phrase in boundary, f"claim boundary lost: {phrase}")
    require(
        report.get("required_next_capture")
        == "promote the versioned seam-repair transactions as the physical event carrier, construct a certified OPH refinement family/map, and freeze a held-out comparison before evaluating that family",
        "required next physical capture changed",
    )

    print(json.dumps({
        "status": "VERIFIED_RER_CAUSET_LIKENESS_RECEIPT",
        "report_sha256": digest,
        "source_events": 24,
        "local_domain_events": 2304,
        "matched_intervals": 736,
        "intervals_in_exploratory_4d_band": 0,
        "physical_promotion_allowed": False,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
