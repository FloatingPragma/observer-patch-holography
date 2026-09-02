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

    compiler = report.get("constructive_source_grammar_3p1_control", {})
    require(
        compiler.get("status")
        == "ATTAINED_GEOMETRY_SEEDED_3P1_CAUSESET_TO_PROVENANCE_COMPILER_CONTROL",
        "geometry-seeded causet compiler status changed",
    )
    require(
        compiler.get("GEOMETRY_SEEDED_3P1_CAUSESET_TO_PROVENANCE_COMPILER_CONTROL_RECEIPT")
        is True,
        "geometry-seeded causet compiler control is not attained",
    )
    require(
        compiler.get("compiler_acceptance_independent_of_seeded_order_embedding_and_count_volume_controls")
        is True,
        "compiler acceptance was not separated from seeded diagnostics",
    )
    require(
        compiler.get("GEOMETRY_SEEDED_FINITE_ORDER_EMBEDDING_CONTROL_RECEIPT")
        is True
        and compiler.get("GEOMETRY_SEEDED_COUNT_VOLUME_CALIBRATION_CONTROL_RECEIPT")
        is True,
        "a separately evaluated seeded control failed",
    )
    require(
        compiler.get("all_seeded_order_embeddings_attained") is True
        and compiler.get("causet_faithful_embedding_seeded_by_construction") is True
        and compiler.get("GEOMETRY_SEEDED_CAUSESET_FAITHFUL_EMBEDDING_CONTROL_RECEIPT")
        is True,
        "the seeded order-plus-density faithful-embedding control failed",
    )
    faithful_definition = compiler.get("causet_faithful_embedding_control_definition", {})
    require(
        faithful_definition.get("aggregate_predicate")
        == "all_seeded_order_embeddings_attained AND GEOMETRY_SEEDED_COUNT_VOLUME_CALIBRATION_CONTROL_RECEIPT"
        and faithful_definition.get("geometry_density_and_poisson_process_seeded_upstream")
        is True
        and faithful_definition.get("OPH_source_derived") is False,
        "causal-set faithful-embedding definition or nonpromotion boundary changed",
    )
    require(
        compiler.get("source_provenance_grammar_uses_no_declared_parent_lists") is True,
        "declared-parent-list-free grammar boundary changed",
    )
    require(
        compiler.get("SOURCE_GRAMMAR_GEOMETRY_SEEDED_3P1_FAITHFUL_REFINEMENT_CONTROL_RECEIPT")
        is False
        and compiler.get("OPH_SOURCE_LOG_REFINEMENT_RECEIPT") is False,
        "seeded controls were unexpectedly promoted to source-log refinement",
    )
    require(
        compiler.get("regulator_specific_event_material") is True
        and compiler.get("shared_semantic_commit_material_preserved") is False
        and compiler.get("only_carrier_and_transitive_order_are_induced") is True,
        "regulator-specific semantic-material boundary changed",
    )
    require(
        compiler.get("aggregate_dimension_diagnostic_role")
        == "DIAGNOSTIC_ONLY_NOT_COMPILER_ACCEPTANCE_CRITERION",
        "dimension diagnostic was promoted into compiler acceptance",
    )
    require(
        compiler.get("aggregate_dimension_diagnostic")
        == [
            {
                "target_poisson_mean_count": 64,
                "mean_myrheim_meyer_dimension_estimate": 4.028237157921849,
                "absolute_error_from_seeded_dimension_four": 0.028237157921848777,
            },
            {
                "target_poisson_mean_count": 128,
                "mean_myrheim_meyer_dimension_estimate": 3.845947760942414,
                "absolute_error_from_seeded_dimension_four": 0.1540522390575858,
            },
            {
                "target_poisson_mean_count": 256,
                "mean_myrheim_meyer_dimension_estimate": 3.9034561962046683,
                "absolute_error_from_seeded_dimension_four": 0.09654380379533167,
            },
        ],
        "seeded Myrheim-Meyer diagnostics changed",
    )
    families = compiler.get("families", [])
    require(len(families) == 4 and all(len(family.get("levels", [])) == 3 for family in families),
            "expected all twelve geometry-seeded control levels")
    require(
        [[level.get("event_count") for level in family.get("levels", [])]
         for family in families]
        == [[70, 136, 281], [62, 123, 251], [60, 130, 254], [51, 120, 244]],
        "seeded compiler populations changed",
    )
    mutation_names = {
        "one_causal_support_deletion_breaks_order_identity",
        "duplicate_writer_injection_rejects_without_accepted_order",
        "last_writer_mutation_breaks_pre_state_parent_clause",
        "read_value_mutation_breaks_parent_post_value_continuity",
        "version_mutation_breaks_version_continuity",
        "missing_writer_mutation_rejects_unresolved_read",
        "sequence_mutation_rejects_nonprior_writer",
        "self_read_mutation_rejects_self_parent",
    }
    require(
        all(
            family.get("OPH_SOURCE_LOG_REFINEMENT_RECEIPT") is False
            and family.get("shared_semantic_commit_material_preserved") is False
            and all(
                level.get("overall_acceptance") is True
                and level.get("provenance_direct_edges_equal_geometric_cover_relation") is True
                and level.get("provenance_closure_equals_full_geometric_causal_order") is True
                and level.get("declared_parent_lists_absent") is True
                and level.get("coordinate_map_injective") is True
                and level.get("coordinates_inside_fixed_minkowski_alexandrov_interval") is True
                and level.get("finite_order_embedding_control") is True
                and level.get("seeded_order_embedding_by_construction") is True
                and level.get("count_within_declared_three_sigma_poisson_interval") is True
                and set(level.get("negative_controls", {})) == mutation_names
                and all(level.get("negative_controls", {}).values())
                for level in family.get("levels", [])
            )
            for family in families
        ),
        "a seeded compiler reconstruction or mutation control failed",
    )
    require(
        compiler.get("OPH_SOURCE_SELECTS_THIS_EVENT_POPULATION_OR_READ_PATTERN")
        is False
        and compiler.get("OPH_CAUSAL_SET_SIMILARITY_RECEIPT") is False
        and compiler.get("OPH_CAUSAL_3P1_MANIFOLD_DERIVATION_RECEIPT") is False
        and compiler.get("physical_promotion_allowed") is False,
        "seeded compiler control was physically promoted",
    )

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
        "geometry_seeded_3p1_compiler_control": True,
        "physical_promotion_allowed": False,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
