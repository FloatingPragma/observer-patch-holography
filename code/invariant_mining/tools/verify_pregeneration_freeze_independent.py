#!/usr/bin/env python3
"""Independently verify the issue-647 source-only pre-generation freeze.

This verifier intentionally imports no builder and no project helper.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
from typing import Any

STATUS = (
    "GENERATION_ENABLED__"
    "COMPARISON_SEALED"
)
EXPECTED_FEATURE_IDS = {
    "conditional_maximum_randomness",
    "diagonal_z6_descent",
    "direct_n_bounded_counterfamily_nonidentifiability",
    "family_band_response",
    "finite_local_domain_boundary",
    "finite_port_action",
    "observer_agreement_and_readback",
    "oriented_twelve_port_carrier",
}
EXPECTED_SLOT_IDS = {
    "a5_angular_rules",
    "cross_channel_source_identities",
    "direct_n_capacity_closure",
    "family_trace_mixing_invariants",
    "observer_overlap_cross_spectra",
    "wz_scale_free_response",
    "z6_charge_line_congruences",
}
EXPECTED_GRAMMAR_IDS = {
    "angular_selection_cross_level_correlation",
    "discrete_congruence_parity_ordering_multiplicity",
    "exact_zero_forbidden_transition",
    "homogeneous_scale_free_combination",
    "normalized_residue_asymmetry_cross_spectrum",
    "refinement_scaling_exponent",
    "rg_invariant_combination",
    "shared_response_identity",
    "trace_determinant_character_index",
}
EXPECTED_NUISANCE_TYPES = {
    "alternative_source_admissible_completion",
    "continuous_physical_parameter",
    "coordinate_equivalence_nuisance",
    "discrete_branch_selector",
    "measurement_calibration_nuisance",
}
EXPECTED_NUISANCE_IDS = {
    "basis_and_coordinate_choice",
    "common_dimensionful_scale",
    "detector_transfer_calibration",
    "physical_sector_choice",
    "source_admissible_completion",
}
EXPECTED_EXPOSURE_CLASSES = [
    "TARGET_CLEAN_PROSPECTIVE_HOLDOUT",
    "TARGET_ISOLATED_BLIND_POSTDICTION",
    "EXPOSED_RETROSPECTIVE_COMPARISON",
    "EXPLORATORY_DISCOVERY_DATA_CATALOG_ONLY",
]
EXPECTED_EXPOSURE_IDS = {
    "cmb_angular_products",
    "codata_fundamental_constants",
    "gravitational_wave_catalogs",
    "neutrino_oscillation_global_fits",
    "pdg_particle_listings",
    "precision_qed_ratios",
}
EXPECTED_OUTPUT_NAMES = {
    "pregeneration_freeze.json",
    "source_projection.json",
}
OPTIONAL_OUTPUT_NAMES = {
    "candidate_registry.json",
}


class VerificationError(ValueError):
    """A frozen invariant-mining control failed independent verification."""


def require(condition: bool, code: str, detail: str) -> None:
    if not condition:
        raise VerificationError(f"{code}: {detail}")


def canonical_bytes(value: Any) -> bytes:
    return (
        json.dumps(value, sort_keys=True, indent=2, ensure_ascii=False) + "\n"
    ).encode("utf-8")


def digest(payload: bytes) -> str:
    return "sha256:" + hashlib.sha256(payload).hexdigest()


def load(path: Path) -> dict[str, Any]:
    require(path.is_file(), "DOCUMENT_MISSING", str(path))
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise VerificationError(f"DOCUMENT_INVALID: {path}: {error}") from error
    require(isinstance(value, dict), "DOCUMENT_NOT_OBJECT", str(path))
    return value


def repo_path(value: Any) -> str:
    require(isinstance(value, str) and value != "", "PATH_INVALID", repr(value))
    require("\\" not in value, "PATH_NOT_POSIX", value)
    parsed = PurePosixPath(value)
    require(
        not parsed.is_absolute()
        and parsed.as_posix() == value
        and "." not in parsed.parts
        and ".." not in parsed.parts,
        "PATH_NOT_CANONICAL",
        value,
    )
    return value


def file_bytes(repo_root: Path, relative: str) -> bytes:
    relative = repo_path(relative)
    unresolved = repo_root / relative
    require(not unresolved.is_symlink(), "SYMLINK_REFUSED", relative)
    root = repo_root.resolve()
    resolved = unresolved.resolve()
    require(root in resolved.parents, "PATH_ESCAPE", relative)
    require(resolved.is_file(), "PINNED_FILE_MISSING", relative)
    payload = resolved.read_bytes()
    require(len(payload) > 0, "PINNED_FILE_EMPTY", relative)
    return payload


def pin(repo_root: Path, relative: str) -> dict[str, Any]:
    payload = file_bytes(repo_root, relative)
    return {"path": relative, "bytes": len(payload), "sha256": digest(payload)}


def rows_by_id(rows: Any, key: str, label: str) -> dict[str, dict[str, Any]]:
    require(isinstance(rows, list), "REGISTRY_ROWS_INVALID", label)
    result: dict[str, dict[str, Any]] = {}
    for row in rows:
        require(isinstance(row, dict), "REGISTRY_ROW_INVALID", label)
        identifier = row.get(key)
        require(
            isinstance(identifier, str) and identifier,
            "REGISTRY_ID_INVALID",
            f"{label}.{key}",
        )
        require(identifier not in result, "REGISTRY_ID_DUPLICATE", identifier)
        result[identifier] = row
    return result


def all_keys(value: Any) -> set[str]:
    result: set[str] = set()
    if isinstance(value, dict):
        for key, child in value.items():
            result.add(str(key))
            result.update(all_keys(child))
    elif isinstance(value, list):
        for child in value:
            result.update(all_keys(child))
    return result


def verify(
    package_root: Path,
    repo_root: Path,
) -> None:
    policy = load(package_root / "policy" / "pregeneration_policy.json")
    source = load(package_root / "data" / "source_feature_registry.json")
    producers = load(package_root / "data" / "producer_slots.json")
    grammar = load(package_root / "data" / "observable_grammar.json")
    nuisances = load(package_root / "data" / "nuisance_registry.json")
    ranking = load(package_root / "data" / "ranking_policy.json")
    exposure = load(package_root / "data" / "exposed_data_registry.json")
    projection = load(package_root / "outputs" / "source_projection.json")
    freeze = load(package_root / "outputs" / "pregeneration_freeze.json")
    document_schema = load(
        package_root / "schemas" / "pregeneration_documents.schema.json"
    )
    output_schema = load(
        package_root / "schemas" / "pregeneration_outputs.schema.json"
    )

    require(
        document_schema.get("$schema")
        == "https://json-schema.org/draft/2020-12/schema",
        "SCHEMA_DRAFT_DRIFT",
        "pregeneration_documents",
    )
    require(
        output_schema.get("$schema")
        == "https://json-schema.org/draft/2020-12/schema",
        "SCHEMA_DRAFT_DRIFT",
        "pregeneration_outputs",
    )

    require(
        policy.get("schema") == "oph.invariant_mining.pregeneration_policy.v3",
        "POLICY_SCHEMA_DRIFT",
        str(policy.get("schema")),
    )
    require(policy.get("status") == STATUS, "POLICY_STATUS_DRIFT", str(policy.get("status")))
    require(
        policy.get("registry_finalization_complete") is True,
        "REGISTRY_FINALIZATION_FLAG_DRIFT",
        str(policy.get("registry_finalization_complete")),
    )
    for key in ("candidate_generator_enabled", "candidate_evaluator_enabled"):
        require(policy.get(key) is True, "ENABLEMENT_FLAG_DRIFT", key)
    for key in ("public_data_access_enabled", "target_payloads_registered"):
        require(policy.get(key) is False, "EXECUTION_BOUNDARY_OPEN", key)
    require(
        policy.get("declared_candidate_output")
        == "code/invariant_mining/outputs/candidate_registry.json",
        "DECLARED_CANDIDATE_OUTPUT_DRIFT",
        str(policy.get("declared_candidate_output")),
    )
    require(policy.get("candidate_count") == 0, "CANDIDATE_COUNT_NONZERO", str(policy.get("candidate_count")))
    budget = policy.get("campaign_comparison_budget")
    require(isinstance(budget, dict), "COMPARISON_BUDGET_MISSING", str(budget))
    require(
        budget.get("maximum_physical_comparisons") == 1
        and budget.get("comparisons_consumed") == 0
        and budget.get("terminate_after_first_physical_comparison") is True,
        "COMPARISON_BUDGET_DRIFT",
        str(budget),
    )
    require(
        set(policy.get("required_exposure_ids", [])) == EXPECTED_EXPOSURE_IDS,
        "POLICY_EXPOSURE_SET_DRIFT",
        "required_exposure_ids",
    )
    require(
        exposure.get("schema")
        == "oph.invariant_mining.exposed_data_registry.v1",
        "EXPOSURE_SCHEMA_DRIFT",
        str(exposure.get("schema")),
    )
    require(
        exposure.get("freeze_status") == "FROZEN_BEFORE_CANDIDATE_GENERATION",
        "EXPOSURE_NOT_FROZEN",
        str(exposure.get("freeze_status")),
    )
    exposure_rows = rows_by_id(
        exposure.get("surfaces"), "exposure_id", "exposed-data surfaces"
    )
    require(
        set(exposure_rows) == EXPECTED_EXPOSURE_IDS,
        "EXPOSURE_SURFACE_OMISSION",
        str(sorted(exposure_rows)),
    )
    require(
        exposure.get("exposure_classes") == EXPECTED_EXPOSURE_CLASSES,
        "EXPOSURE_CLASS_VOCABULARY_DRIFT",
        str(exposure.get("exposure_classes")),
    )
    producer_slot_ids = set(
        rows_by_id(producers.get("slots"), "slot_id", "producer slots")
    )
    covered_by_exposure: set[str] = set()
    for exposure_id, row in exposure_rows.items():
        require(
            row.get("quarantine_evidence") is None,
            "EXPOSURE_QUARANTINE_PRECLAIMED",
            exposure_id,
        )
        require(
            row.get("default_exposure_class") in EXPECTED_EXPOSURE_CLASSES,
            "EXPOSURE_CLASS_UNKNOWN",
            exposure_id,
        )
        applies = row.get("applies_to_slot_ids")
        require(
            isinstance(applies, list)
            and applies
            and set(applies) <= producer_slot_ids,
            "EXPOSURE_SLOT_UNKNOWN",
            exposure_id,
        )
        covered_by_exposure.update(applies)
    direct_slot_id = "direct_n_capacity_closure"
    exposure_gap = producer_slot_ids - covered_by_exposure - {direct_slot_id}
    require(
        not exposure_gap,
        "EXPOSURE_COVERAGE_GAP",
        str(sorted(exposure_gap)),
    )
    completion_rows = [
        row
        for row in (nuisances.get("unresolved_directions") or [])
        if isinstance(row, dict)
        and row.get("nuisance_id") == "source_admissible_completion"
    ]
    require(
        len(completion_rows) == 1,
        "COMPLETION_NUISANCE_MISSING",
        str(len(completion_rows)),
    )
    require(
        set(completion_rows[0].get("applies_to_slot_ids", []))
        == producer_slot_ids,
        "COMPLETION_NUISANCE_COVERAGE_GAP",
        str(sorted(producer_slot_ids - set(completion_rows[0].get("applies_to_slot_ids", [])))),
    )

    require(
        set(policy.get("required_feature_ids", [])) == EXPECTED_FEATURE_IDS,
        "POLICY_FEATURE_SET_DRIFT",
        "required_feature_ids",
    )
    require(
        set(policy.get("required_producer_slot_ids", [])) == EXPECTED_SLOT_IDS,
        "POLICY_SLOT_SET_DRIFT",
        "required_producer_slot_ids",
    )
    require(
        set(policy.get("required_grammar_class_ids", [])) == EXPECTED_GRAMMAR_IDS,
        "POLICY_GRAMMAR_SET_DRIFT",
        "required_grammar_class_ids",
    )
    require(
        set(policy.get("required_nuisance_type_ids", []))
        == EXPECTED_NUISANCE_TYPES,
        "POLICY_NUISANCE_TYPE_SET_DRIFT",
        "required_nuisance_type_ids",
    )
    require(
        set(policy.get("required_nuisance_ids", [])) == EXPECTED_NUISANCE_IDS,
        "POLICY_NUISANCE_SET_DRIFT",
        "required_nuisance_ids",
    )

    require(
        source.get("schema")
        == "oph.invariant_mining.source_feature_registry.v1",
        "SOURCE_REGISTRY_SCHEMA_DRIFT",
        str(source.get("schema")),
    )
    require(
        source.get("completeness", {}).get("status")
        == "FINAL_FOR_FROZEN_CAMPAIGN_SCOPE",
        "SOURCE_REGISTRY_FALSE_COMPLETENESS",
        str(source.get("completeness")),
    )
    feature_rows = rows_by_id(source.get("features"), "feature_id", "features")
    require(set(feature_rows) == EXPECTED_FEATURE_IDS, "SOURCE_FEATURE_OMISSION", str(sorted(feature_rows)))

    require(
        producers.get("schema") == "oph.invariant_mining.producer_slots.v1",
        "PRODUCER_REGISTRY_SCHEMA_DRIFT",
        str(producers.get("schema")),
    )
    require(
        producers.get("execution_state")
        == "ALL_PRODUCERS_DISABLED_PENDING_ENABLEMENT_REVIEW",
        "PRODUCER_EXECUTION_STATE_OPEN",
        str(producers.get("execution_state")),
    )
    slot_rows = rows_by_id(producers.get("slots"), "slot_id", "slots")
    require(set(slot_rows) == EXPECTED_SLOT_IDS, "PRODUCER_SLOT_OMISSION", str(sorted(slot_rows)))

    require(
        grammar.get("schema") == "oph.invariant_mining.observable_grammar.v1",
        "GRAMMAR_SCHEMA_DRIFT",
        str(grammar.get("schema")),
    )
    require(
        grammar.get("freeze_status") == "FROZEN_BEFORE_CANDIDATE_GENERATION",
        "GRAMMAR_NOT_FROZEN",
        str(grammar.get("freeze_status")),
    )
    grammar_ids = grammar.get("class_ids")
    require(
        isinstance(grammar_ids, list)
        and set(grammar_ids) == EXPECTED_GRAMMAR_IDS
        and len(grammar_ids) == len(set(grammar_ids)),
        "GRAMMAR_CLASS_DRIFT",
        str(grammar_ids),
    )
    budget = grammar.get("complexity_budget")
    require(isinstance(budget, dict), "GRAMMAR_BUDGET_MISSING", "complexity_budget")
    for key in (
        "maximum_ast_depth",
        "maximum_ast_nodes",
        "maximum_distinct_registered_terminals",
        "maximum_generated_candidates",
    ):
        require(
            isinstance(budget.get(key), int) and budget[key] > 0,
            "GRAMMAR_BUDGET_INVALID",
            key,
        )
    require(
        grammar.get("equivalence_relation", {}).get(
            "physical_parameters_or_branch_selectors_never_quotiented"
        )
        is True,
        "PHYSICAL_NUISANCE_QUOTIENTED",
        "grammar equivalence relation",
    )

    require(
        nuisances.get("schema") == "oph.invariant_mining.nuisance_registry.v1",
        "NUISANCE_SCHEMA_DRIFT",
        str(nuisances.get("schema")),
    )
    require(
        nuisances.get("freeze_status") == "FROZEN_BEFORE_CANDIDATE_GENERATION",
        "NUISANCE_REGISTRY_NOT_FROZEN",
        str(nuisances.get("freeze_status")),
    )
    taxonomy = rows_by_id(nuisances.get("taxonomy"), "type_id", "taxonomy")
    unresolved = rows_by_id(
        nuisances.get("unresolved_directions"),
        "nuisance_id",
        "unresolved_directions",
    )
    require(set(taxonomy) == EXPECTED_NUISANCE_TYPES, "NUISANCE_TYPE_OMISSION", str(sorted(taxonomy)))
    require(set(unresolved) == EXPECTED_NUISANCE_IDS, "NUISANCE_OMISSION", str(sorted(unresolved)))
    require(
        taxonomy["coordinate_equivalence_nuisance"].get("quotient_allowed")
        is True,
        "COORDINATE_QUOTIENT_DISABLED",
        "coordinate_equivalence_nuisance",
    )
    for identifier, row in taxonomy.items():
        if identifier != "coordinate_equivalence_nuisance":
            require(
                row.get("quotient_allowed") is False,
                "PHYSICAL_NUISANCE_QUOTIENTED",
                identifier,
            )

    require(
        ranking.get("schema") == "oph.invariant_mining.ranking_policy.v1",
        "RANKING_SCHEMA_DRIFT",
        str(ranking.get("schema")),
    )
    require(
        ranking.get("freeze_status") == "FROZEN_BEFORE_CANDIDATE_GENERATION",
        "RANKING_NOT_FROZEN",
        str(ranking.get("freeze_status")),
    )
    require(
        ranking.get("score_rule") == "integer_sum_of_applicable_frozen_weights",
        "RANKING_RULE_DRIFT",
        str(ranking.get("score_rule")),
    )
    weights = ranking.get("weights")
    require(
        isinstance(weights, dict)
        and weights
        and all(isinstance(value, int) for value in weights.values()),
        "RANKING_WEIGHTS_INVALID",
        str(weights),
    )
    require(
        ranking.get("tie_breaker")
        == [
            "lower_expression_complexity",
            "lexicographically_smaller_canonical_expression",
            "lexicographically_smaller_candidate_id",
        ],
        "RANKING_TIE_BREAKER_DRIFT",
        str(ranking.get("tie_breaker")),
    )

    source_paths: set[str] = set()
    forbidden_prefixes = tuple(policy.get("forbidden_repository_prefixes", []))
    for identifier, row in feature_rows.items():
        require(
            row.get("candidate_generation_eligible") is False,
            "FEATURE_ENABLED_EARLY",
            identifier,
        )
        require(
            row.get("target_ancestry")
            == "SOURCE_ONLY_NO_PUBLIC_COMPARISON_INPUT",
            "TARGET_ANCESTRY_DRIFT",
            identifier,
        )
        artifacts = row.get("artifacts")
        require(isinstance(artifacts, list) and artifacts, "FEATURE_ARTIFACTS_MISSING", identifier)
        for artifact in artifacts:
            require(
                isinstance(artifact, dict)
                and set(artifact) == {"path", "kind", "role"},
                "FEATURE_ARTIFACT_SHAPE_DRIFT",
                identifier,
            )
            relative = repo_path(artifact["path"])
            require(
                not relative.startswith(forbidden_prefixes),
                "COMPARISON_SURFACE_IMPORTED",
                relative,
            )
            require(relative not in source_paths, "SOURCE_ARTIFACT_DUPLICATE", relative)
            source_paths.add(relative)

    for identifier, row in slot_rows.items():
        require(
            row.get("execution_status") == "REGISTERED_DISABLED",
            "PRODUCER_ENABLED_EARLY",
            identifier,
        )
        required_features = row.get("required_feature_ids")
        allowed_classes = row.get("allowed_grammar_class_ids")
        require(
            isinstance(required_features, list)
            and required_features
            and set(required_features) <= EXPECTED_FEATURE_IDS,
            "PRODUCER_FEATURE_UNKNOWN",
            identifier,
        )
        require(
            isinstance(allowed_classes, list)
            and allowed_classes
            and set(allowed_classes) <= EXPECTED_GRAMMAR_IDS,
            "PRODUCER_GRAMMAR_UNKNOWN",
            identifier,
        )

    direct = slot_rows["direct_n_capacity_closure"]
    require(
        direct.get("fallback_eligibility") == "EXCLUDED_NO_FALLBACK_REENTRY",
        "DIRECT_N_REENTRY",
        str(direct.get("fallback_eligibility")),
    )
    require(
        direct.get("execution_status") == "REGISTERED_DISABLED",
        "DIRECT_N_ENABLED",
        str(direct.get("execution_status")),
    )
    require(
        policy.get("direct_n_contract")
        == {
            "slot_id": "direct_n_capacity_closure",
            "required_fallback_eligibility": "EXCLUDED_NO_FALLBACK_REENTRY",
            "required_execution_status": "REGISTERED_DISABLED",
        },
        "DIRECT_N_POLICY_DRIFT",
        str(policy.get("direct_n_contract")),
    )

    forbidden_keys = set(policy.get("forbidden_document_keys", []))
    for label, document in (
        ("source", source),
        ("producers", producers),
        ("grammar", grammar),
        ("nuisances", nuisances),
        ("ranking", ranking),
        ("exposure", exposure),
    ):
        leaked = sorted(all_keys(document) & forbidden_keys)
        require(not leaked, "FORBIDDEN_TARGET_KEY", f"{label}: {leaked}")

    control_paths = policy.get("control_artifacts")
    require(
        isinstance(control_paths, list)
        and control_paths == sorted(control_paths)
        and len(control_paths) == len(set(control_paths)),
        "CONTROL_PATHS_INVALID",
        str(control_paths),
    )
    required_control_paths = {
        "code/invariant_mining/data/exposed_data_registry.json",
        "code/invariant_mining/data/nuisance_registry.json",
        "code/invariant_mining/data/observable_grammar.json",
        "code/invariant_mining/data/producer_slots.json",
        "code/invariant_mining/data/ranking_policy.json",
        "code/invariant_mining/data/source_feature_registry.json",
        "code/invariant_mining/policy/pregeneration_policy.json",
        "code/invariant_mining/tools/build_pregeneration_freeze.py",
        "code/invariant_mining/tools/build_source_projection.py",
        "code/invariant_mining/tools/verify_pregeneration_freeze_independent.py",
    }
    require(
        required_control_paths <= set(control_paths),
        "CONTROL_ARTIFACT_OMISSION",
        str(sorted(required_control_paths - set(control_paths))),
    )

    expected_projection = {
        "schema": "oph.invariant_mining.source_projection.v1",
        "issue": 647,
        "status": STATUS,
        "projection_id": "oph-invariant-source-projection-v1",
        "registry_finalization_complete": True,
        "candidate_generator_enabled": True,
        "candidate_evaluator_enabled": True,
        "candidate_count": 0,
        "control_documents": [pin(repo_root, path) for path in control_paths],
        "source_artifacts": [pin(repo_root, path) for path in sorted(source_paths)],
        "registered_source_paths_sha256": digest(
            canonical_bytes(sorted(source_paths))
        ),
    }
    require(
        projection == expected_projection,
        "SOURCE_PROJECTION_DRIFT",
        "committed projection differs from independent reconstruction",
    )

    document_relatives = [
        "code/invariant_mining/data/exposed_data_registry.json",
        "code/invariant_mining/data/nuisance_registry.json",
        "code/invariant_mining/data/observable_grammar.json",
        "code/invariant_mining/data/producer_slots.json",
        "code/invariant_mining/data/ranking_policy.json",
        "code/invariant_mining/data/source_feature_registry.json",
        "code/invariant_mining/policy/pregeneration_policy.json",
        "code/invariant_mining/schemas/pregeneration_documents.schema.json",
        "code/invariant_mining/schemas/pregeneration_outputs.schema.json",
    ]
    projection_relative = "code/invariant_mining/outputs/source_projection.json"
    expected_freeze: dict[str, Any] = {
        "schema": "oph.invariant_mining.pregeneration_freeze.v1",
        "issue": 647,
        "status": STATUS,
        "freeze_id": "",
        "source_projection": pin(repo_root, projection_relative),
        "document_bindings": [
            pin(repo_root, relative) for relative in document_relatives
        ],
        "direct_n_fallback_contract": {
            "slot_id": "direct_n_capacity_closure",
            "fallback_eligibility": "EXCLUDED_NO_FALLBACK_REENTRY",
            "execution_status": "REGISTERED_DISABLED",
        },
        "execution_boundary": {
            "registry_finalization_complete": True,
            "candidate_generator_enabled": True,
            "candidate_evaluator_enabled": True,
            "public_data_access_enabled": False,
            "target_payloads_registered": False,
            "candidate_count": 0,
            "comparison_scoring_enabled": False,
        },
        "promotion_rule": policy.get("promotion_rule"),
    }
    freeze_digest_payload = dict(expected_freeze)
    freeze_digest_payload.pop("freeze_id")
    expected_freeze["freeze_id"] = digest(canonical_bytes(freeze_digest_payload))
    require(
        freeze == expected_freeze,
        "PREGENERATION_FREEZE_DRIFT",
        "committed freeze differs from independent reconstruction",
    )

    output_dir = package_root / "outputs"
    observed_output_names = {
        path.name for path in output_dir.iterdir() if path.is_file()
    }
    require(
        EXPECTED_OUTPUT_NAMES <= observed_output_names
        and observed_output_names
        <= EXPECTED_OUTPUT_NAMES | OPTIONAL_OUTPUT_NAMES,
        "UNREGISTERED_OUTPUT",
        str(sorted(observed_output_names)),
    )
    candidate_registry_path = output_dir / "candidate_registry.json"
    if candidate_registry_path.is_file():
        registry = load(candidate_registry_path)
        require(
            registry.get("schema")
            == "oph.invariant_mining.candidate_registry.v1",
            "CANDIDATE_REGISTRY_SCHEMA_DRIFT",
            str(registry.get("schema")),
        )
        require(
            registry.get("status") == "GENERATION_COMPLETE__SCORING_SEALED",
            "CANDIDATE_REGISTRY_STATUS_DRIFT",
            str(registry.get("status")),
        )
        slot_states = registry.get("slot_generation_states", {})
        require(
            set(slot_states)
            == set(rows_by_id(producers.get("slots"), "slot_id", "slots"))
            and "REGISTERED_AWAITING_GENERATION" not in slot_states.values(),
            "CANDIDATE_REGISTRY_INCOMPLETE",
            str(slot_states),
        )
        skip_records = registry.get("skip_records", {})
        frozen_skip_conditions = set(
            ranking.get("skip_policy", {}).get("skip_conditions", [])
        )
        for slot_id, state in slot_states.items():
            if state == "NOT_EVALUABLE_TYPED_SKIP":
                record = skip_records.get(slot_id)
                require(
                    isinstance(record, dict)
                    and record.get("skip_condition") in frozen_skip_conditions,
                    "SKIP_RECORD_INVALID",
                    slot_id,
                )
        boundary = registry.get("scoring_boundary", {})
        require(
            boundary.get("physical_scoring_permitted") is False
            and boundary.get("comparison_access_permitted") is False,
            "CANDIDATE_SCORING_UNSEALED",
            str(boundary),
        )
        cleanliness = registry.get("target_cleanliness", {})
        require(
            set(cleanliness)
            == {
                "public_measurement_read",
                "comparison_payload_read",
                "target_value_read",
            }
            and all(value is False for value in cleanliness.values()),
            "CANDIDATE_TARGET_CLEANLINESS_DRIFT",
            str(cleanliness),
        )
        expected_table = []
        for triality in range(3):
            for duality in range(2):
                for six_y in range(6):
                    phase = (2 * triality + 3 * duality + six_y) % 6
                    expected_table.append(
                        {
                            "triality": triality,
                            "duality": duality,
                            "six_y_mod_6": six_y,
                            "kernel_phase_class": phase,
                            "descends": phase == 0,
                        }
                    )
        require(
            registry.get("descent_congruence_table") == expected_table,
            "DESCENT_TABLE_DRIFT",
            "committed table differs from independent arithmetic",
        )
        grammar_budget = grammar["complexity_budget"]
        weights = ranking["weights"]
        forbidden_terminals = set(grammar["forbidden_terminals"])
        slot_ids = set(rows_by_id(producers.get("slots"), "slot_id", "slots"))
        candidates = registry.get("candidates")
        require(
            isinstance(candidates, list)
            and 0
            < len(candidates)
            <= grammar_budget["maximum_generated_candidates"],
            "CANDIDATE_COUNT_INVALID",
            str(registry.get("candidate_count")),
        )
        require(
            registry.get("candidate_count") == len(candidates),
            "CANDIDATE_COUNT_MISMATCH",
            str(registry.get("candidate_count")),
        )
        recomputed = []
        for candidate in candidates:
            require(
                candidate.get("slot_id") in slot_ids,
                "CANDIDATE_SLOT_UNKNOWN",
                str(candidate.get("candidate_id")),
            )
            require(
                candidate.get("grammar_class")
                in set(grammar.get("class_ids", [])),
                "CANDIDATE_CLASS_UNKNOWN",
                str(candidate.get("candidate_id")),
            )
            expression = candidate.get("expression", {})
            require(
                expression.get("ast_depth", 99)
                <= grammar_budget["maximum_ast_depth"]
                and expression.get("ast_nodes", 99)
                <= grammar_budget["maximum_ast_nodes"],
                "CANDIDATE_BUDGET_EXCEEDED",
                str(candidate.get("candidate_id")),
            )
            require(
                not (
                    set(expression.get("registered_terminals", []))
                    & forbidden_terminals
                ),
                "CANDIDATE_FORBIDDEN_TERMINAL",
                str(candidate.get("candidate_id")),
            )
            claims = candidate.get("weight_claims", {})
            score = 0
            for claim_key, weight_key in (
                ("exact_global_certificate", "exact_global_certificate"),
                ("independent_recomputation", "independent_recomputation"),
                ("physicalization_complete", "physicalization_complete"),
                (
                    "baseline_freedom_counterexample",
                    "baseline_freedom_counterexample",
                ),
                (
                    "all_registered_completions_covered",
                    "all_registered_completions_covered",
                ),
                (
                    "all_continuous_parameters_covered",
                    "all_continuous_parameters_covered",
                ),
                ("conditional_branch", "conditional_branch_penalty"),
                ("open_physical_map", "open_physical_map_penalty"),
            ):
                if claims.get(claim_key):
                    score += weights[weight_key]
            score += (
                weights["expression_complexity_unit_penalty"]
                * expression.get("complexity_units", 0)
            )
            require(
                candidate.get("score") == score,
                "CANDIDATE_SCORE_DRIFT",
                f"{candidate.get('candidate_id')}: {candidate.get('score')} vs {score}",
            )
            recomputed.append(candidate)
        expected_order = sorted(
            recomputed,
            key=lambda candidate: (
                -candidate["score"],
                candidate["expression"]["complexity_units"],
                candidate["expression"]["form"],
                candidate["candidate_id"],
            ),
        )
        require(
            [candidate["candidate_id"] for candidate in expected_order]
            == [candidate["candidate_id"] for candidate in candidates]
            and all(
                candidate.get("rank") == position
                for position, candidate in enumerate(candidates, start=1)
            ),
            "CANDIDATE_ORDER_DRIFT",
            "committed order differs from frozen ranking recomputation",
        )
        from fractions import Fraction

        def q5_mul(x, y):
            return (x[0] * y[0] + 5 * x[1] * y[1], x[0] * y[1] + x[1] * y[0])

        def q5_div(x, y):
            norm = y[0] * y[0] - 5 * y[1] * y[1]
            numerator = q5_mul(x, (y[0], -y[1]))
            return (numerator[0] / norm, numerator[1] / norm)

        ratio = q5_div((Fraction(5), Fraction(1)), (Fraction(5), Fraction(-1)))
        phi = (Fraction(1, 2), Fraction(1, 2))
        require(
            ratio == q5_mul(phi, phi)
            and ratio == (Fraction(3, 2), Fraction(1, 2)),
            "GALOIS_RATIO_DRIFT",
            str(ratio),
        )

        ports = (
            ["north"]
            + [f"upper_{i}" for i in range(5)]
            + [f"lower_{i}" for i in range(5)]
            + ["south"]
        )
        edge_list = []
        for i in range(5):
            edge_list.extend(
                [
                    ("north", f"upper_{i}"),
                    (f"upper_{i}", f"upper_{(i + 1) % 5}"),
                    ("south", f"lower_{i}"),
                    (f"lower_{i}", f"lower_{(i + 1) % 5}"),
                    (f"upper_{i}", f"lower_{i}"),
                    (f"upper_{i}", f"lower_{(i + 1) % 5}"),
                ]
            )
        position = {port: index for index, port in enumerate(ports)}
        size = len(ports)
        adjacency = [[0] * size for _ in range(size)]
        for left, right in edge_list:
            adjacency[position[left]][position[right]] = 1
            adjacency[position[right]][position[left]] = 1

        def matmul(x, y):
            return [
                [
                    sum(x[row][k] * y[k][column] for k in range(size))
                    for column in range(size)
                ]
                for row in range(size)
            ]

        def shift(matrix, value):
            return [
                [
                    matrix[row][column] + (value if row == column else 0)
                    for column in range(size)
                ]
                for row in range(size)
            ]

        a_squared = matmul(adjacency, adjacency)
        product = matmul(
            matmul(shift(adjacency, -5), shift(a_squared, -5)),
            shift(adjacency, 1),
        )
        require(
            all(value == 0 for row in product for value in row)
            and sum(adjacency[i][i] for i in range(size)) == 0
            and sum(a_squared[i][i] for i in range(size)) == 60,
            "ADJACENCY_BAND_DRIFT",
            "annihilating polynomial or traces failed",
        )

        phi_pair = (Fraction(1, 2), Fraction(1, 2))
        psi_pair = (Fraction(1, 2), Fraction(-1, 2))
        class_sizes = [1, 15, 20, 12, 12]
        chi3 = [(Fraction(3), Fraction(0)), (Fraction(-1), Fraction(0)), (Fraction(0), Fraction(0)), phi_pair, psi_pair]
        chi3p = [(Fraction(3), Fraction(0)), (Fraction(-1), Fraction(0)), (Fraction(0), Fraction(0)), psi_pair, phi_pair]
        total = (Fraction(0), Fraction(0))
        for weight, x, y in zip(class_sizes, chi3, chi3p):
            term = q5_mul((Fraction(weight), Fraction(0)), q5_mul(x, y))
            total = (total[0] + term[0], total[1] + term[1])
        cross = q5_div(total, (Fraction(60), Fraction(0)))
        require(
            cross == (Fraction(0), Fraction(0)),
            "A5_PAIRING_DRIFT",
            str(cross),
        )

        by_id = {
            candidate["candidate_id"]: candidate for candidate in candidates
        }
        ratio_row = by_id.get("wz-galois-cost-ratio-phi-squared")
        require(
            ratio_row is not None
            and ratio_row["relation_certificate"]["ratio"] == "3/2+1/2*sqrt5",
            "CANDIDATE_RATIO_DRIFT",
            str(ratio_row and ratio_row["relation_certificate"]),
        )
        pairing_row = by_id.get("a5-cross-band-pairing-zero")
        require(
            pairing_row is not None
            and pairing_row["relation_certificate"]["invariant_multiplicity"]
            == "0+0*sqrt5",
            "CANDIDATE_PAIRING_DRIFT",
            str(pairing_row and pairing_row["relation_certificate"]),
        )

        registry_keys = all_keys(registry)
        leaked = sorted(
            registry_keys & set(policy.get("forbidden_document_keys", []))
        )
        require(not leaked, "CANDIDATE_TARGET_KEY", str(leaked))

    forbidden_output_names = set(policy.get("forbidden_output_names", []))
    package_files = {
        path.relative_to(package_root).as_posix()
        for path in package_root.rglob("*")
        if path.is_file()
        and "__pycache__" not in path.parts
        and ".pytest_cache" not in path.parts
    }
    require(
        not {Path(path).name for path in package_files} & forbidden_output_names,
        "CANDIDATE_OR_SCORE_OUTPUT_PRESENT",
        str(sorted(package_files)),
    )
    require(
        not any(
            part in {"candidates", "comparison_payloads", "scores"}
            for path in package_files
            for part in PurePosixPath(path).parts
        ),
        "FORBIDDEN_OUTPUT_DIRECTORY",
        "candidate, comparison, or score directory exists",
    )


def main() -> None:
    script = Path(__file__).resolve()
    default_package = script.parents[1]
    default_repo = default_package.parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package-root", type=Path, default=default_package)
    parser.add_argument("--repo-root", type=Path, default=default_repo)
    args = parser.parse_args()
    verify(args.package_root.resolve(), args.repo_root.resolve())
    print("GENERATION_LOCK_VALID")


if __name__ == "__main__":
    main()
