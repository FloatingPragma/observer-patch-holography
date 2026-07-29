#!/usr/bin/env python3
"""Build the exact source projection for the issue-647 pre-generation lock."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
from typing import Any

STATUS = (
    "PREGENERATION_POLICY_LOCKED__"
    "GENERATOR_DISABLED_PENDING_REGISTRY_FINALIZATION"
)
SCHEMA = "oph.invariant_mining.source_projection.v1"


class ProjectionError(ValueError):
    """The pre-generation source projection is incomplete or unsafe."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ProjectionError(message)


def canonical_bytes(value: Any) -> bytes:
    return (
        json.dumps(value, sort_keys=True, indent=2, ensure_ascii=False) + "\n"
    ).encode("utf-8")


def sha256_bytes(payload: bytes) -> str:
    return "sha256:" + hashlib.sha256(payload).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    require(isinstance(value, dict), f"JSON document must be an object: {path}")
    return value


def require_repo_path(value: str) -> str:
    require(value != "" and "\\" not in value, f"invalid repository path: {value!r}")
    parsed = PurePosixPath(value)
    require(
        not parsed.is_absolute()
        and parsed.as_posix() == value
        and "." not in parsed.parts
        and ".." not in parsed.parts,
        f"path is not canonical repository-relative POSIX syntax: {value}",
    )
    return value


def resolve_file(repo_root: Path, relative: str) -> Path:
    require_repo_path(relative)
    unresolved = repo_root / relative
    require(not unresolved.is_symlink(), f"symbolic-link pin refused: {relative}")
    resolved_root = repo_root.resolve()
    resolved = unresolved.resolve()
    require(
        resolved_root in resolved.parents,
        f"repository path escapes root: {relative}",
    )
    require(resolved.is_file(), f"registered file is missing: {relative}")
    return resolved


def pin_file(repo_root: Path, relative: str) -> dict[str, Any]:
    payload = resolve_file(repo_root, relative).read_bytes()
    require(len(payload) > 0, f"registered file is empty: {relative}")
    return {
        "path": relative,
        "bytes": len(payload),
        "sha256": sha256_bytes(payload),
    }


def recursive_keys(value: Any) -> set[str]:
    keys: set[str] = set()
    if isinstance(value, dict):
        for key, child in value.items():
            keys.add(str(key))
            keys.update(recursive_keys(child))
    elif isinstance(value, list):
        for child in value:
            keys.update(recursive_keys(child))
    return keys


def exact_ids(rows: Any, key: str, *, label: str) -> list[str]:
    require(isinstance(rows, list), f"{label} must be a list")
    values = [row.get(key) for row in rows if isinstance(row, dict)]
    require(
        len(values) == len(rows)
        and all(isinstance(value, str) and value for value in values),
        f"{label} has a missing or invalid {key}",
    )
    require(len(set(values)) == len(values), f"{label} contains duplicate {key} values")
    return sorted(values)


def validate_documents(
    *,
    policy: dict[str, Any],
    source_registry: dict[str, Any],
    producer_registry: dict[str, Any],
    grammar: dict[str, Any],
    nuisances: dict[str, Any],
    ranking: dict[str, Any],
) -> set[str]:
    require(policy.get("schema") == "oph.invariant_mining.pregeneration_policy.v1", "policy schema drift")
    require(policy.get("status") == STATUS, "policy status drift")
    require(policy.get("registry_finalization_complete") is False, "registry may not be finalized in this lock")
    require(policy.get("candidate_generator_enabled") is False, "candidate generator must remain disabled")
    require(policy.get("candidate_evaluator_enabled") is False, "candidate evaluator must remain disabled")
    require(policy.get("public_data_access_enabled") is False, "public-data access must remain disabled")
    require(policy.get("target_payloads_registered") is False, "target payloads may not be registered")
    require(policy.get("candidate_count") == 0, "pre-generation candidate count must be zero")

    require(
        source_registry.get("schema")
        == "oph.invariant_mining.source_feature_registry.v1",
        "source-feature registry schema drift",
    )
    require(
        source_registry.get("completeness", {}).get("status")
        == "BOUNDED_SEED_REGISTRY__NOT_FINAL",
        "source-feature registry must remain explicitly non-final",
    )
    feature_rows = source_registry.get("features")
    feature_ids = exact_ids(feature_rows, "feature_id", label="source features")
    require(feature_ids == sorted(policy["required_feature_ids"]), "required source-feature set drift")

    require(
        producer_registry.get("schema")
        == "oph.invariant_mining.producer_slots.v1",
        "producer-slot registry schema drift",
    )
    require(
        producer_registry.get("execution_state")
        == "ALL_PRODUCERS_DISABLED_PENDING_REGISTRY_FINALIZATION",
        "producer execution-state drift",
    )
    slot_rows = producer_registry.get("slots")
    slot_ids = exact_ids(slot_rows, "slot_id", label="producer slots")
    require(slot_ids == sorted(policy["required_producer_slot_ids"]), "required producer-slot set drift")

    require(
        grammar.get("schema") == "oph.invariant_mining.observable_grammar.v1",
        "observable grammar schema drift",
    )
    require(
        grammar.get("freeze_status") == "FROZEN_BEFORE_CANDIDATE_GENERATION",
        "observable grammar is not frozen",
    )
    grammar_ids = grammar.get("class_ids")
    require(
        isinstance(grammar_ids, list)
        and sorted(grammar_ids) == sorted(policy["required_grammar_class_ids"])
        and len(set(grammar_ids)) == len(grammar_ids),
        "required observable grammar class set drift",
    )
    budget = grammar.get("complexity_budget", {})
    require(
        all(
            isinstance(budget.get(key), int) and budget[key] > 0
            for key in (
                "maximum_ast_depth",
                "maximum_ast_nodes",
                "maximum_distinct_registered_terminals",
                "maximum_generated_candidates",
            )
        ),
        "observable grammar complexity budget must be positive and integral",
    )

    require(
        nuisances.get("schema") == "oph.invariant_mining.nuisance_registry.v1",
        "nuisance registry schema drift",
    )
    require(
        nuisances.get("freeze_status") == "FROZEN_BEFORE_CANDIDATE_GENERATION",
        "nuisance registry is not frozen",
    )
    taxonomy = nuisances.get("taxonomy")
    nuisance_type_ids = exact_ids(taxonomy, "type_id", label="nuisance taxonomy")
    require(
        nuisance_type_ids == sorted(policy["required_nuisance_type_ids"]),
        "required nuisance-type set drift",
    )
    unresolved = nuisances.get("unresolved_directions")
    nuisance_ids = exact_ids(unresolved, "nuisance_id", label="unresolved nuisances")
    require(nuisance_ids == sorted(policy["required_nuisance_ids"]), "required unresolved-nuisance set drift")
    taxonomy_by_id = {row["type_id"]: row for row in taxonomy}
    require(
        taxonomy_by_id["coordinate_equivalence_nuisance"]["quotient_allowed"] is True,
        "coordinate equivalence must be the sole quotientable nuisance type",
    )
    require(
        all(
            row["quotient_allowed"] is False
            for type_id, row in taxonomy_by_id.items()
            if type_id != "coordinate_equivalence_nuisance"
        ),
        "a physical nuisance type was incorrectly made quotientable",
    )

    require(
        ranking.get("schema") == "oph.invariant_mining.ranking_policy.v1",
        "ranking-policy schema drift",
    )
    require(
        ranking.get("freeze_status") == "FROZEN_BEFORE_CANDIDATE_GENERATION",
        "ranking policy is not frozen",
    )
    require(
        ranking.get("score_rule") == "integer_sum_of_applicable_frozen_weights",
        "ranking score rule drift",
    )
    weights = ranking.get("weights")
    require(
        isinstance(weights, dict)
        and weights
        and all(isinstance(value, int) for value in weights.values()),
        "ranking weights must be a nonempty integer map",
    )
    require(
        ranking.get("tie_breaker")
        == [
            "lower_expression_complexity",
            "lexicographically_smaller_canonical_expression",
            "lexicographically_smaller_candidate_id",
        ],
        "ranking tie-breaker drift",
    )

    feature_set = set(feature_ids)
    grammar_set = set(grammar_ids)
    for row in feature_rows:
        require(row.get("candidate_generation_eligible") is False, f"feature became generation-eligible: {row['feature_id']}")
        require(
            row.get("target_ancestry") == "SOURCE_ONLY_NO_PUBLIC_COMPARISON_INPUT",
            f"source ancestry drift: {row['feature_id']}",
        )
        require(
            isinstance(row.get("artifacts"), list) and row["artifacts"],
            f"source feature has no artifacts: {row['feature_id']}",
        )

    for row in slot_rows:
        require(row.get("execution_status") == "REGISTERED_DISABLED", f"producer slot became executable: {row['slot_id']}")
        required_features = row.get("required_feature_ids")
        allowed_classes = row.get("allowed_grammar_class_ids")
        require(
            isinstance(required_features, list)
            and required_features
            and set(required_features) <= feature_set,
            f"producer slot has unknown source features: {row['slot_id']}",
        )
        require(
            isinstance(allowed_classes, list)
            and allowed_classes
            and set(allowed_classes) <= grammar_set,
            f"producer slot has unknown grammar classes: {row['slot_id']}",
        )

    direct_contract = policy["direct_n_contract"]
    direct_rows = [row for row in slot_rows if row["slot_id"] == direct_contract["slot_id"]]
    require(len(direct_rows) == 1, "direct-N slot must occur exactly once")
    direct_row = direct_rows[0]
    require(
        direct_row.get("fallback_eligibility")
        == direct_contract["required_fallback_eligibility"],
        "direct-N fallback re-entry is forbidden",
    )
    require(
        direct_row.get("execution_status")
        == direct_contract["required_execution_status"],
        "direct-N producer must remain disabled",
    )

    forbidden_keys = set(policy.get("forbidden_document_keys", []))
    for name, document in (
        ("source-feature registry", source_registry),
        ("producer-slot registry", producer_registry),
        ("observable grammar", grammar),
        ("nuisance registry", nuisances),
        ("ranking policy", ranking),
    ):
        found = sorted(recursive_keys(document) & forbidden_keys)
        require(not found, f"{name} contains forbidden target keys: {found}")

    source_paths: set[str] = set()
    forbidden_prefixes = tuple(policy.get("forbidden_repository_prefixes", []))
    for row in feature_rows:
        for artifact in row["artifacts"]:
            require(isinstance(artifact, dict), f"invalid artifact row in {row['feature_id']}")
            require(set(artifact) == {"path", "kind", "role"}, f"artifact shape drift in {row['feature_id']}")
            relative = require_repo_path(str(artifact["path"]))
            require(
                not relative.startswith(forbidden_prefixes),
                f"registered source enters a forbidden comparison surface: {relative}",
            )
            require(relative not in source_paths, f"duplicate registered source artifact: {relative}")
            source_paths.add(relative)

    return source_paths


def build_projection(package_root: Path, repo_root: Path) -> dict[str, Any]:
    policy_path = package_root / "policy" / "pregeneration_policy.json"
    source_path = package_root / "data" / "source_feature_registry.json"
    producer_path = package_root / "data" / "producer_slots.json"
    grammar_path = package_root / "data" / "observable_grammar.json"
    nuisance_path = package_root / "data" / "nuisance_registry.json"
    ranking_path = package_root / "data" / "ranking_policy.json"

    policy = load_json(policy_path)
    source_registry = load_json(source_path)
    producer_registry = load_json(producer_path)
    grammar = load_json(grammar_path)
    nuisances = load_json(nuisance_path)
    ranking = load_json(ranking_path)
    source_paths = validate_documents(
        policy=policy,
        source_registry=source_registry,
        producer_registry=producer_registry,
        grammar=grammar,
        nuisances=nuisances,
        ranking=ranking,
    )

    control_paths = policy.get("control_artifacts")
    require(
        isinstance(control_paths, list)
        and control_paths == sorted(control_paths)
        and len(control_paths) == len(set(control_paths)),
        "control artifacts must be a sorted unique list",
    )
    required_controls = {
        "code/invariant_mining/data/source_feature_registry.json",
        "code/invariant_mining/data/producer_slots.json",
        "code/invariant_mining/data/observable_grammar.json",
        "code/invariant_mining/data/nuisance_registry.json",
        "code/invariant_mining/data/ranking_policy.json",
        "code/invariant_mining/policy/pregeneration_policy.json",
        "code/invariant_mining/tools/build_source_projection.py",
        "code/invariant_mining/tools/build_pregeneration_freeze.py",
        "code/invariant_mining/tools/verify_pregeneration_freeze_independent.py",
    }
    require(required_controls <= set(control_paths), "mandatory control artifact omitted")

    control_pins = [pin_file(repo_root, relative) for relative in control_paths]
    source_pins = [pin_file(repo_root, relative) for relative in sorted(source_paths)]
    path_digest = sha256_bytes(canonical_bytes(sorted(source_paths)))
    return {
        "schema": SCHEMA,
        "issue": 647,
        "status": STATUS,
        "projection_id": "oph-invariant-source-projection-v1",
        "registry_finalization_complete": False,
        "candidate_generator_enabled": False,
        "candidate_evaluator_enabled": False,
        "candidate_count": 0,
        "control_documents": control_pins,
        "source_artifacts": source_pins,
        "registered_source_paths_sha256": path_digest,
    }


def main() -> None:
    script = Path(__file__).resolve()
    default_package = script.parents[1]
    default_repo = default_package.parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package-root", type=Path, default=default_package)
    parser.add_argument("--repo-root", type=Path, default=default_repo)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    package_root = args.package_root.resolve()
    repo_root = args.repo_root.resolve()
    output = (args.output or package_root / "outputs" / "source_projection.json").resolve()
    payload = canonical_bytes(build_projection(package_root, repo_root))
    if args.check:
        require(output.is_file(), f"source projection is missing: {output}")
        require(output.read_bytes() == payload, "source projection differs from exact rebuild")
        print("SOURCE_PROJECTION_VALID")
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes(payload)
    print(output)


if __name__ == "__main__":
    main()
