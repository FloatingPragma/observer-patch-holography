#!/usr/bin/env python3
"""Build the content-addressed issue-647 pre-generation campaign freeze."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

STATUS = (
    "REGISTRY_FINALIZED__"
    "GENERATOR_DISABLED_PENDING_ENABLEMENT_REVIEW"
)
SCHEMA = "oph.invariant_mining.pregeneration_freeze.v1"


class FreezeError(ValueError):
    """The pre-generation freeze cannot be constructed safely."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise FreezeError(message)


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


def pin(path: Path, relative: str) -> dict[str, Any]:
    require(path.is_file() and not path.is_symlink(), f"freeze input missing or symbolic: {relative}")
    payload = path.read_bytes()
    require(payload, f"freeze input is empty: {relative}")
    return {
        "path": relative,
        "bytes": len(payload),
        "sha256": sha256_bytes(payload),
    }


def verify_pin(pin_row: dict[str, Any], repo_root: Path) -> None:
    relative = pin_row.get("path")
    require(isinstance(relative, str) and relative and "\\" not in relative, "invalid source-projection path")
    path = repo_root / relative
    require(path.is_file() and not path.is_symlink(), f"source-projection input missing: {relative}")
    payload = path.read_bytes()
    require(pin_row.get("bytes") == len(payload), f"source-projection byte-count drift: {relative}")
    require(pin_row.get("sha256") == sha256_bytes(payload), f"source-projection hash drift: {relative}")


def build_freeze(package_root: Path, repo_root: Path) -> dict[str, Any]:
    policy_path = package_root / "policy" / "pregeneration_policy.json"
    projection_path = package_root / "outputs" / "source_projection.json"
    policy = load_json(policy_path)
    projection = load_json(projection_path)

    require(policy.get("status") == STATUS, "policy status drift")
    require(projection.get("schema") == "oph.invariant_mining.source_projection.v1", "source-projection schema drift")
    require(projection.get("status") == STATUS, "source-projection status drift")
    require(
        projection.get("registry_finalization_complete") is True,
        "source projection lost the registry finalization flag",
    )
    for key in (
        "candidate_generator_enabled",
        "candidate_evaluator_enabled",
    ):
        require(projection.get(key) is False, f"source projection changed execution boundary: {key}")
    require(projection.get("candidate_count") == 0, "source projection contains candidates")

    all_pins = projection.get("control_documents", []) + projection.get("source_artifacts", [])
    require(all_pins and all(isinstance(row, dict) for row in all_pins), "source projection has no typed pins")
    pinned_paths = [row.get("path") for row in all_pins]
    require(len(pinned_paths) == len(set(pinned_paths)), "source projection contains duplicate paths")
    for row in all_pins:
        verify_pin(row, repo_root)

    control_paths = policy.get("control_artifacts")
    projected_controls = [row["path"] for row in projection["control_documents"]]
    require(projected_controls == control_paths, "source projection does not exactly bind control-artifact policy")

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
    bindings = [pin(repo_root / relative, relative) for relative in document_relatives]
    projection_relative = "code/invariant_mining/outputs/source_projection.json"
    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "issue": 647,
        "status": STATUS,
        "freeze_id": "",
        "source_projection": pin(projection_path, projection_relative),
        "document_bindings": bindings,
        "direct_n_fallback_contract": {
            "slot_id": policy["direct_n_contract"]["slot_id"],
            "fallback_eligibility": policy["direct_n_contract"][
                "required_fallback_eligibility"
            ],
            "execution_status": policy["direct_n_contract"][
                "required_execution_status"
            ],
        },
        "execution_boundary": {
            "registry_finalization_complete": True,
            "candidate_generator_enabled": False,
            "candidate_evaluator_enabled": False,
            "public_data_access_enabled": False,
            "target_payloads_registered": False,
            "candidate_count": 0,
            "comparison_scoring_enabled": False,
        },
        "promotion_rule": policy["promotion_rule"],
    }
    digest_payload = dict(payload)
    digest_payload.pop("freeze_id")
    payload["freeze_id"] = sha256_bytes(canonical_bytes(digest_payload))
    return payload


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
    output = (args.output or package_root / "outputs" / "pregeneration_freeze.json").resolve()
    payload = canonical_bytes(build_freeze(package_root, repo_root))
    if args.check:
        require(output.is_file(), f"pre-generation freeze is missing: {output}")
        require(output.read_bytes() == payload, "pre-generation freeze differs from exact rebuild")
        print("PREGENERATION_FREEZE_VALID")
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes(payload)
    print(output)


if __name__ == "__main__":
    main()
