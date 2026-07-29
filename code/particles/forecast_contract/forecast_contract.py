#!/usr/bin/env python3
"""Issue-639 forecast-governance packet.

This module validates the draft contract format, records the provisional
candidate inventory, checks declared vocabulary for forbidden inputs, and
verifies file pins required by a future freeze.  It deliberately does not
claim comparison quarantine, execute a forecast, score a verdict, or enforce
single-use unsealing.  Those controls require an eligible source-visible
candidate, independent custody, two executable implementations, and a durable
scoring record.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path
from pathlib import PurePosixPath
from typing import Any

from jsonschema import Draft202012Validator

PACKET_DIR = Path(__file__).resolve().parent
REPO_ROOT = PACKET_DIR.parents[2]
SCHEMA_PATH = PACKET_DIR / "data" / "prediction_contract_schema_v2.json"
INVENTORY_PATH = PACKET_DIR / "data" / "candidate_inventory_v2.json"
POLICY_PATH = PACKET_DIR / "data" / "forecast_freeze_policy_v1.json"
OUT_PATH = PACKET_DIR / "outputs" / "forecast_contract_state.json"
FZ_REGISTER_PATH = REPO_ROOT / "claims" / "frozen_prediction_register.json"

STATE_SCHEMA = "oph.forecast_contract_state.v2"
CONTRACT_SCHEMA = "oph.prediction_contract.v2"
DISCRIMINATOR_SCHEMA = "oph.forecast_discriminator_receipt.v1"
INVENTORY_SCHEMA = "oph.forecast_candidate_inventory.v2"
OMISSION_AUDIT_SCHEMA = "oph.forecast_inventory_omission_audit.v1"
ANCESTRY_EVIDENCE_SCHEMA = "oph.forecast_target_ancestry_evidence.v1"
CUSTODY_SCHEMA = "oph.forecast_comparison_custody.v1"
POLICY_SCHEMA = "oph.forecast_freeze_policy.v1"
FROZEN_SELECTION_ALGORITHM = (
    "candidate_class_then_priority_then_candidate_id_first_eligible_v1"
)
CONTRACT_DIGEST_SCOPE = (
    "canonical_contract_without_freeze_canonical_payload_sha256"
)
INVENTORY_DIGEST_SCOPE = (
    "canonical_inventory_without_completeness_"
    "adversarial_omission_audit"
)
REQUIRED_PROMOTION_CONDITIONS = {
    "executable generator replay receipt",
    "checker implementation independence receipt",
    "comparison quarantine and access-control receipt",
    "durable single-use independent scoring receipt",
}
ISSUE = 639


class ContractError(ValueError):
    def __init__(self, code: str, message: str) -> None:
        super().__init__(f"{code}: {message}")
        self.code = code


def require(condition: bool, code: str, message: str) -> None:
    if not condition:
        raise ContractError(code, message)


def sha256_bytes(data: bytes) -> str:
    return "sha256:" + hashlib.sha256(data).hexdigest()


def _is_tagged_sha256(value: Any) -> bool:
    if not isinstance(value, str) or not value.startswith("sha256:"):
        return False
    digest = value.removeprefix("sha256:")
    return len(digest) == 64 and all(
        character in "0123456789abcdef" for character in digest
    )


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _require_repo_relative_path(relative: str) -> PurePosixPath:
    require(
        relative != ""
        and "\\" not in relative
        and not relative.startswith("/"),
        "PIN_PATH_INVALID",
        f"pin path must use nonempty repository-relative POSIX syntax: {relative}",
    )
    path = PurePosixPath(relative)
    require(
        path.as_posix() == relative
        and "." not in path.parts
        and ".." not in path.parts,
        "PIN_PATH_INVALID",
        f"pin path is not canonical repository-relative POSIX syntax: {relative}",
    )
    return path


def _resolve_repo_file(repo_root: Path, relative: str) -> Path:
    _require_repo_relative_path(relative)
    root = repo_root.resolve()
    unresolved = root / relative
    require(
        not unresolved.is_symlink(),
        "PIN_SYMLINK_REFUSED",
        f"pinned file may not be a symbolic link: {relative}",
    )
    candidate = unresolved.resolve()
    require(
        candidate == root or root in candidate.parents,
        "PIN_PATH_ESCAPE",
        f"pin path escapes the repository: {relative}",
    )
    require(
        candidate.is_file(),
        "PIN_FILE_MISSING",
        f"pinned file does not exist: {relative}",
    )
    return candidate


def _run_git(
    repo_root: Path,
    arguments: list[str],
    *,
    code: str,
    message: str,
) -> subprocess.CompletedProcess[bytes]:
    result = subprocess.run(
        ["git", "-C", str(repo_root), *arguments],
        check=False,
        capture_output=True,
    )
    require(result.returncode == 0, code, message)
    return result


def verify_repository_commit(repo_root: Path, commit_sha: str) -> None:
    """Require a real commit in the current repository ancestry."""

    resolved = _run_git(
        repo_root,
        ["rev-parse", "--verify", f"{commit_sha}^{{commit}}"],
        code="FREEZE_COMMIT_NOT_FOUND",
        message=f"freeze repository state is not a real commit: {commit_sha}",
    ).stdout.decode("ascii").strip()
    require(
        resolved == commit_sha,
        "FREEZE_COMMIT_RESOLUTION_DRIFT",
        f"freeze commit resolves to {resolved}, not {commit_sha}",
    )
    _run_git(
        repo_root,
        ["merge-base", "--is-ancestor", commit_sha, "HEAD"],
        code="FREEZE_COMMIT_NOT_ANCESTOR",
        message="freeze repository state is not an ancestor of repository HEAD",
    )


def verify_file_pin(
    pin: dict[str, Any],
    repo_root: Path,
    *,
    path_key: str,
) -> None:
    """Resolve a repository-relative path and verify its exact byte hash."""

    path = _resolve_repo_file(repo_root, str(pin[path_key]))
    observed = sha256_bytes(path.read_bytes())
    require(
        observed == pin["sha256"],
        "PIN_HASH_MISMATCH",
        f"{pin[path_key]}: expected {pin['sha256']}, observed {observed}",
    )


def _require_artifact_pin_shape(
    pin: Any,
    *,
    code: str,
    message: str,
) -> dict[str, Any]:
    require(
        isinstance(pin, dict)
        and set(pin) == {"artifact", "sha256", "role"}
        and isinstance(pin.get("artifact"), str)
        and bool(pin["artifact"])
        and isinstance(pin.get("role"), str)
        and bool(pin["role"])
        and _is_tagged_sha256(pin.get("sha256")),
        code,
        message,
    )
    return pin


def verify_file_pin_at_commit(
    pin: dict[str, Any],
    repo_root: Path,
    commit_sha: str,
    *,
    path_key: str,
) -> None:
    """Verify that a pin names regular bytes in the frozen repository state."""

    relative = str(pin[path_key])
    _require_repo_relative_path(relative)
    tree_row = _run_git(
        repo_root,
        ["ls-tree", commit_sha, "--", relative],
        code="COMMIT_PIN_MISSING",
        message=f"pinned file is absent from freeze commit: {relative}",
    ).stdout.decode("utf-8").strip()
    require(
        bool(tree_row),
        "COMMIT_PIN_MISSING",
        f"pinned file is absent from freeze commit: {relative}",
    )
    mode = tree_row.split(maxsplit=1)[0]
    require(
        mode in {"100644", "100755"},
        "COMMIT_PIN_NOT_REGULAR",
        f"pinned commit object is not a regular file: {relative}",
    )
    data = _run_git(
        repo_root,
        ["show", f"{commit_sha}:{relative}"],
        code="COMMIT_PIN_MISSING",
        message=f"could not read pinned bytes from freeze commit: {relative}",
    ).stdout
    observed = sha256_bytes(data)
    require(
        observed == pin["sha256"],
        "COMMIT_PIN_HASH_MISMATCH",
        f"{relative}: freeze commit has {observed}, pin declares {pin['sha256']}",
    )


def validate_contract(contract: dict[str, Any]) -> None:
    """Apply the complete JSON Schema and cross-field draft checks."""

    schema = load_json(SCHEMA_PATH)
    errors = sorted(
        Draft202012Validator(schema).iter_errors(contract),
        key=lambda error: tuple(str(part) for part in error.absolute_path),
    )
    if errors:
        first = errors[0]
        location = "/".join(str(part) for part in first.absolute_path) or "/"
        raise ContractError(
            "CONTRACT_SCHEMA_INVALID",
            f"{location}: {first.message}",
        )

    require(
        contract["generator"]["module"]
        != contract["independent_checker"]["module"],
        "CHECKER_NOT_INDEPENDENT",
        "generator and checker module paths must differ",
    )
    require(
        contract["generator"]["sha256"]
        != contract["independent_checker"]["sha256"],
        "CHECKER_NOT_INDEPENDENT",
        "generator and checker byte hashes must differ",
    )
    require(
        set(contract["promotion_conditions"])
        == REQUIRED_PROMOTION_CONDITIONS,
        "PROMOTION_CONDITIONS_INCOMPLETE",
        "contract must retain every unimplemented scoring-readiness condition",
    )


def canonical_contract_payload(contract: dict[str, Any]) -> dict[str, Any]:
    """Return the non-self-referential payload covered by the contract digest.

    The repository-state commit is part of this payload.  The digest field
    itself is omitted.  The contract file need not claim to live in the commit
    it names; that commit freezes the pinned producer state, while an external
    record freezes this canonical contract payload.
    """

    payload = json.loads(json.dumps(contract))
    freeze = payload.get("freeze")
    require(
        isinstance(freeze, dict),
        "CONTRACT_FREEZE_BLOCK_INVALID",
        "contract has no freeze object",
    )
    freeze.pop("canonical_payload_sha256", None)
    return payload


def canonical_contract_payload_sha256(contract: dict[str, Any]) -> str:
    return sha256_bytes(
        canonical_json(canonical_contract_payload(contract)).encode("utf-8")
    )


def canonical_inventory_audit_payload(
    inventory: dict[str, Any],
) -> dict[str, Any]:
    """Return the non-circular inventory payload reviewed by an omission audit."""

    payload = json.loads(json.dumps(inventory))
    completeness = payload.get("completeness")
    require(
        isinstance(completeness, dict),
        "INVENTORY_COMPLETENESS_INVALID",
        "candidate inventory has no completeness object",
    )
    completeness.pop("adversarial_omission_audit", None)
    return payload


def canonical_inventory_audit_payload_sha256(
    inventory: dict[str, Any],
) -> str:
    return sha256_bytes(
        canonical_json(canonical_inventory_audit_payload(inventory)).encode(
            "utf-8"
        )
    )


def canonical_candidate_row_sha256(candidate: dict[str, Any]) -> str:
    return sha256_bytes(canonical_json(candidate).encode("utf-8"))


def canonical_pin_list_sha256(pins: list[dict[str, Any]]) -> str:
    return sha256_bytes(canonical_json(pins).encode("utf-8"))


def check_target_ancestry(
    candidate: dict[str, Any], forbidden_classes: dict[str, Any]
) -> dict[str, Any]:
    """Scan all declared candidate text for forbidden vocabulary.

    Negation words do not suppress a hit.  The result is a bounded lexical
    warning over declared text.  It is not a transitive provenance proof.
    """

    entries: list[dict[str, str]] = []
    for index, entry in enumerate(candidate.get("allowed_ancestry", [])):
        entries.append({"field": f"allowed_ancestry/{index}", "text": str(entry)})
    for field in ("blocking_condition", "description", "public_knowledge_caveat"):
        if field in candidate:
            entries.append({"field": field, "text": str(candidate[field])})

    hits: list[dict[str, str]] = []
    for entry in entries:
        text = entry["text"].lower()
        for name, row in forbidden_classes.items():
            for fragment in row["match_fragments"]:
                if fragment in text:
                    hits.append(
                        {
                            "field": entry["field"],
                            "entry": entry["text"],
                            "forbidden_class": name,
                            "fragment": fragment,
                        }
                    )
    return {
        "check_type": "bounded_declared_vocabulary_warning",
        "semantic_input_closure_proved": False,
        "ancestry_hits": hits,
        "declared_vocabulary_hit_free": not hits,
    }


def seal_comparison(payload_bytes: bytes, storage_note: str) -> dict[str, Any]:
    """Return an integrity digest without claiming secrecy or quarantine."""

    return {
        "payload_sha256": sha256_bytes(payload_bytes),
        "byte_count": len(payload_bytes),
        "storage_note": storage_note,
        "seal_scope": "integrity_only",
    }


def verify_frozen_contract(
    contract: dict[str, Any],
    external_digest_input: str | None,
) -> None:
    """Check the embedded digest against a separately supplied value.

    The supplied value is validated input.  Passing it does not demonstrate
    durable external recording or custody.
    """

    embedded = contract["freeze"]["canonical_payload_sha256"]
    current = canonical_contract_payload_sha256(contract)
    require(
        current == embedded,
        "CONTRACT_PAYLOAD_DIGEST_MISMATCH",
        "the canonical contract payload differs from its embedded digest",
    )
    require(
        external_digest_input is not None,
        "EXTERNAL_FREEZE_DIGEST_INPUT_REQUIRED",
        "freeze validation requires a separately supplied contract digest",
    )
    require(
        current == external_digest_input,
        "EXTERNAL_FREEZE_DIGEST_INPUT_MISMATCH",
        "the canonical contract payload differs from the external digest input",
    )


def _eligible_candidate_ids(inventory: dict[str, Any]) -> list[str]:
    candidates = inventory["candidates"]
    return [
        candidate_id
        for _, _, candidate_id in sorted(
            (
                (
                    candidate["candidate_class"],
                    candidate["selection_priority"],
                    candidate_id,
                )
                for candidate_id, candidate in candidates.items()
                if candidate["eligibility"] == "ELIGIBLE_SOURCE_VISIBLE"
            )
        )
    ]


def _validate_anchored_policy(
    contract: dict[str, Any],
    *,
    repo_root: Path,
    commit_sha: str,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Load the immutable issue policy and the policy-required live register."""

    policy_pin = contract["freeze_policy"]
    verify_file_pin(policy_pin, repo_root, path_key="artifact")
    verify_file_pin_at_commit(
        policy_pin,
        repo_root,
        commit_sha,
        path_key="artifact",
    )
    canonical_policy_raw = POLICY_PATH.read_bytes()
    require(
        policy_pin["sha256"] == sha256_bytes(canonical_policy_raw),
        "FREEZE_POLICY_NOT_CANONICAL",
        "contract policy pin differs from the validator's issue-639 policy",
    )
    try:
        pinned_policy = load_json(
            _resolve_repo_file(repo_root, policy_pin["artifact"])
        )
        canonical_policy = json.loads(canonical_policy_raw.decode("utf-8"))
    except (json.JSONDecodeError, UnicodeDecodeError) as error:
        raise ContractError(
            "FREEZE_POLICY_INVALID",
            "freeze policy is not valid UTF-8 JSON",
        ) from error
    require(
        pinned_policy == canonical_policy
        and canonical_policy.get("schema") == POLICY_SCHEMA
        and canonical_policy.get("issue") == ISSUE,
        "FREEZE_POLICY_NOT_CANONICAL",
        "pinned freeze policy differs from the canonical issue-639 policy",
    )

    register_pin = contract["live_fz_register"]
    require(
        register_pin["artifact"]
        == canonical_policy["live_register_projection"]["artifact"],
        "LIVE_REGISTER_PATH_MISMATCH",
        "contract does not pin the policy-required live register path",
    )
    verify_file_pin(register_pin, repo_root, path_key="artifact")
    verify_file_pin_at_commit(
        register_pin,
        repo_root,
        commit_sha,
        path_key="artifact",
    )
    canonical_register_raw = FZ_REGISTER_PATH.read_bytes()
    require(
        register_pin["sha256"] == sha256_bytes(canonical_register_raw),
        "LIVE_REGISTER_NOT_CANONICAL",
        "contract live-register pin differs from the validator's registry",
    )
    try:
        live_register = load_json(
            _resolve_repo_file(repo_root, register_pin["artifact"])
        )
        canonical_live_register = json.loads(
            canonical_register_raw.decode("utf-8")
        )
    except (json.JSONDecodeError, UnicodeDecodeError) as error:
        raise ContractError(
            "LIVE_REGISTER_INVALID",
            "live frozen-prediction register is not valid UTF-8 JSON",
        ) from error
    require(
        live_register == canonical_live_register,
        "LIVE_REGISTER_NOT_CANONICAL",
        "pinned live register differs from the validator's canonical registry",
    )
    require(
        live_register.get("schema")
        == canonical_policy["live_register_projection"]["schema"],
        "LIVE_REGISTER_SCHEMA",
        "live register does not match the policy-required schema",
    )
    for collection in canonical_policy["live_register_projection"][
        "collections"
    ]:
        require(
            isinstance(live_register.get(collection), list),
            "LIVE_REGISTER_COLLECTION_MISSING",
            f"live register lacks policy-required collection: {collection}",
        )
    build_fz_crosswalk(live_register)
    return canonical_policy, live_register


def _validate_frozen_inventory(
    inventory: dict[str, Any],
    *,
    repo_root: Path,
    commit_sha: str,
    policy: dict[str, Any],
    live_register: dict[str, Any],
    policy_sha256: str,
) -> list[str]:
    require(
        inventory.get("schema") == INVENTORY_SCHEMA,
        "INVENTORY_SCHEMA_INVALID",
        "candidate inventory has the wrong schema",
    )
    require(
        inventory.get("issue") == ISSUE,
        "INVENTORY_ISSUE_MISMATCH",
        f"candidate inventory must be owned by issue {ISSUE}",
    )
    require(
        inventory.get("inventory_scope") == policy["inventory_scope"],
        "INVENTORY_SCOPE_POLICY_MISMATCH",
        "candidate inventory scope differs from the anchored issue policy",
    )
    require(
        inventory.get("forbidden_input_classes")
        == policy["forbidden_input_classes"],
        "FORBIDDEN_INPUT_POLICY_MISMATCH",
        "candidate inventory weakens or changes the anchored forbidden policy",
    )
    completeness = inventory.get("completeness", {})
    require(
        completeness.get("status") == "FROZEN_EXHAUSTIVE"
        and completeness.get("closure_criterion_frozen") is True,
        "INVENTORY_NOT_FROZEN_EXHAUSTIVE",
        "a contract freeze requires a frozen exhaustive inventory",
    )
    closure_criterion_id = completeness.get("closure_criterion_id")
    closure_criterion = completeness.get("closure_criterion")
    require(
        isinstance(closure_criterion_id, str)
        and bool(closure_criterion_id.strip())
        and isinstance(closure_criterion, str)
        and bool(closure_criterion.strip()),
        "INVENTORY_CLOSURE_CRITERION_MISSING",
        "frozen exhaustiveness requires a named, nonempty closure criterion",
    )
    require(
        closure_criterion_id
        == policy["closure_policy"]["closure_criterion_id"]
        and closure_criterion
        == policy["closure_policy"]["closure_criterion"],
        "INVENTORY_CLOSURE_POLICY_MISMATCH",
        "inventory closure criterion differs from the anchored issue policy",
    )
    audit_pin = completeness.get("adversarial_omission_audit")
    require(
        isinstance(audit_pin, dict)
        and set(audit_pin) == {"artifact", "sha256", "role"}
        and isinstance(audit_pin.get("artifact"), str)
        and isinstance(audit_pin.get("role"), str)
        and _is_tagged_sha256(audit_pin.get("sha256")),
        "INVENTORY_OMISSION_AUDIT_MISSING",
        "frozen exhaustiveness requires a hash-pinned omission audit",
    )
    verify_file_pin(audit_pin, repo_root, path_key="artifact")
    verify_file_pin_at_commit(
        audit_pin,
        repo_root,
        commit_sha,
        path_key="artifact",
    )
    try:
        omission_audit = load_json(
            _resolve_repo_file(repo_root, audit_pin["artifact"])
        )
    except (json.JSONDecodeError, UnicodeDecodeError) as error:
        raise ContractError(
            "INVENTORY_OMISSION_AUDIT_INVALID",
            "omission audit is not valid UTF-8 JSON",
        ) from error
    require(
        set(omission_audit)
        == {
            "schema",
            "issue",
            "inventory_payload_sha256",
            "inventory_digest_scope",
            "scope_id",
            "closure_criterion_id",
            "reviewed_candidate_ids",
            "reviewed_static_surface_ids",
            "reviewed_live_register_surface_ids",
            "freeze_policy_sha256",
            "verdict",
            "adversarial_review_complete",
            "reviewer_id",
            "reviewer_attestation",
        }
        and omission_audit.get("schema") == OMISSION_AUDIT_SCHEMA
        and omission_audit.get("issue") == ISSUE
        and omission_audit.get("inventory_payload_sha256")
        == canonical_inventory_audit_payload_sha256(inventory)
        and omission_audit.get("inventory_digest_scope")
        == INVENTORY_DIGEST_SCOPE
        and omission_audit.get("scope_id")
        == inventory.get("inventory_scope", {}).get("scope_id")
        and omission_audit.get("closure_criterion_id")
        == closure_criterion_id
        and omission_audit.get("reviewed_candidate_ids")
        == sorted(inventory.get("candidates", {}))
        and omission_audit.get("reviewed_static_surface_ids")
        == sorted(
            row["surface"]
            for row in inventory.get("known_surface_crosswalk", [])
        )
        and omission_audit.get("reviewed_live_register_surface_ids")
        == sorted(
            row["surface"] for row in build_fz_crosswalk(live_register)
        )
        and omission_audit.get("freeze_policy_sha256") == policy_sha256
        and omission_audit.get("verdict")
        == "NO_OMITTED_IN_SCOPE_CANDIDATE"
        and omission_audit.get("adversarial_review_complete") is True,
        "INVENTORY_OMISSION_AUDIT_NOT_POSITIVE",
        "omission audit does not certify the frozen scope and criterion",
    )
    require(
        isinstance(omission_audit.get("reviewer_id"), str)
        and bool(omission_audit["reviewer_id"].strip())
        and isinstance(omission_audit.get("reviewer_attestation"), str)
        and len(omission_audit["reviewer_attestation"].strip()) >= 16,
        "INVENTORY_OMISSION_AUDIT_ATTESTATION_MISSING",
        "omission audit lacks reviewer identity or attestation",
    )

    selection = inventory.get("selection_rule", {})
    require(
        selection.get("status") == "FROZEN_DETERMINISTIC"
        and selection.get("algorithm")
        == policy["selection_algorithm"]
        == FROZEN_SELECTION_ALGORITHM,
        "SELECTION_RULE_NOT_FROZEN",
        "candidate selection must use the frozen deterministic algorithm",
    )
    candidates = inventory.get("candidates")
    forbidden_classes = inventory.get("forbidden_input_classes")
    require(
        isinstance(forbidden_classes, dict) and bool(forbidden_classes),
        "FORBIDDEN_INPUT_INVENTORY_INVALID",
        "frozen inventory has no forbidden-input class map",
    )
    require(
        isinstance(candidates, dict) and bool(candidates),
        "CANDIDATE_INVENTORY_EMPTY",
        "frozen candidate inventory is empty",
    )
    require(
        set(candidates) == set(policy["required_candidates"]),
        "REQUIRED_CANDIDATE_SET_MISMATCH",
        "candidate inventory omits or adds rows outside the anchored policy",
    )
    for candidate_id, candidate in candidates.items():
        require(
            isinstance(candidate_id, str)
            and isinstance(candidate, dict)
            and candidate.get("candidate_id") == candidate_id,
            "CANDIDATE_ID_DRIFT",
            f"candidate key differs from candidate_id: {candidate_id}",
        )
        require(
            isinstance(candidate.get("candidate_class"), int)
            and not isinstance(candidate.get("candidate_class"), bool),
            "CANDIDATE_CLASS_INVALID",
            f"candidate class is not an integer: {candidate_id}",
        )
        require(
            1 <= candidate["candidate_class"] <= 4,
            "CANDIDATE_CLASS_INVALID",
            f"candidate class is outside the frozen range: {candidate_id}",
        )
        require(
            isinstance(candidate.get("selection_priority"), int)
            and not isinstance(candidate.get("selection_priority"), bool)
            and candidate["selection_priority"] >= 1,
            "CANDIDATE_PRIORITY_INVALID",
            f"candidate selection priority is invalid: {candidate_id}",
        )
        require(
            {
                "candidate_class": candidate["candidate_class"],
                "selection_priority": candidate["selection_priority"],
            }
            == policy["required_candidates"][candidate_id],
            "CANDIDATE_CLASS_PRIORITY_POLICY_MISMATCH",
            f"candidate class or priority differs from policy: {candidate_id}",
        )
        if candidate_id in policy["permanently_ineligible_candidates"]:
            require(
                candidate.get("eligibility")
                == policy["permanently_ineligible_candidates"][candidate_id],
                "PERMANENTLY_INELIGIBLE_CANDIDATE_PROMOTED",
                f"policy forbids forecast promotion of candidate: {candidate_id}",
            )
        require(
            isinstance(candidate.get("eligibility"), str),
            "CANDIDATE_ELIGIBILITY_INVALID",
            f"candidate eligibility is not a string: {candidate_id}",
        )
    selected_candidates = [
        candidate
        for candidate in candidates.values()
        if candidate["eligibility"] == "ELIGIBLE_SOURCE_VISIBLE"
    ]
    for candidate in selected_candidates:
        lexical = check_target_ancestry(
            candidate,
            forbidden_classes,
        )
        require(
            lexical["declared_vocabulary_hit_free"],
            "SELECTED_CANDIDATE_DECLARED_ANCESTRY_HIT",
            "an eligible candidate contains forbidden declared vocabulary",
        )
    static_rows = inventory.get("known_surface_crosswalk", [])
    require(
        static_rows == policy["required_static_crosswalk_rows"],
        "STATIC_CROSSWALK_POLICY_MISMATCH",
        "inventory static crosswalk content differs from anchored policy",
    )
    static_surfaces = [
        row["surface"]
        for row in static_rows
        if isinstance(row, dict) and isinstance(row.get("surface"), str)
    ]
    require(
        len(static_surfaces) == len(set(static_surfaces))
        and set(static_surfaces)
        == set(policy["required_static_crosswalk_surfaces"]),
        "REQUIRED_STATIC_CROSSWALK_MISMATCH",
        "inventory static crosswalk differs from the anchored policy",
    )
    eligible = _eligible_candidate_ids(inventory)
    require(
        bool(eligible),
        "NO_ELIGIBLE_CANDIDATE",
        "the frozen inventory contains no eligible source-visible candidate",
    )
    return eligible


def _validate_discriminator(
    contract: dict[str, Any],
    inventory: dict[str, Any],
    repo_root: Path,
    commit_sha: str,
) -> None:
    discriminator_pin = contract["discriminator_receipt"]
    verify_file_pin(discriminator_pin, repo_root, path_key="artifact")
    try:
        receipt = load_json(
            _resolve_repo_file(repo_root, discriminator_pin["artifact"])
        )
    except (json.JSONDecodeError, UnicodeDecodeError) as error:
        raise ContractError(
            "DISCRIMINATOR_INVALID",
            "discriminator receipt is not valid UTF-8 JSON",
        ) from error
    expected_keys = {
        "schema",
        "candidate_id",
        "freeze_policy_sha256",
        "inventory_sha256",
        "live_fz_register_sha256",
        "candidate_row_sha256",
        "allowed_ancestry_artifacts_sha256",
        "reviewed_source_ids",
        "evidence_pins",
        "verdict",
        "source_visible",
        "declared_ancestry_review_complete",
        "semantic_input_closure_proved",
        "review_scope",
    }
    require(
        set(receipt) == expected_keys
        and receipt.get("schema") == DISCRIMINATOR_SCHEMA,
        "DISCRIMINATOR_SCHEMA",
        "discriminator receipt does not match its exact schema",
    )
    require(
        receipt.get("candidate_id") == contract["candidate_id"],
        "DISCRIMINATOR_CANDIDATE_MISMATCH",
        "discriminator receipt names another candidate",
    )
    candidate = inventory["candidates"][contract["candidate_id"]]
    expected = candidate.get("allowed_ancestry_artifacts")
    require(
        isinstance(expected, list) and expected,
        "CANDIDATE_ANCESTRY_NOT_PINNED",
        "candidate inventory lacks exact allowed-ancestry artifact pins",
    )
    expected_source_ids = sorted(pin["artifact"] for pin in expected)
    require(
        contract["allowed_ancestry"] == expected,
        "CONTRACT_ANCESTRY_DRIFT",
        "contract ancestry differs from the promoted inventory row",
    )
    require(
        receipt.get("freeze_policy_sha256")
        == contract["freeze_policy"]["sha256"]
        and receipt.get("inventory_sha256")
        == contract["candidate_inventory"]["sha256"]
        and receipt.get("live_fz_register_sha256")
        == contract["live_fz_register"]["sha256"]
        and receipt.get("candidate_row_sha256")
        == canonical_candidate_row_sha256(candidate)
        and receipt.get("allowed_ancestry_artifacts_sha256")
        == canonical_pin_list_sha256(expected)
        and receipt.get("reviewed_source_ids") == expected_source_ids,
        "DISCRIMINATOR_BINDING_MISMATCH",
        "discriminator does not bind the selected inventory and ancestry",
    )
    require(
        receipt.get("verdict")
        == "ELIGIBLE_SOURCE_VISIBLE_WITH_BOUNDED_DECLARED_ANCESTRY_REVIEW"
        and receipt.get("source_visible") is True
        and receipt.get("declared_ancestry_review_complete") is True
        and receipt.get("semantic_input_closure_proved") is False
        and receipt.get("review_scope")
        == "bounded_declared_artifacts_not_transitive_semantic_closure",
        "DISCRIMINATOR_NOT_POSITIVE",
        "discriminator does not establish bounded source-visible eligibility",
    )
    evidence_pins = receipt.get("evidence_pins")
    require(
        isinstance(evidence_pins, list) and bool(evidence_pins),
        "DISCRIMINATOR_EVIDENCE_MISSING",
        "discriminator has no hash-pinned review evidence",
    )
    for raw_pin in evidence_pins:
        evidence_pin = _require_artifact_pin_shape(
            raw_pin,
            code="DISCRIMINATOR_EVIDENCE_PIN_INVALID",
            message="discriminator evidence pin has the wrong shape",
        )
        verify_file_pin(evidence_pin, repo_root, path_key="artifact")
        verify_file_pin_at_commit(
            evidence_pin,
            repo_root,
            commit_sha,
            path_key="artifact",
        )
        try:
            evidence = load_json(
                _resolve_repo_file(repo_root, evidence_pin["artifact"])
            )
        except (json.JSONDecodeError, UnicodeDecodeError) as error:
            raise ContractError(
                "DISCRIMINATOR_EVIDENCE_INVALID",
                "discriminator evidence is not valid UTF-8 JSON",
            ) from error
        require(
            set(evidence)
            == {
                "schema",
                "candidate_id",
                "reviewed_source_ids",
                "reviewed_source_pins_sha256",
                "forbidden_input_classes_sha256",
                "review_method",
                "verdict",
                "semantic_input_closure_proved",
                "reviewer_id",
                "reviewer_attestation",
            }
            and evidence.get("schema") == ANCESTRY_EVIDENCE_SCHEMA
            and evidence.get("candidate_id") == contract["candidate_id"]
            and evidence.get("reviewed_source_ids") == expected_source_ids
            and evidence.get("reviewed_source_pins_sha256")
            == canonical_pin_list_sha256(expected)
            and evidence.get("forbidden_input_classes_sha256")
            == sha256_bytes(
                canonical_json(
                    inventory["forbidden_input_classes"]
                ).encode("utf-8")
            )
            and evidence.get("review_method")
            == "bounded_declared_artifact_ancestry_review_v1"
            and evidence.get("verdict")
            == "NO_DECLARED_FORBIDDEN_INPUT_FOUND"
            and evidence.get("semantic_input_closure_proved") is False,
            "DISCRIMINATOR_EVIDENCE_NOT_POSITIVE",
            "discriminator evidence does not bind a positive bounded review",
        )
        require(
            isinstance(evidence.get("reviewer_id"), str)
            and bool(evidence["reviewer_id"].strip())
            and isinstance(evidence.get("reviewer_attestation"), str)
            and len(evidence["reviewer_attestation"].strip()) >= 16,
            "DISCRIMINATOR_EVIDENCE_ATTESTATION_MISSING",
            "discriminator evidence lacks reviewer identity or attestation",
        )


def _validate_custody(
    contract: dict[str, Any],
    repo_root: Path,
) -> None:
    custody_pin = contract["sealed_comparison"]["custody_record"]
    try:
        custody = load_json(
            _resolve_repo_file(repo_root, custody_pin["artifact"])
        )
    except (json.JSONDecodeError, UnicodeDecodeError) as error:
        raise ContractError(
            "CUSTODY_RECORD_INVALID",
            "comparison custody record is not valid UTF-8 JSON",
        ) from error
    require(
        set(custody)
        == {
            "schema",
            "contract_id",
            "payload_sha256",
            "byte_count",
            "independent_custodian_id",
            "independence_attestation",
            "attestation_id",
            "storage_locator",
            "access_history",
            "access_history_complete",
            "producer_access_before_freeze",
            "comparison_disclosed_before_freeze",
        }
        and custody.get("schema") == CUSTODY_SCHEMA,
        "CUSTODY_RECORD_SCHEMA",
        "comparison custody record does not match its exact schema",
    )
    sealed = contract["sealed_comparison"]
    require(
        custody.get("contract_id") == contract["contract_id"]
        and custody.get("payload_sha256") == sealed["payload_sha256"]
        and custody.get("byte_count") == sealed["byte_count"],
        "CUSTODY_COMPARISON_BINDING_MISMATCH",
        "custody record does not bind the sealed comparison",
    )
    require(
        custody.get("access_history") == []
        and custody.get("access_history_complete") is True
        and custody.get("producer_access_before_freeze") is False
        and custody.get("comparison_disclosed_before_freeze") is False,
        "CUSTODY_PRE_FREEZE_ACCESS_NOT_DENIED",
        "custody record does not deny and exhaustively record pre-freeze access",
    )
    custodian_id = custody.get("independent_custodian_id")
    require(
        isinstance(custodian_id, str)
        and bool(custodian_id.strip())
        and custodian_id
        not in {
            contract["generator"]["module"],
            contract["independent_checker"]["module"],
            "producer",
        }
        and isinstance(custody.get("independence_attestation"), str)
        and len(custody["independence_attestation"].strip()) >= 16
        and isinstance(custody.get("attestation_id"), str)
        and bool(custody["attestation_id"].strip())
        and isinstance(custody.get("storage_locator"), str)
        and bool(custody["storage_locator"].strip()),
        "CUSTODY_INDEPENDENCE_ATTESTATION_MISSING",
        "custody record lacks an independent identity and attestation",
    )


def validate_contract_for_freeze(
    contract: dict[str, Any],
    *,
    repo_root: Path = REPO_ROOT,
    inventory: dict[str, Any] | None = None,
    external_contract_digest_input: str | None = None,
) -> None:
    """Validate declared pre-freeze controls, without claiming scoring readiness."""

    validate_contract(contract)
    verify_frozen_contract(contract, external_contract_digest_input)
    commit_sha = contract["freeze"]["repository_state_commit_sha"]
    verify_repository_commit(repo_root, commit_sha)
    policy, live_register = _validate_anchored_policy(
        contract,
        repo_root=repo_root,
        commit_sha=commit_sha,
    )

    inventory_pin = contract["candidate_inventory"]
    verify_file_pin(inventory_pin, repo_root, path_key="artifact")
    verify_file_pin_at_commit(
        inventory_pin,
        repo_root,
        commit_sha,
        path_key="artifact",
    )
    try:
        pinned_inventory = load_json(
            _resolve_repo_file(repo_root, inventory_pin["artifact"])
        )
    except (json.JSONDecodeError, UnicodeDecodeError) as error:
        raise ContractError(
            "INVENTORY_INVALID",
            "candidate inventory is not valid UTF-8 JSON",
        ) from error
    if inventory is not None:
        require(
            inventory == pinned_inventory,
            "INVENTORY_ARGUMENT_DRIFT",
            "supplied inventory differs from the hash-pinned inventory",
        )
    inventory = pinned_inventory
    eligible = _validate_frozen_inventory(
        inventory,
        repo_root=repo_root,
        commit_sha=commit_sha,
        policy=policy,
        live_register=live_register,
        policy_sha256=contract["freeze_policy"]["sha256"],
    )
    candidate_id = contract["candidate_id"]
    require(
        candidate_id in inventory["candidates"],
        "CANDIDATE_NOT_REGISTERED",
        f"candidate is absent from the inventory: {candidate_id}",
    )
    candidate = inventory["candidates"][candidate_id]
    require(
        candidate["candidate_id"] == candidate_id,
        "CANDIDATE_ID_DRIFT",
        "candidate key and candidate_id field differ",
    )
    require(
        candidate["eligibility"] == "ELIGIBLE_SOURCE_VISIBLE",
        "CANDIDATE_NOT_ELIGIBLE",
        f"candidate status is {candidate['eligibility']}",
    )
    require(
        candidate_id == eligible[0],
        "CANDIDATE_SELECTION_DRIFT",
        f"frozen selection requires {eligible[0]}, not {candidate_id}",
    )

    expected_forbidden = set(inventory["forbidden_input_classes"])
    require(
        set(contract["forbidden_input_classes"]) == expected_forbidden,
        "FORBIDDEN_INPUT_SET_DRIFT",
        "contract must carry the full frozen forbidden-input class set",
    )

    pins = [
        (pin, "artifact") for pin in contract["allowed_ancestry"]
    ] + [
        (contract["generator"], "module"),
        (contract["independent_checker"], "module"),
        (contract["sealed_comparison"]["custody_record"], "artifact"),
        (contract["discriminator_receipt"], "artifact"),
    ]
    for pin, path_key in pins:
        verify_file_pin(pin, repo_root, path_key=path_key)
        verify_file_pin_at_commit(
            pin,
            repo_root,
            commit_sha,
            path_key=path_key,
        )
    _validate_custody(contract, repo_root)
    _validate_discriminator(
        contract,
        inventory,
        repo_root,
        commit_sha,
    )


def unseal_and_score(*_args: Any, **_kwargs: Any) -> dict[str, Any]:
    """Fail closed until an executable independent scorer is implemented."""

    raise ContractError(
        "SCORER_NOT_IMPLEMENTED",
        "this packet contains no generator execution, independent scoring, "
        "comparison quarantine, or durable single-use record",
    )


def build_fz_crosswalk(register: dict[str, Any]) -> list[dict[str, Any]]:
    """Project the live v3 prospective and retrospective register surfaces."""

    require(
        register.get("schema") == "oph.frozen_prediction_register.v3",
        "FZ_REGISTER_SCHEMA",
        "forecast crosswalk requires the live frozen-register v3 schema",
    )
    projected: list[dict[str, Any]] = []
    seen: set[str] = set()
    for row in register.get("rows", []):
        surface = row["id"]
        require(
            surface not in seen,
            "FZ_REGISTER_DUPLICATE_ID",
            f"duplicate live prospective register row: {surface}",
        )
        seen.add(surface)
        projected.append(
            {
                "surface": surface,
                "classification": (
                    "PROSPECTIVE_REGISTER_ROW__"
                    + str(row["status"]).upper()
                ),
                "forecast_use": (
                    "Use only under the live registered comparison protocol: "
                    + str(row["comparison_protocol"])
                ),
                "registry_kind": "prospective",
                "registry_status": row["status"],
                "owning_issue": row.get("owning_issue"),
            }
        )
    for row in register.get("retrospective_results", []):
        surface = row["id"]
        require(
            surface not in seen,
            "FZ_REGISTER_DUPLICATE_ID",
            f"duplicate live retrospective register row: {surface}",
        )
        seen.add(surface)
        projected.append(
            {
                "surface": surface,
                "classification": (
                    "RETROSPECTIVE_RESULT__"
                    + str(row["status"]).upper()
                ),
                "forecast_use": str(row["evidential_boundary"]),
                "registry_kind": "retrospective",
                "registry_status": row["status"],
                "owning_issue": row.get("owning_issue"),
                "former_ladder_reservation": row.get(
                    "former_ladder_reservation"
                ),
            }
        )
    return projected


def build_state() -> dict[str, Any]:
    """Emit the bounded state of the draft governance packet."""

    schema_raw = SCHEMA_PATH.read_bytes()
    policy_raw = POLICY_PATH.read_bytes()
    policy = json.loads(policy_raw.decode("utf-8"))
    inventory_raw = INVENTORY_PATH.read_bytes()
    inventory = json.loads(inventory_raw.decode("utf-8"))
    fz_register_raw = FZ_REGISTER_PATH.read_bytes()
    fz_register = json.loads(fz_register_raw.decode("utf-8"))

    require(
        policy.get("schema") == POLICY_SCHEMA and policy.get("issue") == ISSUE,
        "FREEZE_POLICY_INVALID",
        "canonical issue-639 policy has the wrong schema or issue",
    )
    require(
        inventory.get("inventory_scope") == policy["inventory_scope"],
        "INVENTORY_SCOPE_POLICY_MISMATCH",
        "provisional inventory scope differs from the anchored issue policy",
    )
    require(
        inventory.get("forbidden_input_classes")
        == policy["forbidden_input_classes"],
        "FORBIDDEN_INPUT_POLICY_MISMATCH",
        "provisional inventory forbidden map differs from anchored policy",
    )
    require(
        inventory.get("completeness", {}).get("closure_criterion_id")
        == policy["closure_policy"]["closure_criterion_id"]
        and inventory.get("completeness", {}).get("closure_criterion")
        == policy["closure_policy"]["closure_criterion"],
        "INVENTORY_CLOSURE_POLICY_MISMATCH",
        "provisional inventory closure criterion differs from anchored policy",
    )
    require(
        inventory.get("selection_rule", {}).get("algorithm")
        == policy["selection_algorithm"]
        == FROZEN_SELECTION_ALGORITHM,
        "SELECTION_RULE_POLICY_MISMATCH",
        "provisional inventory selection algorithm differs from anchored policy",
    )
    candidates = inventory.get("candidates", {})
    require(
        set(candidates) == set(policy["required_candidates"]),
        "REQUIRED_CANDIDATE_SET_MISMATCH",
        "provisional inventory candidate set differs from anchored policy",
    )
    for candidate_id, required_order in policy["required_candidates"].items():
        require(
            {
                "candidate_class": candidates[candidate_id].get(
                    "candidate_class"
                ),
                "selection_priority": candidates[candidate_id].get(
                    "selection_priority"
                ),
            }
            == required_order,
            "CANDIDATE_CLASS_PRIORITY_POLICY_MISMATCH",
            f"provisional class or priority differs from policy: {candidate_id}",
        )
        if candidate_id in policy["permanently_ineligible_candidates"]:
            require(
                candidates[candidate_id].get("eligibility")
                == policy["permanently_ineligible_candidates"][candidate_id],
                "PERMANENTLY_INELIGIBLE_CANDIDATE_PROMOTED",
                f"provisional inventory promotes excluded row: {candidate_id}",
            )
    require(
        inventory.get("known_surface_crosswalk")
        == policy["required_static_crosswalk_rows"],
        "STATIC_CROSSWALK_POLICY_MISMATCH",
        "provisional static crosswalk content differs from anchored policy",
    )
    static_surfaces = [
        row.get("surface")
        for row in inventory.get("known_surface_crosswalk", [])
        if isinstance(row, dict)
    ]
    require(
        len(static_surfaces) == len(set(static_surfaces))
        and set(static_surfaces)
        == set(policy["required_static_crosswalk_surfaces"]),
        "REQUIRED_STATIC_CROSSWALK_MISMATCH",
        "provisional static crosswalk differs from anchored policy",
    )

    forbidden = inventory["forbidden_input_classes"]
    ancestry_verdicts: dict[str, Any] = {}
    eligible: list[str] = []
    potential_upgrade_issues: set[int] = set()
    for name, candidate in sorted(inventory["candidates"].items()):
        require(
            candidate["candidate_id"] == name,
            "CANDIDATE_ID_DRIFT",
            f"inventory key differs from candidate_id: {name}",
        )
        verdict = check_target_ancestry(candidate, forbidden)
        ancestry_verdicts[name] = {
            "candidate_class": candidate["candidate_class"],
            "selection_priority": candidate["selection_priority"],
            "eligibility": candidate["eligibility"],
            "blocking_condition": candidate["blocking_condition"],
            "potential_upgrade_issues": candidate.get(
                "potential_upgrade_issues", []
            ),
            **verdict,
        }
        potential_upgrade_issues.update(
            candidate.get("potential_upgrade_issues", [])
        )
        if candidate["eligibility"] == "ELIGIBLE_SOURCE_VISIBLE":
            eligible.append(name)
    if eligible:
        eligible = _eligible_candidate_ids(inventory)

    payload = {
        "schema": STATE_SCHEMA,
        "issue": ISSUE,
        "pins": {
            "contract_schema_sha256": sha256_bytes(schema_raw),
            "forecast_freeze_policy_sha256": sha256_bytes(policy_raw),
            "candidate_inventory_sha256": sha256_bytes(inventory_raw),
            "frozen_prediction_register_sha256": sha256_bytes(fz_register_raw),
            "module_sha256": sha256_bytes(Path(__file__).read_bytes()),
        },
        "inventory_scope": inventory["inventory_scope"],
        "inventory_completeness": inventory["completeness"],
        "known_surface_crosswalk": (
            build_fz_crosswalk(fz_register)
            + inventory["known_surface_crosswalk"]
        ),
        "forbidden_input_classes": forbidden,
        "selection_rule": inventory["selection_rule"],
        "candidate_ledger": ancestry_verdicts,
        "eligible_candidates": eligible,
        "contract_freeze_status": (
            "DRAFT_GOVERNANCE_PACKET__NO_ELIGIBLE_CANDIDATE"
            if not eligible
            else "ELIGIBLE_ROW_REQUIRES_FULL_FREEZE_VALIDATION"
        ),
        "controls": {
            "full_json_schema_validation": "IMPLEMENTED",
            "cross_field_checker_separation": (
                "IMPLEMENTED_PATH_AND_BYTE_INEQUALITY_ONLY"
            ),
            "freeze_time_file_pin_resolution": "IMPLEMENTED",
            "freeze_time_discriminator_verification": (
                "IMPLEMENTED_EXACT_BINDINGS_AND_PINNED_BOUNDED_EVIDENCE"
            ),
            "freeze_time_comparison_custody_declaration_validation": (
                "IMPLEMENTED_EXACT_SCHEMA_BINDING_AND_PRE_FREEZE_DENIAL"
            ),
            "freeze_time_inventory_exhaustiveness_gate": (
                "IMPLEMENTED_REQUIRES_FROZEN_CRITERION_AND_PINNED_POSITIVE_"
                "OMISSION_AUDIT"
            ),
            "immutable_issue_freeze_policy_anchor": (
                "IMPLEMENTED_EXACT_SCOPE_FORBIDDEN_MAP_CANDIDATES_"
                "NONPROMOTABLE_STATUSES_CROSSWALK_CONTENT_ORDER_AND_"
                "LIVE_REGISTER_PROJECTION"
            ),
            "freeze_time_deterministic_first_eligible_selection": "IMPLEMENTED",
            "freeze_time_repository_commit_blob_verification": "IMPLEMENTED",
            "non_self_referential_canonical_contract_digest": "IMPLEMENTED",
            "external_contract_digest_input_validation": (
                "IMPLEMENTED_CALLER_SUPPLIED_INPUT_NOT_DURABLE_CUSTODY"
            ),
            "durable_external_contract_digest_custody": "NOT_IMPLEMENTED",
            "declared_vocabulary_scan": (
                "IMPLEMENTED_BOUNDED_WARNING_NOT_SEMANTIC_CLOSURE"
            ),
            "comparison_integrity_digest": "IMPLEMENTED",
            "comparison_access_quarantine": "NOT_IMPLEMENTED",
            "executable_generator_validation": "NOT_IMPLEMENTED",
            "checker_true_independence_validation": "NOT_IMPLEMENTED",
            "executable_independent_scorer": "NOT_IMPLEMENTED",
            "durable_single_use_unsealing_record": "NOT_IMPLEMENTED",
            "candidate_inventory_exhaustiveness": "NOT_ESTABLISHED",
            "freeze_validator_scope": (
                "DECLARED_PRE_FREEZE_CONTROLS_ONLY__NOT_SCORING_READY"
            ),
        },
        "potential_upgrade_routes": {
            "issues": sorted(potential_upgrade_issues),
            "dependency_status": (
                "candidate-specific routes only; none is a hard issue "
                "dependency before selection"
            ),
        },
        "claim_boundary": (
            "Draft issue-639 governance packet with deterministic state replay, "
            "complete JSON Schema enforcement for the draft contract, bounded "
            "declared-vocabulary warnings, and fail-closed freeze-time gates "
            "for exhaustive inventory evidence, deterministic selection, real "
            "repository-commit blobs, canonical contract payload integrity, "
            "an immutable issue-policy anchor, a caller-supplied digest check, "
            "exact custody declarations, source pins, and bounded declared-"
            "ancestry eligibility. The external "
            "digest input is not durable custody. The candidate inventory is "
            "provisional and not proved exhaustive. Generator execution, "
            "checker independence, comparison quarantine, scoring, and durable "
            "single-use recording are not implemented. This validator covers "
            "declared pre-freeze controls and does not establish scoring "
            "readiness. No "
            "contract, prediction output, independent score, or persistent "
            "single-use unsealing record is emitted."
        ),
    }
    payload["state_sha256"] = sha256_bytes(
        canonical_json(
            {k: v for k, v in payload.items() if k != "state_sha256"}
        ).encode("utf-8")
    )
    return payload


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args(argv)
    payload = build_state()
    if args.verify:
        stored = json.loads(OUT_PATH.read_text(encoding="utf-8"))
        if stored != payload:
            print("CONTRACT_STATE_DRIFT", file=sys.stderr)
            return 1
        print("CONTRACT_STATE_VERIFIED")
        return 0
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(
        json.dumps(payload, sort_keys=True, indent=1) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"status": payload["contract_freeze_status"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
