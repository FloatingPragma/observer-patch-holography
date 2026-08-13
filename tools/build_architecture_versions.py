#!/usr/bin/env python3
"""Validate and render the immutable V3 architecture-version register."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

import strict_json

ROOT = Path(__file__).resolve().parents[1]
REGISTER_PATH = ROOT / "tracking" / "architecture_versions.json"
SURFACE_PATH = ROOT / "docs" / "ARCHITECTURE_VERSION_REGISTER.md"
PREMISE_PATH = ROOT / "tracking" / "premise_register.json"

SCHEMA = "oph.architecture_version_register.v3"
LEGACY_ROOT_SCHEMA = "oph.architecture_version_register.v2"
AV0_FIRST_APPEARANCE_REVISION = "3ab5bc2064235a740bb5574ea165564e43046bca"
AV0_FIRST_APPEARANCE_RECORD_SHA256 = (
    "6c1c826255c88011d606c34e2317833de44573885d5e7bc1abcf7fb0189b1291"
)
ISSUE = 741
VERSION_RE = re.compile(r"^AV-(0|[1-9][0-9]*)$")
DIGEST_RE = re.compile(r"^[0-9a-f]{64}$")
GIT_BLOB_RE = re.compile(r"^[0-9a-f]{40}$")
GIT_REVISION_RE = re.compile(r"^[0-9a-f]{40}$")
STATUS_ENUM = {"exploratory_uninhabited", "inhabited_conditional", "retired"}
PENDING_TIP_PROMOTION_STATUS = (
    "At its origin commit this version is an unanchored, promotion-ineligible "
    "tip. After its complete record receives an immutable origin anchor, each "
    "row still requires replayed evidence and an independent version-matching "
    "audit before promotion."
)
PREMISE_REGISTER_RELATIVE = "tracking/premise_register.json"

TOP_KEYS = {
    "schema",
    "issue",
    "current_version",
    "version_anchors",
    "policy",
    "versions",
}
ANCHOR_KEYS = {"id", "origin_revision", "record_sha256"}
VERSION_KEYS = {
    "id",
    "created_on",
    "status",
    "snapshot_sha256",
    "predecessor_register_revision",
    "basis",
    "normative_files",
    "protocol_decisions",
    "common_world_status",
    "promotion_status",
    "invalidation_triggers",
    "replay_surfaces",
}
BASIS_KEYS = {
    "core_axioms",
    "closure_principle",
    "proposed_fundamental_numerical_parameters",
}
FILE_KEYS = {"path", "role", "sha256", "git_blob_sha1"}
DECISION_KEYS = {"id", "decision", "premises"}
LEGACY_SNAPSHOT_KEYS = (
    "id",
    "created_on",
    "basis",
    "normative_files",
    "protocol_decisions",
    "invalidation_triggers",
    "replay_surfaces",
)
SNAPSHOT_KEYS = (
    "id",
    "created_on",
    "status",
    "predecessor_register_revision",
    "basis",
    "normative_files",
    "protocol_decisions",
    "common_world_status",
    "promotion_status",
    "invalidation_triggers",
    "replay_surfaces",
)


def fail(message: str) -> None:
    raise SystemExit(f"architecture versions: {message}")


def _display_path(path: Path) -> str:
    try:
        return path.relative_to(ROOT).as_posix()
    except ValueError:
        return str(path)


def load_json(path: Path) -> dict:
    try:
        return strict_json.load(path)
    except FileNotFoundError:
        fail(f"missing {_display_path(path)}")
    except (json.JSONDecodeError, strict_json.DuplicateKeyError) as error:
        fail(f"invalid JSON in {_display_path(path)}: {error}")
    raise AssertionError("unreachable")


def _nonempty(where: str, value: object) -> str:
    if not isinstance(value, str) or not value.strip():
        fail(f"{where} must be a nonempty string")
    return value


def _string_list(where: str, value: object, *, allow_empty: bool = False) -> list[str]:
    if not isinstance(value, list) or (not value and not allow_empty):
        fail(
            f"{where} must be a {'possibly empty ' if allow_empty else 'nonempty '}list"
        )
    if any(not isinstance(item, str) or not item.strip() for item in value):
        fail(f"{where} entries must be nonempty strings")
    if len(value) != len(set(value)):
        fail(f"{where} must be duplicate-free")
    return value


def _git_blob_sha1(payload: bytes) -> str:
    """Return Git's SHA-1 object id for one blob without writing an object."""

    header = f"blob {len(payload)}\0".encode("ascii")
    return hashlib.sha1(header + payload).hexdigest()


def _snapshot_payload(version: dict, *, legacy: bool = False) -> bytes:
    """Canonical immutable portion of one architecture version."""

    keys = LEGACY_SNAPSHOT_KEYS if legacy else SNAPSHOT_KEYS
    projection = {key: version[key] for key in keys}
    return json.dumps(
        projection,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")


def snapshot_sha256(version: dict) -> str:
    return hashlib.sha256(_snapshot_payload(version)).hexdigest()


def legacy_snapshot_sha256(version: dict) -> str:
    return hashlib.sha256(_snapshot_payload(version, legacy=True)).hexdigest()


def _git_output(*args: str) -> bytes:
    try:
        result = subprocess.run(
            ["git", *args],
            cwd=ROOT,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except OSError as error:
        fail(f"cannot execute git: {error}")
    if result.returncode != 0:
        detail = result.stderr.decode("utf-8", errors="replace").strip()
        fail(f"git {' '.join(args)} failed: {detail or 'object is unavailable'}")
    return result.stdout


def _historical_blob(version_id: str, path: str, object_id: str) -> bytes:
    try:
        return _git_output("cat-file", "blob", object_id)
    except SystemExit as error:
        fail(
            f"{version_id}: historical Git blob {object_id} for {path} "
            f"is unavailable ({error})"
        )
    raise AssertionError("unreachable")


def _pinned_register(revision: str) -> dict:
    object_type = _git_output("cat-file", "-t", revision).strip()
    if object_type != b"commit":
        fail(f"predecessor revision {revision} is not a commit")
    raw = _git_output("show", f"{revision}:tracking/architecture_versions.json")
    try:
        result = strict_json.loads(raw.decode("utf-8"))
    except (
        UnicodeDecodeError,
        json.JSONDecodeError,
        strict_json.DuplicateKeyError,
    ) as error:
        fail(f"predecessor revision {revision} has an invalid register: {error}")
    if not isinstance(result, dict):
        fail(f"predecessor revision {revision} has a malformed register")
    return result


def _optional_pinned_register(revision: str) -> dict | None:
    result = subprocess.run(
        ["git", "show", f"{revision}:tracking/architecture_versions.json"],
        cwd=ROOT,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if result.returncode != 0:
        return None
    try:
        payload = strict_json.loads(result.stdout.decode("utf-8"))
    except (
        UnicodeDecodeError,
        json.JSONDecodeError,
        strict_json.DuplicateKeyError,
    ) as error:
        fail(f"revision {revision} has an invalid architecture register: {error}")
    if not isinstance(payload, dict):
        fail(f"revision {revision} has a malformed architecture register")
    return payload


def _parents(revision: str) -> list[str]:
    return _git_output("show", "-s", "--format=%P", revision).decode().strip().split()


def _git_is_ancestor(ancestor: str, descendant: str) -> bool:
    result = subprocess.run(
        ["git", "merge-base", "--is-ancestor", ancestor, descendant],
        cwd=ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if result.returncode == 0:
        return True
    if result.returncode == 1:
        return False
    detail = result.stderr.decode("utf-8", errors="replace").strip()
    fail(
        "git merge-base --is-ancestor "
        f"{ancestor} {descendant} failed: {detail or 'object unavailable'}"
    )
    raise AssertionError("unreachable")


def _git_is_first_parent_ancestor(ancestor: str, descendant: str) -> bool:
    revisions = {
        item
        for item in _git_output("rev-list", "--first-parent", descendant)
        .decode()
        .splitlines()
        if item
    }
    return ancestor in revisions


def _single_parent(revision: str, version_id: str) -> str:
    parents = _parents(revision)
    if len(parents) != 1:
        fail(f"{version_id}: origin revision must be a single-parent commit")
    return parents[0]


def _validate_committed_history_guard(
    versions: list[dict], anchors: list[dict]
) -> None:
    """Reject version rewrites, anchor shrinkage, and delete/reintroduce gaps."""

    head = _git_output("rev-parse", "HEAD").decode().strip()
    revisions = (
        _git_output(
            "log",
            "--first-parent",
            "--full-history",
            "--format=%H",
            head,
            "--",
            "tracking/architecture_versions.json",
        )
        .decode()
        .strip()
        .splitlines()
    )
    current_ids = [version["id"] for version in versions]
    current_by_id = {version["id"]: version for version in versions}
    current_anchors = {anchor["id"]: anchor for anchor in anchors}
    saw_schema = False
    left_schema_span = False
    for revision in revisions:
        historical = _optional_pinned_register(revision)
        if historical is None or historical.get("schema") != SCHEMA:
            if saw_schema:
                left_schema_span = True
            continue
        if left_schema_span:
            fail(
                "architecture register cannot disappear and later reappear in "
                "Git history"
            )
        saw_schema = True
        historical_versions = historical.get("versions")
        historical_anchors = historical.get("version_anchors")
        if not isinstance(historical_versions, list) or not isinstance(
            historical_anchors, list
        ):
            fail(f"{revision}: historical architecture register is malformed")
        historical_ids = [
            version.get("id") if isinstance(version, dict) else None
            for version in historical_versions
        ]
        if current_ids[: len(historical_ids)] != historical_ids:
            fail("architecture version history is not append-only")
        for historical_version in historical_versions:
            version_id = historical_version["id"]
            if current_by_id[version_id] != historical_version:
                fail(f"{version_id}: record rewrites committed architecture history")
        for historical_anchor in historical_anchors:
            version_id = historical_anchor.get("id")
            if current_anchors.get(version_id) != historical_anchor:
                fail(f"{version_id}: origin anchor rewrites committed history")


def _validate_committed_head_anchors(anchors: list[dict]) -> None:
    """An origin anchor is operative only after its declaration is committed.

    The origin revision necessarily predates the anchor declaration.  Accepting
    a worktree-only anchor would collapse that two-commit protocol and let
    generated consumers treat an uncommitted declaration as durable custody.
    """

    head = _git_output("rev-parse", "HEAD").decode().strip()
    committed = _optional_pinned_register(head)
    committed_anchors = (
        committed.get("version_anchors")
        if isinstance(committed, dict) and committed.get("schema") == SCHEMA
        else []
    )
    if not isinstance(committed_anchors, list):
        fail("committed HEAD architecture anchors are malformed")
    # The canonical checkout is presently performing the one-time v2 -> v3
    # migration.  AV-0's declaration is independently hard-bound above to its
    # audited first appearance and complete-record digest, so allow exactly
    # that bootstrap declaration until the first v3 register commit lands.
    # This exception cannot admit a successor or an arbitrary fixture anchor.
    if (
        not isinstance(committed, dict)
        or committed.get("schema") != SCHEMA
    ) and ROOT == REGISTER_PATH.parents[1] and anchors == [
        {
            "id": "AV-0",
            "origin_revision": AV0_FIRST_APPEARANCE_REVISION,
            "record_sha256": AV0_FIRST_APPEARANCE_RECORD_SHA256,
        }
    ]:
        return
    if len(anchors) > len(committed_anchors) or anchors != committed_anchors[: len(anchors)]:
        fail("every operative version anchor must already be committed at HEAD")


def replay_surface_bindings(
    data: dict, version: dict, *, is_current: bool
) -> list[dict[str, str]]:
    """Resolve replay surfaces at the version origin, or live for a pending tip."""

    anchors = {
        anchor["id"]: anchor
        for anchor in data.get("version_anchors", [])
        if isinstance(anchor, dict) and isinstance(anchor.get("id"), str)
    }
    version_id = version["id"]
    anchor = anchors.get(version_id)
    bindings: list[dict[str, str]] = []
    for path in version["replay_surfaces"]:
        if path.startswith("/") or ".." in path.split("/"):
            fail(f"{version_id}: replay surfaces must be safe repo-relative paths")
        if anchor is not None:
            revision = anchor["origin_revision"]
            payload = _git_output("show", f"{revision}:{path}")
            blob = _git_output("rev-parse", f"{revision}:{path}").decode().strip()
            if not GIT_BLOB_RE.fullmatch(blob):
                fail(f"{version_id}: replay surface is not a Git blob: {path}")
            object_type = _git_output("cat-file", "-t", blob).strip()
            if object_type != b"blob":
                fail(f"{version_id}: replay surface is not a Git blob: {path}")
            custody = revision
        elif is_current:
            resolved = ROOT / path
            if not resolved.is_file():
                fail(f"{version_id}: live replay surface missing: {path}")
            payload = resolved.read_bytes()
            blob = _git_blob_sha1(payload)
            custody = "pending_live_tip"
        else:
            fail(f"{version_id}: historical replay surface lacks an origin anchor")
        bindings.append(
            {
                "path": path,
                "custody": custody,
                "sha256": hashlib.sha256(payload).hexdigest(),
                "git_blob_sha1": blob,
            }
        )
    return bindings


def _verify_revision_file_bindings(
    version_id: str, revision: str, version: dict
) -> None:
    files = version.get("normative_files")
    if not isinstance(files, list):
        fail(f"{version_id}: predecessor revision carries malformed normative files")
    for item in files:
        if not isinstance(item, dict):
            fail(f"{version_id}: predecessor revision carries a malformed file row")
        path = item.get("path")
        expected_blob = item.get("git_blob_sha1")
        if not isinstance(path, str) or not isinstance(expected_blob, str):
            fail(f"{version_id}: predecessor revision carries malformed file custody")
        actual_blob = (
            _git_output("rev-parse", f"{revision}:{path}")
            .decode("ascii", errors="replace")
            .strip()
        )
        if actual_blob != expected_blob:
            fail(
                f"{version_id}: predecessor revision binds {path} to Git blob "
                f"{actual_blob}, not registered blob {expected_blob}"
            )


def _premise_ids(version_id: str, payload: bytes) -> set[str]:
    try:
        data = strict_json.loads(payload.decode("utf-8"))
    except (
        UnicodeDecodeError,
        json.JSONDecodeError,
        strict_json.DuplicateKeyError,
    ) as error:
        fail(f"{version_id}: pinned premise register is invalid JSON: {error}")
    rows = data.get("rows") if isinstance(data, dict) else None
    if not isinstance(rows, list) or not rows:
        fail(f"{version_id}: pinned premise register must contain rows")
    ids: list[str] = []
    for row in rows:
        premise_id = row.get("id") if isinstance(row, dict) else None
        if not isinstance(premise_id, str) or not premise_id:
            fail(f"{version_id}: pinned premise register has a malformed id")
        ids.append(premise_id)
    if len(ids) != len(set(ids)):
        fail(f"{version_id}: pinned premise register has duplicate ids")
    return set(ids)


def validate(data: dict) -> list[dict]:
    if not isinstance(data, dict) or set(data) != TOP_KEYS:
        fail(f"top-level keys must equal {sorted(TOP_KEYS)}")
    if data["schema"] != SCHEMA or data["issue"] != ISSUE:
        fail(f"schema and issue must equal {SCHEMA} and {ISSUE}")
    _nonempty("policy", data["policy"])

    versions = data["versions"]
    if not isinstance(versions, list) or not versions:
        fail("versions must be a nonempty list")
    expected_ids = [f"AV-{index}" for index in range(len(versions))]
    if [
        entry.get("id") for entry in versions if isinstance(entry, dict)
    ] != expected_ids:
        fail("version ids must be contiguous AV-0, AV-1, ... in order")
    if data["current_version"] != expected_ids[-1]:
        fail("current_version must name the final version")

    anchors = data["version_anchors"]
    if not isinstance(anchors, list) or not anchors:
        fail("version_anchors must be a nonempty contiguous prefix")
    anchor_ids = [anchor.get("id") if isinstance(anchor, dict) else None for anchor in anchors]
    if anchor_ids != expected_ids[: len(anchors)]:
        fail("version_anchors must cover a contiguous prefix starting at AV-0")
    if len(anchors) > len(versions):
        fail("version_anchors cannot name nonexistent versions")
    if len(versions) - len(anchors) > 1:
        fail("at most one final architecture version may await its origin anchor")
    head_revision = _git_output("rev-parse", "HEAD").decode().strip()
    for index, anchor in enumerate(anchors):
        version_id = expected_ids[index]
        if not isinstance(anchor, dict) or set(anchor) != ANCHOR_KEYS:
            fail(f"{version_id}: malformed version anchor")
        revision = anchor["origin_revision"]
        if not isinstance(revision, str) or not GIT_REVISION_RE.fullmatch(revision):
            fail(f"{version_id}: origin_revision must be a full Git commit")
        if not _git_is_ancestor(revision, head_revision):
            fail(f"{version_id}: origin_revision must be an ancestor of HEAD")
        if not _git_is_first_parent_ancestor(revision, head_revision):
            fail(
                f"{version_id}: origin_revision must lie on HEAD's first-parent "
                "custody history"
            )
        if version_id == "AV-0" and ROOT == REGISTER_PATH.parents[1]:
            if (
                revision != AV0_FIRST_APPEARANCE_REVISION
                or anchor["record_sha256"] != AV0_FIRST_APPEARANCE_RECORD_SHA256
            ):
                fail(
                    "AV-0: origin anchor must equal its audited first appearance "
                    f"{AV0_FIRST_APPEARANCE_REVISION} / "
                    f"{AV0_FIRST_APPEARANCE_RECORD_SHA256}"
                )
        origin = _pinned_register(revision)
        if origin.get("schema") not in {LEGACY_ROOT_SCHEMA, SCHEMA}:
            fail(f"{version_id}: origin revision has an incompatible schema")
        origin_versions = origin.get("versions")
        if (
            origin.get("current_version") != version_id
            or not isinstance(origin_versions, list)
            or len(origin_versions) != index + 1
            or not isinstance(origin_versions[index], dict)
        ):
            fail(f"{version_id}: origin revision must end exactly at {version_id}")
        origin_record = origin_versions[index]
        origin_payload = json.dumps(
            origin_record,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
        digest = hashlib.sha256(origin_payload).hexdigest()
        if anchor["record_sha256"] != digest:
            fail(f"{version_id}: anchor record_sha256 drifted")
        current_payload = json.dumps(
            versions[index],
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
        if current_payload != origin_payload:
            fail(f"{version_id}: record differs from its origin anchor")
        parents = _parents(revision)
        if version_id == "AV-0":
            if len(parents) > 1:
                fail("AV-0: origin revision cannot be a merge commit")
            parent_register = (
                _optional_pinned_register(parents[0]) if parents else None
            )
            if parent_register is not None:
                fail("AV-0: origin parent already contains an architecture register")
        else:
            if origin_record.get("status") == "retired":
                fail(f"{version_id}: origin cannot introduce a retired current tip")
            if (
                origin_record.get("promotion_status")
                != PENDING_TIP_PROMOTION_STATUS
            ):
                fail(
                    f"{version_id}: origin record must preserve the canonical "
                    "pre-anchor promotion-ineligible state"
                )
            origin_parent = _single_parent(revision, version_id)
            if origin_parent != versions[index]["predecessor_register_revision"]:
                fail(
                    f"{version_id}: origin parent must equal "
                    "predecessor_register_revision"
                )
            parent_register = _optional_pinned_register(origin_parent)
            parent_versions = (
                parent_register.get("versions")
                if isinstance(parent_register, dict)
                else None
            )
            if (
                parent_register is None
                or parent_register.get("schema") != SCHEMA
                or not isinstance(parent_versions, list)
                or len(parent_versions) != index
            ):
                fail(
                    f"{version_id}: origin revision is not the version's first "
                    "append-only appearance"
                )
        _verify_revision_file_bindings(version_id, revision, origin_record)

    for version_index, version in enumerate(versions):
        version_id = version["id"] if isinstance(version, dict) else "<malformed>"
        if not isinstance(version, dict) or set(version) != VERSION_KEYS:
            fail(f"{version_id}: keys must equal {sorted(VERSION_KEYS)}")
        if not VERSION_RE.fullmatch(version_id):
            fail(f"{version_id}: malformed version id")
        _nonempty(f"{version_id}.created_on", version["created_on"])
        if version["status"] not in STATUS_ENUM:
            fail(f"{version_id}: unknown status {version['status']!r}")
        is_current = version_index == len(versions) - 1
        if is_current and version["status"] == "retired":
            fail(f"{version_id}: current version cannot be retired")
        is_anchored = version_id in set(anchor_ids)
        if (
            is_current
            and not is_anchored
            and version["promotion_status"] != PENDING_TIP_PROMOTION_STATUS
        ):
            fail(f"{version_id}: an unanchored current tip is promotion-ineligible")

        recorded_snapshot = version["snapshot_sha256"]
        actual_snapshot = (
            legacy_snapshot_sha256(version)
            if version_index == 0
            else snapshot_sha256(version)
        )
        if not isinstance(recorded_snapshot, str) or not DIGEST_RE.fullmatch(
            recorded_snapshot
        ):
            fail(f"{version_id}: invalid snapshot_sha256")
        if recorded_snapshot != actual_snapshot:
            fail(f"{version_id}: immutable snapshot digest mismatch: {actual_snapshot}")

        predecessor_revision = version["predecessor_register_revision"]
        if version_index == 0:
            if predecessor_revision is not None:
                fail("AV-0: predecessor_register_revision must be null")
        else:
            if not isinstance(
                predecessor_revision, str
            ) or not GIT_REVISION_RE.fullmatch(predecessor_revision):
                fail(
                    f"{version_id}: predecessor_register_revision must be a "
                    "full 40-character Git revision"
                )
            predecessor_id = f"AV-{version_index - 1}"
            pinned = _pinned_register(predecessor_revision)
            if pinned.get("schema") != SCHEMA:
                fail(f"{version_id}: predecessor revision register must use {SCHEMA}")
            if pinned.get("current_version") != predecessor_id:
                fail(
                    f"{version_id}: predecessor revision must have "
                    f"current_version {predecessor_id}"
                )
            pinned_versions = pinned.get("versions")
            if (
                not isinstance(pinned_versions, list)
                or len(pinned_versions) != version_index
            ):
                fail(
                    f"{version_id}: predecessor revision must contain exactly "
                    f"AV-0 through {predecessor_id}"
                )
            pinned_predecessor = pinned_versions[-1]
            predecessor = versions[version_index - 1]
            if not isinstance(pinned_predecessor, dict):
                fail(f"{version_id}: predecessor revision carries a malformed version")
            pinned_snapshot = pinned_predecessor.get("snapshot_sha256")
            try:
                pinned_actual_snapshot = (
                    legacy_snapshot_sha256(pinned_predecessor)
                    if version_index - 1 == 0
                    else snapshot_sha256(pinned_predecessor)
                )
            except (KeyError, TypeError) as error:
                fail(
                    f"{version_id}: predecessor revision carries a malformed "
                    f"{predecessor_id} snapshot ({error})"
                )
            if pinned_snapshot != pinned_actual_snapshot:
                fail(
                    f"{version_id}: predecessor revision carries an invalid "
                    f"{predecessor_id} snapshot digest"
                )
            if _snapshot_payload(
                pinned_predecessor, legacy=version_index - 1 == 0
            ) != _snapshot_payload(predecessor, legacy=version_index - 1 == 0):
                fail(
                    f"{version_id}: {predecessor_id} differs from its pinned "
                    "predecessor revision; architecture history is append-only"
                )
            _verify_revision_file_bindings(
                predecessor_id, predecessor_revision, pinned_predecessor
            )

        basis = version["basis"]
        if not isinstance(basis, dict) or set(basis) != BASIS_KEYS:
            fail(f"{version_id}.basis keys must equal {sorted(BASIS_KEYS)}")
        if basis["core_axioms"] != ["A1", "A2", "A3"]:
            fail(f"{version_id}: core_axioms must be exactly A1, A2, A3")
        if basis["proposed_fundamental_numerical_parameters"] != ["P", "N"]:
            fail(f"{version_id}: proposed fundamental parameters must be P and N")
        _nonempty(f"{version_id}.closure_principle", basis["closure_principle"])

        files = version["normative_files"]
        if not isinstance(files, list) or not files:
            fail(f"{version_id}: normative_files must be nonempty")
        paths: list[str] = []
        normative_payloads: dict[str, bytes] = {}
        for item in files:
            if not isinstance(item, dict) or set(item) != FILE_KEYS:
                fail(f"{version_id}: malformed normative file entry")
            path = _nonempty(f"{version_id}.normative_files.path", item["path"])
            if path.startswith("/") or ".." in path.split("/"):
                fail(f"{version_id}: normative paths must be repo-relative")
            if path in paths:
                fail(f"{version_id}: duplicate normative path {path}")
            paths.append(path)
            _nonempty(f"{version_id}.{path}.role", item["role"])
            digest = item["sha256"]
            if not isinstance(digest, str) or not DIGEST_RE.fullmatch(digest):
                fail(f"{version_id}: invalid sha256 for {path}")
            git_blob = item["git_blob_sha1"]
            if not isinstance(git_blob, str) or not GIT_BLOB_RE.fullmatch(git_blob):
                fail(f"{version_id}: invalid git_blob_sha1 for {path}")
            if is_current:
                resolved = ROOT / path
                if not resolved.is_file():
                    fail(f"{version_id}: missing current normative file {path}")
                payload = resolved.read_bytes()
            else:
                payload = _historical_blob(version_id, path, git_blob)
            actual_blob = _git_blob_sha1(payload)
            if git_blob != actual_blob:
                fail(f"{version_id}: Git blob mismatch for {path}: {actual_blob}")
            actual = hashlib.sha256(payload).hexdigest()
            if digest != actual:
                fail(f"{version_id}: digest mismatch for {path}: {actual}")
            normative_payloads[path] = payload

        premise_payload = normative_payloads.get(PREMISE_REGISTER_RELATIVE)
        if premise_payload is None:
            fail(
                f"{version_id}: normative_files must include "
                f"{PREMISE_REGISTER_RELATIVE}"
            )
        premise_ids = _premise_ids(version_id, premise_payload)

        decisions = version["protocol_decisions"]
        if not isinstance(decisions, list) or not decisions:
            fail(f"{version_id}: protocol_decisions must be nonempty")
        for index, decision in enumerate(decisions, start=1):
            if not isinstance(decision, dict) or set(decision) != DECISION_KEYS:
                fail(f"{version_id}: malformed protocol decision")
            if decision["id"] != f"PD-{index:02d}":
                fail(f"{version_id}: protocol decision ids must be contiguous")
            _nonempty(f"{version_id}.{decision['id']}.decision", decision["decision"])
            premises = _string_list(
                f"{version_id}.{decision['id']}.premises",
                decision["premises"],
                allow_empty=True,
            )
            unknown = set(premises) - premise_ids
            if unknown:
                fail(
                    f"{version_id}.{decision['id']}: unknown premises {sorted(unknown)}"
                )

        _nonempty(f"{version_id}.common_world_status", version["common_world_status"])
        _nonempty(f"{version_id}.promotion_status", version["promotion_status"])
        _string_list(
            f"{version_id}.invalidation_triggers", version["invalidation_triggers"]
        )
        _string_list(
            f"{version_id}.replay_surfaces", version["replay_surfaces"]
        )
        replay_surface_bindings(data, version, is_current=is_current)
    _validate_committed_history_guard(versions, anchors)
    _validate_committed_head_anchors(anchors)
    return versions


def render(data: dict, versions: list[dict]) -> str:
    lines = [
        "# Architecture version register",
        "",
        "Generated by `tools/build_architecture_versions.py` from "
        "`tracking/architecture_versions.json`; edit the JSON and regenerate. "
        "Issue [#741](https://github.com/FloatingPragma/observer-patch-holography/issues/741) "
        "established the bootstrap; this register is the durable version and "
        "promotion-custody surface.",
        "",
        data["policy"],
        "",
        "An architecture version uses a two-commit custody protocol: first commit "
        "the complete promotion-ineligible version record, then declare its origin "
        "anchor in a later commit. An anchor is not operative while it exists only "
        "in the worktree. The sole bootstrap exception is the hard-pinned AV-0 "
        "v2-to-v3 migration, which already binds its audited first appearance.",
        "",
        "",
        "Version origin anchors:",
        "",
        *[
            f"- `{anchor['id']}`: origin `{anchor['origin_revision']}`, complete-record SHA-256 `{anchor['record_sha256']}`."
            for anchor in data["version_anchors"]
        ],
    ]
    for version in versions:
        basis = version["basis"]
        lines.extend(
            [
                "",
                f"## {version['id']} ({version['status']})",
                "",
                f"Created: `{version['created_on']}`.",
                f"Immutable snapshot SHA-256: `{version['snapshot_sha256']}`.",
                (
                    "Predecessor register revision: none (root version)."
                    if version["predecessor_register_revision"] is None
                    else "Predecessor register revision: "
                    f"`{version['predecessor_register_revision']}`."
                ),
                (
                    "Normative-byte custody: current live files, with their Git "
                    "blob ids precomputed for historical replay after commit."
                    if version["id"] == data["current_version"]
                    else "Normative-byte custody: immutable Git blobs; live files "
                    "may have changed or disappeared."
                ),
                "",
                f"Basis: {', '.join(basis['core_axioms'])}; closure principle: "
                f"{basis['closure_principle']}; proposed fundamental numerical "
                f"parameters: {', '.join(basis['proposed_fundamental_numerical_parameters'])}.",
                "",
                "### Normative files",
                "",
                "| Path | Role | SHA-256 | Git blob SHA-1 |",
                "| --- | --- | --- | --- |",
            ]
        )
        for item in version["normative_files"]:
            lines.append(
                f"| `{item['path']}` | {item['role']} | `{item['sha256']}` | "
                f"`{item['git_blob_sha1']}` |"
            )
        lines.extend(["", "### Protocol decisions", ""])
        for decision in version["protocol_decisions"]:
            premises = ", ".join(decision["premises"]) or "none"
            lines.append(
                f"- **{decision['id']}** {decision['decision']} Premises: {premises}."
            )
        lines.extend(
            [
                "",
                f"Common-world status: {version['common_world_status']}.",
                "",
                f"Promotion status: {version['promotion_status']}",
                "",
                "Invalidation triggers:",
                "",
                *[f"- {item}" for item in version["invalidation_triggers"]],
                "",
                "Replay surfaces:",
                "",
                *[
                    f"- `{item['path']}` at `{item['custody']}`: SHA-256 "
                    f"`{item['sha256']}`, Git blob `{item['git_blob_sha1']}`."
                    for item in replay_surface_bindings(
                        data,
                        version,
                        is_current=version["id"] == data["current_version"],
                    )
                ],
            ]
        )
    lines.append("")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)
    data = load_json(REGISTER_PATH)
    versions = validate(data)
    rendered = render(data, versions).encode("utf-8")
    if args.check:
        committed = SURFACE_PATH.read_bytes() if SURFACE_PATH.is_file() else b""
        if committed != rendered:
            print(
                "architecture versions: docs/ARCHITECTURE_VERSION_REGISTER.md "
                "is stale; run python tools/build_architecture_versions.py",
                file=sys.stderr,
            )
            return 1
        print(
            f"architecture versions: {len(versions)} version(s), current {data['current_version']}"
        )
        return 0
    SURFACE_PATH.write_bytes(rendered)
    print(f"architecture versions: wrote {SURFACE_PATH.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
