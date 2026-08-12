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

SCHEMA = "oph.architecture_version_register.v2"
ISSUE = 741
VERSION_RE = re.compile(r"^AV-(0|[1-9][0-9]*)$")
DIGEST_RE = re.compile(r"^[0-9a-f]{64}$")
GIT_BLOB_RE = re.compile(r"^[0-9a-f]{40}$")
GIT_REVISION_RE = re.compile(r"^[0-9a-f]{40}$")
STATUS_ENUM = {"exploratory_uninhabited", "inhabited_conditional", "retired"}
PREMISE_REGISTER_RELATIVE = "tracking/premise_register.json"

TOP_KEYS = {"schema", "issue", "current_version", "policy", "versions"}
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
SNAPSHOT_KEYS = (
    "id",
    "created_on",
    "basis",
    "normative_files",
    "protocol_decisions",
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


def _snapshot_payload(version: dict) -> bytes:
    """Canonical immutable portion of one architecture version."""

    projection = {key: version[key] for key in SNAPSHOT_KEYS}
    return json.dumps(
        projection,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")


def snapshot_sha256(version: dict) -> str:
    return hashlib.sha256(_snapshot_payload(version)).hexdigest()


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

        recorded_snapshot = version["snapshot_sha256"]
        actual_snapshot = snapshot_sha256(version)
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
                pinned_actual_snapshot = snapshot_sha256(pinned_predecessor)
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
            if _snapshot_payload(pinned_predecessor) != _snapshot_payload(predecessor):
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
        replay = _string_list(
            f"{version_id}.replay_surfaces", version["replay_surfaces"]
        )
        for path in replay:
            if not (ROOT / path).is_file():
                fail(f"{version_id}: replay surface missing: {path}")
    return versions


def render(data: dict, versions: list[dict]) -> str:
    lines = [
        "# Architecture version register",
        "",
        "Generated by `tools/build_architecture_versions.py` from "
        "`tracking/architecture_versions.json`; edit the JSON and regenerate. "
        "Issue [#741](https://github.com/FloatingPragma/observer-patch-holography/issues/741) "
        "owns version and promotion custody.",
        "",
        data["policy"],
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
                *[f"- `{item}`" for item in version["replay_surfaces"]],
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
