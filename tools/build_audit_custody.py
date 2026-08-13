#!/usr/bin/env python3
"""Validate and render the versioned V3 promotion-audit custody index.

Audit evidence is checked at the Git revision that was actually reviewed, not
against mutable live paths.  This makes an audit record append-only while
still allowing the reviewed surfaces to evolve in later architecture
versions.  Observation-ledger promotion pointers are validated by the
architecture replay builder, which consumes this register.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

import strict_json

ROOT = Path(__file__).resolve().parents[1]
REGISTER_PATH = ROOT / "tracking" / "audit_custody.json"
SURFACE_PATH = ROOT / "docs" / "AUDIT_CUSTODY.md"

SCHEMA = "oph.audit_custody.v2"
ISSUE = 738
AUDIT_ID_RE = re.compile(r"^AUD-[A-Z0-9-]+$")
COMMIT_RE = re.compile(r"^[0-9a-f]{40}$")
DIGEST_RE = re.compile(r"^[0-9a-f]{64}$")
BLOB_RE = re.compile(r"^[0-9a-f]{40}$")
ROW_RE = re.compile(r"^OL-[A-N][1-9]$")
DATE_RE = re.compile(r"^20[0-9]{2}-[01][0-9]-[0-3][0-9]$")

TOP_KEYS = {"schema", "issue", "record_anchors", "policy", "records"}
ANCHOR_KEYS = {"id", "origin_revision", "record_sha256"}
RECORD_KEYS = {
    "id",
    "provenance_kind",
    "started_on",
    "completed_on",
    "baseline_commit",
    "reviewed_commit",
    "repair_commit",
    "audit_class",
    "scope",
    "reviewers",
    "reviewed_rows",
    "promoted_rows",
    "findings",
    "artifact_pins",
    "qualifies_for",
    "does_not_qualify_for",
    "limitations",
    "source_record_pin",
}
PROVENANCE_KINDS = {"legacy_migration", "native"}
REVIEWER_KEYS = {"name", "task", "model", "role"}
MODEL_DISCLOSURE = "Codex GPT-5 family; exact deployment identifier was not exposed"
LEGACY_SOURCE_PINS = {
    "AUD-V3-2026-08-12": {
        "revision": "9c1dd24728218966705b82a900ec2ba966d57394",
        "path": "tracking/audit_register.json",
        "bytes": 9441,
        "sha256": "44b8ac5ec087c28bef1eca0f7df6a74ace891b5411967f241d1fa186cdcc509c",
        "git_blob_sha1": "2b623d3f1d9f3465826111b96f6251476e41341e",
    }
}
FINDING_KEYS = {"id", "severity", "disposition", "summary", "owner_issues"}
PIN_KEYS = {"revision", "path", "bytes", "sha256", "git_blob_sha1"}
SOURCE_PIN_KEYS = {"revision", "path", "bytes", "sha256", "git_blob_sha1"}
SEVERITIES = {"critical", "high", "medium", "low", "boundary"}
DISPOSITIONS = {"fixed", "tracked_open", "verified_boundary", "no_finding"}


def fail(message: str) -> None:
    raise SystemExit(f"audit custody: {message}")


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


def _git_bytes(*args: str) -> bytes:
    result = subprocess.run(
        ["git", *args],
        cwd=ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if result.returncode != 0:
        detail = result.stderr.decode("utf-8", errors="replace").strip()
        fail(f"git {' '.join(args)} failed: {detail or 'object unavailable'}")
    return result.stdout


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
        for item in _git_bytes("rev-list", "--first-parent", descendant)
        .decode()
        .splitlines()
        if item
    }
    return ancestor in revisions


def _nonempty(where: str, value: object) -> str:
    if not isinstance(value, str) or not value.strip():
        fail(f"{where} must be a nonempty string")
    return value


def _strings(where: str, value: object, *, allow_empty: bool = False) -> list[str]:
    if not isinstance(value, list) or (not value and not allow_empty):
        fail(f"{where} must be a {'possibly empty ' if allow_empty else 'nonempty '}list")
    if any(not isinstance(item, str) or not item.strip() for item in value):
        fail(f"{where} entries must be nonempty strings")
    if len(value) != len(set(value)):
        fail(f"{where} must be duplicate-free")
    return value


def _historical_json(revision: str, path: str) -> dict:
    raw = _git_bytes("show", f"{revision}:{path}")
    try:
        payload = strict_json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError, strict_json.DuplicateKeyError) as error:
        fail(f"{revision}:{path} is not strict JSON: {error}")
    if not isinstance(payload, dict):
        fail(f"{revision}:{path} must contain a JSON object")
    return payload


def _optional_historical_json(revision: str, path: str) -> dict | None:
    result = subprocess.run(
        ["git", "show", f"{revision}:{path}"],
        cwd=ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if result.returncode != 0:
        return None
    try:
        payload = strict_json.loads(result.stdout.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError, strict_json.DuplicateKeyError) as error:
        fail(f"{revision}:{path} is not strict JSON: {error}")
    if not isinstance(payload, dict):
        fail(f"{revision}:{path} must contain a JSON object")
    return payload


def _git_parents(revision: str) -> list[str]:
    return _git_bytes("show", "-s", "--format=%P", revision).decode().strip().split()


def _validate_committed_history_guard(
    data: dict, record_ids: list[str], anchor_ids: list[str]
) -> None:
    """Reject record rewrites and anchor shrinkage visible in Git history."""

    head = _git_bytes("rev-parse", "HEAD").decode().strip()
    revisions = (
        _git_bytes(
            "log",
            "--first-parent",
            "--full-history",
            "--format=%H",
            head,
            "--",
            "tracking/audit_custody.json",
        )
        .decode()
        .strip()
        .splitlines()
    )
    current_by_id = {record["id"]: record for record in data["records"]}
    current_anchors = {
        anchor["id"]: anchor
        for anchor in data["record_anchors"]
        if isinstance(anchor, dict) and isinstance(anchor.get("id"), str)
    }
    saw_schema = False
    left_schema_span = False
    for revision in revisions:
        historical = _optional_historical_json(
            revision, "tracking/audit_custody.json"
        )
        if historical is None or historical.get("schema") != SCHEMA:
            if saw_schema:
                left_schema_span = True
            continue
        if left_schema_span:
            fail("audit custody cannot disappear and later reappear in Git history")
        saw_schema = True
        historical_records = historical.get("records")
        historical_anchors = historical.get("record_anchors")
        if not isinstance(historical_records, list) or not isinstance(
            historical_anchors, list
        ):
            fail(f"{revision}: historical audit custody is malformed")
        historical_ids = [
            record.get("id") if isinstance(record, dict) else None
            for record in historical_records
        ]
        if record_ids[: len(historical_ids)] != historical_ids:
            fail("audit record history is not append-only")
        for historical_record in historical_records:
            audit_id = historical_record["id"]
            if _canonical_record_bytes(current_by_id[audit_id]) != _canonical_record_bytes(
                historical_record
            ):
                fail(f"{audit_id}: record rewrites committed audit history")
        for historical_anchor in historical_anchors:
            audit_id = historical_anchor.get("id")
            if current_anchors.get(audit_id) != historical_anchor:
                fail(f"{audit_id}: origin anchor rewrites committed history")


def _validate_committed_head_anchors(anchors: list[dict]) -> None:
    """Reject worktree-only native anchors as non-operative declarations."""

    head = _git_bytes("rev-parse", "HEAD").decode().strip()
    committed = _optional_historical_json(head, "tracking/audit_custody.json")
    committed_anchors = (
        committed.get("record_anchors")
        if isinstance(committed, dict) and committed.get("schema") == SCHEMA
        else []
    )
    if not isinstance(committed_anchors, list):
        fail("committed HEAD audit anchors are malformed")
    if len(anchors) > len(committed_anchors) or anchors != committed_anchors[: len(anchors)]:
        fail("every operative audit anchor must already be committed at HEAD")


def _canonical_record_bytes(record: dict) -> bytes:
    return json.dumps(
        record,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")


def reviewed_row_index(record: dict) -> dict[str, dict]:
    """Return the exact observation rows stored at an audit's repair commit."""

    payload = _historical_json(
        record["repair_commit"], "tracking/observation_ledger.json"
    )
    rows = payload.get("rows")
    if not isinstance(rows, list):
        fail(f"{record['id']}: pinned observation ledger has no rows")
    return {
        row["id"]: row
        for row in rows
        if isinstance(row, dict) and isinstance(row.get("id"), str)
    }


def validate(data: dict) -> list[dict]:
    if not isinstance(data, dict) or set(data) != TOP_KEYS:
        fail(f"top-level keys must equal {sorted(TOP_KEYS)}")
    if data["schema"] != SCHEMA or data["issue"] != ISSUE:
        fail(f"schema and issue must equal {SCHEMA} and {ISSUE}")
    _nonempty("policy", data["policy"])
    records = data["records"]
    if not isinstance(records, list) or not records:
        fail("records must be a nonempty list")

    seen_ids: set[str] = set()
    native_ids: list[str] = []
    seen_native = False
    for index, record in enumerate(records):
        where = f"records[{index}]"
        if not isinstance(record, dict) or set(record) != RECORD_KEYS:
            fail(f"{where}: keys must equal {sorted(RECORD_KEYS)}")
        audit_id = record["id"]
        if not isinstance(audit_id, str) or not AUDIT_ID_RE.fullmatch(audit_id):
            fail(f"{where}: malformed audit id")
        if audit_id in seen_ids:
            fail(f"duplicate audit id {audit_id}")
        seen_ids.add(audit_id)
        provenance = record["provenance_kind"]
        if provenance not in PROVENANCE_KINDS:
            fail(f"{audit_id}: unknown provenance_kind")
        if provenance == "native":
            seen_native = True
            native_ids.append(audit_id)
        elif seen_native:
            fail("legacy-migration records must precede native records")

    anchors = data["record_anchors"]
    if not isinstance(anchors, list):
        fail("record_anchors must be a list")
    anchor_ids = [
        anchor.get("id") if isinstance(anchor, dict) else None for anchor in anchors
    ]
    if anchor_ids != native_ids[: len(anchors)] or len(anchors) > len(native_ids):
        fail("record_anchors must cover a contiguous prefix of native records")
    if len(native_ids) - len(anchors) > 1:
        fail("at most one final native audit may await its origin anchor")
    record_ids = [record["id"] for record in records]
    for anchor in anchors:
        audit_id = anchor.get("id") if isinstance(anchor, dict) else "<malformed>"
        if not isinstance(anchor, dict) or set(anchor) != ANCHOR_KEYS:
            fail(f"{audit_id}: malformed record anchor")
        revision = anchor["origin_revision"]
        digest = anchor["record_sha256"]
        if not isinstance(revision, str) or not COMMIT_RE.fullmatch(revision):
            fail(f"{audit_id}: origin_revision must be a full Git commit")
        if not isinstance(digest, str) or not DIGEST_RE.fullmatch(digest):
            fail(f"{audit_id}: malformed origin record SHA-256")
        if _git_bytes("cat-file", "-t", revision).strip() != b"commit":
            fail(f"{audit_id}: origin_revision is not a commit")
        if not _git_is_ancestor(revision, _git_bytes("rev-parse", "HEAD").decode().strip()):
            fail(f"{audit_id}: origin_revision must be an ancestor of HEAD")
        if not _git_is_first_parent_ancestor(
            revision, _git_bytes("rev-parse", "HEAD").decode().strip()
        ):
            fail(
                f"{audit_id}: origin_revision must lie on HEAD's first-parent "
                "custody history"
            )
        origin = _historical_json(revision, "tracking/audit_custody.json")
        if origin.get("schema") != SCHEMA:
            fail(f"{audit_id}: origin register must use {SCHEMA}")
        origin_records = origin.get("records")
        record_index = record_ids.index(audit_id)
        if (
            not isinstance(origin_records, list)
            or len(origin_records) != record_index + 1
            or [
                item.get("id") if isinstance(item, dict) else None
                for item in origin_records
            ]
            != record_ids[: record_index + 1]
        ):
            fail(f"{audit_id}: origin revision must end exactly at this audit")
        origin_record = origin_records[-1]
        origin_payload = _canonical_record_bytes(origin_record)
        if hashlib.sha256(origin_payload).hexdigest() != digest:
            fail(f"{audit_id}: anchor record_sha256 drifted")
        if origin_payload != _canonical_record_bytes(records[record_index]):
            fail(f"{audit_id}: record differs from its origin anchor")
        if not _git_is_ancestor(records[record_index]["repair_commit"], revision):
            fail(f"{audit_id}: repair_commit must be an ancestor of origin_revision")
        parents = _git_parents(revision)
        if len(parents) != 1:
            fail(f"{audit_id}: native origin must be a single-parent commit")
        parent = _optional_historical_json(
            parents[0], "tracking/audit_custody.json"
        )
        expected_prior_ids = record_ids[:record_index]
        if parent is None:
            if expected_prior_ids:
                fail(f"{audit_id}: origin parent omits prior audit history")
        else:
            parent_records = parent.get("records")
            parent_ids = (
                [item.get("id") for item in parent_records]
                if isinstance(parent_records, list)
                else None
            )
            if parent.get("schema") != SCHEMA or parent_ids != expected_prior_ids:
                fail(
                    f"{audit_id}: origin_revision is not the record's first "
                    "append-only appearance"
                )

    anchored_native_ids = set(anchor_ids)
    validated_records: list[dict] = []
    for index, record in enumerate(records):
        audit_id = record["id"]
        where = audit_id

        for key in ("started_on", "completed_on"):
            if not isinstance(record[key], str) or not DATE_RE.fullmatch(record[key]):
                fail(f"{where}: malformed {key}")
        for key in ("baseline_commit", "reviewed_commit", "repair_commit"):
            revision = record[key]
            if not isinstance(revision, str) or not COMMIT_RE.fullmatch(revision):
                fail(f"{where}: malformed {key}")
            if _git_bytes("cat-file", "-t", revision).strip() != b"commit":
                fail(f"{where}: {key} is not a commit")
        if not _git_is_ancestor(record["baseline_commit"], record["reviewed_commit"]):
            fail(f"{where}: baseline_commit must be an ancestor of reviewed_commit")
        if not _git_is_ancestor(record["reviewed_commit"], record["repair_commit"]):
            fail(f"{where}: reviewed_commit must be an ancestor of repair_commit")
        for key in ("audit_class", "scope"):
            _nonempty(f"{where}.{key}", record[key])

        provenance = record["provenance_kind"]
        source_pin = record["source_record_pin"]
        source_record: dict | None = None
        if provenance == "legacy_migration":
            if not isinstance(source_pin, dict) or set(source_pin) != SOURCE_PIN_KEYS:
                fail(f"{where}: malformed legacy source_record_pin")
            expected_source_pin = LEGACY_SOURCE_PINS.get(audit_id)
            if expected_source_pin is None or source_pin != expected_source_pin:
                fail(
                    f"{where}: legacy source_record_pin must equal its audited "
                    "immutable migration source"
                )
            source_revision = source_pin["revision"]
            source_path = source_pin["path"]
            if (
                not isinstance(source_revision, str)
                or not COMMIT_RE.fullmatch(source_revision)
                or source_path != "tracking/audit_register.json"
            ):
                fail(
                    f"{where}: legacy source_record_pin must bind the historical "
                    "audit register"
                )
            source_bytes = _git_bytes("show", f"{source_revision}:{source_path}")
            source_blob = _git_bytes(
                "rev-parse", f"{source_revision}:{source_path}"
            ).decode().strip()
            source_digest = hashlib.sha256(source_bytes).hexdigest()
            if (
                source_pin["bytes"],
                source_pin["sha256"],
                source_pin["git_blob_sha1"],
            ) != (len(source_bytes), source_digest, source_blob):
                fail(f"{where}: source audit-record pin drift")
            source_data = _historical_json(source_revision, source_path)
            source_records = source_data.get("records")
            source_record = (
                next(
                    (
                        item
                        for item in source_records
                        if isinstance(item, dict) and item.get("id") == audit_id
                    ),
                    None,
                )
                if isinstance(source_records, list)
                else None
            )
            if source_record is None:
                fail(f"{where}: source audit record is absent")
            immutable_pairs = {
                "started_on": "started_on",
                "completed_on": "completed_on",
                "baseline_commit": "baseline_commit",
                "reviewed_commit": "audited_head",
                "repair_commit": "repair_commit",
                "audit_class": "audit_class",
                "scope": "scope",
                "reviewed_rows": "reviewed_observation_rows",
                "promoted_rows": "attained_rows_reviewed",
                "findings": "findings",
                "qualifies_for": "qualifies_for",
                "does_not_qualify_for": "does_not_qualify_for",
                "limitations": "limitations",
            }
            for current_key, source_key in immutable_pairs.items():
                if record[current_key] != source_record.get(source_key):
                    fail(
                        f"{where}: {current_key} rewrites the pinned source audit record"
                    )
        elif source_pin is not None:
            fail(f"{where}: native audits must set source_record_pin to null")

        reviewers = record["reviewers"]
        if not isinstance(reviewers, list) or not reviewers:
            fail(f"{where}: reviewers must be nonempty")
        tasks: set[str] = set()
        for reviewer in reviewers:
            if not isinstance(reviewer, dict) or set(reviewer) != REVIEWER_KEYS:
                fail(f"{where}: malformed reviewer")
            for key in REVIEWER_KEYS:
                _nonempty(f"{where}.reviewer.{key}", reviewer[key])
            if provenance == "legacy_migration" and reviewer["model"] != MODEL_DISCLOSURE:
                fail(
                    f"{where}: reviewer model must use the fixed disclosure for "
                    "legacy audits whose exact deployment id was not recorded"
                )
            if reviewer["task"] in tasks:
                fail(f"{where}: reviewer tasks must be distinct")
            tasks.add(reviewer["task"])
        if source_record is not None:
            source_reviewers = source_record.get("reviewers")
            if not isinstance(source_reviewers, list) or len(source_reviewers) != len(reviewers):
                fail(f"{where}: reviewer inventory rewrites the pinned source audit record")
            for current, source in zip(reviewers, source_reviewers, strict=True):
                if any(
                    current.get(key) != source.get(key)
                    for key in ("name", "task", "role")
                ):
                    fail(
                        f"{where}: reviewer identity rewrites the pinned source audit record"
                    )

        reviewed = _strings(f"{where}.reviewed_rows", record["reviewed_rows"])
        promoted = _strings(f"{where}.promoted_rows", record["promoted_rows"], allow_empty=True)
        if any(not ROW_RE.fullmatch(row_id) for row_id in reviewed + promoted):
            fail(f"{where}: malformed observation row id")
        if set(promoted) - set(reviewed):
            fail(f"{where}: promoted_rows must be a subset of reviewed_rows")
        if promoted and "attained_status_review" not in record["qualifies_for"]:
            fail(f"{where}: promoted rows require attained_status_review qualification")
        if (
            "attained_status_review" in record["qualifies_for"]
            and "attained_status_review" in record["does_not_qualify_for"]
        ):
            fail(f"{where}: qualification and exclusion contradict")

        pins = record["artifact_pins"]
        if not isinstance(pins, list) or not pins:
            fail(f"{where}: artifact_pins must be nonempty")
        pin_keys: set[tuple[str, str]] = set()
        ledger_payload: dict | None = None
        for pin in pins:
            if not isinstance(pin, dict) or set(pin) != PIN_KEYS:
                fail(f"{where}: malformed artifact pin")
            revision = pin["revision"]
            path = pin["path"]
            if not isinstance(revision, str) or not COMMIT_RE.fullmatch(revision):
                fail(f"{where}: malformed pin revision")
            if (
                not isinstance(path, str)
                or not path
                or path.startswith("/")
                or ".." in path.split("/")
            ):
                fail(f"{where}: artifact paths must be safe repo-relative paths")
            key = (revision, path)
            if key in pin_keys:
                fail(f"{where}: duplicate artifact pin {revision}:{path}")
            pin_keys.add(key)
            payload = _git_bytes("show", f"{revision}:{path}")
            actual_blob = _git_bytes("rev-parse", f"{revision}:{path}").decode().strip()
            actual_digest = hashlib.sha256(payload).hexdigest()
            if not isinstance(pin["bytes"], int) or isinstance(pin["bytes"], bool):
                fail(f"{where}: malformed byte count for {path}")
            if not isinstance(pin["sha256"], str) or not DIGEST_RE.fullmatch(pin["sha256"]):
                fail(f"{where}: malformed SHA-256 for {path}")
            if not isinstance(pin["git_blob_sha1"], str) or not BLOB_RE.fullmatch(pin["git_blob_sha1"]):
                fail(f"{where}: malformed Git blob id for {path}")
            if (pin["bytes"], pin["sha256"], pin["git_blob_sha1"]) != (
                len(payload),
                actual_digest,
                actual_blob,
            ):
                fail(f"{where}: historical artifact pin drift for {revision}:{path}")
            if path == "tracking/observation_ledger.json" and revision == record["repair_commit"]:
                ledger_payload = _historical_json(revision, path)
        if ledger_payload is None:
            fail(f"{where}: repair-commit observation-ledger pin is mandatory")
        historical_rows = ledger_payload.get("rows")
        if not isinstance(historical_rows, list):
            fail(f"{where}: pinned observation ledger has no rows")
        by_id = reviewed_row_index(record)
        if set(reviewed) - set(by_id):
            fail(f"{where}: reviewed row absent from pinned observation ledger")
        attained = {row_id for row_id, row in by_id.items() if row.get("status") == "attained"}
        if set(promoted) != attained.intersection(reviewed):
            fail(f"{where}: promoted_rows must exactly match attained reviewed rows")
        pin_paths = {pin["path"] for pin in pins}
        pin_revision_paths = {(pin["revision"], pin["path"]) for pin in pins}
        if source_record is not None:
            source_evidence = source_record.get("evidence")
            if not isinstance(source_evidence, list):
                fail(f"{where}: pinned source audit has malformed evidence")
            source_evidence_paths = {
                item["path"]
                if isinstance(item, dict) and isinstance(item.get("path"), str)
                else item
                for item in source_evidence
            }
            if not source_evidence_paths.issubset(pin_paths):
                fail(f"{where}: artifact pins omit source audit evidence")
            pins_by_path = {pin["path"]: pin for pin in pins}
            for evidence in source_evidence:
                if not isinstance(evidence, dict):
                    continue
                source_sha = evidence.get("sha256")
                if isinstance(source_sha, str):
                    source_sha = source_sha.removeprefix("sha256:")
                    if pins_by_path[evidence["path"]]["sha256"] != source_sha:
                        fail(f"{where}: artifact pin rewrites source audit evidence")
                if (
                    isinstance(evidence.get("bytes"), int)
                    and pins_by_path[evidence["path"]]["bytes"] != evidence["bytes"]
                ):
                    fail(f"{where}: artifact pin rewrites source audit evidence")
        required_row_evidence = {
            (record["repair_commit"], evidence_path)
            for row_id in promoted
            for evidence_path in by_id[row_id].get("evidence", [])
        }
        if not required_row_evidence.issubset(pin_revision_paths):
            fail(
                f"{where}: promoted-row evidence is not pinned at the repair commit"
            )

        findings = record["findings"]
        if not isinstance(findings, list) or not findings:
            fail(f"{where}: findings must be nonempty")
        finding_ids: set[str] = set()
        for finding in findings:
            if not isinstance(finding, dict) or set(finding) != FINDING_KEYS:
                fail(f"{where}: malformed finding")
            finding_id = _nonempty(f"{where}.finding.id", finding["id"])
            if finding_id in finding_ids:
                fail(f"{where}: duplicate finding id {finding_id}")
            finding_ids.add(finding_id)
            if finding["severity"] not in SEVERITIES:
                fail(f"{where}.{finding_id}: unknown severity")
            if finding["disposition"] not in DISPOSITIONS:
                fail(f"{where}.{finding_id}: unknown disposition")
            _nonempty(f"{where}.{finding_id}.summary", finding["summary"])
            owners = finding["owner_issues"]
            if (
                not isinstance(owners, list)
                or any(not isinstance(owner, int) or isinstance(owner, bool) for owner in owners)
                or len(owners) != len(set(owners))
            ):
                fail(f"{where}.{finding_id}: owner_issues must be unique integers")

        for key in ("qualifies_for", "does_not_qualify_for", "limitations"):
            _strings(f"{where}.{key}", record[key])
        effective = copy.deepcopy(record)
        if provenance == "legacy_migration":
            effective["_origin_state"] = "legacy_source_pinned"
        elif audit_id in anchored_native_ids:
            effective["_origin_state"] = "native_origin_anchored"
        else:
            effective["_origin_state"] = "pending_origin_anchor"
            effective["_declared_promoted_rows"] = list(record["promoted_rows"])
            effective["_declared_qualifies_for"] = list(record["qualifies_for"])
            effective["promoted_rows"] = []
            effective["qualifies_for"] = [
                item
                for item in record["qualifies_for"]
                if item != "attained_status_review"
            ]
        validated_records.append(effective)
    _validate_committed_history_guard(data, record_ids, anchor_ids)
    _validate_committed_head_anchors(anchors)
    return validated_records


def render(data: dict, records: list[dict]) -> str:
    lines = [
        "# Promotion audit custody",
        "",
        "Generated by `tools/build_audit_custody.py` from `tracking/audit_custody.json`. "
        "Standing custody belongs to [issue #738](https://github.com/FloatingPragma/observer-patch-holography/issues/738).",
        "",
        data["policy"],
        "",
        "## Native two-commit protocol",
        "",
        "1. Append the complete native record with `provenance_kind` set to `native`, "
        "`source_record_pin` set to `null`, and no anchor for its id; commit it. The "
        "validator renders that record as pending and removes its promotion qualification "
        "from every validated consumer.",
        "2. In a later commit, append its `record_anchors` entry using the first commit as "
        "`origin_revision` and the canonical full-record SHA-256 as `record_sha256`. Do not "
        "alter the record. The anchor declaration must itself be committed at `HEAD`; a "
        "worktree-only anchor cannot qualify. Only then can its declared promotion scope qualify.",
    ]
    anchors_by_id = {anchor["id"]: anchor for anchor in data["record_anchors"]}
    for record in records:
        origin_state = record["_origin_state"]
        if origin_state == "legacy_source_pinned":
            source_pin = record["source_record_pin"]
            origin_text = (
                "Legacy migration, qualified by the exact historical source record "
                f"at `{source_pin['revision']}:{source_pin['path']}`."
            )
        elif origin_state == "native_origin_anchored":
            anchor = anchors_by_id[record["id"]]
            origin_text = (
                "Native record, qualified by its immutable origin at "
                f"`{anchor['origin_revision']}:tracking/audit_custody.json` "
                f"(record SHA-256 `{anchor['record_sha256']}`)."
            )
        else:
            declared = record.get("_declared_promoted_rows", [])
            origin_text = (
                "Native record pending its origin anchor. It cannot qualify a promotion; "
                "its declared post-anchor promotion scope is "
                f"{', '.join(declared) or 'none'}."
            )
        lines.extend(
            [
                "",
                f"## {record['id']}",
                "",
                f"Window: `{record['started_on']}` through `{record['completed_on']}`. "
                f"Baseline `{record['baseline_commit']}`; reviewed `{record['reviewed_commit']}`; "
                f"repair `{record['repair_commit']}`.",
                "",
                origin_text,
                "",
                record["audit_class"],
                "",
                record["scope"],
                "",
                "Reviewers:",
                "",
            ]
        )
        for reviewer in record["reviewers"]:
            lines.append(
                f"- **{reviewer['name']}** (`{reviewer['task']}`; {reviewer['model']}): "
                f"{reviewer['role']}"
            )
        lines.extend(
            [
                "",
                f"Reviewed rows: {', '.join(record['reviewed_rows'])}.",
                "",
                f"Promotion-qualified rows: {', '.join(record['promoted_rows']) or 'none'}.",
                "",
                "| Finding | Severity | Disposition | Owners | Summary |",
                "| --- | --- | --- | --- | --- |",
            ]
        )
        for finding in record["findings"]:
            owners = ", ".join(f"#{owner}" for owner in finding["owner_issues"]) or "none"
            lines.append(
                f"| {finding['id']} | {finding['severity']} | {finding['disposition']} | "
                f"{owners} | {finding['summary']} |"
            )
        lines.extend(["", "Historical artifact pins:", ""])
        for pin in record["artifact_pins"]:
            lines.append(
                f"- `{pin['revision']}:{pin['path']}`: {pin['bytes']} bytes, "
                f"SHA-256 `{pin['sha256']}`, Git blob `{pin['git_blob_sha1']}`."
            )
        lines.extend(["", "Limitations:", "", *[f"- {item}" for item in record["limitations"]]])
    lines.append("")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)
    data = load_json(REGISTER_PATH)
    records = validate(data)
    rendered = render(data, records).encode("utf-8")
    if args.check:
        if not SURFACE_PATH.is_file() or SURFACE_PATH.read_bytes() != rendered:
            print("audit custody: generated surface is stale", file=sys.stderr)
            return 1
        print(f"audit custody: {len(records)} immutable record(s)")
        return 0
    SURFACE_PATH.write_bytes(rendered)
    print(f"audit custody: wrote {SURFACE_PATH.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
