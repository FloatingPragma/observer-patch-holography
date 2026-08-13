#!/usr/bin/env python3
"""Validate append-only, architecture-bound frozen-prediction lineage custody.

The frozen-prediction register predates the V3 architecture register.  Its
rows therefore remain historical inputs, never silently becoming predictions
of an AV-n.  New predictive custody is represented by an append-only freeze
event.  The event is committed once, then anchored by a following commit to
the exact first-appearance revision.  Only an anchored event on the current,
anchored architecture can qualify an attained predictive observation row.
"""

from __future__ import annotations

import hashlib
import json
import re
import subprocess
from pathlib import Path

import strict_json

ROOT = Path(__file__).resolve().parents[1]
LINEAGE_PATH = ROOT / "claims" / "frozen_prediction_architecture_lineages.json"
FROZEN_PATH = ROOT / "claims" / "frozen_prediction_register.json"
ARCHITECTURE_PATH = ROOT / "tracking" / "architecture_versions.json"

SCHEMA = "oph.frozen_prediction_architecture_lineages.v2"
ISSUE = 741
FROZEN_BOOTSTRAP_REVISION = "68663e9e52a3931c322676a127dd0af144a01de3"
FROZEN_BOOTSTRAP_GIT_BLOB_SHA1 = "731ff48c6005a0d2b3930c1d33a6f00c2dd3e345"
# Filled from the canonical v2 static projection below. A schema migration is
# required to change policy, bootstrap custody, baseline rows, or OL bindings.
STATIC_V2_SHA256 = "f950ea41ab51317704e5fad5fa6072dc36669e7fdaaec6b1e28d8fbe39498aa5"
GIT_REVISION_RE = re.compile(r"^[0-9a-f]{40}$")
GIT_BLOB_RE = re.compile(r"^[0-9a-f]{40}$")
EVENT_ID_RE = re.compile(r"^FZE-([0-9]{3})$")

TOP_KEYS = {
    "schema",
    "issue",
    "frozen_register_revision",
    "frozen_register_git_blob_sha1",
    "policy",
    "rows",
    "predictive_observation_bindings",
    "events",
    "event_anchors",
}
ROW_KEYS = {
    "id",
    "registered_content_sha256",
    "source_architecture_version",
    "lineage_status",
    "current_version_eligibility",
    "transition_policy",
}
BINDING_KEYS = {
    "observation_row_id",
    "candidate_baseline_lineages",
    "historical_only_lineages",
    "no_candidate_reason",
    "qualification_rule",
}
EVENT_KEYS = {
    "id",
    "source_architecture_version",
    "source_architecture_snapshot_sha256",
    "observation_row_id",
    "observation_contract_sha256",
    "supersedes_lineage_ids",
    "frozen_utc",
    "target",
    "comparison_protocol",
    "kill_band",
    "precomparison_statement",
    "target_payload_sha256",
    "predecessor_register_revision",
}
ANCHOR_KEYS = {
    "id",
    "origin_revision",
    "event_sha256",
    "register_git_blob_sha1",
}
QUALIFICATION_RULE = "anchored_current_architecture_event_covering_full_row"
PREDICTION_EVENT_POINTER_KEYS = {"id", "target_payload_sha256"}

LEGACY_POLICY = (
    "Preserve historical custody; preregister a new target in an append-only "
    "event bound to one current architecture version before predictive use."
)
PENDING_POLICY = (
    "Bind the eventual freeze through an append-only event to the then-current "
    "AV-n; any later architecture transition invalidates that event by default."
)
VOID_POLICY = (
    "Preserve the void historical record; no architecture transition or event "
    "can revive it. A scientifically distinct successor must be a new event."
)


def fail(message: str) -> None:
    raise SystemExit(f"prediction lineage custody: {message}")


def load_json(path: Path) -> dict:
    try:
        return strict_json.load(path)
    except FileNotFoundError:
        fail(f"missing {path}")
    except (json.JSONDecodeError, strict_json.DuplicateKeyError) as error:
        fail(f"invalid JSON in {path}: {error}")
    raise AssertionError("unreachable")


def _git(root: Path, *args: str, check: bool = True) -> bytes:
    result = subprocess.run(
        ["git", *args],
        cwd=root,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if check and result.returncode != 0:
        fail(
            f"git {' '.join(args)} failed: "
            f"{result.stderr.decode('utf-8', errors='replace').strip()}"
        )
    return result.stdout


def _git_text(root: Path, *args: str, check: bool = True) -> str:
    return _git(root, *args, check=check).decode("utf-8").strip()


def _is_ancestor(root: Path, ancestor: str, descendant: str = "HEAD") -> bool:
    result = subprocess.run(
        ["git", "merge-base", "--is-ancestor", ancestor, descendant],
        cwd=root,
        check=False,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return result.returncode == 0


def _is_first_parent_ancestor(root: Path, revision: str) -> bool:
    return revision in set(
        _git_text(root, "rev-list", "--first-parent", "HEAD").splitlines()
    )


def canonical_sha256(value: object) -> str:
    payload = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def event_payload(event: dict) -> dict:
    return {
        "source_architecture_version": event["source_architecture_version"],
        "source_architecture_snapshot_sha256": event[
            "source_architecture_snapshot_sha256"
        ],
        "observation_row_id": event["observation_row_id"],
        "observation_contract_sha256": event["observation_contract_sha256"],
        "supersedes_lineage_ids": event["supersedes_lineage_ids"],
        "frozen_utc": event["frozen_utc"],
        "target": event["target"],
        "comparison_protocol": event["comparison_protocol"],
        "kill_band": event["kill_band"],
        "precomparison_statement": event["precomparison_statement"],
    }


def event_sha256(event: dict) -> str:
    return canonical_sha256(event)


def static_projection(data: dict) -> dict:
    return {
        key: data[key]
        for key in TOP_KEYS
        if key not in {"events", "event_anchors"}
    }


def observation_contract(row: dict) -> dict:
    """Fields whose scientific meaning must already be fixed at freeze time."""

    return {
        "id": row["id"],
        "target": row["target"],
        "rung": row["rung"],
        "lane_issue": row["lane_issue"],
        "architecture_version": row["architecture_version"],
        "premises": row["premises"],
        "open_premises": row["open_premises"],
    }


def _historical_json(root: Path, revision: str, relative_path: str) -> dict | None:
    raw = _git(root, "show", f"{revision}:{relative_path}", check=False)
    if not raw:
        return None
    try:
        value = strict_json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError, strict_json.DuplicateKeyError) as error:
        fail(f"invalid historical {relative_path} at {revision}: {error}")
    if not isinstance(value, dict):
        fail(f"historical {relative_path} at {revision} is not an object")
    return value


def _baseline_row(frozen_row: dict) -> dict:
    row_id = frozen_row.get("id")
    if frozen_row.get("status") == "superseded_void":
        status = "superseded_void"
        eligibility = "permanently_ineligible"
        policy = VOID_POLICY
    elif frozen_row.get("frozen_utc") is not None:
        status = "legacy_unversioned"
        eligibility = "ineligible_new_freeze_required"
        policy = LEGACY_POLICY
    else:
        status = "unbound_pending_freeze"
        eligibility = "ineligible_until_version_bound_freeze"
        policy = PENDING_POLICY
    return {
        "id": row_id,
        "registered_content_sha256": frozen_row.get("content_sha256"),
        "source_architecture_version": None,
        "lineage_status": status,
        "current_version_eligibility": eligibility,
        "transition_policy": policy,
    }


def _architecture_index(architecture: dict) -> tuple[dict[str, dict], set[str], str]:
    versions = architecture.get("versions") if isinstance(architecture, dict) else None
    anchors = architecture.get("version_anchors") if isinstance(architecture, dict) else None
    current = architecture.get("current_version") if isinstance(architecture, dict) else None
    if not isinstance(versions, list) or not isinstance(anchors, list) or not isinstance(current, str):
        fail("architecture register is malformed")
    by_id = {
        row.get("id"): row
        for row in versions
        if isinstance(row, dict) and isinstance(row.get("id"), str)
    }
    anchored = {
        row.get("id")
        for row in anchors
        if isinstance(row, dict) and isinstance(row.get("id"), str)
    }
    if current not in by_id:
        fail("architecture current_version is not registered")
    return by_id, anchored, current


def _validate_binding_rows(
    bindings: object,
    baseline_by_id: dict[str, dict],
    observation_rows: list[dict],
) -> tuple[list[dict], dict[str, dict]]:
    if not isinstance(bindings, list):
        fail("predictive_observation_bindings must be a list")
    predictive_ids = [
        row["id"] for row in observation_rows if row.get("rung") == "predictive"
    ]
    binding_ids = [
        row.get("observation_row_id") if isinstance(row, dict) else None
        for row in bindings
    ]
    if binding_ids != predictive_ids:
        fail(
            "predictive observation bindings must exactly match predictive ledger "
            "rows in ledger order"
        )
    by_id: dict[str, dict] = {}
    for row in bindings:
        row_id = row.get("observation_row_id") if isinstance(row, dict) else "<malformed>"
        if not isinstance(row, dict) or set(row) != BINDING_KEYS:
            fail(f"{row_id}: binding keys must equal {sorted(BINDING_KEYS)}")
        candidates = row["candidate_baseline_lineages"]
        historical = row["historical_only_lineages"]
        reason = row["no_candidate_reason"]
        if not isinstance(candidates, list) or not isinstance(historical, list):
            fail(f"{row_id}: lineage binding lists must be lists")
        if len(candidates) != len(set(candidates)) or len(historical) != len(set(historical)):
            fail(f"{row_id}: lineage binding lists must be duplicate-free")
        if set(candidates) & set(historical):
            fail(f"{row_id}: candidate and historical lineage lists must be disjoint")
        for lineage_id in (*candidates, *historical):
            if lineage_id not in baseline_by_id:
                fail(f"{row_id}: unknown baseline lineage {lineage_id}")
        for lineage_id in candidates:
            if baseline_by_id[lineage_id]["lineage_status"] != "unbound_pending_freeze":
                fail(f"{row_id}: candidate {lineage_id} is not pending a new freeze")
        for lineage_id in historical:
            if baseline_by_id[lineage_id]["lineage_status"] not in {
                "legacy_unversioned",
                "superseded_void",
            }:
                fail(
                    f"{row_id}: historical-only lineage {lineage_id} is neither "
                    "legacy frozen nor permanently void"
                )
        if candidates:
            if reason is not None:
                fail(f"{row_id}: a row with candidates must set no_candidate_reason to null")
        elif not isinstance(reason, str) or not reason.strip():
            fail(f"{row_id}: a row without candidates needs an explicit reason")
        if row["qualification_rule"] != QUALIFICATION_RULE:
            fail(f"{row_id}: qualification_rule must equal {QUALIFICATION_RULE}")
        by_id[row_id] = row
    return bindings, by_id


def _validate_event_origin(
    event: dict,
    anchor: dict,
    event_index: int,
    *,
    root: Path,
    lineage_relative_path: str,
    architecture_relative_path: str,
    observation_relative_path: str,
) -> None:
    event_id = event["id"]
    origin = anchor["origin_revision"]
    if not _is_first_parent_ancestor(root, origin):
        fail(f"{event_id}: origin revision is not on the HEAD first-parent chain")
    actual_blob = _git_text(root, "rev-parse", f"{origin}:{lineage_relative_path}")
    if actual_blob != anchor["register_git_blob_sha1"]:
        fail(f"{event_id}: origin register Git blob drifted")
    origin_data = _historical_json(root, origin, lineage_relative_path)
    if origin_data is None or origin_data.get("schema") != SCHEMA:
        fail(f"{event_id}: origin revision has no schema-v2 lineage register")
    origin_events = origin_data.get("events")
    origin_anchors = origin_data.get("event_anchors")
    if not isinstance(origin_events, list) or not isinstance(origin_anchors, list):
        fail(f"{event_id}: origin register has malformed event arrays")
    if len(origin_events) != event_index + 1 or origin_events[event_index] != event:
        fail(f"{event_id}: origin must be the exact first append of this event")
    if any(item.get("id") == event_id for item in origin_anchors if isinstance(item, dict)):
        fail(f"{event_id}: event and its origin anchor must be separate commits")
    parents = _git_text(root, "rev-list", "--parents", "-n", "1", origin).split()
    if len(parents) != 2:
        fail(f"{event_id}: origin revision must have exactly one parent")
    parent = parents[1]
    if parent != event["predecessor_register_revision"]:
        fail(f"{event_id}: predecessor_register_revision is not origin parent")
    parent_data = _historical_json(root, parent, lineage_relative_path)
    if parent_data is None or parent_data.get("schema") != SCHEMA:
        fail(f"{event_id}: predecessor has no schema-v2 lineage register")
    if parent_data.get("events") != origin_events[:event_index]:
        fail(f"{event_id}: predecessor event prefix drifted")
    if parent_data.get("event_anchors") != origin_anchors:
        fail(f"{event_id}: origin commit may append only the event")

    origin_architecture = _historical_json(root, origin, architecture_relative_path)
    if origin_architecture is None:
        fail(f"{event_id}: origin revision has no architecture register")
    versions, anchored, current = _architecture_index(origin_architecture)
    source = event["source_architecture_version"]
    if source != current or source not in anchored:
        fail(f"{event_id}: source AV was not current and anchored at event origin")
    if versions[source].get("snapshot_sha256") != event["source_architecture_snapshot_sha256"]:
        fail(f"{event_id}: source architecture snapshot drifted at event origin")

    origin_observations = _historical_json(root, origin, observation_relative_path)
    origin_rows = origin_observations.get("rows") if origin_observations else None
    if not isinstance(origin_rows, list):
        fail(f"{event_id}: origin revision has no observation ledger rows")
    row = next(
        (
            item
            for item in origin_rows
            if isinstance(item, dict)
            and item.get("id") == event["observation_row_id"]
        ),
        None,
    )
    if row is None:
        fail(f"{event_id}: origin observation row is missing")
    if canonical_sha256(observation_contract(row)) != event["observation_contract_sha256"]:
        fail(f"{event_id}: origin observation contract differs from the freeze event")
    if row.get("status") == "attained":
        fail(f"{event_id}: origin observation row was already attained")
    if row.get("prediction_event") is not None:
        fail(f"{event_id}: origin observation row was not pointer-free")
    audit_pointers = origin_observations.get("audit_pointers")
    if not isinstance(audit_pointers, dict):
        fail(f"{event_id}: origin observation audit pointers are malformed")
    if event["observation_row_id"] in audit_pointers:
        fail(f"{event_id}: origin observation row already had a promotion audit pointer")

    parent_observations = _historical_json(root, parent, observation_relative_path)
    parent_rows = parent_observations.get("rows") if parent_observations else None
    if not isinstance(parent_rows, list):
        fail(f"{event_id}: predecessor has no observation ledger rows")
    parent_row = next(
        (
            item
            for item in parent_rows
            if isinstance(item, dict)
            and item.get("id") == event["observation_row_id"]
        ),
        None,
    )
    if parent_row != row:
        fail(f"{event_id}: origin commit modified the bound observation row")
    parent_audits = parent_observations.get("audit_pointers")
    if not isinstance(parent_audits, dict):
        fail(f"{event_id}: predecessor observation audit pointers are malformed")
    if parent_audits.get(event["observation_row_id"]) != audit_pointers.get(
        event["observation_row_id"]
    ):
        fail(f"{event_id}: origin commit modified the row's audit pointer")

    # Reject two-step laundering as well: downgrading an attained row in one
    # commit and appending a nominal "precomparison" freeze in the next does
    # not restore predictive eligibility for the same scientific contract.
    observation_revisions = _git_text(
        root,
        "log",
        "--first-parent",
        "--full-history",
        "--format=%H",
        parent,
        "--",
        observation_relative_path,
    ).splitlines()
    for revision in observation_revisions:
        historical = _historical_json(root, revision, observation_relative_path)
        historical_rows = historical.get("rows") if historical else None
        if not isinstance(historical_rows, list):
            continue
        historical_row = next(
            (
                item
                for item in historical_rows
                if isinstance(item, dict)
                and item.get("id") == event["observation_row_id"]
            ),
            None,
        )
        if historical_row is None:
            continue
        try:
            historical_contract = observation_contract(historical_row)
        except KeyError:
            # Rows predating the V3 architecture/premise schema cannot be the
            # exact scientific contract frozen by a V3 event.
            continue
        if canonical_sha256(historical_contract) != event["observation_contract_sha256"]:
            continue
        historical_audits = historical.get("audit_pointers")
        had_audit = isinstance(historical_audits, dict) and event[
            "observation_row_id"
        ] in historical_audits
        if historical_row.get("status") == "attained" or had_audit:
            fail(
                f"{event_id}: the same observation contract was already attained "
                "or promotion-audited before its freeze"
            )


def _validate_committed_history(
    data: dict,
    *,
    root: Path,
    lineage_relative_path: str,
) -> None:
    revisions = _git_text(
        root,
        "log",
        "--first-parent",
        "--full-history",
        "--reverse",
        "--format=%H",
        "HEAD",
        "--",
        lineage_relative_path,
    ).splitlines()
    seen_v2 = False
    gap_after_v2 = False
    current_events = data["events"]
    current_anchors = data["event_anchors"]
    for revision in revisions:
        historical = _historical_json(root, revision, lineage_relative_path)
        if historical is None or historical.get("schema") != SCHEMA:
            if seen_v2:
                gap_after_v2 = True
            continue
        if gap_after_v2:
            fail("schema-v2 lineage register disappeared and later reappeared")
        seen_v2 = True
        events = historical.get("events")
        anchors = historical.get("event_anchors")
        if not isinstance(events, list) or not isinstance(anchors, list):
            fail(f"historical lineage arrays malformed at {revision}")
        static = static_projection(historical)
        current_static = static_projection(data)
        if static != current_static:
            fail(f"lineage bootstrap, policy, baseline, or bindings drifted at {revision}")
        if events != current_events[: len(events)]:
            fail(f"historical lineage event was deleted or rewritten at {revision}")
        if anchors != current_anchors[: len(anchors)]:
            fail(f"historical lineage anchor was deleted or rewritten at {revision}")
    if seen_v2 and gap_after_v2:
        fail("schema-v2 lineage register disappeared before the current state")


def validate(
    data: dict,
    frozen: dict,
    architecture: dict,
    observation_rows: list[dict],
    *,
    root: Path = ROOT,
    lineage_relative_path: str = "claims/frozen_prediction_architecture_lineages.json",
    architecture_relative_path: str = "tracking/architecture_versions.json",
    observation_relative_path: str = "tracking/observation_ledger.json",
) -> dict:
    """Return validated baseline, bindings, events, and current eligibility."""

    if not isinstance(data, dict) or set(data) != TOP_KEYS:
        fail(f"top-level keys must equal {sorted(TOP_KEYS)}")
    if data["schema"] != SCHEMA or data["issue"] != ISSUE:
        fail(f"schema and issue must equal {SCHEMA} and {ISSUE}")
    if not isinstance(data["policy"], str) or not data["policy"].strip():
        fail("policy must be nonempty")
    frozen_revision = data["frozen_register_revision"]
    frozen_blob = data["frozen_register_git_blob_sha1"]
    if frozen_revision != FROZEN_BOOTSTRAP_REVISION:
        fail("frozen_register_revision differs from the audited bootstrap revision")
    if frozen_blob != FROZEN_BOOTSTRAP_GIT_BLOB_SHA1:
        fail("frozen_register_git_blob_sha1 differs from the audited bootstrap blob")
    if not _is_first_parent_ancestor(root, frozen_revision):
        fail("frozen-register bootstrap is not on the HEAD first-parent chain")
    actual_blob = _git_text(
        root, "rev-parse", f"{frozen_revision}:claims/frozen_prediction_register.json"
    )
    if actual_blob != frozen_blob:
        fail("frozen-register bootstrap Git blob drifted")
    historical_frozen = _historical_json(
        root, frozen_revision, "claims/frozen_prediction_register.json"
    )
    historical_rows = historical_frozen.get("rows") if historical_frozen else None
    current_rows = frozen.get("rows") if isinstance(frozen, dict) else None
    if not isinstance(historical_rows, list) or not isinstance(current_rows, list):
        fail("frozen prediction registers must carry row lists")
    if current_rows != historical_rows:
        fail("frozen prediction register differs from the lineage bootstrap payload")
    expected_rows = [_baseline_row(row) for row in historical_rows]
    if data["rows"] != expected_rows:
        fail("baseline lineage rows drifted from their deterministic frozen-register mapping")
    baseline_by_id = {row["id"]: row for row in expected_rows}

    bindings, binding_by_id = _validate_binding_rows(
        data["predictive_observation_bindings"], baseline_by_id, observation_rows
    )
    if canonical_sha256(static_projection(data)) != STATIC_V2_SHA256:
        fail("complete static v2 lineage surface differs from its canonical hash")
    versions, anchored_architectures, current_architecture = _architecture_index(architecture)

    events = data["events"]
    anchors = data["event_anchors"]
    if not isinstance(events, list) or not isinstance(anchors, list):
        fail("events and event_anchors must be lists")
    event_ids = [event.get("id") if isinstance(event, dict) else None for event in events]
    expected_event_ids = [f"FZE-{index:03d}" for index in range(1, len(events) + 1)]
    if event_ids != expected_event_ids:
        fail("event ids must be contiguous FZE-001..FZE-n in append order")
    anchor_ids = [anchor.get("id") if isinstance(anchor, dict) else None for anchor in anchors]
    if anchor_ids != event_ids[: len(anchors)]:
        fail("event anchors must be an exact prefix of events")
    if len(events) - len(anchors) > 1:
        fail("at most one final event may await its following origin anchor commit")

    # An anchor present only in the worktree is not custody. It becomes
    # operative only after a distinct following commit records it. Later
    # uncommitted events may coexist with already committed anchors.
    head_data = _historical_json(root, "HEAD", lineage_relative_path)
    head_anchors = (
        head_data.get("event_anchors")
        if isinstance(head_data, dict) and head_data.get("schema") == SCHEMA
        else []
    )
    if not isinstance(head_anchors, list):
        fail("committed HEAD lineage anchors are malformed")
    if anchors != head_anchors[: len(anchors)]:
        fail("every operative event anchor must already be committed at HEAD")

    event_by_id: dict[str, dict] = {}
    current_contract_match: dict[str, bool] = {}
    row_by_id = {row["id"]: row for row in observation_rows}
    qualifying_by_row: dict[str, list[dict]] = {row["id"]: [] for row in observation_rows if row.get("rung") == "predictive"}
    for index, event in enumerate(events):
        event_id = event.get("id") if isinstance(event, dict) else "<malformed>"
        if not isinstance(event, dict) or set(event) != EVENT_KEYS:
            fail(f"{event_id}: event keys must equal {sorted(EVENT_KEYS)}")
        if not EVENT_ID_RE.fullmatch(event_id):
            fail(f"{event_id}: malformed event id")
        for field in (
            "frozen_utc",
            "target",
            "comparison_protocol",
            "kill_band",
            "precomparison_statement",
        ):
            if not isinstance(event[field], str) or not event[field].strip():
                fail(f"{event_id}: {field} must be nonempty")
        observation_id = event["observation_row_id"]
        supersedes = event["supersedes_lineage_ids"]
        if observation_id not in binding_by_id:
            fail(f"{event_id}: {observation_id} is not a predictive observation row")
        row = row_by_id[observation_id]
        current_contract_match[event_id] = (
            event["observation_contract_sha256"]
            == canonical_sha256(observation_contract(row))
        )
        if not isinstance(supersedes, list) or len(supersedes) != len(set(supersedes)):
            fail(f"{event_id}: supersedes_lineage_ids must be a duplicate-free list")
        allowed_supersedes = set(
            binding_by_id[observation_id]["candidate_baseline_lineages"]
        ) | set(binding_by_id[observation_id]["historical_only_lineages"])
        if not set(supersedes) <= allowed_supersedes:
            fail(
                f"{event_id}: superseded lineages must be declared candidate or "
                f"historical context for {observation_id}"
            )
        if event["target_payload_sha256"] != canonical_sha256(event_payload(event)):
            fail(f"{event_id}: target payload hash drifted")
        source = event["source_architecture_version"]
        if source not in versions or source not in anchored_architectures:
            fail(f"{event_id}: source architecture must be registered and anchored")
        if versions[source].get("snapshot_sha256") != event["source_architecture_snapshot_sha256"]:
            fail(f"{event_id}: source architecture snapshot hash drifted")
        predecessor = event["predecessor_register_revision"]
        if not isinstance(predecessor, str) or not GIT_REVISION_RE.fullmatch(predecessor):
            fail(f"{event_id}: predecessor_register_revision is malformed")
        event_by_id[event_id] = event

        if index < len(anchors):
            anchor = anchors[index]
            if not isinstance(anchor, dict) or set(anchor) != ANCHOR_KEYS:
                fail(f"{event_id}: anchor keys must equal {sorted(ANCHOR_KEYS)}")
            if anchor["id"] != event_id:
                fail(f"{event_id}: anchor id mismatch")
            if not isinstance(anchor["origin_revision"], str) or not GIT_REVISION_RE.fullmatch(anchor["origin_revision"]):
                fail(f"{event_id}: origin_revision is malformed")
            if not isinstance(anchor["register_git_blob_sha1"], str) or not GIT_BLOB_RE.fullmatch(anchor["register_git_blob_sha1"]):
                fail(f"{event_id}: register_git_blob_sha1 is malformed")
            if anchor["event_sha256"] != event_sha256(event):
                fail(f"{event_id}: event anchor hash drifted")
            _validate_event_origin(
                event,
                anchor,
                index,
                root=root,
                lineage_relative_path=lineage_relative_path,
                architecture_relative_path=architecture_relative_path,
                observation_relative_path=observation_relative_path,
            )
            if source == current_architecture and current_contract_match[event_id]:
                qualifying_by_row[observation_id].append(event)

    for row_id, binding in binding_by_id.items():
        row = row_by_id[row_id]
        pointer = row.get("prediction_event")
        if pointer is None:
            continue
        if not isinstance(pointer, dict) or set(pointer) != PREDICTION_EVENT_POINTER_KEYS:
            fail(f"{row_id}: prediction_event pointer is malformed")
        event = event_by_id.get(pointer["id"])
        if event is None:
            fail(f"{row_id}: prediction_event names an unknown event")
        if event["observation_row_id"] != row_id:
            fail(f"{row_id}: prediction_event belongs to a different observation row")
        if event["target_payload_sha256"] != pointer["target_payload_sha256"]:
            fail(f"{row_id}: prediction_event target hash drifted")
        if not current_contract_match[event["id"]]:
            fail(f"{row_id}: prediction_event contract differs from the current row")

    _validate_committed_history(
        data, root=root, lineage_relative_path=lineage_relative_path
    )
    return {
        "baseline_rows": expected_rows,
        "bindings": bindings,
        "binding_by_id": binding_by_id,
        "events": events,
        "event_by_id": event_by_id,
        "anchored_event_ids": anchor_ids,
        "qualifying_events_by_row": qualifying_by_row,
        "current_architecture_version": current_architecture,
    }


def require_predictive_promotions(
    observation_rows: list[dict], state: dict
) -> None:
    current = state["current_architecture_version"]
    for row in observation_rows:
        if row.get("rung") != "predictive" or row.get("status") != "attained":
            continue
        row_id = row["id"]
        if row.get("architecture_version") != current:
            fail(f"{row_id}: attained predictive row is not recorded on the current AV-n")
        pointer = row.get("prediction_event")
        if not isinstance(pointer, dict):
            fail(f"{row_id}: attained predictive row has no exact event pointer")
        qualifying = [
            event
            for event in state["qualifying_events_by_row"].get(row_id, [])
            if event["id"] == pointer.get("id")
            and event["target_payload_sha256"]
            == pointer.get("target_payload_sha256")
        ]
        if len(qualifying) != 1:
            fail(
                f"{row_id}: attained predictive row lacks exactly one anchored, "
                "current-AV precomparison freeze event matching its audited pointer"
            )
