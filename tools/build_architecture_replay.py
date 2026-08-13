#!/usr/bin/env python3
"""Build the architecture-transition replay and promotion index for V3.

The generator joins architecture versions, observation promotions, historical
audit custody, and frozen-prediction lineages.  It fails closed when a
promoted row lacks a version-matching audit or when a frozen target could be
silently inherited by a different architecture version.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import build_architecture_versions as architecture
import build_audit_custody as audit_custody
import build_observation_ledger as observation
import prediction_lineage_custody
import strict_json

ROOT = Path(__file__).resolve().parents[1]
LINEAGE_PATH = ROOT / "claims" / "frozen_prediction_architecture_lineages.json"
FROZEN_PATH = ROOT / "claims" / "frozen_prediction_register.json"
OUTPUT_PATH = ROOT / "tracking" / "architecture_replay_index.json"
SURFACE_PATH = ROOT / "docs" / "ARCHITECTURE_REPLAY_INDEX.md"

SCHEMA = "oph.architecture_replay_index.v1"
ISSUE = 741


def fail(message: str) -> None:
    raise SystemExit(f"architecture replay: {message}")


def load_json(path: Path) -> dict:
    try:
        return strict_json.load(path)
    except FileNotFoundError:
        fail(f"missing {path.relative_to(ROOT)}")
    except (json.JSONDecodeError, strict_json.DuplicateKeyError) as error:
        fail(f"invalid JSON in {path.relative_to(ROOT)}: {error}")
    raise AssertionError("unreachable")


def _normative_map(version: dict) -> dict[str, dict]:
    return {item["path"]: item for item in version["normative_files"]}


def _decision_map(version: dict) -> dict[str, dict]:
    return {item["id"]: item for item in version["protocol_decisions"]}


def _pinned_payload(version: dict, path: str, *, is_current: bool) -> bytes:
    item = _normative_map(version).get(path)
    if item is None:
        return b""
    if is_current:
        return (ROOT / path).read_bytes()
    return architecture._historical_blob(version["id"], path, item["git_blob_sha1"])


def _premise_rows(payload: bytes) -> dict[str, dict]:
    if not payload:
        return {}
    try:
        data = strict_json.loads(payload.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError, strict_json.DuplicateKeyError) as error:
        fail(f"cannot compare pinned premise register: {error}")
    rows = data.get("rows") if isinstance(data, dict) else None
    if not isinstance(rows, list):
        fail("pinned premise register has no rows")
    return {
        row["id"]: row
        for row in rows
        if isinstance(row, dict) and isinstance(row.get("id"), str)
    }


def semantic_diff(
    predecessor: dict,
    successor: dict,
    observation_rows: list[dict],
    lineage_rows: list[dict],
    *,
    successor_is_current: bool,
) -> dict:
    predecessor_files = _normative_map(predecessor)
    successor_files = _normative_map(successor)
    file_paths = sorted(set(predecessor_files) | set(successor_files))
    changed_files = [
        path
        for path in file_paths
        if predecessor_files.get(path) != successor_files.get(path)
    ]
    predecessor_decisions = _decision_map(predecessor)
    successor_decisions = _decision_map(successor)
    decision_ids = sorted(set(predecessor_decisions) | set(successor_decisions))
    changed_decisions = [
        decision_id
        for decision_id in decision_ids
        if predecessor_decisions.get(decision_id) != successor_decisions.get(decision_id)
    ]
    basis_changed = predecessor["basis"] != successor["basis"]
    trigger_policy_changed = (
        predecessor["invalidation_triggers"] != successor["invalidation_triggers"]
    )
    replay_policy_changed = predecessor["replay_surfaces"] != successor["replay_surfaces"]

    premise_path = "tracking/premise_register.json"
    predecessor_premises = _premise_rows(
        _pinned_payload(predecessor, premise_path, is_current=False)
    )
    successor_premises = _premise_rows(
        _pinned_payload(successor, premise_path, is_current=successor_is_current)
    )
    premise_ids = sorted(set(predecessor_premises) | set(successor_premises))
    changed_premises = [
        premise_id
        for premise_id in premise_ids
        if predecessor_premises.get(premise_id) != successor_premises.get(premise_id)
    ]

    non_premise_file_changed = any(path != premise_path for path in changed_files)
    global_change = bool(
        basis_changed
        or changed_decisions
        or trigger_policy_changed
        or replay_policy_changed
        or non_premise_file_changed
        or changed_premises
    )
    if global_change:
        affected_rows = [row["id"] for row in observation_rows]
        reason = (
            "changed_premise_without_versioned_consumer_snapshot"
            if changed_premises
            and not (
                basis_changed
                or changed_decisions
                or trigger_policy_changed
                or replay_policy_changed
                or non_premise_file_changed
            )
            else "global_architecture_or_promotion_semantics_changed"
        )
    else:
        affected_rows = []
        reason = "no_observation_semantics_changed"

    affected_lineages = [
        row["id"]
        for row in lineage_rows
        if row["lineage_status"] == "version_bound"
        and row["source_architecture_version"] == predecessor["id"]
    ]
    return {
        "from_version": predecessor["id"],
        "to_version": successor["id"],
        "basis_changed": basis_changed,
        "changed_normative_files": changed_files,
        "changed_protocol_decisions": changed_decisions,
        "changed_premises": changed_premises,
        "invalidation_policy_changed": trigger_policy_changed,
        "replay_surface_policy_changed": replay_policy_changed,
        "affected_observation_rows": affected_rows,
        "observation_invalidation_reason": reason,
        "affected_frozen_prediction_lineages": affected_lineages,
        "prediction_transition_rule": "Every affected version-bound lineage is ineligible on the successor until a new target is frozen before comparison data are inspected.",
    }


def _expanded_audit(record: dict) -> dict:
    return {
        "id": record["id"],
        "origin_state": record["_origin_state"],
        "reviewed_commit": record["reviewed_commit"],
        "repair_commit": record["repair_commit"],
        "reviewers": [
            {"name": reviewer["name"], "model": reviewer["model"]}
            for reviewer in record["reviewers"]
        ],
        "finding_ids": [finding["id"] for finding in record["findings"]],
        "artifact_hashes": [
            {
                "revision": pin["revision"],
                "path": pin["path"],
                "sha256": pin["sha256"],
                "git_blob_sha1": pin["git_blob_sha1"],
            }
            for pin in record["artifact_pins"]
        ],
    }


def build() -> dict:
    architecture_data = architecture.load_json(architecture.REGISTER_PATH)
    versions = architecture.validate(architecture_data)
    current_version = architecture_data["current_version"]
    anchored_architecture_ids = {
        anchor["id"] for anchor in architecture_data["version_anchors"]
    }

    ledger = observation.load_ledger(observation.LEDGER_PATH)
    observation_rows = observation.validate(ledger)
    audit_data = audit_custody.load_json(audit_custody.REGISTER_PATH)
    audit_records = audit_custody.validate(audit_data)
    audits_by_id = {record["id"]: record for record in audit_records}

    frozen = load_json(FROZEN_PATH)
    lineage_data = load_json(LINEAGE_PATH)
    lineage_state = prediction_lineage_custody.validate(
        lineage_data,
        frozen,
        architecture_data,
        observation_rows,
        root=ROOT,
    )
    prediction_lineage_custody.require_predictive_promotions(
        observation_rows, lineage_state
    )
    anchored_event_ids = set(lineage_state["anchored_event_ids"])
    lineage_rows = list(lineage_state["baseline_rows"])
    lineage_rows.extend(
        {
            "id": event["id"],
            "registered_content_sha256": event["target_payload_sha256"],
            "source_architecture_version": event["source_architecture_version"],
            "lineage_status": (
                "version_bound"
                if event["id"] in anchored_event_ids
                else "event_pending_origin_anchor"
            ),
            "current_version_eligibility": (
                "eligible_on_source_version_only"
                if event["id"] in anchored_event_ids
                else "ineligible_until_origin_anchor"
            ),
            "transition_policy": (
                "A successor architecture requires a new precomparison event; "
                "the anchored event remains historical on its source AV-n."
            ),
            "observation_rows": [event["observation_row_id"]],
            "supersedes_lineage_ids": event["supersedes_lineage_ids"],
        }
        for event in lineage_state["events"]
    )

    transitions = [
        semantic_diff(
            versions[index - 1],
            versions[index],
            observation_rows,
            lineage_rows,
            successor_is_current=index == len(versions) - 1,
        )
        for index in range(1, len(versions))
    ]

    observation_index: list[dict] = []
    for row in observation_rows:
        audit_ids = ledger["audit_pointers"].get(row["id"], [])
        is_current = row["architecture_version"] == current_version
        is_anchored = row["architecture_version"] in anchored_architecture_ids
        if row["status"] == "attained" and is_current and is_anchored:
            promotion_state = "qualified_on_current_version"
            action = "none"
        elif row["status"] == "attained" and is_current:
            promotion_state = "invalidated_unanchored_architecture_tip"
            action = "anchor the complete AV-n record, then replay evidence and independent audit"
        elif row["status"] == "attained":
            promotion_state = "invalidated_on_current_version"
            action = "replay evidence and independent audit on the current AV-n before promotion"
        else:
            promotion_state = "not_promoted"
            action = "complete the row contract before promotion"
        prediction_state = "not_predictive"
        binding = None
        qualifying_events: list[dict] = []
        if row["rung"] == "predictive":
            binding = lineage_state["binding_by_id"][row["id"]]
            qualifying_events = lineage_state["qualifying_events_by_row"][row["id"]]
            if not is_current:
                prediction_state = "invalidated_by_architecture_transition"
            elif row["status"] == "attained":
                prediction_state = "qualified_current_architecture_freeze"
            elif qualifying_events:
                prediction_state = "current_freeze_available_row_not_attained"
            elif binding["candidate_baseline_lineages"]:
                prediction_state = "candidate_requires_version_bound_event"
            else:
                prediction_state = "no_current_candidate"
        observation_index.append(
            {
                "row_id": row["id"],
                "rung": row["rung"],
                "status": row["status"],
                "recorded_architecture_version": row["architecture_version"],
                "current_architecture_version": current_version,
                "promotion_state": promotion_state,
                "prediction_lineage_state": prediction_state,
                "candidate_baseline_lineages": (
                    binding["candidate_baseline_lineages"] if binding else []
                ),
                "historical_only_lineages": (
                    binding["historical_only_lineages"] if binding else []
                ),
                "qualifying_freeze_events": [event["id"] for event in qualifying_events],
                "required_action": action,
                "audits": [_expanded_audit(audits_by_id[audit_id]) for audit_id in audit_ids],
            }
        )

    frozen_index: list[dict] = []
    for row in lineage_rows:
        source = row["source_architecture_version"]
        if (
            row["lineage_status"] == "version_bound"
            and source == current_version
            and source in anchored_architecture_ids
        ):
            effective = "eligible_on_current_version_only"
        elif row["lineage_status"] == "version_bound":
            effective = "invalidated_by_architecture_transition"
        else:
            effective = row["current_version_eligibility"]
        frozen_index.append(
            {
                **row,
                "current_architecture_version": current_version,
                "effective_eligibility": effective,
            }
        )

    return {
        "schema": SCHEMA,
        "issue": ISSUE,
        "current_architecture_version": current_version,
        "current_architecture_snapshot_sha256": versions[-1]["snapshot_sha256"],
        "policy": "Architecture transitions invalidate affected observation promotions and version-bound predictions by default. Until every AV-n pins its own observation-consumer snapshot, any changed premise conservatively invalidates every observation row. Reinstatement requires explicit replay on the successor and a version-matching audit; a prediction requires a new pre-comparison freeze. Historical results remain visible under their original version.",
        "transitions": transitions,
        "observation_rows": observation_index,
        "frozen_prediction_lineages": frozen_index,
    }


def render(data: dict) -> str:
    lines = [
        "# Architecture replay and promotion index",
        "",
        "Generated by `tools/build_architecture_replay.py`. Issue [#741](https://github.com/FloatingPragma/observer-patch-holography/issues/741) established the replay bootstrap; this generated register is the durable replay surface, while standing audit custody remains under [#738](https://github.com/FloatingPragma/observer-patch-holography/issues/738).",
        "",
        data["policy"],
        "",
        f"Current architecture: `{data['current_architecture_version']}` (`{data['current_architecture_snapshot_sha256']}`).",
        "",
        "## Transitions",
        "",
    ]
    if not data["transitions"]:
        lines.append("No successor transition has been registered yet.")
        lines.append("")
    for transition in data["transitions"]:
        lines.extend(
            [
                f"### {transition['from_version']} to {transition['to_version']}",
                "",
                f"Changed normative files: {', '.join(transition['changed_normative_files']) or 'none'}.",
                f"Changed protocol decisions: {', '.join(transition['changed_protocol_decisions']) or 'none'}.",
                f"Changed premises: {', '.join(transition['changed_premises']) or 'none'}.",
                f"Affected observation rows: {', '.join(transition['affected_observation_rows']) or 'none'}.",
                f"Affected frozen lineages: {', '.join(transition['affected_frozen_prediction_lineages']) or 'none'}.",
                "",
            ]
        )
    lines.extend(
        [
            "## Observation promotions",
            "",
            "| Row | Rung | Status | Recorded AV | Promotion | Prediction custody | Audit | Required action |",
            "| --- | --- | --- | --- | --- | --- | --- | --- |",
        ]
    )
    for row in data["observation_rows"]:
        audits = ", ".join(audit["id"] for audit in row["audits"]) or "none"
        lines.append(
            f"| {row['row_id']} | {row['rung']} | {row['status']} | "
            f"{row['recorded_architecture_version']} | {row['promotion_state']} | "
            f"{row['prediction_lineage_state']} | {audits} | {row['required_action']} |"
        )
    lines.extend(
        [
            "",
            "## Frozen-prediction lineages",
            "",
            "| Row | Source AV | Lineage | Current eligibility |",
            "| --- | --- | --- | --- |",
        ]
    )
    for row in data["frozen_prediction_lineages"]:
        lines.append(
            f"| {row['id']} | {row['source_architecture_version'] or 'none'} | "
            f"{row['lineage_status']} | {row['effective_eligibility']} |"
        )
    lines.append("")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)
    data = build()
    payload = json.dumps(data, indent=2, ensure_ascii=False).encode("utf-8") + b"\n"
    surface = render(data).encode("utf-8")
    if args.check:
        if not OUTPUT_PATH.is_file() or OUTPUT_PATH.read_bytes() != payload:
            print("architecture replay: machine index is stale", file=sys.stderr)
            return 1
        if not SURFACE_PATH.is_file() or SURFACE_PATH.read_bytes() != surface:
            print("architecture replay: rendered surface is stale", file=sys.stderr)
            return 1
        print(
            f"architecture replay: {len(data['observation_rows'])} rows, "
            f"{len(data['frozen_prediction_lineages'])} frozen lineages"
        )
        return 0
    OUTPUT_PATH.write_bytes(payload)
    SURFACE_PATH.write_bytes(surface)
    print(f"architecture replay: wrote {OUTPUT_PATH.relative_to(ROOT)} and {SURFACE_PATH.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
