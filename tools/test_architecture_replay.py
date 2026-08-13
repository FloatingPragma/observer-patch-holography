"""Cross-surface and transition tests for the architecture replay index."""

from __future__ import annotations

import copy
import hashlib
import json
import shutil
import subprocess
from pathlib import Path

import pytest

import build_architecture_replay as replay
import build_audit_custody as audit_custody
import build_architecture_versions as architecture
import build_observation_ledger as observation
import prediction_lineage_custody


def _git(root: Path, *args: str) -> str:
    return subprocess.run(
        ["git", *args],
        cwd=root,
        check=True,
        stdout=subprocess.PIPE,
        text=True,
    ).stdout.strip()


def test_committed_replay_outputs_are_current() -> None:
    data = replay.build()
    assert replay.OUTPUT_PATH.read_text(encoding="utf-8") == replay.json.dumps(
        data, indent=2, ensure_ascii=False
    ) + "\n"
    assert replay.SURFACE_PATH.read_text(encoding="utf-8") == replay.render(data)


def test_protocol_change_invalidates_every_observation_and_bound_prediction() -> None:
    architecture_data = architecture.load_json(architecture.REGISTER_PATH)
    predecessor = copy.deepcopy(architecture_data["versions"][0])
    successor = copy.deepcopy(predecessor)
    successor["id"] = "AV-1"
    successor["protocol_decisions"][0]["decision"] += " Changed semantics."
    ledger = observation.load_ledger(observation.LEDGER_PATH)
    rows = observation.validate(ledger)
    lineages = [
        {
            "id": "FZ-TEST",
            "lineage_status": "version_bound",
            "source_architecture_version": "AV-0",
        }
    ]
    diff = replay.semantic_diff(
        predecessor,
        successor,
        rows,
        lineages,
        successor_is_current=True,
    )
    assert diff["changed_protocol_decisions"] == ["PD-01"]
    assert diff["affected_observation_rows"] == [row["id"] for row in rows]
    assert diff["affected_frozen_prediction_lineages"] == ["FZ-TEST"]


def test_unknown_normative_semantics_file_change_fails_closed_globally() -> None:
    architecture_data = architecture.load_json(architecture.REGISTER_PATH)
    predecessor = copy.deepcopy(architecture_data["versions"][0])
    successor = copy.deepcopy(predecessor)
    successor["id"] = "AV-1"
    successor["normative_files"].append(
        {
            "path": "new/physical_semantics.json",
            "role": "future normative semantics",
            "sha256": "0" * 64,
            "git_blob_sha1": "0" * 40,
        }
    )
    rows = observation.validate(observation.load_ledger(observation.LEDGER_PATH))
    diff = replay.semantic_diff(
        predecessor,
        successor,
        rows,
        [],
        successor_is_current=True,
    )
    assert diff["affected_observation_rows"] == [row["id"] for row in rows]
    assert diff["observation_invalidation_reason"] == (
        "global_architecture_or_promotion_semantics_changed"
    )


def test_changed_premise_invalidates_every_row_without_versioned_consumer_snapshot(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    architecture_data = architecture.load_json(architecture.REGISTER_PATH)
    predecessor = copy.deepcopy(architecture_data["versions"][0])
    successor = copy.deepcopy(predecessor)
    successor["id"] = "AV-1"
    before = json.dumps(
        {"rows": [{"id": "PR-01", "statement": "before"}]}
    ).encode()
    after = json.dumps(
        {"rows": [{"id": "PR-01", "statement": "after"}]}
    ).encode()

    def pinned_payload(version: dict, path: str, *, is_current: bool) -> bytes:
        assert path == "tracking/premise_register.json"
        return after if version["id"] == "AV-1" else before

    monkeypatch.setattr(replay, "_pinned_payload", pinned_payload)
    rows = observation.validate(observation.load_ledger(observation.LEDGER_PATH))
    diff = replay.semantic_diff(
        predecessor,
        successor,
        rows,
        [],
        successor_is_current=True,
    )
    assert diff["changed_premises"] == ["PR-01"]
    assert diff["affected_observation_rows"] == [row["id"] for row in rows]
    assert diff["observation_invalidation_reason"] == (
        "changed_premise_without_versioned_consumer_snapshot"
    )


def test_real_temp_git_av0_to_av1_build_invalidates_old_promotions(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Exercise the complete builders against a real successor Git manifest."""

    source_root = replay.ROOT
    repo = tmp_path / "repo"
    subprocess.run(
        ["git", "clone", "-q", "--shared", str(source_root), str(repo)],
        check=True,
    )

    current_files = [
        "tracking/architecture_versions.json",
        "tracking/audit_custody.json",
        "tracking/observation_ledger.json",
        "claims/frozen_prediction_register.json",
        "claims/frozen_prediction_architecture_lineages.json",
    ]
    ledger = json.loads(
        (source_root / "tracking/observation_ledger.json").read_text(encoding="utf-8")
    )
    current_files.extend(
        evidence
        for row in ledger["rows"]
        for evidence in row["evidence"]
        if evidence not in current_files
    )
    for relative in current_files:
        source = source_root / relative
        target = repo / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
    _git(repo, "config", "user.name", "Replay Test")
    _git(repo, "config", "user.email", "replay-test@example.invalid")
    _git(repo, "add", "-A")
    # The source checkout may already have every copied custody surface committed.
    # Keep the fixture independent of unrelated source-worktree dirtiness while
    # still creating the explicit test baseline commit used below.
    _git(repo, "commit", "-q", "--allow-empty", "-m", "materialize AV-0 custody")
    predecessor_revision = _git(repo, "rev-parse", "HEAD")

    register_path = repo / "tracking/architecture_versions.json"
    architecture_data = json.loads(register_path.read_text(encoding="utf-8"))
    premise_path = repo / "tracking/premise_register.json"
    premise_data = json.loads(premise_path.read_text(encoding="utf-8"))
    premise_data["rows"][0]["statement"] += " AV-1 integration-test change."
    premise_path.write_text(
        json.dumps(premise_data, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    premise_payload = premise_path.read_bytes()

    successor = copy.deepcopy(architecture_data["versions"][0])
    successor["id"] = "AV-1"
    successor["created_on"] = "2026-08-13"
    successor["predecessor_register_revision"] = predecessor_revision
    successor["status"] = "exploratory_uninhabited"
    successor["promotion_status"] = architecture.PENDING_TIP_PROMOTION_STATUS
    premise_pin = next(
        item
        for item in successor["normative_files"]
        if item["path"] == "tracking/premise_register.json"
    )
    premise_pin["sha256"] = hashlib.sha256(premise_payload).hexdigest()
    premise_pin["git_blob_sha1"] = architecture._git_blob_sha1(premise_payload)
    successor["snapshot_sha256"] = architecture.snapshot_sha256(successor)
    architecture_data["versions"].append(successor)
    architecture_data["current_version"] = "AV-1"
    register_path.write_text(
        json.dumps(architecture_data, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    monkeypatch.setattr(architecture, "ROOT", repo)
    monkeypatch.setattr(architecture, "REGISTER_PATH", register_path)
    monkeypatch.setattr(
        architecture, "PREMISE_PATH", repo / "tracking/premise_register.json"
    )
    monkeypatch.setattr(audit_custody, "ROOT", repo)
    monkeypatch.setattr(
        audit_custody, "REGISTER_PATH", repo / "tracking/audit_custody.json"
    )
    monkeypatch.setattr(observation, "ROOT", repo)
    monkeypatch.setattr(
        observation, "LEDGER_PATH", repo / "tracking/observation_ledger.json"
    )
    monkeypatch.setattr(
        observation, "PREMISE_REGISTER_PATH", repo / "tracking/premise_register.json"
    )
    monkeypatch.setattr(observation, "ARCHITECTURE_REGISTER_PATH", register_path)
    monkeypatch.setattr(
        observation, "AUDIT_CUSTODY_PATH", repo / "tracking/audit_custody.json"
    )
    monkeypatch.setattr(
        observation,
        "PREDICTION_LINEAGE_PATH",
        repo / "claims/frozen_prediction_architecture_lineages.json",
    )
    monkeypatch.setattr(
        observation,
        "FROZEN_PREDICTION_PATH",
        repo / "claims/frozen_prediction_register.json",
    )
    monkeypatch.setattr(replay, "ROOT", repo)
    monkeypatch.setattr(
        replay,
        "LINEAGE_PATH",
        repo / "claims/frozen_prediction_architecture_lineages.json",
    )
    monkeypatch.setattr(
        replay, "FROZEN_PATH", repo / "claims/frozen_prediction_register.json"
    )

    built = replay.build()
    assert built["current_architecture_version"] == "AV-1"
    assert len(built["transitions"]) == 1
    transition = built["transitions"][0]
    assert transition["changed_premises"] == ["PR-01"]
    assert transition["affected_observation_rows"] == [
        row["id"] for row in ledger["rows"]
    ]
    attained = {
        row["row_id"]: row["promotion_state"]
        for row in built["observation_rows"]
        if row["status"] == "attained"
    }
    assert attained == {
        "OL-C1": "invalidated_on_current_version",
        "OL-E2": "invalidated_on_current_version",
        "OL-E3": "invalidated_on_current_version",
        "OL-G3": "invalidated_on_current_version",
    }


def test_legacy_frozen_target_cannot_be_mislabeled_version_bound() -> None:
    frozen = replay.load_json(replay.FROZEN_PATH)
    lineages = replay.load_json(replay.LINEAGE_PATH)
    lineages["rows"][0]["lineage_status"] = "version_bound"
    lineages["rows"][0]["source_architecture_version"] = "AV-0"
    lineages["rows"][0]["current_version_eligibility"] = (
        "eligible_on_source_version_only"
    )
    rows = observation.validate(observation.load_ledger(observation.LEDGER_PATH))
    with pytest.raises(SystemExit, match="baseline lineage rows drifted"):
        prediction_lineage_custody.validate(
            lineages,
            frozen,
            architecture.load_json(architecture.REGISTER_PATH),
            rows,
        )


def test_frozen_content_binding_tamper_fails() -> None:
    frozen = replay.load_json(replay.FROZEN_PATH)
    lineages = replay.load_json(replay.LINEAGE_PATH)
    lineages["rows"][-1]["registered_content_sha256"] = "0" * 64
    rows = observation.validate(observation.load_ledger(observation.LEDGER_PATH))
    with pytest.raises(SystemExit, match="baseline lineage rows drifted"):
        prediction_lineage_custody.validate(
            lineages,
            frozen,
            architecture.load_json(architecture.REGISTER_PATH),
            rows,
        )


def test_frozen_row_rewrite_behind_old_declared_hash_fails() -> None:
    frozen = replay.load_json(replay.FROZEN_PATH)
    lineages = replay.load_json(replay.LINEAGE_PATH)
    row = next(item for item in frozen["rows"] if item["id"] == "FZ-10")
    original_hash = row["content_sha256"]
    row["content"] = "A stronger post-hoc target"
    row["kill_band"] = "A weaker post-hoc decision rule"
    assert row["content_sha256"] == original_hash
    rows = observation.validate(observation.load_ledger(observation.LEDGER_PATH))
    with pytest.raises(SystemExit, match="lineage bootstrap payload"):
        prediction_lineage_custody.validate(
            lineages,
            frozen,
            architecture.load_json(architecture.REGISTER_PATH),
            rows,
        )
