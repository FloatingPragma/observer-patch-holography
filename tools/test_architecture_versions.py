"""Tamper tests for the architecture-version register."""

from __future__ import annotations

import copy
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

import build_architecture_versions as architecture


def _data() -> dict:
    return copy.deepcopy(architecture.load_json(architecture.REGISTER_PATH))


def _git(root: Path, *args: str) -> str:
    result = subprocess.run(
        ["git", *args],
        cwd=root,
        check=True,
        stdout=subprocess.PIPE,
        text=True,
    )
    return result.stdout.strip()


def _git_bytes(root: Path, *args: str) -> bytes:
    return subprocess.run(
        ["git", *args],
        cwd=root,
        check=True,
        stdout=subprocess.PIPE,
    ).stdout


def _materialize_av0_commit(root: Path, data: dict) -> str:
    paths = {item["path"] for item in data["versions"][0]["normative_files"]} | set(
        data["versions"][0]["replay_surfaces"]
    )
    for relative in paths:
        target = root / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes((architecture.ROOT / relative).read_bytes())
    register = root / "tracking/architecture_versions.json"
    register.parent.mkdir(parents=True, exist_ok=True)
    # AV-0 first appeared under the legacy v2 schema, before origin anchors
    # existed.  Reproduce that history so the later v3 anchor is genuinely a
    # separate custody migration rather than a self-referential same-commit pin.
    legacy = copy.deepcopy(data)
    legacy["schema"] = architecture.LEGACY_ROOT_SCHEMA
    legacy.pop("version_anchors")
    register.write_text(json.dumps(legacy, indent=2) + "\n", encoding="utf-8")
    _git(root, "init", "-q")
    _git(root, "config", "user.name", "Architecture Test")
    _git(root, "config", "user.email", "architecture-test@example.invalid")
    _git(root, "add", "-A")
    _git(root, "commit", "-q", "-m", "AV-0 snapshot")
    return _git(root, "rev-parse", "HEAD")


def _av1_transition(root: Path) -> dict:
    data = _data()
    origin_revision = _materialize_av0_commit(root, data)
    data["version_anchors"] = [
        {
            "id": "AV-0",
            "origin_revision": origin_revision,
            "record_sha256": hashlib.sha256(
                json.dumps(
                    data["versions"][0],
                    ensure_ascii=False,
                    sort_keys=True,
                    separators=(",", ":"),
                ).encode()
            ).hexdigest(),
        }
    ]
    register = root / "tracking/architecture_versions.json"
    register.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    _git(root, "add", "tracking/architecture_versions.json")
    _git(root, "commit", "-q", "-m", "migrate and anchor AV-0 custody")
    predecessor_revision = _git(root, "rev-parse", "HEAD")
    successor = copy.deepcopy(data["versions"][0])
    successor["id"] = "AV-1"
    successor["created_on"] = "2026-08-13"
    successor["predecessor_register_revision"] = predecessor_revision
    successor["promotion_status"] = architecture.PENDING_TIP_PROMOTION_STATUS
    changed = root / successor["normative_files"][0]["path"]
    changed.write_bytes(changed.read_bytes() + b"\n# AV-1 test mutation\n")
    changed_payload = changed.read_bytes()
    successor["normative_files"][0]["sha256"] = hashlib.sha256(
        changed_payload
    ).hexdigest()
    successor["normative_files"][0]["git_blob_sha1"] = architecture._git_blob_sha1(
        changed_payload
    )
    successor["protocol_decisions"][0]["decision"] += " AV-1 fixture change."
    successor["snapshot_sha256"] = architecture.snapshot_sha256(successor)
    data["versions"].append(successor)
    data["current_version"] = "AV-1"
    return data


def test_committed_register_and_surface_are_current() -> None:
    data = _data()
    versions = architecture.validate(data)
    assert data["current_version"] == "AV-0"
    assert (
        architecture.render(data, versions).encode()
        == architecture.SURFACE_PATH.read_bytes()
    )


def test_normative_byte_mutation_fails(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    data = _data()
    data["versions"][0]["normative_files"][0]["path"] = "claims/axiom_registry.yaml"
    mutated = tmp_path / "claims" / "axiom_registry.yaml"
    mutated.parent.mkdir(parents=True)
    mutated.write_bytes(
        (architecture.ROOT / "claims/axiom_registry.yaml").read_bytes() + b"\n"
    )
    data["versions"][0]["normative_files"][0]["sha256"] = hashlib.sha256(
        (architecture.ROOT / "claims/axiom_registry.yaml").read_bytes()
    ).hexdigest()
    data["versions"][0]["normative_files"][0]["git_blob_sha1"] = (
        architecture._git_blob_sha1(mutated.read_bytes())
    )
    data["versions"][0]["snapshot_sha256"] = architecture.legacy_snapshot_sha256(
        data["versions"][0]
    )
    root_revision = _materialize_av0_commit(tmp_path, data)
    data["version_anchors"] = [
        {
            "id": "AV-0",
            "origin_revision": root_revision,
            "record_sha256": hashlib.sha256(
                json.dumps(
                    data["versions"][0],
                    ensure_ascii=False,
                    sort_keys=True,
                    separators=(",", ":"),
                ).encode()
            ).hexdigest(),
        }
    ]
    monkeypatch.setattr(architecture, "ROOT", tmp_path)
    with pytest.raises(SystemExit, match="not registered blob"):
        architecture.validate(data)


def test_current_git_blob_pin_mutation_fails() -> None:
    data = _data()
    data["versions"][0]["normative_files"][0]["git_blob_sha1"] = "0" * 40
    data["versions"][0]["snapshot_sha256"] = architecture.snapshot_sha256(
        data["versions"][0]
    )
    with pytest.raises(SystemExit, match="origin anchor"):
        architecture.validate(data)


def test_common_world_status_is_part_of_immutable_snapshot() -> None:
    data = _data()
    data["versions"][0]["common_world_status"] = "proved_inhabited"
    with pytest.raises(SystemExit, match="origin anchor"):
        architecture.validate(data)


def test_promotion_status_is_part_of_immutable_snapshot() -> None:
    data = _data()
    data["versions"][0]["promotion_status"] = "Every row promoted."
    with pytest.raises(SystemExit, match="origin anchor"):
        architecture.validate(data)


def test_av0_cannot_be_rewritten_with_a_recomputed_snapshot() -> None:
    data = _data()
    data["versions"][0]["protocol_decisions"][0]["decision"] += " Rewritten."
    data["versions"][0]["snapshot_sha256"] = architecture.snapshot_sha256(
        data["versions"][0]
    )
    with pytest.raises(SystemExit, match="origin anchor"):
        architecture.validate(data)


def test_av0_origin_cannot_be_reanchored_to_a_later_commit() -> None:
    data = _data()
    data["version_anchors"][0]["origin_revision"] = (
        "68663e9e52a3931c322676a127dd0af144a01de3"
    )
    with pytest.raises(SystemExit, match="audited first appearance"):
        architecture.validate(data)


def test_unknown_premise_fails() -> None:
    data = _data()
    data["versions"][0]["protocol_decisions"][1]["premises"].append("PR-99")
    data["versions"][0]["snapshot_sha256"] = architecture.snapshot_sha256(
        data["versions"][0]
    )
    with pytest.raises(SystemExit, match="origin anchor"):
        architecture.validate(data)


def test_noncontiguous_version_id_fails() -> None:
    data = _data()
    successor = copy.deepcopy(data["versions"][0])
    successor["id"] = "AV-2"
    data["versions"].append(successor)
    data["current_version"] = "AV-2"
    with pytest.raises(SystemExit, match="contiguous"):
        architecture.validate(data)


def test_replay_surface_path_traversal_fails() -> None:
    data = _data()
    data["versions"][0]["replay_surfaces"][0] = "../outside.json"
    data["versions"][0]["snapshot_sha256"] = architecture.legacy_snapshot_sha256(
        data["versions"][0]
    )
    with pytest.raises(SystemExit, match="origin anchor"):
        # The immutable origin binding rejects the rewrite before the path
        # guard in production; the guard remains independently exercised for
        # a legitimate future unanchored tip below.
        architecture.validate(data)


def test_unanchored_successor_replay_surface_path_traversal_fails(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    data = _av1_transition(tmp_path)
    data["versions"][1]["replay_surfaces"][0] = "../outside.json"
    data["versions"][1]["snapshot_sha256"] = architecture.snapshot_sha256(
        data["versions"][1]
    )
    monkeypatch.setattr(architecture, "ROOT", tmp_path)
    with pytest.raises(SystemExit, match="safe repo-relative"):
        architecture.validate(data)


def test_anchored_replay_surface_tree_is_not_accepted_as_blob(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    data = _data()
    origin = _materialize_av0_commit(tmp_path, data)
    data["version_anchors"] = [
        {
            "id": "AV-0",
            "origin_revision": origin,
            "record_sha256": hashlib.sha256(
                json.dumps(
                    data["versions"][0],
                    ensure_ascii=False,
                    sort_keys=True,
                    separators=(",", ":"),
                ).encode()
            ).hexdigest(),
        }
    ]
    candidate = copy.deepcopy(data["versions"][0])
    candidate["replay_surfaces"] = ["tracking"]
    monkeypatch.setattr(architecture, "ROOT", tmp_path)
    with pytest.raises(SystemExit, match="not a Git blob"):
        architecture.replay_surface_bindings(data, candidate, is_current=True)


def test_architecture_register_cannot_disappear_then_reappear(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    data = _data()
    origin = _materialize_av0_commit(tmp_path, data)
    data["version_anchors"] = [
        {
            "id": "AV-0",
            "origin_revision": origin,
            "record_sha256": hashlib.sha256(
                json.dumps(
                    data["versions"][0],
                    ensure_ascii=False,
                    sort_keys=True,
                    separators=(",", ":"),
                ).encode()
            ).hexdigest(),
        }
    ]
    register = tmp_path / "tracking/architecture_versions.json"
    register.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    _git(tmp_path, "add", "tracking/architecture_versions.json")
    _git(tmp_path, "commit", "-q", "-m", "anchor architecture custody")
    register.unlink()
    _git(tmp_path, "add", "-A")
    _git(tmp_path, "commit", "-q", "-m", "erase architecture custody")
    register.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    _git(tmp_path, "add", "tracking/architecture_versions.json")
    _git(tmp_path, "commit", "-q", "-m", "reintroduce architecture custody")
    monkeypatch.setattr(architecture, "ROOT", tmp_path)
    with pytest.raises(SystemExit, match="disappear and later reappear"):
        architecture.validate(data)


def test_av0_to_av1_keeps_historical_bytes_verifiable(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    data = _av1_transition(tmp_path)
    old_blob = data["versions"][0]["normative_files"][0]["git_blob_sha1"]
    old_bytes = _git_bytes(tmp_path, "cat-file", "blob", old_blob)
    live_bytes = (
        tmp_path / data["versions"][1]["normative_files"][0]["path"]
    ).read_bytes()
    assert old_bytes != live_bytes
    assert (
        hashlib.sha256(old_bytes).hexdigest()
        == data["versions"][0]["normative_files"][0]["sha256"]
    )
    monkeypatch.setattr(architecture, "ROOT", tmp_path)
    versions = architecture.validate(data)
    assert [version["id"] for version in versions] == ["AV-0", "AV-1"]


def test_successor_rejects_rewritten_predecessor_snapshot(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    data = _av1_transition(tmp_path)
    data["versions"][0]["protocol_decisions"][0]["decision"] += " Rewritten."
    data["versions"][0]["snapshot_sha256"] = architecture.snapshot_sha256(
        data["versions"][0]
    )
    monkeypatch.setattr(architecture, "ROOT", tmp_path)
    with pytest.raises(SystemExit, match="origin anchor"):
        architecture.validate(data)


def test_unanchored_av1_tip_cannot_claim_promotion_after_self_rehash(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    data = _av1_transition(tmp_path)
    data["versions"][1]["status"] = "inhabited_conditional"
    data["versions"][1]["promotion_status"] = "Every row promoted."
    data["versions"][1]["snapshot_sha256"] = architecture.snapshot_sha256(
        data["versions"][1]
    )
    monkeypatch.setattr(architecture, "ROOT", tmp_path)
    with pytest.raises(SystemExit, match="unanchored current tip"):
        architecture.validate(data)


def test_inhabited_successor_can_be_anchored_before_row_audits(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    data = _av1_transition(tmp_path)
    data["versions"][1]["status"] = "inhabited_conditional"
    data["versions"][1]["common_world_status"] = (
        "A supplied conditional inhabitant exists; issue #740 still owns the "
        "full common-world compatibility proof."
    )
    data["versions"][1]["snapshot_sha256"] = architecture.snapshot_sha256(
        data["versions"][1]
    )
    register = tmp_path / "tracking/architecture_versions.json"
    register.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    monkeypatch.setattr(architecture, "ROOT", tmp_path)
    assert architecture.validate(data)[-1]["status"] == "inhabited_conditional"

    _git(tmp_path, "add", "-A")
    _git(tmp_path, "commit", "-q", "-m", "append inhabited AV-1 tip")
    origin = _git(tmp_path, "rev-parse", "HEAD")
    data["version_anchors"].append(
        {
            "id": "AV-1",
            "origin_revision": origin,
            "record_sha256": hashlib.sha256(
                json.dumps(
                    data["versions"][1],
                    ensure_ascii=False,
                    sort_keys=True,
                    separators=(",", ":"),
                ).encode()
            ).hexdigest(),
        }
    )
    register.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    _git(tmp_path, "add", "tracking/architecture_versions.json")
    _git(tmp_path, "commit", "-q", "-m", "commit AV-1 origin anchor")
    assert architecture.validate(data)[-1]["status"] == "inhabited_conditional"


def test_worktree_only_successor_anchor_is_not_operational(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    data = _av1_transition(tmp_path)
    register = tmp_path / "tracking/architecture_versions.json"
    register.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    _git(tmp_path, "add", "-A")
    _git(tmp_path, "commit", "-q", "-m", "append AV-1 origin")
    origin = _git(tmp_path, "rev-parse", "HEAD")
    data["version_anchors"].append(
        {
            "id": "AV-1",
            "origin_revision": origin,
            "record_sha256": hashlib.sha256(
                json.dumps(
                    data["versions"][1],
                    ensure_ascii=False,
                    sort_keys=True,
                    separators=(",", ":"),
                ).encode()
            ).hexdigest(),
        }
    )
    monkeypatch.setattr(architecture, "ROOT", tmp_path)
    with pytest.raises(SystemExit, match="must already be committed at HEAD"):
        architecture.validate(data)


def test_anchored_successor_cannot_launder_invalid_origin_promotion(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    data = _av1_transition(tmp_path)
    data["versions"][1]["promotion_status"] = "Every row promoted before custody."
    data["versions"][1]["snapshot_sha256"] = architecture.snapshot_sha256(
        data["versions"][1]
    )
    register = tmp_path / "tracking/architecture_versions.json"
    register.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    _git(tmp_path, "add", "-A")
    _git(tmp_path, "commit", "-q", "-m", "invalid AV-1 origin")
    origin = _git(tmp_path, "rev-parse", "HEAD")
    data["version_anchors"].append(
        {
            "id": "AV-1",
            "origin_revision": origin,
            "record_sha256": hashlib.sha256(
                json.dumps(
                    data["versions"][1],
                    ensure_ascii=False,
                    sort_keys=True,
                    separators=(",", ":"),
                ).encode()
            ).hexdigest(),
        }
    )
    monkeypatch.setattr(architecture, "ROOT", tmp_path)
    with pytest.raises(SystemExit, match="pre-anchor promotion-ineligible"):
        architecture.validate(data)


def test_duplicate_json_key_rejected(tmp_path: Path) -> None:
    duplicate = tmp_path / "duplicate.json"
    duplicate.write_text('{"schema":"one","schema":"two"}', encoding="utf-8")
    with pytest.raises(SystemExit, match="duplicate JSON key"):
        architecture.load_json(duplicate)
