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
    register.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    _git(root, "init", "-q")
    _git(root, "config", "user.name", "Architecture Test")
    _git(root, "config", "user.email", "architecture-test@example.invalid")
    _git(root, "add", "-A")
    _git(root, "commit", "-q", "-m", "AV-0 snapshot")
    return _git(root, "rev-parse", "HEAD")


def _av1_transition(root: Path) -> dict:
    data = _data()
    predecessor_revision = _materialize_av0_commit(root, data)
    successor = copy.deepcopy(data["versions"][0])
    successor["id"] = "AV-1"
    successor["created_on"] = "2026-08-13"
    successor["predecessor_register_revision"] = predecessor_revision
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
    mutated = tmp_path / "axiom_registry.yaml"
    mutated.write_bytes(
        (architecture.ROOT / "claims/axiom_registry.yaml").read_bytes() + b"\n"
    )
    data["versions"][0]["normative_files"][0]["path"] = "axiom_registry.yaml"
    data["versions"][0]["normative_files"][0]["sha256"] = hashlib.sha256(
        (architecture.ROOT / "claims/axiom_registry.yaml").read_bytes()
    ).hexdigest()
    data["versions"][0]["normative_files"][0]["git_blob_sha1"] = (
        architecture._git_blob_sha1(mutated.read_bytes())
    )
    monkeypatch.setattr(architecture, "ROOT", tmp_path)
    premise = tmp_path / "tracking/premise_register.json"
    premise.parent.mkdir()
    premise.write_bytes(architecture.PREMISE_PATH.read_bytes())
    monkeypatch.setattr(architecture, "PREMISE_PATH", premise)
    for replay in data["versions"][0]["replay_surfaces"]:
        target = tmp_path / replay
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text("{}\n", encoding="utf-8")
    with pytest.raises(SystemExit, match="digest mismatch"):
        architecture.validate(data)


def test_current_git_blob_pin_mutation_fails() -> None:
    data = _data()
    data["versions"][0]["normative_files"][0]["git_blob_sha1"] = "0" * 40
    data["versions"][0]["snapshot_sha256"] = architecture.snapshot_sha256(
        data["versions"][0]
    )
    with pytest.raises(SystemExit, match="Git blob mismatch"):
        architecture.validate(data)


def test_unknown_premise_fails() -> None:
    data = _data()
    data["versions"][0]["protocol_decisions"][1]["premises"].append("PR-99")
    data["versions"][0]["snapshot_sha256"] = architecture.snapshot_sha256(
        data["versions"][0]
    )
    with pytest.raises(SystemExit, match="unknown premises"):
        architecture.validate(data)


def test_noncontiguous_version_id_fails() -> None:
    data = _data()
    successor = copy.deepcopy(data["versions"][0])
    successor["id"] = "AV-2"
    data["versions"].append(successor)
    data["current_version"] = "AV-2"
    with pytest.raises(SystemExit, match="contiguous"):
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
    with pytest.raises(SystemExit, match="append-only"):
        architecture.validate(data)


def test_duplicate_json_key_rejected(tmp_path: Path) -> None:
    duplicate = tmp_path / "duplicate.json"
    duplicate.write_text('{"schema":"one","schema":"two"}', encoding="utf-8")
    with pytest.raises(SystemExit, match="duplicate JSON key"):
        architecture.load_json(duplicate)
