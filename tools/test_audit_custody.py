"""Tamper tests for historical promotion-audit custody."""

from __future__ import annotations

import copy
import hashlib
import json
import subprocess
from pathlib import Path

import pytest

import build_audit_custody as custody


def _data() -> dict:
    return copy.deepcopy(custody.load_json(custody.REGISTER_PATH))


def _git(root: Path, *args: str) -> str:
    return subprocess.run(
        ["git", *args],
        cwd=root,
        check=True,
        stdout=subprocess.PIPE,
        text=True,
    ).stdout.strip()


def _pin(root: Path, revision: str, path: str) -> dict:
    payload = subprocess.run(
        ["git", "show", f"{revision}:{path}"],
        cwd=root,
        check=True,
        stdout=subprocess.PIPE,
    ).stdout
    return {
        "revision": revision,
        "path": path,
        "bytes": len(payload),
        "sha256": hashlib.sha256(payload).hexdigest(),
        "git_blob_sha1": _git(root, "rev-parse", f"{revision}:{path}"),
    }


def _native_origin(root: Path) -> tuple[dict, str]:
    (root / "tracking").mkdir(parents=True)
    _git(root, "init", "-q")
    _git(root, "config", "user.name", "Audit Test")
    _git(root, "config", "user.email", "audit-test@example.invalid")
    (root / "evidence.txt").write_text("stale draft\n", encoding="utf-8")
    _git(root, "add", "evidence.txt")
    _git(root, "commit", "-q", "-m", "baseline evidence")
    baseline = _git(root, "rev-parse", "HEAD")

    (root / "evidence.txt").write_text("independent replay\n", encoding="utf-8")
    (root / "tracking/observation_ledger.json").write_text(
        json.dumps(
            {
                "rows": [
                    {
                        "id": "OL-A1",
                        "status": "attained",
                        "evidence": ["evidence.txt"],
                    }
                ]
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    _git(root, "add", "-A")
    _git(root, "commit", "-q", "-m", "repair evidence")
    repair = _git(root, "rev-parse", "HEAD")
    record = {
        "id": "AUD-NATIVE-TEST",
        "provenance_kind": "native",
        "started_on": "2026-08-13",
        "completed_on": "2026-08-13",
        "baseline_commit": baseline,
        "reviewed_commit": repair,
        "repair_commit": repair,
        "audit_class": "independent fixture audit",
        "scope": "One exact attained-row payload and its evidence.",
        "reviewers": [
            {
                "name": "Fixture Reviewer",
                "task": "native_origin_test",
                "model": "Test model with exact fixture identity",
                "role": "independent replay",
            }
        ],
        "reviewed_rows": ["OL-A1"],
        "promoted_rows": ["OL-A1"],
        "findings": [
            {
                "id": "NATIVE-01",
                "severity": "boundary",
                "disposition": "verified_boundary",
                "summary": "The fixture evidence is exact.",
                "owner_issues": [738],
            }
        ],
        "artifact_pins": [
            _pin(root, repair, "tracking/observation_ledger.json"),
            _pin(root, repair, "evidence.txt"),
        ],
        "qualifies_for": ["attained_status_review"],
        "does_not_qualify_for": ["external_scientific_replication"],
        "limitations": ["Synthetic test fixture only."],
        "source_record_pin": None,
    }
    data = {
        "schema": custody.SCHEMA,
        "issue": custody.ISSUE,
        "record_anchors": [],
        "policy": "Native fixture custody.",
        "records": [record],
    }
    (root / "tracking/audit_custody.json").write_text(
        json.dumps(data, indent=2) + "\n", encoding="utf-8"
    )
    _git(root, "add", "tracking/audit_custody.json")
    _git(root, "commit", "-q", "-m", "append pending native audit")
    return data, _git(root, "rev-parse", "HEAD")


def test_committed_register_and_surface_are_current() -> None:
    data = _data()
    records = custody.validate(data)
    assert custody.SURFACE_PATH.read_bytes() == custody.render(data, records).encode()


def test_historical_digest_tamper_fails() -> None:
    data = _data()
    data["records"][0]["artifact_pins"][0]["sha256"] = "0" * 64
    with pytest.raises(SystemExit, match="historical artifact pin drift"):
        custody.validate(data)


def test_reviewer_model_is_mandatory() -> None:
    data = _data()
    del data["records"][0]["reviewers"][0]["model"]
    with pytest.raises(SystemExit, match="malformed reviewer"):
        custody.validate(data)


def test_reviewer_model_disclosure_cannot_be_rewritten() -> None:
    data = _data()
    data["records"][0]["reviewers"][0]["model"] = "A more impressive model"
    with pytest.raises(SystemExit, match="fixed disclosure"):
        custody.validate(data)


def test_audit_scope_cannot_be_rewritten_after_custody() -> None:
    data = _data()
    data["records"][0]["scope"] = "Everything passed."
    with pytest.raises(SystemExit, match="rewrites the pinned source audit record"):
        custody.validate(data)


def test_legacy_source_cannot_be_rebound_with_coherent_new_pin() -> None:
    data = _data()
    data["records"][0]["source_record_pin"] = {
        "revision": "3ab5bc2064235a740bb5574ea165564e43046bca",
        "path": "tracking/audit_register.json",
        "bytes": 9427,
        "sha256": "0" * 64,
        "git_blob_sha1": "0" * 40,
    }
    with pytest.raises(SystemExit, match="audited immutable migration source"):
        custody.validate(data)


def test_audit_cannot_remove_attained_qualification() -> None:
    data = _data()
    data["records"][0]["qualifies_for"].remove("attained_status_review")
    with pytest.raises(SystemExit, match="rewrites the pinned source audit record"):
        custody.validate(data)


def test_promoted_row_evidence_pin_is_mandatory() -> None:
    data = _data()
    data["records"][0]["artifact_pins"] = [
        pin
        for pin in data["records"][0]["artifact_pins"]
        if pin["path"] != "Lean/EventAlgebra/FiniteBuschGleason.lean"
    ]
    with pytest.raises(SystemExit, match="promoted-row evidence"):
        custody.validate(data)


def test_promoted_row_evidence_must_be_pinned_at_repair_commit(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    data, _ = _native_origin(tmp_path)
    repair = data["records"][0]["repair_commit"]
    parent = _git(tmp_path, "rev-parse", f"{repair}^")
    data["records"][0]["artifact_pins"] = [
        _pin(tmp_path, parent, "evidence.txt")
        if pin["path"] == "evidence.txt"
        else pin
        for pin in data["records"][0]["artifact_pins"]
    ]
    monkeypatch.setattr(custody, "ROOT", tmp_path)
    with pytest.raises(SystemExit, match="pinned at the repair commit"):
        custody.validate(data)


def test_commit_window_must_follow_git_ancestry(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    data, _ = _native_origin(tmp_path)
    repair = data["records"][0]["repair_commit"]
    parent = _git(tmp_path, "rev-parse", f"{repair}^")
    data["records"][0]["baseline_commit"] = repair
    data["records"][0]["reviewed_commit"] = parent
    monkeypatch.setattr(custody, "ROOT", tmp_path)
    with pytest.raises(SystemExit, match="baseline_commit must be an ancestor"):
        custody.validate(data)


def test_promoted_scope_cannot_drop_attained_row() -> None:
    data = _data()
    data["records"][0]["promoted_rows"].pop()
    with pytest.raises(SystemExit, match="rewrites the pinned source audit record"):
        custody.validate(data)


def test_unanchored_native_record_is_pending_and_cannot_qualify(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    data, _ = _native_origin(tmp_path)
    monkeypatch.setattr(custody, "ROOT", tmp_path)
    records = custody.validate(data)
    assert records[0]["_origin_state"] == "pending_origin_anchor"
    assert records[0]["promoted_rows"] == []
    assert "attained_status_review" not in records[0]["qualifies_for"]
    assert records[0]["_declared_promoted_rows"] == ["OL-A1"]


def test_worktree_only_native_anchor_is_not_operational(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    data, origin_revision = _native_origin(tmp_path)
    data["record_anchors"] = [
        {
            "id": "AUD-NATIVE-TEST",
            "origin_revision": origin_revision,
            "record_sha256": hashlib.sha256(
                custody._canonical_record_bytes(data["records"][0])
            ).hexdigest(),
        }
    ]
    monkeypatch.setattr(custody, "ROOT", tmp_path)
    with pytest.raises(SystemExit, match="must already be committed at HEAD"):
        custody.validate(data)


def test_native_record_qualifies_only_after_separate_origin_anchor_commit(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    data, origin_revision = _native_origin(tmp_path)
    data["record_anchors"] = [
        {
            "id": "AUD-NATIVE-TEST",
            "origin_revision": origin_revision,
            "record_sha256": hashlib.sha256(
                custody._canonical_record_bytes(data["records"][0])
            ).hexdigest(),
        }
    ]
    register = tmp_path / "tracking/audit_custody.json"
    register.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    _git(tmp_path, "add", "tracking/audit_custody.json")
    _git(tmp_path, "commit", "-q", "-m", "commit native audit anchor")
    monkeypatch.setattr(custody, "ROOT", tmp_path)
    records = custody.validate(data)
    assert records[0]["_origin_state"] == "native_origin_anchored"
    assert records[0]["promoted_rows"] == ["OL-A1"]
    assert "attained_status_review" in records[0]["qualifies_for"]


def test_anchored_native_record_cannot_be_rewritten(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    data, origin_revision = _native_origin(tmp_path)
    data["record_anchors"] = [
        {
            "id": "AUD-NATIVE-TEST",
            "origin_revision": origin_revision,
            "record_sha256": hashlib.sha256(
                custody._canonical_record_bytes(data["records"][0])
            ).hexdigest(),
        }
    ]
    register = tmp_path / "tracking/audit_custody.json"
    register.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    _git(tmp_path, "add", "tracking/audit_custody.json")
    _git(tmp_path, "commit", "-q", "-m", "commit native audit anchor")
    data["records"][0]["scope"] = "Rewritten after origin."
    monkeypatch.setattr(custody, "ROOT", tmp_path)
    with pytest.raises(SystemExit, match="differs from its origin anchor"):
        custody.validate(data)


def test_native_origin_must_be_on_first_parent_history(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    data, origin_revision = _native_origin(tmp_path)
    main_tip = origin_revision
    parent = _git(tmp_path, "rev-parse", f"{origin_revision}^")
    _git(tmp_path, "checkout", "-q", "-b", "side-origin", parent)
    register_path = tmp_path / "tracking/audit_custody.json"
    register_path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    _git(tmp_path, "add", "tracking/audit_custody.json")
    _git(tmp_path, "commit", "-q", "-m", "side-branch audit origin")
    side_origin = _git(tmp_path, "rev-parse", "HEAD")
    _git(tmp_path, "checkout", "-q", "--detach", main_tip)
    _git(
        tmp_path,
        "merge",
        "--no-ff",
        "-q",
        "-m",
        "merge side origin without making it first-parent custody",
        side_origin,
    )
    data["record_anchors"] = [
        {
            "id": "AUD-NATIVE-TEST",
            "origin_revision": side_origin,
            "record_sha256": hashlib.sha256(
                custody._canonical_record_bytes(data["records"][0])
            ).hexdigest(),
        }
    ]
    monkeypatch.setattr(custody, "ROOT", tmp_path)
    with pytest.raises(SystemExit, match="ancestor of HEAD|first-parent"):
        custody.validate(data)


def test_committed_anchor_cannot_be_erased_then_reintroduced(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    data, origin_revision = _native_origin(tmp_path)
    data["record_anchors"] = [
        {
            "id": "AUD-NATIVE-TEST",
            "origin_revision": origin_revision,
            "record_sha256": hashlib.sha256(
                custody._canonical_record_bytes(data["records"][0])
            ).hexdigest(),
        }
    ]
    register_path = tmp_path / "tracking/audit_custody.json"
    register_path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    _git(tmp_path, "add", "tracking/audit_custody.json")
    _git(tmp_path, "commit", "-q", "-m", "anchor native audit")

    register_path.unlink()
    _git(tmp_path, "add", "-A")
    _git(tmp_path, "commit", "-q", "-m", "erase custody")

    rewritten = copy.deepcopy(data)
    rewritten["record_anchors"] = []
    rewritten["records"][0]["scope"] = "Rewritten after erasing history."
    register_path.write_text(
        json.dumps(rewritten, indent=2) + "\n", encoding="utf-8"
    )
    _git(tmp_path, "add", "tracking/audit_custody.json")
    _git(tmp_path, "commit", "-q", "-m", "reintroduce rewritten audit")
    rewritten_origin = _git(tmp_path, "rev-parse", "HEAD")
    rewritten["record_anchors"] = [
        {
            "id": "AUD-NATIVE-TEST",
            "origin_revision": rewritten_origin,
            "record_sha256": hashlib.sha256(
                custody._canonical_record_bytes(rewritten["records"][0])
            ).hexdigest(),
        }
    ]
    monkeypatch.setattr(custody, "ROOT", tmp_path)
    with pytest.raises(SystemExit, match="disappear and later reappear"):
        custody.validate(rewritten)


def test_duplicate_key_is_rejected(tmp_path) -> None:
    path = tmp_path / "duplicate.json"
    path.write_text('{"schema":"one","schema":"two"}', encoding="utf-8")
    with pytest.raises(SystemExit, match="duplicate JSON key"):
        custody.load_json(path)
