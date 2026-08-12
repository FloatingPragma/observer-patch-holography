"""Drift and tamper rejection tests for the V3 emergent-instrument register.

The committed register carries INS-01 in SPECIFIED status. The FROZEN, RUNNING,
and verdict transitions are covered abstractly through fixture rows, so a
promotion of INS-01 lands as a register edit with no tool change.
"""

from __future__ import annotations

import copy
import hashlib
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

import build_instrument_register as register_tool

FIXTURE_RULE = (
    "Frozen fixture rule: the arm passes on the reference reading, with"
    " verdict REPLICATED, FAILED, or INCONCLUSIVE."
)
FIXTURE_FREEZE = [
    {
        "path": "oph-physics-sim/receipts/fixture_freeze.json",
        "sha256": "0" * 64,
    }
]
FIXTURE_RECEIPTS = [
    {
        "path": "oph-physics-sim/receipts/fixture_verdict.json",
        "sha256": "1" * 64,
    }
]


def _register() -> dict:
    return copy.deepcopy(register_tool.load_json(register_tool.REGISTER_PATH))


def _fixture_row(
    status: str,
    *,
    row_id: str = "INS-02",
    frozen: bool = True,
    receipts: bool = False,
) -> dict:
    return {
        "id": row_id,
        "title": "Fixture transition instrument",
        "owning_issue": 737,
        "ledger_row": "OL-A1",
        "status": status,
        "spec_pointer": "plan/OL_A1_SIGNATURE_REPLICATION_SPEC.md (oph-meta planning workspace)",
        "decision_rule": FIXTURE_RULE,
        "seeds_policy": "One declared seed and five declared replicate ids.",
        "controls": ["Fixture control: matched-density arm."],
        "freeze_artifacts": copy.deepcopy(FIXTURE_FREEZE) if frozen else [],
        "frozen_utc": "2026-08-11T00:00:00Z" if frozen else None,
        "verdict_receipts": copy.deepcopy(FIXTURE_RECEIPTS) if receipts else [],
    }


def test_committed_register_validates() -> None:
    register = _register()
    rows = register_tool.validate(register)
    assert rows[0]["id"] == "INS-01"
    assert rows[0]["ledger_row"] == "OL-A1"
    assert rows[0]["status"] == "FAILED"
    assert len(rows[0]["freeze_artifacts"]) == 3
    assert rows[0]["frozen_utc"] is not None
    assert len(rows[0]["verdict_receipts"]) == 1


def test_rebuild_parity_with_committed_surface() -> None:
    register = _register()
    rows = register_tool.validate(register)
    rendered = register_tool.render(register, rows).encode("utf-8")
    committed = register_tool.SURFACE_PATH.read_bytes()
    assert rendered == committed


def test_check_mode_passes_on_committed_artifacts() -> None:
    result = subprocess.run(
        [
            sys.executable,
            str(Path(__file__).parent / "build_instrument_register.py"),
            "--check",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr


def test_unknown_status_rejected() -> None:
    register = _register()
    register["rows"][0]["status"] = "PENDING"
    with pytest.raises(SystemExit, match="status must be one of"):
        register_tool.validate(register)


def test_specified_row_with_freeze_artifacts_rejected() -> None:
    register = _register()
    row = _fixture_row("SPECIFIED", frozen=False)
    register["rows"].append(row)
    row["freeze_artifacts"] = copy.deepcopy(FIXTURE_FREEZE)
    with pytest.raises(SystemExit, match="unfrozen row carries no freeze artifacts"):
        register_tool.validate(register)


def test_specified_row_with_freeze_time_rejected() -> None:
    register = _register()
    row = _fixture_row("SPECIFIED", frozen=False)
    register["rows"].append(row)
    row["frozen_utc"] = "2026-08-11T00:00:00Z"
    with pytest.raises(SystemExit, match="SPECIFIED row carries no freeze time"):
        register_tool.validate(register)


def test_specified_row_with_verdict_receipts_rejected() -> None:
    register = _register()
    row = _fixture_row("SPECIFIED", frozen=False)
    register["rows"].append(row)
    row["verdict_receipts"] = copy.deepcopy(FIXTURE_RECEIPTS)
    with pytest.raises(SystemExit, match="SPECIFIED row carries no verdict receipts"):
        register_tool.validate(register)


def test_future_freeze_time_rejected() -> None:
    register = _register()
    row = _fixture_row("FROZEN")
    row["frozen_utc"] = "2126-01-01T00:00:00Z"
    register["rows"].append(row)
    with pytest.raises(SystemExit, match="cannot be in the future"):
        register_tool.validate(register)


def test_malformed_freeze_time_rejected() -> None:
    register = _register()
    row = _fixture_row("FROZEN")
    row["frozen_utc"] = "2026-08-11 00:00:00"
    register["rows"].append(row)
    with pytest.raises(SystemExit, match="ending in Z"):
        register_tool.validate(register)


def test_frozen_fixture_row_validates() -> None:
    register = _register()
    register["rows"].append(_fixture_row("FROZEN"))
    rows = register_tool.validate(register)
    assert rows[1]["status"] == "FROZEN"


def test_running_fixture_row_validates() -> None:
    register = _register()
    register["rows"].append(_fixture_row("RUNNING"))
    rows = register_tool.validate(register)
    assert rows[1]["status"] == "RUNNING"


def test_replicated_fixture_row_validates() -> None:
    register = _register()
    register["rows"].append(_fixture_row("REPLICATED", receipts=True))
    rows = register_tool.validate(register)
    assert rows[1]["verdict_receipts"]


def test_failed_fixture_row_validates() -> None:
    register = _register()
    register["rows"].append(_fixture_row("FAILED", receipts=True))
    rows = register_tool.validate(register)
    assert rows[1]["status"] == "FAILED"


def test_inconclusive_fixture_row_validates() -> None:
    register = _register()
    register["rows"].append(_fixture_row("INCONCLUSIVE", receipts=True))
    rows = register_tool.validate(register)
    assert rows[1]["status"] == "INCONCLUSIVE"


def test_frozen_row_without_artifacts_rejected() -> None:
    register = _register()
    row = _fixture_row("FROZEN")
    row["freeze_artifacts"] = []
    register["rows"].append(row)
    with pytest.raises(SystemExit, match="requires freeze artifacts"):
        register_tool.validate(register)


def test_frozen_row_with_receipts_rejected() -> None:
    register = _register()
    row = _fixture_row("FROZEN", receipts=True)
    register["rows"].append(row)
    with pytest.raises(SystemExit, match="FROZEN row carries no verdict receipts"):
        register_tool.validate(register)


def test_verdict_row_without_receipts_rejected() -> None:
    register = _register()
    register["rows"].append(_fixture_row("REPLICATED", receipts=False))
    with pytest.raises(SystemExit, match="requires verdict receipts"):
        register_tool.validate(register)


def test_verdict_row_without_freeze_rejected() -> None:
    register = _register()
    row = _fixture_row("FAILED", frozen=False, receipts=True)
    register["rows"].append(row)
    with pytest.raises(SystemExit, match="requires a freeze time"):
        register_tool.validate(register)


def test_frozen_rule_missing_verdict_label_rejected() -> None:
    register = _register()
    row = _fixture_row("FROZEN")
    row["decision_rule"] = "Frozen fixture rule with REPLICATED and FAILED only."
    register["rows"].append(row)
    with pytest.raises(SystemExit, match="INCONCLUSIVE"):
        register_tool.validate(register)


def test_void_row_with_receipts_rejected() -> None:
    register = _register()
    row = _fixture_row("VOID", receipts=True)
    register["rows"].append(row)
    with pytest.raises(SystemExit, match="VOID row carries no verdict receipts"):
        register_tool.validate(register)


def test_void_row_coherence() -> None:
    register = _register()
    register["rows"].append(_fixture_row("VOID", frozen=True))
    rows = register_tool.validate(register)
    assert rows[1]["status"] == "VOID"

    register = _register()
    register["rows"].append(_fixture_row("VOID", frozen=False))
    rows = register_tool.validate(register)
    assert rows[1]["frozen_utc"] is None

    register = _register()
    row = _fixture_row("VOID", frozen=True)
    row["freeze_artifacts"] = []
    register["rows"].append(row)
    with pytest.raises(SystemExit, match="keeps its freeze artifacts"):
        register_tool.validate(register)


def test_bad_sha256_rejected() -> None:
    register = _register()
    row = _fixture_row("FROZEN")
    row["freeze_artifacts"][0]["sha256"] = "DEADBEEF"
    register["rows"].append(row)
    with pytest.raises(SystemExit, match="lowercase SHA-256"):
        register_tool.validate(register)


def test_absolute_artifact_path_rejected() -> None:
    register = _register()
    row = _fixture_row("FROZEN")
    row["freeze_artifacts"][0]["path"] = "/etc/fixture_freeze.json"
    register["rows"].append(row)
    with pytest.raises(SystemExit, match="must be relative"):
        register_tool.validate(register)


def test_duplicate_artifact_path_rejected() -> None:
    register = _register()
    row = _fixture_row("FROZEN")
    row["freeze_artifacts"].append(copy.deepcopy(FIXTURE_FREEZE[0]))
    register["rows"].append(row)
    with pytest.raises(SystemExit, match="duplicate artifact path"):
        register_tool.validate(register)


def test_in_repo_artifact_hash_binding() -> None:
    pinned = "claims/emergent_instrument_register.json"
    digest = hashlib.sha256(
        (register_tool.ROOT / pinned).read_bytes()
    ).hexdigest()

    register = _register()
    row = _fixture_row("FROZEN")
    row["freeze_artifacts"] = [{"path": pinned, "sha256": digest}]
    register["rows"].append(row)
    rows = register_tool.validate(register)
    assert rows[1]["freeze_artifacts"][0]["path"] == pinned

    register = _register()
    row = _fixture_row("FROZEN")
    row["freeze_artifacts"] = [{"path": pinned, "sha256": "2" * 64}]
    register["rows"].append(row)
    with pytest.raises(SystemExit, match="hash mismatch"):
        register_tool.validate(register)


def test_unknown_ledger_row_rejected() -> None:
    register = _register()
    register["rows"][0]["ledger_row"] = "OL-I9"
    with pytest.raises(SystemExit, match="is not on the observation ledger"):
        register_tool.validate(register)


def test_malformed_ledger_row_rejected() -> None:
    register = _register()
    register["rows"][0]["ledger_row"] = "OLA1"
    with pytest.raises(SystemExit, match="ledger_row must match"):
        register_tool.validate(register)


def test_lane_out_of_range_rejected() -> None:
    register = _register()
    register["rows"][0]["owning_issue"] = 607
    with pytest.raises(SystemExit, match="owning_issue must lie in"):
        register_tool.validate(register)


def test_first_row_must_be_ins01() -> None:
    register = _register()
    register["rows"][0]["id"] = "INS-02"
    with pytest.raises(SystemExit, match="first row must be INS-01"):
        register_tool.validate(register)


def test_nonascending_ids_rejected() -> None:
    register = _register()
    register["rows"].append(_fixture_row("FROZEN", row_id="INS-01"))
    with pytest.raises(SystemExit, match="ascend strictly"):
        register_tool.validate(register)


def test_extra_row_key_rejected() -> None:
    register = _register()
    register["rows"][0]["comment"] = "off-schema"
    with pytest.raises(SystemExit, match="keys mismatch"):
        register_tool.validate(register)


def test_banned_dash_rejected() -> None:
    register = _register()
    register["rows"][0]["title"] = "Signature replication — instrument"
    with pytest.raises(SystemExit, match="banned dash character"):
        register_tool.validate(register)


def test_wrong_schema_rejected() -> None:
    register = _register()
    register["schema"] = "oph.instrument_register.v0"
    with pytest.raises(SystemExit, match="schema must equal"):
        register_tool.validate(register)


def test_stale_surface_fails_check(tmp_path, monkeypatch) -> None:
    stale = tmp_path / "INSTRUMENT_REGISTER_V3.md"
    stale.write_bytes(b"stale\n")
    monkeypatch.setattr(register_tool, "SURFACE_PATH", stale)
    assert register_tool.main(["--check"]) == 1


def test_missing_in_repo_artifact_path_rejected():
    register = _register()
    row = _fixture_row("FROZEN")
    row["freeze_artifacts"] = [
        {"path": "claims/does_not_exist_receipt.json", "sha256": "4" * 64}
    ]
    register["rows"].append(row)
    with pytest.raises(SystemExit, match="does not exist"):
        register_tool.validate(register)


def test_external_pointer_path_accepted():
    register = _register()
    row = _fixture_row("FROZEN")
    row["freeze_artifacts"] = [
        {"path": "oph-physics-sim/docs/some_prereg.md", "sha256": "4" * 64}
    ]
    register["rows"].append(row)
    register_tool.validate(register)
