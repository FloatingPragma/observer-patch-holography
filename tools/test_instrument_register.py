"""Drift and tamper rejection tests for the V3 emergent-instrument register."""

from __future__ import annotations

import copy
import hashlib
import json
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


def _local_artifact(path: str) -> dict[str, str]:
    return {
        "path": path,
        "sha256": hashlib.sha256((register_tool.ROOT / path).read_bytes()).hexdigest(),
    }


FIXTURE_FREEZE = [_local_artifact("docs/STYLE_GUIDE.md")]
FIXTURE_RECEIPTS = [_local_artifact("docs/AXIOM_REFERENCE.md")]


def _register() -> dict:
    return copy.deepcopy(register_tool.load_json(register_tool.REGISTER_PATH))


def _fixture_row(
    status: str,
    *,
    row_id: str = "INS-03",
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
        "custody_repository": None,
        "freeze_artifacts": copy.deepcopy(FIXTURE_FREEZE) if frozen else [],
        "frozen_utc": "2026-08-11T00:00:00Z" if frozen else None,
        "lineage_predecessor": "INS-02",
        "promotion_eligibility": {
            "state": "INELIGIBLE",
            "architecture_version": None,
            "audit_ids": [],
        },
        "limitations": "Fixture-only validation row with no scientific claim.",
        "verdict_receipts": copy.deepcopy(FIXTURE_RECEIPTS) if receipts else [],
    }


def _replicated_ins02(register: dict) -> dict:
    row = register["rows"][1]
    row["status"] = "REPLICATED"
    row["freeze_artifacts"] = copy.deepcopy(FIXTURE_FREEZE)
    row["frozen_utc"] = "2026-08-11T00:00:00Z"
    row["verdict_receipts"] = copy.deepcopy(FIXTURE_RECEIPTS)
    row["promotion_eligibility"] = {
        "state": "ELIGIBLE",
        "architecture_version": "AV-0",
        "audit_ids": ["AUD-FIXTURE"],
    }
    return row


def _ledger_rows(*, ol_a1_status: str = "owed", architecture: str = "AV-0") -> dict:
    rows = copy.deepcopy(register_tool.load_ledger_rows())
    rows["OL-A1"]["status"] = ol_a1_status
    rows["OL-A1"]["architecture_version"] = architecture
    return rows


def test_committed_register_validates() -> None:
    register = _register()
    rows = register_tool.validate(register)
    assert rows[0]["id"] == "INS-01"
    assert rows[0]["ledger_row"] == "OL-A1"
    assert rows[0]["status"] == "FAILED"
    assert len(rows[0]["freeze_artifacts"]) == 3
    assert rows[0]["frozen_utc"] is not None
    assert len(rows[0]["verdict_receipts"]) == 17
    assert rows[0]["custody_repository"] == {
        "url": register_tool.SIMULATOR_REPOSITORY_URL,
        "commit": register_tool.SIMULATOR_COMMIT,
    }
    limitations = rows[0]["limitations"]
    assert "no raw feature matrices or fit captures" in limitations
    assert "not an independent observable recomputation" in limitations
    assert "reachable from the configured GitHub remote's main branch" in limitations
    assert rows[1]["id"] == "INS-02"
    assert rows[1]["status"] == "SPECIFIED"
    assert rows[1]["freeze_artifacts"] == []
    assert rows[1]["verdict_receipts"] == []
    assert rows[1]["frozen_utc"] is None
    assert "not a preregistration or authorization to run" in rows[1]["limitations"]
    assert rows[1]["lineage_predecessor"] == "INS-01"
    assert rows[1]["promotion_eligibility"] == {
        "state": "INELIGIBLE",
        "architecture_version": None,
        "audit_ids": [],
    }
    assert register["ledger_controls"][0]["controlling_instrument"] == "INS-01"


def test_ins02_design_contract_is_explicit_and_unfrozen() -> None:
    register = _register()
    row = register["rows"][1]
    design = (register_tool.ROOT / row["spec_pointer"]).read_text(encoding="utf-8")
    design_flat = " ".join(design.split())

    assert row["status"] == "SPECIFIED"
    assert row["freeze_artifacts"] == []
    assert row["frozen_utc"] is None
    assert row["verdict_receipts"] == []
    assert "16,384; 65,536; 131,072; 262,144" in design
    assert "absolute support size" in design
    assert "one value fixed at the future preregistration" in design
    assert "target-blind prospective calculation" in design
    assert "rather than assuming that eight" in design
    assert "paired seed block" in design_flat
    assert "randomize execution order" in design_flat
    assert "same fixed estimator" in design_flat
    assert "**destruction:**" in design
    assert "**sham equivalence:**" in design
    assert "**synthetic sensitivity:**" in design
    assert "prospective precision and power" in row["seeds_policy"]
    assert "paired blocks" in row["seeds_policy"]
    assert "randomized within each block" in row["seeds_policy"]
    assert "quantitative ancestry-destruction" in row["decision_rule"]
    assert "through one fixed estimator" in row["decision_rule"]


def test_rebuild_parity_with_committed_surface() -> None:
    register = _register()
    rows = register_tool.validate(register)
    rendered = register_tool.render(register, rows).encode("utf-8")
    committed = register_tool.SURFACE_PATH.read_bytes()
    assert rendered == committed


def test_intro_distinguishes_negative_control_from_positive_eligibility() -> None:
    register = _register()
    rendered = register_tool.render(register, register_tool.validate(register))
    assert "completed decisive instrument" in rendered
    assert "positive qualification additionally requires" in rendered
    assert "completed, eligible instrument" not in rendered


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
    assert rows[-1]["status"] == "FROZEN"


def test_running_fixture_row_validates() -> None:
    register = _register()
    register["rows"].append(_fixture_row("RUNNING"))
    rows = register_tool.validate(register)
    assert rows[-1]["status"] == "RUNNING"


def test_replicated_fixture_row_validates() -> None:
    register = _register()
    register["rows"].append(_fixture_row("REPLICATED", receipts=True))
    rows = register_tool.validate(register)
    assert rows[-1]["verdict_receipts"]


def test_failed_fixture_row_validates() -> None:
    register = _register()
    register["rows"].append(_fixture_row("FAILED", receipts=True))
    rows = register_tool.validate(register)
    assert rows[-1]["status"] == "FAILED"


def test_inconclusive_fixture_row_validates() -> None:
    register = _register()
    register["rows"].append(_fixture_row("INCONCLUSIVE", receipts=True))
    rows = register_tool.validate(register)
    assert rows[-1]["status"] == "INCONCLUSIVE"


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
    assert rows[-1]["status"] == "VOID"

    register = _register()
    register["rows"].append(_fixture_row("VOID", frozen=False))
    rows = register_tool.validate(register)
    assert rows[-1]["frozen_utc"] is None

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
    digest = hashlib.sha256((register_tool.ROOT / pinned).read_bytes()).hexdigest()

    register = _register()
    row = _fixture_row("FROZEN")
    row["freeze_artifacts"] = [{"path": pinned, "sha256": digest}]
    register["rows"].append(row)
    rows = register_tool.validate(register)
    assert rows[-1]["freeze_artifacts"][0]["path"] == pinned

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


def test_duplicate_json_key_rejected(tmp_path: Path) -> None:
    duplicate = tmp_path / "duplicate.json"
    duplicate.write_text('{"schema":"first","schema":"second"}\n', encoding="utf-8")
    with pytest.raises(SystemExit, match="duplicate JSON key 'schema'"):
        register_tool.load_json(duplicate)


def test_nested_duplicate_json_key_rejected(tmp_path: Path) -> None:
    duplicate = tmp_path / "nested-duplicate.json"
    duplicate.write_text(
        json.dumps({"row": {"id": "INS-01"}}).replace(
            '"id": "INS-01"', '"id": "INS-01", "id": "INS-02"'
        )
        + "\n",
        encoding="utf-8",
    )
    with pytest.raises(SystemExit, match="duplicate JSON key 'id'"):
        register_tool.load_json(duplicate)


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


def test_external_pointer_without_custody_rejected():
    register = _register()
    row = _fixture_row("FROZEN")
    row["freeze_artifacts"] = [
        {"path": "oph-physics-sim/docs/some_prereg.md", "sha256": "4" * 64}
    ]
    register["rows"].append(row)
    with pytest.raises(
        SystemExit, match="external artifacts require custody_repository"
    ):
        register_tool.validate(register)


def test_ins01_abbreviated_commit_rejected():
    register = _register()
    register["rows"][0]["custody_repository"]["commit"] = "42aa966"
    with pytest.raises(SystemExit, match="full 40-character SHA-1"):
        register_tool.validate(register)


def test_ins01_wrong_repository_rejected():
    if not (register_tool.ROOT.parent / "oph-physics-sim" / ".git").exists():
        pytest.skip("sibling simulator custody checkout is absent")
    register = _register()
    register["rows"][0]["custody_repository"]["url"] = (
        "https://github.com/example/oph-physics-sim"
    )
    with pytest.raises(SystemExit, match="sibling simulator origin"):
        register_tool.validate(register)


def test_ins01_manifest_is_mandatory():
    register = _register()
    register["rows"][0]["verdict_receipts"].pop(0)
    with pytest.raises(SystemExit, match="inventory is incomplete"):
        register_tool.validate(register)


def test_ins01_all_fifteen_cell_receipts_are_mandatory():
    register = _register()
    register["rows"][0]["verdict_receipts"].pop()
    with pytest.raises(SystemExit, match="inventory is incomplete"):
        register_tool.validate(register)


def test_ins01_sibling_bytes_are_hash_bound():
    if not (register_tool.ROOT.parent / "oph-physics-sim" / ".git").exists():
        pytest.skip("sibling simulator custody checkout is absent")
    register = _register()
    register["rows"][0]["verdict_receipts"][-1]["sha256"] = "4" * 64
    with pytest.raises(SystemExit, match="external artifact hash mismatch"):
        register_tool.validate(register)


def test_ledger_id_range_includes_n() -> None:
    register = _register()
    register["rows"] = [register["rows"][0]]
    register["rows"][0]["ledger_row"] = "OL-N1"
    register["ledger_controls"] = [
        {
            "ledger_row": "OL-N1",
            "controlling_instrument": "INS-01",
            "supersession_policy": "Fixture supersession policy.",
        }
    ]
    rows = register_tool.validate(
        register,
        ledger_rows={
            "OL-N1": {
                "id": "OL-N1",
                "status": "owed",
                "architecture_version": "AV-0",
            }
        },
    )
    assert rows[0]["ledger_row"] == "OL-N1"


def test_ledger_id_range_rejects_o() -> None:
    register = _register()
    register["rows"][0]["ledger_row"] = "OL-O1"
    with pytest.raises(SystemExit, match="ledger_row must match"):
        register_tool.validate(register)


def test_specified_successor_cannot_control() -> None:
    register = _register()
    register["ledger_controls"][0]["controlling_instrument"] = "INS-02"
    with pytest.raises(SystemExit, match="must carry a decisive completed verdict"):
        register_tool.validate(register)


def test_inconclusive_successor_cannot_erase_controlling_failure() -> None:
    register = _register()
    row = register["rows"][1]
    row["status"] = "INCONCLUSIVE"
    row["freeze_artifacts"] = copy.deepcopy(FIXTURE_FREEZE)
    row["frozen_utc"] = "2026-08-11T00:00:00Z"
    row["verdict_receipts"] = copy.deepcopy(FIXTURE_RECEIPTS)
    register["ledger_controls"][0]["controlling_instrument"] = "INS-02"
    with pytest.raises(SystemExit, match="must carry a decisive completed verdict"):
        register_tool.validate(register)


def test_ineligible_positive_successor_cannot_control() -> None:
    register = _register()
    row = register["rows"][1]
    row["status"] = "REPLICATED"
    row["freeze_artifacts"] = copy.deepcopy(FIXTURE_FREEZE)
    row["frozen_utc"] = "2026-08-11T00:00:00Z"
    row["verdict_receipts"] = copy.deepcopy(FIXTURE_RECEIPTS)
    register["ledger_controls"][0]["controlling_instrument"] = "INS-02"
    with pytest.raises(SystemExit, match="REPLICATED verdict must be ELIGIBLE"):
        register_tool.validate(register)


def test_completed_failed_successor_may_control_but_keeps_row_owed() -> None:
    register = _register()
    row = register["rows"][1]
    row["status"] = "FAILED"
    row["freeze_artifacts"] = copy.deepcopy(FIXTURE_FREEZE)
    row["frozen_utc"] = "2026-08-11T00:00:00Z"
    row["verdict_receipts"] = copy.deepcopy(FIXTURE_RECEIPTS)
    register["ledger_controls"][0]["controlling_instrument"] = "INS-02"
    rows = register_tool.validate(register)
    assert rows[1]["status"] == "FAILED"

    with pytest.raises(SystemExit, match="requires the ledger row to be owed"):
        register_tool.validate(
            register,
            ledger_rows=_ledger_rows(ol_a1_status="attained"),
        )


def test_controlling_failed_verdict_forces_owed_ledger_status() -> None:
    register = _register()
    with pytest.raises(
        SystemExit,
        match="controlling FAILED verdict requires the ledger row to be owed",
    ):
        register_tool.validate(
            register,
            ledger_rows=_ledger_rows(ol_a1_status="attained"),
        )


def test_eligible_audited_replicated_successor_can_control() -> None:
    register = _register()
    _replicated_ins02(register)
    register["ledger_controls"][0]["controlling_instrument"] = "INS-02"
    rows = register_tool.validate(
        register,
        ledger_rows=_ledger_rows(ol_a1_status="attained"),
        promotion_gates={
            "anchored_versions": {"AV-0"},
            "current_version": "AV-0",
            "independent_audit_promotions": {"AUD-FIXTURE": {"OL-A1"}},
            "independently_audited_instruments": {
                "AUD-FIXTURE": {
                    "INS-02": register_tool.auditable_instrument_sha256(
                        register["rows"][1]
                    )
                }
            },
        },
    )
    assert rows[1]["promotion_eligibility"]["state"] == "ELIGIBLE"


def test_eligible_successor_rejects_unanchored_architecture() -> None:
    register = _register()
    _replicated_ins02(register)
    with pytest.raises(SystemExit, match="architecture version is not origin-anchored"):
        register_tool.validate(
            register,
            ledger_rows=_ledger_rows(ol_a1_status="attained"),
            promotion_gates={
                "anchored_versions": set(),
                "current_version": "AV-0",
                "independent_audit_promotions": {"AUD-FIXTURE": {"OL-A1"}},
                "independently_audited_instruments": {
                    "AUD-FIXTURE": {
                        "INS-02": register_tool.auditable_instrument_sha256(
                            register["rows"][1]
                        )
                    }
                },
            },
        )


def test_eligible_successor_rejects_nonqualifying_audit() -> None:
    register = _register()
    _replicated_ins02(register)
    with pytest.raises(SystemExit, match="independent promotion audit"):
        register_tool.validate(
            register,
            ledger_rows=_ledger_rows(ol_a1_status="attained"),
            promotion_gates={
                "anchored_versions": {"AV-0"},
                "current_version": "AV-0",
                "independent_audit_promotions": {"AUD-FIXTURE": set()},
                "independently_audited_instruments": {
                    "AUD-FIXTURE": {
                        "INS-02": register_tool.auditable_instrument_sha256(
                            register["rows"][1]
                        )
                    }
                },
            },
        )


def test_eligible_successor_requires_matching_ledger_architecture() -> None:
    register = _register()
    _replicated_ins02(register)
    with pytest.raises(SystemExit, match="does not match the ledger row"):
        register_tool.validate(
            register,
            ledger_rows=_ledger_rows(
                ol_a1_status="attained", architecture="AV-1"
            ),
            promotion_gates={
                "anchored_versions": {"AV-0"},
                "current_version": "AV-0",
                "independent_audit_promotions": {"AUD-FIXTURE": {"OL-A1"}},
                "independently_audited_instruments": {
                    "AUD-FIXTURE": {
                        "INS-02": register_tool.auditable_instrument_sha256(
                            register["rows"][1]
                        )
                    }
                },
            },
        )


def test_eligible_successor_requires_current_architecture() -> None:
    register = _register()
    _replicated_ins02(register)
    with pytest.raises(SystemExit, match="architecture version is not current"):
        register_tool.validate(
            register,
            ledger_rows=_ledger_rows(ol_a1_status="attained"),
            promotion_gates={
                "anchored_versions": {"AV-0", "AV-1"},
                "current_version": "AV-1",
                "independent_audit_promotions": {"AUD-FIXTURE": {"OL-A1"}},
                "independently_audited_instruments": {
                    "AUD-FIXTURE": {
                        "INS-02": register_tool.auditable_instrument_sha256(
                            register["rows"][1]
                        )
                    }
                },
            },
        )


def test_eligible_successor_requires_audit_of_exact_instrument() -> None:
    register = _register()
    _replicated_ins02(register)
    with pytest.raises(SystemExit, match="does not pin this exact replicated instrument"):
        register_tool.validate(
            register,
            ledger_rows=_ledger_rows(ol_a1_status="attained"),
            promotion_gates={
                "anchored_versions": {"AV-0"},
                "current_version": "AV-0",
                "independent_audit_promotions": {"AUD-FIXTURE": {"OL-A1"}},
                "independently_audited_instruments": {
                    "AUD-FIXTURE": {"INS-02": "0" * 64}
                },
            },
        )


def test_replicated_root_cannot_control_without_promotion_gates() -> None:
    register = _register()
    register["rows"] = [register["rows"][0]]
    register["rows"][0]["status"] = "REPLICATED"
    register["ledger_controls"][0]["controlling_instrument"] = "INS-01"
    with pytest.raises(SystemExit, match="REPLICATED verdict must be ELIGIBLE"):
        register_tool.validate(
            register,
            ledger_rows=_ledger_rows(ol_a1_status="attained"),
        )


def test_nonreplicated_successor_cannot_be_eligible() -> None:
    register = _register()
    _replicated_ins02(register)
    register["rows"][1]["status"] = "FAILED"
    with pytest.raises(SystemExit, match="only a REPLICATED instrument can be ELIGIBLE"):
        register_tool.validate(register)


def test_lineage_predecessor_must_name_earlier_same_ledger_instrument() -> None:
    register = _register()
    register["rows"][1]["lineage_predecessor"] = "INS-99"
    with pytest.raises(SystemExit, match="must name an earlier instrument"):
        register_tool.validate(register)


def test_ledger_control_coverage_is_exact() -> None:
    register = _register()
    register["ledger_controls"] = []
    with pytest.raises(SystemExit, match="must be a nonempty list"):
        register_tool.validate(register)
