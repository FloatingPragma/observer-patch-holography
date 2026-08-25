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
    row_id: str = "INS-04",
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
        "limitations": "Fixture-only validation row with no scientific claim.",
        "verdict_receipts": copy.deepcopy(FIXTURE_RECEIPTS) if receipts else [],
    }


def _replicated_ins02(register: dict) -> dict:
    row = register["rows"][1]
    row["status"] = "REPLICATED"
    row["freeze_artifacts"] = copy.deepcopy(FIXTURE_FREEZE)
    row["frozen_utc"] = "2026-08-11T00:00:00Z"
    row["verdict_receipts"] = copy.deepcopy(FIXTURE_RECEIPTS)
    return row


def _completed_ins03(register: dict, status: str) -> dict:
    row = register["rows"][2]
    row["status"] = status
    row["freeze_artifacts"] = copy.deepcopy(FIXTURE_FREEZE)
    row["frozen_utc"] = "2026-08-11T00:00:00Z"
    row["verdict_receipts"] = copy.deepcopy(FIXTURE_RECEIPTS)
    return row


def _ledger_rows(
    *, ol_a1_status: str = "owed", ol_c5_status: str = "partial"
) -> dict:
    rows = copy.deepcopy(register_tool.load_ledger_rows())
    rows["OL-A1"]["status"] = ol_a1_status
    rows["OL-C5"]["status"] = ol_c5_status
    return rows


def test_committed_register_validates() -> None:
    register = _register()
    rows = register_tool.validate(register)
    assert register["schema"] == "oph.emergent_instrument_register.v5"
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
    assert rows[2]["id"] == "INS-03"
    assert rows[2]["ledger_row"] == "OL-C5"
    assert rows[2]["status"] == "SPECIFIED"
    assert rows[2]["custody_repository"] is None
    assert rows[2]["freeze_artifacts"] == []
    assert rows[2]["frozen_utc"] is None
    assert rows[2]["lineage_predecessor"] is None
    assert rows[2]["verdict_receipts"] == []
    assert register["ledger_controls"][0]["controlling_instrument"] == "INS-01"
    assert register["ledger_controls"][0]["replicated_consequence"] == "attain_row"
    assert register["ledger_controls"][0]["failed_consequence"] == "owe_row"
    assert register["ledger_controls"][0]["replicated_preserves_open_premises"] == []
    assert register["ledger_controls"][0]["failed_preserves_open_premises"] == []
    assert register["ledger_controls"][1] == {
        "ledger_row": "OL-C5",
        "controlling_instrument": None,
        "replicated_consequence": "support_only",
        "failed_consequence": "block_attainment",
        "replicated_preserves_open_premises": ["PR-03"],
        "failed_preserves_open_premises": ["PR-03"],
        "supersession_policy": (
            "No completed decisive instrument controls OL-C5. INS-03 is SPECIFIED"
            " design work only and cannot move the ledger row. A later decisive"
            " verdict must be selected explicitly here in the same scientific"
            " update that applies its ledger consequence."
        ),
    }


def test_ins02_design_contract_is_explicit_and_unfrozen() -> None:
    register = _register()
    row = register["rows"][1]

    assert row["status"] == "SPECIFIED"
    assert row["freeze_artifacts"] == []
    assert row["frozen_utc"] is None
    assert row["verdict_receipts"] == []
    assert row["spec_pointer"] == (
        "plan/OL_A1_FACTORIAL_FOLLOWUP_DESIGN.md (oph-meta planning workspace)"
    )
    assert "16,384; 65,536; 131,072; and 262,144" in row["decision_rule"]
    assert "absolute support sizes 48; 96; and 192" in row["decision_rule"]
    assert "observer count fixed at the future preregistration" in row["decision_rule"]
    assert "target-blind freeze" in row["decision_rule"]
    assert "one fixed estimator" in row["decision_rule"]
    assert "prospective precision and power" in row["seeds_policy"]
    assert "paired blocks" in row["seeds_policy"]
    controls = " ".join(row["controls"])
    assert "Degree-preserving ancestry rewiring" in controls
    assert "sham identifier relabeling" in controls
    assert "synthetic ancestry-sensitive statistic" in controls
    assert "randomized within each block" in row["seeds_policy"]
    assert "quantitative ancestry-destruction" in row["decision_rule"]
    assert "through one fixed estimator" in row["decision_rule"]


def test_ins03_design_contract_is_explicit_unfrozen_and_nonpromoting() -> None:
    register = _register()
    row = register["rows"][2]

    assert row["id"] == "INS-03"
    assert row["ledger_row"] == "OL-C5"
    assert row["status"] == "SPECIFIED"
    assert row["spec_pointer"] == (
        "plan/INS03_SOURCE_BOUND_PHASE_INSTRUMENT_DESIGN.md"
        " (oph-meta planning workspace)"
    )
    assert row["custody_repository"] is None
    assert row["freeze_artifacts"] == []
    assert row["frozen_utc"] is None
    assert row["lineage_predecessor"] is None
    assert row["verdict_receipts"] == []

    rule = row["decision_rule"]
    assert "implemented v1 validator is static-only" in rule
    assert "STATIC_COMMITTED_FIXTURE_CONFORMANT" in rule
    assert "PRODUCER_AUTHENTICATION_UNIMPLEMENTED" in rule
    assert "Neither result maps to REPLICATED or FAILED" in rule
    assert "future authenticated v2" in rule
    assert "No numeric placeholder is bound here" in rule
    assert "authorizes no execution" in rule

    seeds = row["seeds_policy"]
    assert "No seed has been drawn" in seeds
    assert "prospective precision calculation" in seeds
    assert "One fresh master seed" in seeds
    assert "no redraw, no substitution, no optional stopping" in seeds

    controls = " ".join(row["controls"])
    assert "positive-semidefinite certificates" in controls
    assert "matrix-unit trace-action table as a diagnostic" in controls
    assert "no promotion authority" in controls
    assert "PRODUCER_AUTHENTICATION_UNIMPLEMENTED" in controls

    limitations = row["limitations"]
    assert "not a preregistration or an authenticated simulator instrument" in limitations
    assert "OL-C5 remains partial" in limitations
    assert "PR-03" in limitations
    assert "public-readback remainder of PR-64" in limitations
    assert "PR-65 remain open" in limitations
    assert "no endpoint in this design supplies" in limitations
    assert "required to remove PR-03" in limitations
    assert "REPLICATED support_only and FAILED block_attainment" in limitations
    assert "neither verdict alone attains or owes the whole ledger row" in limitations
    assert "OL-C5 has no controlling instrument" in limitations
    assert "separate immutable authenticated-v2 preregistration and explicit authorization" in limitations


def test_rebuild_parity_with_committed_surface() -> None:
    register = _register()
    rows = register_tool.validate(register)
    rendered = register_tool.render(register, rows).encode("utf-8")
    committed = register_tool.SURFACE_PATH.read_bytes()
    assert rendered == committed


def test_intro_explains_explicit_ledger_control() -> None:
    register = _register()
    rendered = register_tool.render(register, register_tool.validate(register))
    assert "completed decisive instrument" in rendered
    assert "ledger-control lineage" in rendered
    assert "OL-C5: no completed controlling verdict" in rendered
    assert "REPLICATED consequence `support_only`" in rendered
    assert "FAILED consequence `block_attainment`" in rendered
    assert rendered.count("preserving open premises PR-03") == 2
    assert "This does not create a controlling verdict" in rendered
    assert "none (unfrozen; no artifacts pinned)" in rendered


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
    with pytest.raises(SystemExit, match="owning_issue must equal 737"):
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
            "replicated_consequence": "attain_row",
            "failed_consequence": "owe_row",
            "replicated_preserves_open_premises": [],
            "failed_preserves_open_premises": [],
            "supersession_policy": "Fixture supersession policy.",
        }
    ]
    rows = register_tool.validate(
        register,
        ledger_rows={
            "OL-N1": {
                "id": "OL-N1",
                "status": "owed",
                "open_premises": [],
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


def test_specified_ins03_cannot_control_ol_c5() -> None:
    register = _register()
    register["ledger_controls"][1]["controlling_instrument"] = "INS-03"
    with pytest.raises(SystemExit, match="must carry a decisive completed verdict"):
        register_tool.validate(register)


def test_null_controller_accepts_inconclusive_completed_instrument() -> None:
    register = _register()
    row = register["rows"][2]
    row["status"] = "INCONCLUSIVE"
    row["freeze_artifacts"] = copy.deepcopy(FIXTURE_FREEZE)
    row["frozen_utc"] = "2026-08-11T00:00:00Z"
    row["verdict_receipts"] = copy.deepcopy(FIXTURE_RECEIPTS)
    rows = register_tool.validate(register)
    assert rows[2]["status"] == "INCONCLUSIVE"
    assert register["ledger_controls"][1]["controlling_instrument"] is None


def test_null_controller_rejects_unselected_decisive_instrument() -> None:
    register = _register()
    _completed_ins03(register, "FAILED")
    with pytest.raises(
        SystemExit,
        match="null controller cannot leave a decisive completed instrument unselected",
    ):
        register_tool.validate(register)


def test_null_controller_rejects_attained_ledger_row() -> None:
    register = _register()
    with pytest.raises(
        SystemExit,
        match="null controller cannot govern an attained ledger row",
    ):
        register_tool.validate(
            register,
            ledger_rows=_ledger_rows(ol_c5_status="attained"),
        )


def test_replicated_ins03_can_supply_partial_support_without_attainment() -> None:
    register = _register()
    _completed_ins03(register, "REPLICATED")
    register["ledger_controls"][1]["controlling_instrument"] = "INS-03"
    rows = register_tool.validate(register)
    assert rows[2]["status"] == "REPLICATED"
    assert register["ledger_controls"][1]["replicated_consequence"] == "support_only"


def test_replicated_support_only_instrument_cannot_leave_row_owed() -> None:
    register = _register()
    _completed_ins03(register, "REPLICATED")
    register["ledger_controls"][1]["controlling_instrument"] = "INS-03"
    with pytest.raises(SystemExit, match="support_only consequence requires"):
        register_tool.validate(
            register,
            ledger_rows=_ledger_rows(ol_c5_status="owed"),
        )


def test_replicated_support_only_instrument_cannot_attain_whole_row() -> None:
    register = _register()
    _completed_ins03(register, "REPLICATED")
    register["ledger_controls"][1]["controlling_instrument"] = "INS-03"
    with pytest.raises(SystemExit, match="support_only consequence requires"):
        register_tool.validate(
            register,
            ledger_rows=_ledger_rows(ol_c5_status="attained"),
        )


def test_failed_ins03_blocks_attainment_without_erasing_partial_evidence() -> None:
    register = _register()
    _completed_ins03(register, "FAILED")
    register["ledger_controls"][1]["controlling_instrument"] = "INS-03"
    rows = register_tool.validate(register)
    assert rows[2]["status"] == "FAILED"
    assert register["ledger_controls"][1]["failed_consequence"] == "block_attainment"


def test_failed_block_attainment_instrument_rejects_attained_row() -> None:
    register = _register()
    _completed_ins03(register, "FAILED")
    register["ledger_controls"][1]["controlling_instrument"] = "INS-03"
    with pytest.raises(SystemExit, match="block_attainment consequence requires"):
        register_tool.validate(
            register,
            ledger_rows=_ledger_rows(ol_c5_status="attained"),
        )


def test_failed_block_attainment_instrument_cannot_erase_partial_evidence() -> None:
    register = _register()
    _completed_ins03(register, "FAILED")
    register["ledger_controls"][1]["controlling_instrument"] = "INS-03"
    with pytest.raises(SystemExit, match="block_attainment consequence requires"):
        register_tool.validate(
            register,
            ledger_rows=_ledger_rows(ol_c5_status="owed"),
        )


def test_ins03_consequences_preserve_pr03_open() -> None:
    register = _register()
    _completed_ins03(register, "REPLICATED")
    register["ledger_controls"][1]["controlling_instrument"] = "INS-03"
    ledger_rows = _ledger_rows()
    ledger_rows["OL-C5"]["open_premises"].remove("PR-03")
    with pytest.raises(SystemExit, match="requires open premise PR-03 on OL-C5"):
        register_tool.validate(register, ledger_rows=ledger_rows)


def test_ledger_consequence_enums_fail_closed() -> None:
    register = _register()
    register["ledger_controls"][1]["replicated_consequence"] = "attain_or_support"
    with pytest.raises(SystemExit, match="replicated_consequence must be one of"):
        register_tool.validate(register)


def test_preserved_open_premise_contract_fails_closed() -> None:
    register = _register()
    register["ledger_controls"][1]["replicated_preserves_open_premises"] = "PR-03"
    with pytest.raises(SystemExit, match="must be a list"):
        register_tool.validate(register)

    register = _register()
    register["ledger_controls"][1]["failed_preserves_open_premises"] = [
        "PR-03",
        "PR-03",
    ]
    with pytest.raises(SystemExit, match="must be duplicate-free"):
        register_tool.validate(register)

    register = _register()
    register["ledger_controls"][1]["replicated_preserves_open_premises"] = [
        "premise-03"
    ]
    with pytest.raises(SystemExit, match="entries must match PR-<two digits>"):
        register_tool.validate(register)

    register = _register()
    register["ledger_controls"][1]["failed_consequence"] = "maybe_partial"
    with pytest.raises(SystemExit, match="failed_consequence must be one of"):
        register_tool.validate(register)

    register = _register()
    register["ledger_controls"][1]["replicated_consequence"] = []
    with pytest.raises(SystemExit, match="replicated_consequence must be one of"):
        register_tool.validate(register)


def test_controller_must_be_null_or_an_instrument_id() -> None:
    register = _register()
    register["ledger_controls"][1]["controlling_instrument"] = 3
    with pytest.raises(SystemExit, match="must be null or an instrument id"):
        register_tool.validate(register)


def test_null_controller_ledger_row_must_be_registered() -> None:
    register = _register()
    register["ledger_controls"][1]["ledger_row"] = "OL-Z9"
    with pytest.raises(SystemExit, match="ledger_row must match"):
        register_tool.validate(register)


def test_each_ledger_row_has_exactly_one_lineage_root() -> None:
    register = _register()
    register["rows"][2]["ledger_row"] = "OL-A1"
    with pytest.raises(SystemExit, match="already has a root instrument"):
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

    with pytest.raises(SystemExit, match="requires ledger status owed"):
        register_tool.validate(
            register,
            ledger_rows=_ledger_rows(ol_a1_status="attained"),
        )


def test_controlling_failed_verdict_forces_owed_ledger_status() -> None:
    register = _register()
    with pytest.raises(
        SystemExit,
        match="controlling FAILED verdict with owe_row consequence",
    ):
        register_tool.validate(
            register,
            ledger_rows=_ledger_rows(ol_a1_status="attained"),
        )


def test_replicated_successor_can_control_attained_ledger_row() -> None:
    register = _register()
    _replicated_ins02(register)
    register["ledger_controls"][0]["controlling_instrument"] = "INS-02"
    rows = register_tool.validate(
        register,
        ledger_rows=_ledger_rows(ol_a1_status="attained"),
    )
    assert rows[1]["status"] == "REPLICATED"


def test_replicated_attain_row_successor_rejects_unattained_ledger_row() -> None:
    register = _register()
    _replicated_ins02(register)
    register["ledger_controls"][0]["controlling_instrument"] = "INS-02"
    with pytest.raises(SystemExit, match="attain_row consequence requires"):
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

    register = _register()
    register["ledger_controls"] = [register["ledger_controls"][0]]
    with pytest.raises(SystemExit, match="cover exactly"):
        register_tool.validate(register)
