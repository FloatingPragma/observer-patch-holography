"""Drift and tamper rejection tests for the V3 premise register surface."""

from __future__ import annotations

import copy
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

import build_premise_register as register_tool
import check_claim_registry as claim_tool


def _register() -> dict:
    return copy.deepcopy(register_tool.load_json(register_tool.REGISTER_PATH))


def test_retired_premise_id_cannot_be_reused() -> None:
    register = _register()
    register["rows"][-1]["id"] = "PR-83"
    with pytest.raises(SystemExit, match="retired premise ids must not be reused"):
        register_tool.validate(register)
    assert "PR-83" not in {row[0] for row in register_tool.EXPECTED_ROWS}


def test_committed_register_validates() -> None:
    rows = register_tool.validate(_register())
    assert [row["id"] for row in rows] == [
        expected[0] for expected in register_tool.EXPECTED_ROWS
    ]


def test_rebuild_parity_with_committed_surface() -> None:
    register = _register()
    rows = register_tool.validate(register)
    rendered = register_tool.render(rows, register).encode("utf-8")
    committed = register_tool.SURFACE_PATH.read_bytes()
    assert rendered == committed


def test_register_carries_exactly_the_folded_sections() -> None:
    register = _register()
    assert list(register) == list(register_tool.TOP_KEYS)
    assert register_tool.TOP_KEYS == (
        "schema", "issue", "rows", "model_supply", "sm_coverage", "gravity_rungs"
    )
    for section, schema in register_tool.SECTION_SCHEMAS.items():
        assert register[section]["schema"] == schema


@pytest.mark.parametrize("section", ["model_supply", "sm_coverage", "gravity_rungs"])
def test_missing_section_rejected(section: str) -> None:
    register = _register()
    del register[section]
    with pytest.raises(SystemExit, match="top-level keys must be exactly"):
        register_tool.validate(register)


def test_unknown_top_level_section_rejected() -> None:
    register = _register()
    register["evidence_rows"] = []
    with pytest.raises(SystemExit, match="top-level keys must be exactly"):
        register_tool.validate(register)


@pytest.mark.parametrize("section", ["model_supply", "sm_coverage", "gravity_rungs"])
def test_section_schema_tamper_rejected(section: str) -> None:
    register = _register()
    register[section]["schema"] = "oph.some_other_ledger.v1"
    with pytest.raises(SystemExit, match=f"section {section}: schema must equal"):
        register_tool.validate(register)


def test_folded_section_validators_run_inside_register_validation() -> None:
    register = _register()
    register["sm_coverage"]["rows"][0]["classification"] = "derived"
    with pytest.raises(SystemExit, match="classification"):
        register_tool.validate(register)
    register = _register()
    next(r for r in register["gravity_rungs"]["rungs"] if r["rung"] == 7)[
        "status"
    ] = "axiom_forced"
    with pytest.raises(SystemExit, match="not among its registry"):
        register_tool.validate(register)


def test_page_renders_the_folded_sections_with_anchors() -> None:
    register = _register()
    rendered = register_tool.render(register_tool.validate(register), register)
    introduction = rendered.split("| Row |", 1)[0]
    assert "(#sm-coverage)" in introduction
    assert "(#gravity-rungs)" in introduction
    assert "(../../model/PREMISES.md)" in introduction
    sm_start = rendered.index('<a id="sm-coverage"></a>')
    gravity_start = rendered.index('<a id="gravity-rungs"></a>')
    assert rendered.index("## What the dispositions mean") < sm_start < gravity_start
    assert rendered.count("## The Standard Model term-and-sector correspondence table") == 1
    assert rendered.count("## The gravity premise-elimination ladder") == 1
    assert rendered[sm_start:gravity_start].count("| SML-") == len(
        register["sm_coverage"]["rows"]
    )
    assert len(register["gravity_rungs"]["rungs"]) == 11
    for rung in register["gravity_rungs"]["rungs"]:
        assert f"- **Rung {rung['rung']}**:" in rendered[gravity_start:]


def test_intro_keeps_register_scientific_and_explains_evidence_roles() -> None:
    register = _register()
    rendered = register_tool.render(register_tool.validate(register), register)
    assert "Established under the historical issue #727 (closed)" in rendered
    assert "Live owner (historical lanes)" in rendered
    assert "Current implementation" in rendered
    assert "maintained as a scientific register" in rendered
    assert "implicit reverse-consumer edge" in rendered
    assert "premise-discharge queue" not in rendered


def test_intro_accounts_for_valid_result_specific_assumptions() -> None:
    """A valid no-PR consumer still has hypotheses outside the PR inventory."""
    root = Path(__file__).resolve().parents[1]
    registry = claim_tool.load_json(root / "claims" / "claim_registry.yaml")
    claim = next(
        row for row in registry["claims"]
        if row["claim_id"] == "OPH-DE-FIXED-CAPACITY-WLAW"
    )
    assert claim["premise_dependencies"]["classification"] == "explicit_non_consumer"
    claim_tool.check_premise_dependencies(claim, claim_tool.premise_ids(root))
    assert "fixed_or_monotone_capacity_branch" in claim["assumptions"]
    assert set(claim["assumptions"]) <= claim_tool.dictionary_tokens(root)

    register = _register()
    introduction = register_tool.render(
        register_tool.validate(register), register
    ).split("| Row |", 1)[0]
    assert "canonical PR-numbered premises" in introduction
    assert "result-specific assumptions" in introduction.lower()
    assert "claims/claim_registry.yaml" in introduction
    assert "claims/assumption_dictionary.md" in introduction
    assert "Every named premise the V3 composition lanes consume appears here" not in introduction
    assert "A premise used by a lane theorem and absent from this register" not in introduction


def test_check_mode_passes_on_committed_artifacts() -> None:
    result = subprocess.run(
        [
            sys.executable,
            str(Path(__file__).parent / "build_premise_register.py"),
            "--check",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr


def test_disposition_tamper_rejected() -> None:
    register = _register()
    assert register["rows"][2]["id"] == "PR-03"
    register["rows"][2]["disposition"] = "axiomatize"
    with pytest.raises(SystemExit, match="disposition must equal"):
        register_tool.validate(register)


def test_unknown_type_rejected() -> None:
    register = _register()
    register["rows"][0]["type"] = "aesthetic_preference"
    with pytest.raises(SystemExit, match="not in the enum"):
        register_tool.validate(register)


def test_unknown_disposition_rejected() -> None:
    register = _register()
    register["rows"][0]["disposition"] = "postpone"
    with pytest.raises(SystemExit, match="not in the enum"):
        register_tool.validate(register)


def test_missing_evidence_path_rejected() -> None:
    register = _register()
    register["rows"][0]["evidence"].append("docs/NO_SUCH_EVIDENCE_FILE.md")
    with pytest.raises(SystemExit, match="evidence path missing"):
        register_tool.validate(register)


def test_lane_out_of_range_rejected() -> None:
    register = _register()
    register["rows"][0]["consuming_lanes"] = [727]
    with pytest.raises(SystemExit, match="not a registered historical lane"):
        register_tool.validate(register)


def test_retired_lane_738_rejected() -> None:
    register = _register()
    register["rows"][0]["consuming_lanes"] = [738]
    with pytest.raises(SystemExit, match="not a registered historical lane"):
        register_tool.validate(register)


def test_evidence_role_path_parity_is_exact() -> None:
    register = _register()
    row = register["rows"][0]
    row["evidence_roles"].pop(row["evidence"][0])
    with pytest.raises(SystemExit, match="exact evidence-path parity"):
        register_tool.validate(register)


def test_unknown_evidence_role_rejected() -> None:
    register = _register()
    row = register["rows"][0]
    row["evidence_roles"][row["evidence"][0]] = "related"
    with pytest.raises(SystemExit, match="must be one of"):
        register_tool.validate(register)


def test_pr53_has_no_retired_custody_lane() -> None:
    register = _register()
    row = next(row for row in register["rows"] if row["id"] == "PR-53")
    assert row["consuming_lanes"] == [733, 736, 740]
    assert 738 not in row["consuming_lanes"]


def test_renamed_row_rejected() -> None:
    register = _register()
    register["rows"][4]["name"] = "path reference"
    with pytest.raises(SystemExit, match="name must equal"):
        register_tool.validate(register)


def test_dropped_row_rejected() -> None:
    register = _register()
    register["rows"].pop()
    with pytest.raises(SystemExit, match="exactly 84 entries"):
        register_tool.validate(register)


def test_stale_surface_fails_check(tmp_path, monkeypatch) -> None:
    stale = tmp_path / "PREMISE_REGISTER_V3.md"
    stale.write_bytes(b"stale\n")
    monkeypatch.setattr(register_tool, "SURFACE_PATH", stale)
    monkeypatch.setattr(sys, "argv", ["build_premise_register.py", "--check"])
    assert register_tool.main() == 1


def test_duplicate_json_key_rejected(tmp_path: Path) -> None:
    duplicate = tmp_path / "duplicate.json"
    duplicate.write_text('{"schema":"one","schema":"two"}', encoding="utf-8")
    with pytest.raises(SystemExit, match="duplicate JSON key"):
        register_tool.load_json(duplicate)


def test_every_row_carries_a_valid_discharge_record() -> None:
    register = _register()
    rows = register_tool.validate(register)
    statuses = {row["discharge"]["status"] for row in rows}
    assert statuses <= set(register_tool.DISCHARGE_STATUSES)
    by_id = {row["id"]: row["discharge"]["status"] for row in rows}
    assert by_id["PR-01"] == "obstructed"
    assert by_id["PR-65"] == "obstructed"


def test_unknown_discharge_status_rejected() -> None:
    register = _register()
    register["rows"][0]["discharge"]["status"] = "proved_physical"
    with pytest.raises(SystemExit, match="discharge status must be one of"):
        register_tool.validate(register)


@pytest.mark.parametrize("premise", ["PR-21", "PR-22", "PR-81", "PR-82"])
def test_execution_asserting_rows_bind_executable_evidence(premise: str) -> None:
    register = _register()
    row = next(row for row in register["rows"] if row["id"] == premise)
    media = [claim_tool.evidence_medium(path) for path in row["discharge"]["evidence"]]
    assert all(medium is not None for medium in media)
    assert any(medium not in claim_tool.PROSE_EVIDENCE_MEDIA for medium in media)


def test_execution_asserting_discharge_rejects_prose_only_evidence() -> None:
    register = _register()
    row = next(row for row in register["rows"] if row["id"] == "PR-81")
    row["discharge"]["evidence"] = ["model/README.md"]
    with pytest.raises(SystemExit, match="evidence cannot be prose only"):
        register_tool.validate(register)


def test_declared_discharge_may_rest_on_prose_evidence() -> None:
    register = _register()
    row = next(row for row in register["rows"] if row["id"] == "PR-05")
    assert row["discharge"]["status"] == "declared"
    assert row["discharge"]["evidence"] == ["model/CANDIDATES.md"]
    register_tool.validate(register)


def test_execution_asserting_discharge_rejects_unclassified_medium(monkeypatch) -> None:
    register = _register()
    monkeypatch.setattr(register_tool.claim_tool, "evidence_medium", lambda path: None)
    with pytest.raises(SystemExit, match="medium is unclassified"):
        register_tool.validate(register)
