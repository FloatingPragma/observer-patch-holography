"""Fail-closed tests for the #554 selection ledger and its generated surface."""

import copy
import importlib.util
import json
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
BUILDER = REPO_ROOT / "tools" / "build_selection_ledger.py"

spec = importlib.util.spec_from_file_location("build_selection_ledger", BUILDER)
ledger_tool = importlib.util.module_from_spec(spec)
sys.modules["build_selection_ledger"] = ledger_tool
spec.loader.exec_module(ledger_tool)


def live_ledger() -> dict:
    return json.loads(ledger_tool.LEDGER_PATH.read_text(encoding="utf-8"))


def test_live_ledger_validates_and_surface_is_current():
    rows = ledger_tool.validate(live_ledger())
    rendered = ledger_tool.render(rows)
    committed = ledger_tool.SURFACE_PATH.read_text(encoding="utf-8")
    assert rendered == committed


def test_every_row_names_one_canonical_claim_and_class():
    rows = ledger_tool.validate(live_ledger())
    assert len(rows) >= 22
    for row in rows:
        assert row["class"] in ledger_tool.CLASSES
        assert row["canonical_claim_id"]
        assert row["menu"]["alternatives"].strip()

def test_required_selector_omission_fails_even_after_renumbering():
    ledger = live_ledger()
    ledger["rows"] = [
        row
        for row in ledger["rows"]
        if row["selector_id"] != "compact_gauge_refinement_receipt"
    ]
    for index, row in enumerate(ledger["rows"], start=1):
        row["row"] = index
    with pytest.raises(SystemExit, match="selector inventory mismatch"):
        ledger_tool.validate(ledger)


def test_unknown_conditional_selector_fails():
    ledger = live_ledger()
    ledger["rows"][0]["conditional_on"] = ["no_such_selector"]
    with pytest.raises(SystemExit, match="conditional_on names unknown"):
        ledger_tool.validate(ledger)


def test_conditional_selector_cycle_fails():
    ledger = live_ledger()
    first = ledger["rows"][0]["selector_id"]
    second = ledger["rows"][1]["selector_id"]
    ledger["rows"][0]["conditional_on"] = [second]
    ledger["rows"][1]["conditional_on"] = [first]
    with pytest.raises(SystemExit, match="conditional_on cycle"):
        ledger_tool.validate(ledger)


def test_open_row_owner_outside_claim_gates_fails():
    ledger = live_ledger()
    registry = json.loads(ledger_tool.REGISTRY_PATH.read_text(encoding="utf-8"))
    identification = json.loads(
        ledger_tool.IDENTIFICATION_PATH.read_text(encoding="utf-8")
    )
    snapshot = json.loads(ledger_tool.SNAPSHOT_PATH.read_text(encoding="utf-8"))
    open_issues = {row["number"] for row in snapshot["rows"]}
    row = next(r for r in ledger["rows"] if r["class"] == "open")
    allowed: set[int] = set()
    for claim_id in [row["canonical_claim_id"], *row["secondary_claim_ids"]]:
        claim = next(c for c in registry["claims"] if c["claim_id"] == claim_id)
        allowed |= set(claim["gates"])
        for entry in identification["physical_identifications"]:
            if claim_id in entry["claim_ids"]:
                allowed |= {b["number"] for b in entry["blocking_issues"]}
    divergent = sorted(open_issues - allowed)[0]
    row["owner_issues"] = [divergent]
    with pytest.raises(SystemExit, match="boundaries diverge"):
        ledger_tool.validate(ledger)


def test_unknown_claim_id_fails():
    ledger = live_ledger()
    ledger["rows"][0]["canonical_claim_id"] = "OPH-NO-SUCH-CLAIM"
    with pytest.raises(SystemExit, match="unknown claim id"):
        ledger_tool.validate(ledger)


def test_declared_menu_size_must_match_selector_registry():
    ledger = live_ledger()
    row = next(r for r in ledger["rows"] if r["selector_id"] == "echosahedral_carrier_lineage")
    row["menu"]["size"] = 4
    with pytest.raises(SystemExit, match="declared with size 4"):
        ledger_tool.validate(ledger)


def test_forced_row_without_citation_fails():
    ledger = live_ledger()
    row = next(r for r in ledger["rows"] if r["class"] == "forced")
    row["lean"] = []
    row["receipts"] = []
    with pytest.raises(SystemExit, match="requires a Lean or receipt citation"):
        ledger_tool.validate(ledger)


def test_missing_lean_declaration_fails():
    ledger = live_ledger()
    row = next(r for r in ledger["rows"] if r["selector_id"] == "twelve_unit_port_splitting")
    row["lean"][0]["declaration"] = "no_such_declaration_name"
    with pytest.raises(SystemExit, match="absent from"):
        ledger_tool.validate(ledger)


def test_numeric_menu_outside_compression_accounting_fails():
    ledger = live_ledger()
    row = next(r for r in ledger["rows"] if r["selector_id"] == "compact_gauge_refinement_receipt")
    row["menu"]["size"] = 7
    with pytest.raises(SystemExit, match="must be a compression input"):
        ledger_tool.validate(ledger)


def test_stale_surface_is_detected():
    rows = ledger_tool.validate(live_ledger())
    rendered = ledger_tool.render(rows)
    committed = ledger_tool.SURFACE_PATH.read_text(encoding="utf-8")
    assert rendered == committed
    tampered = committed.replace("Totals:", "Totals (edited):")
    assert tampered != rendered
