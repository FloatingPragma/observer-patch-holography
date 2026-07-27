"""Fail-closed tests for the #618 gravity premise ladder."""

import importlib.util
import json
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
BUILDER = REPO_ROOT / "tools" / "build_gravity_ladder.py"

spec = importlib.util.spec_from_file_location("build_gravity_ladder", BUILDER)
ladder_tool = importlib.util.module_from_spec(spec)
sys.modules["build_gravity_ladder"] = ladder_tool
spec.loader.exec_module(ladder_tool)


def live_ladder() -> dict:
    return json.loads(ladder_tool.LADDER_PATH.read_text(encoding="utf-8"))


def test_live_ladder_validates_and_surface_is_current():
    ladder = live_ladder()
    rungs = ladder_tool.validate(ladder)
    rendered = ladder_tool.render(ladder, rungs)
    committed = ladder_tool.SURFACE_PATH.read_text(encoding="utf-8")
    assert rendered == committed
    assert len(rungs) == 11


def test_every_rung_has_countermodel_and_mutation_test():
    rungs = ladder_tool.validate(live_ladder())
    for row in rungs:
        assert row["theorem_or_countermodel"]
        assert row["mutation_tests"]


def test_status_stronger_than_registry_class_fails():
    ladder = live_ladder()
    row = next(r for r in ladder["rungs"] if r["rung"] == 7)
    row["status"] = "axiom_forced"
    with pytest.raises(SystemExit, match="not among its registry"):
        ladder_tool.validate(ladder)


def test_unknown_interface_fails():
    ladder = live_ladder()
    ladder["rungs"][3]["premise_interfaces"] = ["no_such_interface"]
    with pytest.raises(SystemExit, match="unknown premise interface"):
        ladder_tool.validate(ladder)


def test_missing_lean_declaration_fails():
    ladder = live_ladder()
    row = next(r for r in ladder["rungs"] if r["rung"] == 9)
    row["theorem_or_countermodel"][0]["declaration"] = "no_such_declaration"
    with pytest.raises(SystemExit, match="absent from"):
        ladder_tool.validate(ladder)


def test_missing_mutation_test_fails():
    ladder = live_ladder()
    row = next(r for r in ladder["rungs"] if r["rung"] == 5)
    row["mutation_tests"] = []
    with pytest.raises(SystemExit, match="mutation test is required"):
        ladder_tool.validate(ladder)


def test_pending_identification_requires_live_owner():
    ladder = live_ladder()
    row = next(r for r in ladder["rungs"] if r["rung"] == 10)
    row["owner_issues"] = []
    with pytest.raises(SystemExit, match="requires owner issues"):
        ladder_tool.validate(ladder)


def test_unknown_assumption_token_fails():
    ladder = live_ladder()
    row = next(r for r in ladder["rungs"] if r["rung"] == 4)
    row["assumption_tokens"] = ["no_such_token"]
    with pytest.raises(SystemExit, match="not in the dictionary"):
        ladder_tool.validate(ladder)
