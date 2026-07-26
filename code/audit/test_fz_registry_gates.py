"""Fail-closed tests for the #607 frozen-prediction ladder register."""

import importlib.util
import json
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
BUILDER = REPO_ROOT / "tools" / "build_fz_registry.py"

spec = importlib.util.spec_from_file_location("build_fz_registry", BUILDER)
fz_tool = importlib.util.module_from_spec(spec)
sys.modules["build_fz_registry"] = fz_tool
spec.loader.exec_module(fz_tool)


def live_register() -> dict:
    return json.loads(fz_tool.REGISTER_PATH.read_text(encoding="utf-8"))


def test_live_register_validates_and_surface_is_current():
    register = live_register()
    rows = fz_tool.validate(register)
    rendered = fz_tool.render(register, rows)
    committed = fz_tool.SURFACE_PATH.read_text(encoding="utf-8")
    assert rendered == committed


def test_ladder_is_contiguous_fz01_through_fz09():
    rows = fz_tool.validate(live_register())
    assert [row["id"] for row in rows] == [f"FZ-{i:02d}" for i in range(1, 10)]


def test_fz02_hash_is_bound_to_the_live_receipt():
    register = live_register()
    register["rows"][1]["content_sha256"] = "0" * 64
    with pytest.raises(SystemExit, match="does not equal the live"):
        fz_tool.validate(register)


def test_frozen_row_requires_attestation_fields():
    register = live_register()
    register["rows"][0]["attestation"] = None
    with pytest.raises(SystemExit, match="requires attestation"):
        fz_tool.validate(register)


def test_pending_row_requires_a_live_owning_issue():
    register = live_register()
    register["rows"][2]["owning_issue"] = 599
    with pytest.raises(SystemExit, match="not open in the snapshot"):
        fz_tool.validate(register)


def test_pending_row_requires_a_kill_band():
    register = live_register()
    register["rows"][3]["kill_band"] = ""
    with pytest.raises(SystemExit, match="kill_band must be nonempty"):
        fz_tool.validate(register)
