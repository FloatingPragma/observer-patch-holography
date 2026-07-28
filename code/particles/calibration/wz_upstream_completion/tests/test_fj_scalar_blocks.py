"""Regression tests for the FJ Goldstone, Higgs and mixing blocks."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "producers"))

PAYLOAD = json.loads((ROOT / "outputs" / "fj_direct_scalar_blocks.json").read_text(encoding="utf-8"))
VECTOR = json.loads((ROOT / "outputs" / "fj_direct_vector_blocks.json").read_text(encoding="utf-8"))


def test_custodial_control_passes() -> None:
    control = PAYLOAD["controls"]["custodial_goldstone_poles"]
    assert control["passed"], control["difference"]


def test_all_blocks_present() -> None:
    for name, minimum in (("hh", 30), ("G0G0", 25), ("GpGm", 35),
                          ("WpGm", 20), ("ZG0", 18), ("AG0", 5)):
        assert PAYLOAD["blocks"][name]["diagram_count"] >= minimum, name


def test_tadpole_insertions_on_scalar_blocks() -> None:
    for name in ("hh", "G0G0", "GpGm", "WpGm"):
        kinds = {d["kind"] for d in PAYLOAD["blocks"][name]["diagrams"]}
        assert "tadpole_insertion" in kinds, name


def test_mixed_ghost_pairs_in_ww_block() -> None:
    kinds = {d["kind"] for d in VECTOR["blocks"]["WpWm"]["diagrams"]}
    assert any(k.startswith("ghost_bubble_c") and ("cZ" in k or "cA" in k) for k in kinds), kinds


def test_vector_controls_still_green() -> None:
    for name, control in VECTOR["controls"].items():
        assert control.get("passed", True), name


def test_ag0_mixing_nonzero_for_st_replay() -> None:
    assert sp.simplify(sp.sympify(PAYLOAD["controls"]["ag0_mixing_present_for_st_replay"]["pole"])) != 0
