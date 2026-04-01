#!/usr/bin/env python3
"""Validate the exact same-family selected-sheet quark closure artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys

import pytest


ROOT = pathlib.Path(__file__).resolve().parents[2]
SPREAD_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_spread_map.py"
MEAN_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_sector_mean_split.py"
QUADRATIC_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_current_family_quadratic_readout_theorem.py"
EXACT_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_current_family_exact_readout.py"
SELECTOR_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_relative_sheet_selector.py"
PHYSICAL_BRANCH_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_physical_branch_repair_theorem.py"
SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_current_family_selected_sheet_closure.py"
OUTPUT = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_selected_sheet_closure.json"


def test_quark_current_family_selected_sheet_closure() -> None:
    subprocess.run([sys.executable, str(SPREAD_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(MEAN_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(QUADRATIC_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(EXACT_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SELECTOR_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(PHYSICAL_BRANCH_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_quark_current_family_selected_sheet_exact_closure"
    assert payload["proof_status"] == "closed_current_family_selected_sheet_exact_readout"
    assert payload["theorem_scope"] == "current_family_only"
    assert payload["selected_sheet"]["sigma_id"] == "sigma_ref"
    assert set(payload["exact_outputs_gev"]) == {"u", "c", "t", "d", "s", "b"}
    assert payload["exact_fit_residuals_u"] == pytest.approx([0.0, 0.0, 0.0], rel=1.0e-12, abs=1.0e-12)
    assert payload["exact_fit_residuals_d"] == pytest.approx([0.0, 0.0, 0.0], rel=1.0e-12, abs=1.0e-12)
