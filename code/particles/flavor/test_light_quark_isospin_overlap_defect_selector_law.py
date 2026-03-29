#!/usr/bin/env python3
"""Validate the D12 light-quark isospin selector-law artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SOURCE_READBACK_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_diagonal_common_gap_shift_source_readback.py"
SELECTOR_SCRIPT = ROOT / "particles" / "flavor" / "derive_light_quark_isospin_overlap_defect_selector_law.py"
OUTPUT = ROOT / "particles" / "runs" / "flavor" / "light_quark_isospin_overlap_defect_selector_law.json"


def test_light_quark_isospin_selector_law_is_explicit_but_value_open() -> None:
    subprocess.run([sys.executable, str(SOURCE_READBACK_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SELECTOR_SCRIPT)], check=True, cwd=ROOT)

    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_light_quark_isospin_overlap_defect_selector_law"
    assert payload["proof_status"] == "selector_law_closed_value_open"
    assert payload["scope"] == "D12_continuation_only"
    assert payload["recovered_core_promotion_allowed"] is False
    assert payload["selector_scalar_name"] == "Delta_ud_overlap"
    assert payload["next_single_residual_object"] == "Delta_ud_overlap"
    assert payload["selector_equivalence_formula"] == "Delta_ud_overlap = beta_u_diag_B_source - beta_d_diag_B_source"
    assert payload["odd_budget_neutrality_formula"] == "beta_u_diag_B_source + beta_d_diag_B_source = 0"
    assert payload["Delta_ud_overlap"] is None

