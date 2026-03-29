#!/usr/bin/env python3
"""Validate the D12 quark mass-branch and CKM residual artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_d12_mass_branch_and_ckm_residual.py"
OUTPUT = ROOT / "particles" / "runs" / "flavor" / "quark_d12_mass_branch_and_ckm_residual.json"


def test_quark_d12_diagnostic_branch_stays_continuation_only() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_quark_d12_mass_branch_and_ckm_residual"
    assert payload["status"] == "diagnostic_only_d12_continuation"
    assert payload["public_promotion_allowed"] is False
    assert payload["candidate_mass_branch_from_t1_over_5"]["Delta_ud_overlap"] > 0.0
    assert payload["ckm_cp_residual_generator"]["off_diagonal_abs"]["12"] > 0.2

