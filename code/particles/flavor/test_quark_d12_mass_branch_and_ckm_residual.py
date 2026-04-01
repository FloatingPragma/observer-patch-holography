#!/usr/bin/env python3
"""Validate the D12 quark mass-branch and CKM closure artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_d12_mass_branch_and_ckm_residual.py"
OUTPUT = ROOT / "particles" / "runs" / "flavor" / "quark_d12_mass_branch_and_ckm_residual.json"


def test_quark_d12_ckm_transport_closes_on_d12_continuation_branch() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_quark_d12_mass_branch_and_ckm_closure"
    assert payload["status"] == "d12_current_sheet_ckm_cp_transport_closed_wrong_branch_no_go"
    assert payload["public_promotion_allowed"] is False
    assert payload["quark_relative_sheet_selector"]["sigma_id"] == "sigma_ref"
    assert payload["physical_branch_status"] == "current_d12_sheet_is_strict_no_go_for_physical_ckm_shell"
    assert payload["sample_same_family_point"]["Delta_ud_overlap"] > 0.0
    assert payload["remaining_open_objects"][0] == "quark_d12_t1_value_law"
    assert "D12_ud_mass_ray" not in payload["remaining_open_objects"]
    assert payload["comparison_only_best_same_family_point"]["status"] == "comparison_only_not_promotable"
    assert payload["forward_same_label_transport"]["principal_log_exists_uniquely"] is True
    assert payload["same_label_transport_generator"]["generator_invariants"]["theta_12_K"] > 0.0
    assert payload["physical_ckm_comparison_shell"]["undershoot_factors"]["theta_12"] > 1.0
    assert payload["closure_residual"]["fro_norm"] < 1.0e-12
