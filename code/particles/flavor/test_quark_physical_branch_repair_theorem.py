#!/usr/bin/env python3
"""Validate the D12 quark physical-branch repair theorem artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
BRANCH_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_d12_mass_branch_and_ckm_residual.py"
REPAIR_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_physical_branch_repair_theorem.py"
OUTPUT = ROOT / "particles" / "runs" / "flavor" / "quark_physical_branch_repair_theorem.json"


def test_quark_physical_branch_repair_theorem_marks_current_d12_sheet_as_no_go() -> None:
    subprocess.run([sys.executable, str(BRANCH_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(REPAIR_SCRIPT)], check=True, cwd=ROOT)

    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_quark_physical_branch_repair_theorem"
    assert payload["proof_status"] == "current_d12_sheet_is_strict_no_go_for_physical_ckm_shell"
    assert payload["current_d12_sheet"]["quark_relative_sheet_selector"]["sigma_id"] == "sigma_ref"
    assert payload["current_d12_sheet"]["branch_key"] == ["D12", "sigma_ref"]
    assert payload["current_d12_sheet"]["family_branch_key"] == ["D12", None]
    assert payload["current_d12_sheet"]["local_solver_limit"] == "same_label_left_handed_local_orbit_singleton_closed"
    assert payload["current_d12_sheet"]["emitted_same_label_left_reference_sheet"]["canonical_token"] == "D12::same_label_left::reference_sheet"
    assert payload["minimal_branch_shift_repair_theorem"]["must_emit"] == "quark_relative_sheet_selector"
    assert payload["minimal_branch_shift_repair_theorem"]["selector_domain"] == "left_handed_same_label_relative_sheet_orbit_only"
    assert payload["minimal_branch_shift_repair_theorem"]["branch_key_after_repair"] == ["D12", "sigma_ref"]
    assert payload["minimal_branch_shift_repair_theorem"]["selected_value"]["sigma_id"] == "sigma_ref"
    assert payload["minimal_branch_shift_repair_theorem"]["selected_value"]["branch_key"] == ["D12", "sigma_ref"]
    assert payload["relative_sheet_scan"]["status"] == "singleton_closed_by_uniqueness_theorem"
    assert payload["relative_sheet_scan"]["selector_value"]["sigma_id"] == "sigma_ref"
    assert payload["relative_sheet_scan"]["selector_value"]["branch_key"] == ["D12", "sigma_ref"]
    assert payload["relative_sheet_scan"]["available_elements"][0]["sigma_id"] == "sigma_ref"
    assert payload["insufficiency_theorem"]["id"] == "D12_local_selector_closure_still_wrong_branch"
    assert payload["minimal_solver_extension"]["id"] == "sigma_ud_orbit"
    assert "left-handed" in payload["minimal_solver_extension"]["definition"]
    assert payload["minimal_solver_extension"]["currently_emitted_reference_sheet_count"] == 1
    assert payload["minimal_solver_extension"]["status"] == "closed_to_reference_singleton"
    assert payload["minimal_solver_extension"]["next_exact_object_after_orbit_closure"] == "quark_d12_t1_value_law"
    assert payload["physical_shell_mismatch"]["undershoot_factors"]["theta_12"] > 1.0
    assert payload["current_sheet_invariants"]["jarlskog_fraction_of_max_allowed_by_current_angles"] > 0.9
    assert "ray_modulus motion" in payload["theorem_statement"]
