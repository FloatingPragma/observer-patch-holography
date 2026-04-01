#!/usr/bin/env python3
"""Validate the current-family quark exactness audit artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SPREAD_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_spread_map.py"
MEAN_SPLIT_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_sector_mean_split.py"
DESCENT_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_sector_descent.py"
FORWARD_SCRIPT = ROOT / "particles" / "flavor" / "build_forward_yukawas.py"
J_B_EVALUATOR_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_diagonal_B_odd_source_scalar_evaluator.py"
D12_SELECTOR_SCRIPT = ROOT / "particles" / "flavor" / "derive_light_quark_isospin_overlap_defect_selector_law.py"
D12_OVERLAP_LAW_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_d12_overlap_transport_law.py"
ONE_SCALAR_SPECIALIZATION_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_d12_one_scalar_specialization.py"
MASS_RAY_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_d12_mass_ray.py"
T1_VALUE_LAW_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_d12_t1_value_law.py"
MASS_SIDE_UNDERDETERMINATION_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_d12_mass_side_underdetermination_theorem.py"
PHYSICAL_BRANCH_REPAIR_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_physical_branch_repair_theorem.py"
AUDIT_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_current_family_exactness_audit.py"
DESCENT_OUTPUT = ROOT / "particles" / "runs" / "flavor" / "quark_sector_descent.json"
OUTPUT = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_exactness_audit.json"


def test_quark_exactness_audit_identifies_current_family_residual_after_spread_closure() -> None:
    subprocess.run([sys.executable, str(SPREAD_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(MEAN_SPLIT_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(DESCENT_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(FORWARD_SCRIPT), "--input", str(DESCENT_OUTPUT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(J_B_EVALUATOR_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(D12_SELECTOR_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(D12_OVERLAP_LAW_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(ONE_SCALAR_SPECIALIZATION_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(MASS_RAY_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(T1_VALUE_LAW_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(MASS_SIDE_UNDERDETERMINATION_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(PHYSICAL_BRANCH_REPAIR_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(AUDIT_SCRIPT)], check=True, cwd=ROOT)

    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_quark_current_family_exactness_audit"
    assert payload["centered_ray_fit"]["residual_norm_u"] < 0.1
    assert payload["centered_ray_fit"]["residual_norm_d"] < 0.2
    fit = payload["mean_split_audit"]["exact_two_scalar_mean_fit"]
    assert fit["g_u_exact_fit"] > fit["g_d_exact_fit"] > 0.0
    assert payload["spread_emitter_audit"]["sigma_source_kind"] == "theorem_grade_mean_surface_readback"
    assert payload["spread_emitter_audit"]["spread_emitter_status"] == "closed"
    quad = payload["quadratic_residual_audit"]
    assert quad["residual_norm_u_after_best_quadratic_fit"] < payload["centered_ray_fit"]["residual_norm_u"]
    assert quad["residual_norm_d_after_best_quadratic_fit"] < payload["centered_ray_fit"]["residual_norm_d"]
    assert quad["residual_norm_u_after_best_quadratic_fit"] > 0.05
    assert quad["residual_norm_d_after_best_quadratic_fit"] > 0.1
    diag = payload["diagonal_gap_shift_audit"]
    assert all(abs(value) < 1.0e-12 for value in diag["residual_u_after_best_diagonal_shift"])
    assert all(abs(value) < 1.0e-12 for value in diag["residual_d_after_best_diagonal_shift"])
    assert payload["smallest_constructive_missing_object"] == "source_readback_u_log_per_side_and_source_readback_d_log_per_side"
    blind = payload["b_mode_amplitude_blindness_audit"]
    assert abs(blind["sum_B_ord"]) < 1.0e-12
    assert abs(blind["dot_B_ord_Q_ord"]) < 1.0e-12
    assert payload["b_mode_odd_projector_evaluator"]["artifact"] == "oph_quark_diagonal_B_odd_source_scalar_evaluator"
    assert payload["d12_isospin_selector_law"]["artifact"] == "oph_light_quark_isospin_overlap_defect_selector_law"
    assert payload["d12_overlap_transport_law"]["artifact"] == "oph_quark_d12_overlap_transport_law"
    assert payload["d12_overlap_transport_law"]["next_single_residual_object"] == "Delta_ud_overlap"
    assert payload["d12_one_scalar_specialization"]["artifact"] == "oph_quark_d12_one_scalar_specialization"
    assert payload["d12_one_scalar_specialization"]["scalar_name"] == "ray_modulus"
    assert payload["d12_one_scalar_specialization"]["sample_scalar_name"] == "t1_sample"
    assert payload["d12_one_scalar_specialization"]["next_single_residual_object"] == "D12_ud_mass_ray"
    assert payload["d12_mass_ray"]["artifact"] == "oph_quark_d12_mass_ray"
    assert payload["d12_mass_ray"]["emitted_object"]["id"] == "D12_ud_mass_ray"
    assert payload["d12_mass_ray"]["next_exact_missing_object"] == "quark_d12_t1_value_law"
    assert payload["d12_mass_side_underdetermination_theorem"]["artifact"] == "oph_quark_d12_mass_side_underdetermination_theorem"
    assert payload["d12_mass_side_underdetermination_theorem"]["next_exact_missing_object"] == "quark_d12_t1_value_law"
    assert payload["d12_physical_branch_repair_theorem"]["artifact"] == "oph_quark_physical_branch_repair_theorem"
    assert payload["d12_physical_branch_repair_theorem"]["minimal_branch_shift_repair_theorem"]["must_emit"] == "quark_relative_sheet_selector"
    assert payload["d12_physical_branch_repair_theorem"]["minimal_branch_shift_repair_theorem"]["selected_value"]["sigma_id"] == "sigma_ref"
    assert payload["broader_honest_frontier"] == "quark_d12_t1_value_law"
    assert payload["predictive_J_B_source_law_status"] == "missing"
    assert payload["diagnostic_fit_promotion_allowed"] is False
    assert payload["diagnostic_only_tau_best_fit"]["tau_u_best_fit"] is not None
