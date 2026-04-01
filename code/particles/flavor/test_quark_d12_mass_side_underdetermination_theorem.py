#!/usr/bin/env python3
"""Validate the quark D12 mass-side underdetermination theorem artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
MASS_BRANCH_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_d12_mass_branch_and_ckm_residual.py"
OVERLAP_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_d12_overlap_transport_law.py"
ONE_SCALAR_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_d12_one_scalar_specialization.py"
MASS_RAY_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_d12_mass_ray.py"
T1_VALUE_LAW_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_d12_t1_value_law.py"
SCALARIZED_BUNDLE_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_scalarized_continuation_bundle.py"
PHYSICAL_SCRIPT = ROOT / "particles" / "flavor" / "derive_generation_bundle_same_label_physical_invariant_bundle.py"
QUADRATIC_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_quadratic_even_transport_scalar.py"
SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_d12_mass_side_underdetermination_theorem.py"
OUTPUT = ROOT / "particles" / "runs" / "flavor" / "quark_d12_mass_side_underdetermination_theorem.json"


def test_quark_d12_mass_side_theorem_starts_after_emitted_ray() -> None:
    subprocess.run([sys.executable, str(MASS_BRANCH_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(OVERLAP_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(ONE_SCALAR_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(MASS_RAY_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(T1_VALUE_LAW_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(QUADRATIC_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(PHYSICAL_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCALARIZED_BUNDLE_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)

    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_quark_d12_mass_side_underdetermination_theorem"
    assert payload["proof_status"] == "same_family_mass_ray_emitted_intrinsic_scale_open"
    assert payload["emitted_same_family_ray"]["artifact"] == "oph_quark_d12_mass_ray"
    assert payload["remaining_object"]["id"] == "quark_d12_t1_value_law"
    assert payload["next_exact_missing_object"] == "quark_d12_t1_value_law"
    assert "D12_ud_mass_ray" not in payload["honest_remaining_value_laws"]
