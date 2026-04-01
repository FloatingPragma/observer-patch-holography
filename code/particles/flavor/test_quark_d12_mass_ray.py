#!/usr/bin/env python3
"""Validate the emitted D12 quark mass-ray artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
OVERLAP_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_d12_overlap_transport_law.py"
ONE_SCALAR_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_d12_one_scalar_specialization.py"
MASS_BRANCH_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_d12_mass_branch_and_ckm_residual.py"
SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_d12_mass_ray.py"
OUTPUT = ROOT / "particles" / "runs" / "flavor" / "quark_d12_mass_ray.json"


def test_quark_d12_mass_ray_is_emitted_and_leaves_only_t1_value_law_open() -> None:
    subprocess.run([sys.executable, str(MASS_BRANCH_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(OVERLAP_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(ONE_SCALAR_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)

    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_quark_d12_mass_ray"
    assert payload["proof_status"] == "same_family_mass_ray_emitted_modulus_open"
    assert payload["emitted_object"]["id"] == "D12_ud_mass_ray"
    assert payload["emitted_object"]["unresolved_coordinate"] == "ray_modulus"
    assert payload["same_family_ray"]["ray_formulas"]["Delta_ud_overlap"] == "ray_modulus / 5"
    assert payload["next_exact_missing_object"] == "quark_d12_t1_value_law"
    assert payload["t1_value_law_contract"]["identifies"] == "ray_modulus = t1"
    assert payload["intrinsic_scale_law_contract"]["unique_intersection_with"] == "D12_ud_mass_ray"
