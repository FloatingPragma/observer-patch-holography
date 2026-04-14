#!/usr/bin/env python3
"""Validate the charged P-to-affine-anchor reduction theorem."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
FAMILY_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_observable_family.py"
SOURCE_PAIR_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_source_transport_pair.py"
FORWARD_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_forward_transmutation_certificate.py"
REPAIR_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_target_free_repair_value_law.py"
DETERMINANT_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_determinant_line_section_extension.py"
ANCHOR_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_absolute_anchor_section.py"
SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_p_to_affine_anchor_reduction.py"
OUTPUT = ROOT / "particles" / "runs" / "leptons" / "charged_p_to_affine_anchor_reduction.json"


def test_charged_p_to_affine_anchor_reduction() -> None:
    subprocess.run([sys.executable, str(FAMILY_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SOURCE_PAIR_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(FORWARD_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(REPAIR_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(DETERMINANT_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(ANCHOR_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)

    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_charged_p_to_affine_anchor_reduction"
    assert payload["proof_status"] == "closed_bridge_reduction_to_determinant_line_or_physical_Y_e"
    assert payload["reduction_theorem"]["id"] == "charged_P_to_A_ch_reduces_to_determinant_line_bridge"
    assert payload["exact_smallest_bridge_target"]["id"] == "d10_to_charged_determinant_line_bridge"
    assert "theorem_grade_physical_Y_e(P)" in payload["exact_smallest_bridge_target"]["equivalent_landings"]
    assert (
        payload["P_threaded_source_exactness_reformulation"]["id"]
        == "d10_selected_charged_source_exactness_forces_physical_affine_scalar_mu"
    )
    assert payload["induced_objects_once_bridge_closes"]["affine_anchor"]["formula"] == (
        "A_ch(P) = (1/3) * s_det(P) = (1/3) * log(det(Y_e(P)))"
    )
