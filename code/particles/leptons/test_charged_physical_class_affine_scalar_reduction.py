#!/usr/bin/env python3
"""Guard the charged physical-class affine-scalar reduction theorem artifact."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BRIDGE_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_p_to_affine_anchor_reduction.py"
LIFT_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_physical_determinant_line_canonical_uncentered_lift.py"
SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_physical_class_affine_scalar_reduction.py"
OUTPUT = ROOT / "particles" / "runs" / "leptons" / "charged_physical_class_affine_scalar_reduction.json"


def test_charged_physical_class_affine_scalar_reduction() -> None:
    subprocess.run([sys.executable, str(BRIDGE_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(LIFT_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_charged_physical_class_affine_scalar_reduction"
    assert payload["proof_status"] == "closed_selected_surface_reduction_theorem"
    assert payload["theorem_id"] == "charged_physical_class_affine_scalar_reduction"
    assert payload["dependencies"]["bridge_reduction"]["theorem_id"] == (
        "charged_P_to_A_ch_reduces_to_determinant_line_bridge"
    )
    assert payload["dependencies"]["canonical_lift"]["theorem_id"] == (
        "charged_physical_determinant_line_canonical_uncentered_lift"
    )
    assert payload["canonical_consequences_on_fill"]["affine_anchor"] == (
        "A_ch(P) = mu_phys(Y_e(P)) = (1/3) * log(det(Y_e(P)))"
    )
    assert "No public selected-class charged source descent" in payload["scope_exclusions"][0]
