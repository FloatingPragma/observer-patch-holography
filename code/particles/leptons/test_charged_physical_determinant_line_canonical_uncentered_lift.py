#!/usr/bin/env python3
"""Guard the canonical charged uncentered lift theorem on theorem-grade physical Y_e."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DETERMINANT_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_determinant_line_section_extension.py"
ANCHOR_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_absolute_anchor_section.py"
EQUALIZER_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_physical_identity_mode_equalizer.py"
SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_physical_determinant_line_canonical_uncentered_lift.py"
OUTPUT = ROOT / "particles" / "runs" / "leptons" / "charged_physical_determinant_line_canonical_uncentered_lift.json"


def test_charged_physical_determinant_line_canonical_uncentered_lift() -> None:
    subprocess.run([sys.executable, str(DETERMINANT_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(ANCHOR_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(EQUALIZER_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_charged_physical_determinant_line_canonical_uncentered_lift"
    assert payload["proof_status"] == "closed_conditional_theorem_on_theorem_grade_physical_Y_e"
    assert payload["theorem_id"] == "charged_physical_determinant_line_canonical_uncentered_lift"
    assert payload["canonical_definitions"]["descended_scalar"] == "mu_phys(Y_e) = (1/3) * log(det(Y_e))"
    assert payload["canonical_definitions"]["uncentered_lift"] == (
        "C_tilde_e(Y_e) = C_hat_e(Y_e) + mu_phys(Y_e) I"
    )
    assert payload["forced_equalizer_on_declared_surface"]["fiber_rule"] == (
        "delta(r,r') = 0 for refinement representatives r,r' of the same physical Y_e"
    )
    assert payload["canonical_consequences"]["determinant_line"] == "s_det(Y_e) = 3 * mu_phys(Y_e)"
