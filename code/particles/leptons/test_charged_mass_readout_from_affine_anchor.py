#!/usr/bin/env python3
"""Guard the algebraic charged-mass readout theorem from theorem-grade A_ch(P)."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ANCHOR_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_absolute_anchor_section.py"
LIFT_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_physical_determinant_line_canonical_uncentered_lift.py"
REDUCTION_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_physical_class_affine_scalar_reduction.py"
SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_mass_readout_from_affine_anchor.py"
OUTPUT = ROOT / "particles" / "runs" / "leptons" / "charged_mass_readout_from_affine_anchor.json"


def test_charged_mass_readout_from_affine_anchor() -> None:
    subprocess.run([sys.executable, str(ANCHOR_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(LIFT_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(REDUCTION_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_charged_mass_readout_from_affine_anchor"
    assert payload["proof_status"] == "closed_algebraic_readout_theorem"
    assert payload["theorem_id"] == "charged_mass_readout_from_affine_anchor"
    assert payload["exact_readout"]["absolute_scale"] == "g_e(P) = exp(A_ch(P))"
    assert payload["exact_readout"]["charged_masses"] == "m_i(P) = exp(A_ch(P) + ell_i_centered(P))"
    assert payload["issue_repair_interpretation"]["separate_open_upstream_math"] == (
        "d10_to_charged_determinant_line_bridge"
    )
