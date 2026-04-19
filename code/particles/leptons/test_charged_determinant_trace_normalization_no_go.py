#!/usr/bin/env python3
"""Guard the charged determinant trace-normalization no-go artifact."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_determinant_trace_normalization_no_go.py"
OUTPUT = ROOT / "particles" / "runs" / "leptons" / "charged_determinant_trace_normalization_no_go.json"


def test_charged_determinant_trace_normalization_no_go() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_charged_determinant_trace_normalization_no_go"
    assert payload["proof_status"] == "source_character_defined_additive_normalization_open"
    assert payload["unresolved_scalar_defect"]["id"] == "charged_determinant_normalization_defect"
    assert "mu -> mu + kappa" in payload["affine_shift_witness"]["shift_family"]
    assert "3 mu(r) = sum_e M_e^ch log q_e(r)" in payload["exact_remaining_theorem"]["statement"]
