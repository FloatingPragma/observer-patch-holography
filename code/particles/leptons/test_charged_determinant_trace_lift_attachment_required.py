#!/usr/bin/env python3
"""Guard the missing charged determinant trace-lift promotion gate."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
NO_GO_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_determinant_trace_normalization_no_go.py"
SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_determinant_trace_lift_attachment_required.py"
OUTPUT = ROOT / "particles" / "runs" / "leptons" / "charged_determinant_trace_lift_attachment_required.json"


def test_charged_determinant_trace_lift_attachment_required() -> None:
    subprocess.run([sys.executable, str(NO_GO_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_charged_determinant_trace_lift_attachment_required"
    assert payload["status"] == "missing_theorem"
    assert payload["source_only"] is False
    assert payload["promotable"] is False
    assert payload["public_theorem_value"] is None
    assert payload["current_closed_chain"]["A_ch_to_charged_masses"] is True
    assert payload["current_closed_chain"]["P_to_A_ch"] is False
    assert payload["required_identity"] == "3*A_ch(P)=sum_psi M_ch[psi]*log(q_psi(P))"
    assert "charged_determinant_trace_lift_attachment" in payload["missing_for_promotion"]
    assert "NO_TARGET_LEAK_DAG_CHARGED_A_CH" in payload["missing_for_promotion"]
    assert "lepton_current_family_exact_readout" in payload["forbidden_ancestors"]
    assert payload["conditional_readout_if_closed"]["centered_logs"]["ell_e"] == -4.495210108808081
    assert payload["conditional_readout_if_closed"]["formulas"]["electron"] == (
        "m_e(P)=exp(A_ch(P)-4.495210108808081)"
    )
