#!/usr/bin/env python3
"""Guard the charged physical-scalar descent reduction beneath the lift slot."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
TRACE_LIFT_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_uncentered_trace_lift_scaffold.py"
DETERMINANT_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_determinant_line_section_extension.py"
ANCHOR_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_absolute_anchor_section.py"
COCYCLE_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_uncentered_trace_lift_cocycle_reduction.py"
EQUALIZER_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_physical_identity_mode_equalizer.py"
SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_mu_physical_descent_reduction.py"
OUTPUT = ROOT / "particles" / "runs" / "leptons" / "charged_mu_physical_descent_reduction.json"


def test_charged_mu_physical_descent_reduction() -> None:
    subprocess.run([sys.executable, str(TRACE_LIFT_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(DETERMINANT_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(ANCHOR_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(COCYCLE_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(EQUALIZER_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_charged_mu_physical_descent_reduction"
    assert payload["status"] == "exact_descended_scalar_reduction"
    assert payload["larger_missing_object"] == "refinement_stable_uncentered_trace_lift"
    assert payload["exact_smaller_missing_object"] == "charged_physical_affine_scalar_mu"
    assert payload["exact_smaller_forcing_object"]["id"] == "charged_physical_identity_mode_equalizer"
    assert payload["forced_vanishing"]["on_same_physical_Y_e"] == (
        "delta(r,r') = 0 for refinement representatives of the same physical Y_e"
    )
    assert payload["equivalent_presentations_on_fill"]["descended_scalar"] == "mu_phys(Y_e)"
    assert payload["reduction_theorem"]["id"] == "charged_refinement_stable_mu_descends_to_physical_scalar"
    assert (
        payload["P_threaded_source_exactness_reformulation"]["id"]
        == "d10_selected_charged_source_exactness_forces_physical_affine_scalar_mu"
    )
    assert "P -> Y_e(P)" in payload["P_threaded_source_exactness_reformulation"]["statement"]
