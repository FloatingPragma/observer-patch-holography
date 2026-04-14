#!/usr/bin/env python3
"""Validate the current-corpus D10-to-charged affine bridge no-go."""

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
ANCHOR_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_absolute_anchor_section.py"
ROUTE_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_post_promotion_absolute_closure_route.py"
REDUCTION_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_p_to_affine_anchor_reduction.py"
SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_p_to_affine_anchor_bridge_no_go.py"
OUTPUT = ROOT / "particles" / "runs" / "leptons" / "charged_p_to_affine_anchor_bridge_no_go.json"


def test_current_corpus_has_no_d10_to_charged_affine_bridge() -> None:
    subprocess.run([sys.executable, str(FAMILY_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SOURCE_PAIR_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(FORWARD_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(REPAIR_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(ANCHOR_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(ROUTE_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(REDUCTION_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)

    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_charged_p_to_affine_anchor_bridge_no_go"
    assert payload["verdict"] == "no_current_corpus_P_to_A_ch_bridge"
    assert payload["bridge_absence_certificate"]["no_hidden_landing_on_physical_Y_e_from_P"] is True
    assert payload["bridge_absence_certificate"]["no_hidden_landing_on_A_ch_from_P"] is True
    assert payload["exact_missing_object"]["id"] == "d10_to_charged_affine_anchor_bridge"
    assert payload["exact_missing_object"]["exact_reduction_artifact"] == "oph_charged_p_to_affine_anchor_reduction"
    assert payload["exact_missing_object"]["exact_smallest_bridge_target"] == "d10_to_charged_determinant_line_bridge"
    assert "theorem_grade_physical_Y_e(P)" in payload["exact_missing_object"]["admissible_landings"]
    assert payload["future_symbolic_surface_if_bridge_closes"]["g_e(P)"] == "exp(A_ch(P))"
    assert "m_e(P)" in payload["forbidden_current_outputs"]
