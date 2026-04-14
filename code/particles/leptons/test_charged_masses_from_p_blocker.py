#!/usr/bin/env python3
"""Validate the current-corpus blocker for charged masses from P."""

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
DESCENT_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_mu_physical_descent_reduction.py"
NO_GO_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_centered_operator_mu_phys_no_go.py"
ROUTE_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_post_promotion_absolute_closure_route.py"
END_TO_END_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_end_to_end_impossibility_theorem.py"
REDUCTION_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_p_to_affine_anchor_reduction.py"
BRIDGE_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_p_to_affine_anchor_bridge_no_go.py"
SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_masses_from_p_blocker.py"
OUTPUT = ROOT / "particles" / "runs" / "leptons" / "charged_masses_from_p_blocker.json"


def test_current_corpus_does_not_emit_charged_masses_from_p() -> None:
    subprocess.run([sys.executable, str(TRACE_LIFT_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(DETERMINANT_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(ANCHOR_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(COCYCLE_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(EQUALIZER_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(DESCENT_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(NO_GO_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(ROUTE_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(END_TO_END_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(REDUCTION_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(BRIDGE_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)

    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_charged_masses_from_p_blocker"
    assert payload["verdict"] == "charged_masses_from_P_not_emitted_on_current_corpus"
    assert payload["failed_dependencies"][0]["id"] == "oph_generation_bundle_branch_generator_splitting"
    assert payload["failed_dependencies"][1]["id"] == "refinement_stable_uncentered_trace_lift"
    assert payload["failed_dependencies"][2]["id"] == "d10_to_charged_affine_anchor_bridge"
    assert payload["failed_dependencies"][2]["exact_smallest_bridge_target"] == "d10_to_charged_determinant_line_bridge"
    assert (
        payload["failed_dependencies"][2]["P_threaded_reformulation"]
        == "d10_selected_charged_source_exactness_forces_physical_affine_scalar_mu"
    )
    assert payload["bridge_reduction_artifact"]["artifact"] == "oph_charged_p_to_affine_anchor_reduction"
    assert (
        payload["bridge_reduction_artifact"]["P_threaded_reformulation_id"]
        == "d10_selected_charged_source_exactness_forces_physical_affine_scalar_mu"
    )
    assert payload["symbolic_readout_if_dependencies_close"]["g_e(P)"] == "exp(A_ch(P))"
    assert "m_tau(P)" in payload["current_forbidden_outputs"]
