#!/usr/bin/env python3
"""Validate the decisive charged end-to-end non-closure theorem artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
TRACE_LIFT_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_uncentered_trace_lift_scaffold.py"
DETERMINANT_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_determinant_line_section_extension.py"
ANCHOR_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_absolute_anchor_section.py"
COCYCLE_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_uncentered_trace_lift_cocycle_reduction.py"
EQUALIZER_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_physical_identity_mode_equalizer.py"
DESCENT_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_mu_physical_descent_reduction.py"
NO_GO_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_centered_operator_mu_phys_no_go.py"
ROUTE_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_post_promotion_absolute_closure_route.py"
SPLITTING_OBSTRUCTION_SCRIPT = (
    ROOT / "particles" / "flavor" / "derive_generation_bundle_branch_generator_splitting_obstruction.py"
)
SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_end_to_end_impossibility_theorem.py"
OUTPUT = ROOT / "particles" / "runs" / "leptons" / "charged_end_to_end_impossibility_theorem.json"


def test_charged_lane_is_not_end_to_end_closed_on_live_corpus() -> None:
    subprocess.run([sys.executable, str(TRACE_LIFT_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(DETERMINANT_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(ANCHOR_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(COCYCLE_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(EQUALIZER_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(DESCENT_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(NO_GO_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(ROUTE_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SPLITTING_OBSTRUCTION_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_charged_end_to_end_impossibility_theorem"
    assert payload["verdict"] == "no_current_corpus_end_to_end_closure"
    assert payload["closure_now"] is False
    assert payload["charged_public_masses_emitted"] is False
    assert payload["operator_side_no_go"]["theorem_grade_C_hat_e_available_now"] is False
    assert payload["operator_side_no_go"]["exact_missing_theorem"] == "oph_generation_bundle_branch_generator_splitting"
    assert payload["operator_side_no_go"]["exact_missing_clause"] == "compression_descendant_commutator_vanishes_or_is_uniformly_quadratic_small_after_central_split"
    obstruction = payload["operator_side_no_go"]["exact_obstruction_certificate"]["current_attached_data_obstruction"]
    assert abs(float(obstruction["commutator_operator_norm"]) - 0.04861550547372144) < 1.0e-12
    assert payload["minimal_operator_extension"]["id"] == "central_split_quadratic_commutator_transfer"
    assert payload["minimal_operator_extension"]["current_corpus_contains_theorem"] is False
    assert payload["promotion_only_absolute_no_go"]["theorem_id"] == "charged_centered_operator_cannot_emit_mu_phys"
    assert payload["absolute_anchor_no_go"]["exact_missing_object"] == "charged_absolute_anchor_A_ch"
    assert payload["absolute_anchor_no_go"]["theorem_grade_A_ch_available_now"] is False
    assert payload["absolute_anchor_no_go"]["no_replacement_law_on_current_surface"] is True
    assert payload["exact_irreducible_chain"][1]["id"] == "refinement_stable_uncentered_trace_lift"
    assert payload["exact_irreducible_chain"][1]["internal_carrier"] == "scalar_affine_cocycle_primitive"
    assert payload["exact_irreducible_chain"][1]["exact_descended_scalar"] == "charged_physical_affine_scalar_mu"
    assert payload["exact_irreducible_chain"][1]["exact_smaller_forcing_object"] == "charged_physical_identity_mode_equalizer"
    assert payload["induced_after_irreducible_chain"]["charged_absolute_anchor"]["id"] == "charged_absolute_anchor_A_ch"
    assert payload["future_symbolic_forward_surface"]["if_A_ch_exists"]["g_e"] == "exp(A_ch)"
    assert "m_e" in payload["theorem_forbid_emit_now"]
