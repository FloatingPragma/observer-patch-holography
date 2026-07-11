#!/usr/bin/env python3
"""Validate the exact neutrino closure summary artifact."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
WEIGHTED_CYCLE_SCRIPT = ROOT / "particles" / "neutrino" / "derive_neutrino_weighted_cycle_repair.py"
RIGIDITY_SCRIPT = ROOT / "particles" / "neutrino" / "derive_neutrino_bridge_rigidity_theorem.py"
ABSOLUTE_SCRIPT = ROOT / "particles" / "neutrino" / "derive_neutrino_absolute_attachment_theorem.py"
SCRIPT = ROOT / "particles" / "neutrino" / "derive_neutrino_lane_closure_contract.py"
OUTPUT = ROOT / "particles" / "runs" / "neutrino" / "neutrino_lane_closure_contract.json"
CERTIFICATE = ROOT / "particles" / "runs" / "neutrino" / "same_label_scalar_certificate.json"
COCYCLE = ROOT / "particles" / "runs" / "flavor" / "overlap_edge_transport_cocycle.json"
PHASE_SOURCE = ROOT / "particles" / "runs" / "neutrino" / "intrinsic_neutrino_mass_eigenstate_bundle_from_scalar_certificate.json"
ISOTROPIC = ROOT / "particles" / "runs" / "neutrino" / "forward_majorana_matrix.json"
SELECTOR = ROOT / "particles" / "runs" / "neutrino" / "neutrino_transport_load_segment_selector.json"


def test_neutrino_lane_closure_contract_records_emitted_theorem_pair() -> None:
    subprocess.run(
        [
            sys.executable,
            str(WEIGHTED_CYCLE_SCRIPT),
            "--certificate",
            str(CERTIFICATE),
            "--cocycle",
            str(COCYCLE),
            "--phase-source",
            str(PHASE_SOURCE),
            "--isotropic",
            str(ISOTROPIC),
            "--selector",
            str(SELECTOR),
        ],
        check=True,
        cwd=ROOT,
    )
    for script in (RIGIDITY_SCRIPT, ABSOLUTE_SCRIPT, SCRIPT):
        subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT)

    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_neutrino_lane_closure_contract"
    assert payload["proof_status"] == "target_informed_weighted_cycle_candidate_with_compare_only_absolute_attachment"
    assert payload["scale_free_prediction_promotion_allowed"] is False
    assert payload["historical_target_exposure"] is True
    assert payload["current_branch_status"]["no_hidden_discrete_branch"] is False
    rigidity = payload["emitted_bridge_rigidity_theorem"]
    assert rigidity["emitted_formula"] == "sum_gap^2 * prod_qbar * solar_response_over_mstar^-0.5"
    assert rigidity["emitted_value"] is None
    assert rigidity["display_value"] == 0.9994295999075177
    attachment = payload["emitted_absolute_attachment_theorem"]
    assert attachment["status"] == "conditional_absolute_family_blocked_by_compare_only_C_nu"
    assert attachment["public_surface_candidate_allowed"] is False
    assert attachment["B_nu"] == 6.696004159297337
    assert attachment["lambda_nu"] == 1.7237014208357415
    assert payload["public_promotion_allowed"] is False
    assert payload["non_circularity_status"]["missing_source_object"] == (
        "source_emitted_neutrino_C_nu_no_compare_target"
    )
