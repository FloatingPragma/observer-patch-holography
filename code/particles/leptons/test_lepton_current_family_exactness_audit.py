#!/usr/bin/env python3
"""Validate the current-family charged-lepton exactness audit artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
READOUT_SCRIPT = ROOT / "particles" / "leptons" / "derive_lepton_log_spectrum_readout.py"
FORWARD_SCRIPT = ROOT / "particles" / "leptons" / "build_forward_charged_leptons.py"
OBSTRUCTION_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_sector_local_current_support_obstruction_certificate.py"
SUPPORT_EXTENSION_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_sector_local_minimal_source_support_extension_emitter.py"
COMPLETION_LAW_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_sector_local_support_extension_completion_law.py"
PAIR_READBACK_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_sector_local_support_extension_source_scalar_pair_readback.py"
ENDPOINT_RATIO_BREAKER_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_sector_local_support_extension_endpoint_ratio_breaker.py"
ABSOLUTE_SCALE_GAP_IDENTITY_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_absolute_scale_transport_gap_identity.py"
AUDIT_SCRIPT = ROOT / "particles" / "leptons" / "derive_lepton_current_family_exactness_audit.py"
OUTPUT = ROOT / "particles" / "runs" / "leptons" / "lepton_current_family_exactness_audit.json"


def test_lepton_exactness_audit_identifies_common_shift_as_insufficient() -> None:
    subprocess.run([sys.executable, str(READOUT_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(FORWARD_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(OBSTRUCTION_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SUPPORT_EXTENSION_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(COMPLETION_LAW_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(PAIR_READBACK_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(ENDPOINT_RATIO_BREAKER_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(ABSOLUTE_SCALE_GAP_IDENTITY_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(AUDIT_SCRIPT)], check=True, cwd=ROOT)

    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_lepton_current_family_exactness_audit"
    assert payload["centered_hierarchy_audit"]["ratio_invariant_under_common_shift"] is True
    assert payload["centered_hierarchy_audit"]["residual_norm"] > 1.0
    assert payload["centered_hierarchy_audit"]["target_sigma_total_log_per_side"] > payload["centered_hierarchy_audit"]["current_sigma_total_log_per_side"]
    assert payload["centered_hierarchy_audit"]["target_gap_pair"]["rho_gap_ratio"] > payload["centered_hierarchy_audit"]["current_gap_pair"]["rho_gap_ratio"]
    assert payload["active_readout_contract"]["hierarchy_mode"] == "two_scalar_ordered_gap_pair_candidate"
    assert abs(payload["active_readout_contract"]["eta_e_split_log_per_side"]) < 1.0e-12
    assert payload["two_scalar_on_current_sigma_audit"]["residual_norm_after_best_eta_on_current_sigma"] > 1.0
    rel_errors = payload["common_shift_audit"]["relative_errors_after_best_shift"]
    assert max(abs(value) for value in rel_errors) > 0.5
    assert payload["same_support_obstruction_audit"]["same_support_exhausted"] is True
    assert payload["same_support_obstruction_audit"]["sigma_support_gap"] > 0.0
    assert payload["smallest_constructive_missing_object"] == "eta_source_support_extension_log_per_side"
    assert payload["current_support_obstruction_certificate"]["proof_status"] == "current_support_obstruction_certificate_closed"
    assert payload["support_extension_emitter"]["artifact"] == "oph_charged_sector_local_minimal_source_support_extension_emitter"
    assert payload["support_extension_completion_law"]["artifact"] == "oph_charged_sector_local_support_extension_completion_law"
    assert payload["support_extension_endpoint_ratio_breaker"]["artifact"] == "oph_charged_sector_local_support_extension_endpoint_ratio_breaker"
    assert payload["support_extension_source_scalar_pair_readback"]["artifact"] == "oph_charged_sector_local_support_extension_source_scalar_pair_readback"
    assert payload["absolute_scale_gap_identity"]["artifact"] == "oph_charged_absolute_scale_transport_gap_identity"
    assert abs(payload["absolute_scale_gap_identity"]["identity_residual"]) < 1.0e-12
    waiting = payload["exact_waiting_set"]
    assert waiting["mandatory_package_a"]["id"] == "charged_sector_response_pushforward_to_C_hat_e"
    assert waiting["mandatory_package_b"]["id"] == "charged_common_refinement_transport_equalizer"
    assert waiting["optional_package_c"]["id"] == "charged_holonomy_bridge_for_legacy_delta_2_over_9"
