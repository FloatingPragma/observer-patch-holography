#!/usr/bin/env python3
"""Validate the repaired neutrino weighted-cycle branch artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "neutrino" / "derive_neutrino_weighted_cycle_repair.py"
CERTIFICATE = ROOT / "particles" / "runs" / "neutrino" / "same_label_scalar_certificate.json"
COCYCLE = ROOT / "particles" / "runs" / "flavor" / "overlap_edge_transport_cocycle.json"
PHASE_SOURCE = ROOT / "particles" / "runs" / "neutrino" / "intrinsic_neutrino_mass_eigenstate_bundle_from_scalar_certificate.json"
ISOTROPIC = ROOT / "particles" / "runs" / "neutrino" / "forward_majorana_matrix.json"
SELECTOR = ROOT / "particles" / "runs" / "neutrino" / "neutrino_transport_load_segment_selector.json"
FAMILY_KERNEL = ROOT / "particles" / "runs" / "flavor" / "family_transport_kernel.json"


def test_repaired_weighted_cycle_branch_matches_expected_live_numbers() -> None:
    with tempfile.TemporaryDirectory(prefix="oph_neutrino_repair_") as tmpdir:
        out = pathlib.Path(tmpdir) / "repair.json"
        subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
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
                "--family-kernel",
                str(FAMILY_KERNEL),
                "--output",
                str(out),
            ],
            check=True,
            cwd=ROOT,
        )
        payload = json.loads(out.read_text(encoding="utf-8"))
        assert payload["artifact"] == "oph_neutrino_weighted_cycle_repair"
        assert payload["final_verdict"] == "target_informed_template_weighted_cycle_candidate"
        assert payload["theorem_status"] == "exact_linear_algebra_conditional_on_declared_candidate_inputs"
        assert payload["physical_window_status"] == "separate_marginal_box_comparison_only"
        assert payload["historical_target_exposure"] is True
        assert payload["source_only_prediction_eligible"] is False
        assert payload["prediction_promotion_allowed"] is False
        assert payload["public_surface_candidate_allowed"] is False
        assert payload["source_closure_status"]["family_transport_kernel_status"] == "template"
        assert payload["symbolic_absolute_family"]["family_parameter"] == "lambda_nu > 0"
        assert payload["compare_only_atmospheric_anchor"]["adapter_status"] == "hard_separated_compare_only_adapter"
        assert payload["transport_load_selector_status"] == "standard_math_fixed_balanced_segment_selector_closed"
        assert payload["selected_transport_load_selector"] == "balanced_equals_least_distortion_midpoint"
        assert abs(payload["selected_D_nu"] - 1.127883690210334) < 1.0e-15
        assert abs(payload["selected_tau_nu"] - 0.5) < 1.0e-15
        assert abs(payload["weight_exponent"] - 1.395092021318097) < 1.0e-15
        assert abs(payload["diag_loading"] - 1.1045845502912137) < 1.0e-15
        assert abs(payload["pmns_observables"]["theta12_deg"] - 34.225904631810025) < 1.0e-10
        assert abs(payload["pmns_observables"]["theta23_deg"] - 49.72282845058266) < 1.0e-10
        assert abs(payload["pmns_observables"]["theta13_deg"] - 8.686355527700156) < 1.0e-10
        assert abs(payload["dimensionless_ratio_dm21_over_dm32"] - 0.030721110097966534) < 1.0e-15
        assert abs(payload["compare_only_atmospheric_anchor"]["delta_m32_sq_eV2"] - 2.438e-3) < 1.0e-15
        assert payload["pdg_2025_no_3sigma_window"]["within_window"]["theta12_deg"] is True
        assert payload["pdg_2025_no_3sigma_window"]["within_window"]["theta23_deg"] is True
        assert payload["pdg_2025_no_3sigma_window"]["within_window"]["theta13_deg"] is True
