#!/usr/bin/env python3
"""Validate the D10 electroweak transport-kernel boundary artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
FAMILY_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_observable_family.py"
SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_transport_kernel.py"
OUTPUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_transport_kernel.json"


def test_d10_ew_transport_kernel_boundary() -> None:
    subprocess.run([sys.executable, str(FAMILY_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_d10_ew_transport_kernel"
    assert payload["kernel_id"] == "EWTransportKernel_D10"
    assert payload["readout_coherence_clause_id"] == "EWTransportReadoutCoherence_D10"
    assert payload["promotion_gate"]["smaller_exact_missing_clause"] == "EWTransportReadoutCoherence_D10"
    assert payload["promotion_gate"]["strictly_smaller_next_subclause"] == "EWScalarProvenanceEquality_D10"
    assert payload["promotion_gate"]["z_only_surrogate_forbidden"] is True
    assert payload["family_purity_gate"]["no_run_pole_mix"] is True
    assert payload["family_purity_gate"]["common_readout_certified"] is False
    assert payload["source_pair_population_status"] == "candidate_only"
    assert payload["source_pair_two_scalar_family"]["tau_Y_formula"] == "sigma_EW - eta_EW"
    assert payload["source_pair_first_nonzero_trial"]["tau_Y"] is not None
    assert payload["readout_assignments"]["MW_pole"]["origin_kernel_id"] == "EWTransportKernel_D10"
    assert payload["scalar_provenance"]["delta_alpha"]["family_source_id"] == "d10_running_tree"
    assert payload["common_provenance_witness"]["all_equal_origin_kernel_id"] is True
    assert payload["coherence_witness"]["mixed_sources_detected"] is True
