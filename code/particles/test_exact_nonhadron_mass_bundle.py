#!/usr/bin/env python3
"""Validate the exact non-hadron mass bundle."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys
import tempfile

import pytest


ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "particles" / "scripts" / "build_exact_nonhadron_mass_bundle.py"


def test_exact_nonhadron_mass_bundle_is_complete() -> None:
    with tempfile.TemporaryDirectory(prefix="oph_exact_nonhadron_bundle_") as tmpdir:
        tmp = pathlib.Path(tmpdir)
        md = tmp / "EXACT_NONHADRON_MASSES.md"
        js = tmp / "exact_nonhadron_masses.json"
        forward = tmp / "exact_nonhadron_masses_current.json"
        subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--markdown-out",
                str(md),
                "--json-out",
                str(js),
                "--forward-out",
                str(forward),
            ],
            check=True,
            cwd=ROOT,
        )
        payload = json.loads(js.read_text(encoding="utf-8"))
        assert payload["artifact"] == "oph_exact_nonhadron_mass_bundle"
        assert payload["status"] == "public_mass_outputs_with_target_anchored_witnesses_withheld"
        entries = {entry["particle_id"]: entry for entry in payload["entries"]}
        withheld = {entry["particle_id"]: entry for entry in payload["withheld_entries"]}
        assert len(entries) == 4
        assert entries["photon"]["mass_gev"] == pytest.approx(0.0, abs=1.0e-18)
        assert "w_boson" not in entries
        assert "z_boson" not in entries
        assert entries["higgs"]["mass_gev"] == pytest.approx(125.1995304097179, abs=1.0e-12)
        assert "electron" not in entries
        assert "top_quark" not in entries
        assert "tau_neutrino" not in entries
        assert withheld["electron"]["reason"] == "target_anchored_witness_kept_in_exact_fit_audit_not_public_prediction"
        assert withheld["electron"]["public_theorem_value"] is None
        assert withheld["electron"]["source_only"] is False
        assert withheld["electron"]["centered_log"] == -4.495210108808081
        assert withheld["electron"]["formula_if_anchor_exists"] == "m_e(P)=exp(A_ch(P)-4.495210108808081)"
        assert "charged_determinant_trace_lift_attachment" in withheld["electron"]["missing_for_promotion"]
        assert "NO_TARGET_LEAK_DAG_CHARGED_A_CH" in withheld["electron"]["missing_for_promotion"]
        assert withheld["top_quark"]["reason"] == "target_anchored_witness_kept_in_exact_fit_audit_not_public_prediction"
        assert withheld["tau_neutrino"]["reason"] == "compare_only_absolute_or_adapter_surface_kept_out_of_public_prediction_table"
        markdown = md.read_text(encoding="utf-8")
        assert "Public Non-Hadron Mass Outputs" in markdown
        assert "Bottom Quark" not in markdown
        assert "Tau Neutrino" not in markdown
