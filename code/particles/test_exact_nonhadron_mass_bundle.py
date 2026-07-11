#!/usr/bin/env python3
"""Validate the exact non-hadron mass bundle."""

from __future__ import annotations

import importlib.util
import json
import pathlib
import subprocess
import sys
import tempfile

import pytest


ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "particles" / "scripts" / "build_exact_nonhadron_mass_bundle.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("build_exact_nonhadron_mass_bundle", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_neutrino_absolute_promotion_requires_explicit_non_circularity_pass() -> None:
    module = _load_module()
    candidate = {
        "artifact": "oph_neutrino_absolute_attachment_theorem",
        "status": "theorem_grade_emitted",
        "weighted_cycle_base_eligible": True,
        "prediction_promotion_allowed": True,
        "public_surface_candidate_allowed": True,
    }
    assert module._neutrino_absolute_promotable(candidate) is False
    candidate["non_circularity_status"] = {"promotion_allowed": True}
    assert module._neutrino_absolute_promotable(candidate) is True
    assert module._is_public_mass_output({"legacy_particle_id_slot": True, "exact_kind": "theorem_grade"}) is False


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
        assert payload["status"] == (
            "public_mass_outputs_with_classical_carriers_separated_and_target_anchored_witnesses_withheld"
        )
        entries = {entry["particle_id"]: entry for entry in payload["entries"]}
        carriers = {entry["carrier_id"]: entry for entry in payload["classical_carrier_modes"]}
        withheld = {entry["particle_id"]: entry for entry in payload["withheld_entries"]}
        assert len(entries) == 1
        assert not {"photon", "gluon", "graviton"} & set(entries)
        assert set(carriers) == {"photon", "gluon", "graviton"}
        assert all(row["hard_quadratic_mass_parameter_squared"] == 0 for row in carriers.values())
        assert all(row["quantum_particle_gate"]["passed"] is False for row in carriers.values())
        assert all(row["particle_promotion_allowed"] is False for row in carriers.values())
        assert "w_boson" not in entries
        assert "z_boson" not in entries
        assert entries["higgs"]["mass_gev"] == pytest.approx(125.1995304097179, abs=1.0e-12)
        assert "electron" not in entries
        assert "top_quark" not in entries
        assert "tau_neutrino" not in entries
        assert withheld["electron"]["reason"] == "target_anchored_witness_kept_in_exact_fit_audit_not_public_prediction"
        assert withheld["electron"]["public_theorem_value"] is None
        assert withheld["electron"]["source_only"] is False
        assert withheld["electron"]["centered_log"] == -4.495209645475038
        assert withheld["electron"]["formula_if_anchor_exists"] == "m_e(P)=exp(A_ch(P)-4.495209645475038)"
        assert "charged_determinant_trace_lift_attachment" in withheld["electron"]["missing_for_promotion"]
        assert "NO_TARGET_LEAK_DAG_CHARGED_A_CH" in withheld["electron"]["missing_for_promotion"]
        assert withheld["top_quark"]["reason"] == "target_anchored_witness_kept_in_exact_fit_audit_not_public_prediction"
        assert withheld["tau_neutrino"]["reason"] == "target_informed_candidate_rejected_by_correlated_profile"
        markdown = md.read_text(encoding="utf-8")
        assert "Public Non-Hadron Mass Outputs" in markdown
        assert "Separated Classical Carrier Modes" in markdown
        assert "not emitted as `0 GeV` particle predictions" in markdown
        assert "Bottom Quark" not in markdown
        assert "Tau Neutrino" not in markdown
