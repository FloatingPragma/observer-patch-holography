#!/usr/bin/env python3
"""Smoke tests for the Ward-projected spectral-measure contract."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "hadron" / "derive_ward_projected_spectral_measure_contract.py"
OUTPUT = ROOT / "particles" / "runs" / "hadron" / "ward_projected_spectral_measure_contract.json"
SCHEMA = ROOT / "particles" / "hadron" / "ward_projected_spectral_measure.schema.json"


def test_contract_emits_constructive_spectral_measure_target() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    schema = json.loads(SCHEMA.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_ward_projected_spectral_measure_contract"
    assert payload["constructive_next_artifact"] == "oph_qcd_ward_projected_hadronic_spectral_measure"
    assert payload["promotion_allowed"] is False
    assert payload["current_local_scope"] == "closed_out_of_scope_computationally_blocked"
    assert payload["github_issues_closed_out_of_scope"] == [153, 157]
    assert payload["hardware_gate"]["requires_working_oph_hadron_backend"] is True
    assert payload["hardware_gate"]["chrome_workers_useful_for_backend_execution"] is False
    assert "stable_channel_only_backend_export" in payload["forbidden_promotions"]
    assert schema["properties"]["artifact"]["const"] == "oph_qcd_ward_projected_hadronic_spectral_measure"
    required = set(schema["required"])
    assert "finite_volume_levels" in required
    assert "ward_projected_residues" in required
    assert "systematics" in required


def test_contract_names_empirical_companion_without_promoting_it() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))

    companion = payload["empirical_companion"]
    assert companion["artifact"] == "oph_empirical_ward_projected_hadronic_spectral_measure"
    assert companion["row_class"] == "oph_plus_empirical_hadron_closure"
    assert companion["satisfies_constructive_next_artifact"] is False
    # the companion never replaces the production target
    assert payload["constructive_next_artifact"] == (
        "oph_qcd_ward_projected_hadronic_spectral_measure")
    assert payload["promotion_allowed"] is False

    companion_schema = json.loads(
        (ROOT / companion["schema"]).read_text(encoding="utf-8"))
    assert companion_schema["properties"]["artifact"]["const"] == companion["artifact"]
    guard_props = companion_schema["properties"]["guards"]["properties"]
    assert guard_props["promotable_as_oph_source_theorem"]["const"] is False
    assert guard_props["surrogate_hadron_artifact"]["const"] is False
    assert guard_props["satisfies_production_constructive_next_artifact"]["const"] is False


def test_contract_names_local_real_engine_without_promoting_it() -> None:
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    engine = payload["local_real_engine"]
    assert engine["execution_class"] == "real_lattice_diagnostic_toy_scale"
    assert engine["satisfies_constructive_next_artifact"] is False
    assert (ROOT / engine["runner"]).exists()
    assert (ROOT / engine["package"]).is_dir()
