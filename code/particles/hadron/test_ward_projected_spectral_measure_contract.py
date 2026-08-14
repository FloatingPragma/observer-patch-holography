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


def test_contract_emits_scientific_production_boundary() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    schema = json.loads(SCHEMA.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_ward_projected_spectral_measure_contract"
    assert payload["classification"] == "scientific_contract_without_production_data"
    assert payload["production_data_supplied"] is False
    assert payload["promotion_allowed"] is False
    boundary = payload["production_boundary"]
    assert boundary["required_artifact"] == "oph_qcd_ward_projected_hadronic_spectral_measure"
    assert boundary["requires_working_oph_hadron_backend"] is True
    assert boundary["no_go_without_production_payload"] is True
    assert boundary["local_surrogate_promotable"] is False
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
    assert companion["classification"] == "oph_plus_empirical_hadron_closure"
    assert companion["satisfies_production_contract"] is False
    # the companion never replaces the production target
    assert payload["production_boundary"]["required_artifact"] == (
        "oph_qcd_ward_projected_hadronic_spectral_measure")
    assert payload["promotion_allowed"] is False

    companion_schema = json.loads(
        (ROOT / companion["schema"]).read_text(encoding="utf-8"))
    assert companion_schema["properties"]["artifact"]["const"] == companion["artifact"]
    guard_props = companion_schema["properties"]["guards"]["properties"]
    assert guard_props["promotable_as_oph_source_theorem"]["const"] is False
    assert guard_props["surrogate_hadron_artifact"]["const"] is False


def test_contract_names_local_real_engine_without_promoting_it() -> None:
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    engine = payload["local_real_engine"]
    assert engine["classification"] == "real_lattice_diagnostic_toy_scale"
    assert engine["satisfies_production_contract"] is False
    assert (ROOT / engine["runner"]).exists()
    assert (ROOT / engine["package"]).is_dir()


def test_contract_generation_is_byte_deterministic(tmp_path: pathlib.Path) -> None:
    first = tmp_path / "first.json"
    second = tmp_path / "second.json"
    subprocess.run(
        [sys.executable, str(SCRIPT), "--output", str(first)], check=True, cwd=ROOT)
    subprocess.run(
        [sys.executable, str(SCRIPT), "--output", str(second)], check=True, cwd=ROOT)
    assert first.read_bytes() == second.read_bytes()
