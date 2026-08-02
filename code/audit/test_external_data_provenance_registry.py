"""Regression and adversarial checks for the issue-553 provenance registry."""

from __future__ import annotations

import copy
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
sys.path.insert(0, str(TOOLS))

import check_external_data_provenance as provenance  # noqa: E402


REGISTRY_PATH = ROOT / "code/audit/external_data_provenance_registry.json"
HADRON_ARTIFACT = (
    ROOT
    / "code/particles/runs/hadron/empirical_ee_hadronic_spectral_measure.json"
)
HADRON_LOADER = ROOT / "code/particles/hadron/ingest_empirical_ee_hadrons.py"
PLANCK_ARTIFACT = (
    ROOT
    / "code/capacity_readback/planck_posterior/planck_lambda_to_N_propagation.json"
)
PLANCK_LOADER = ROOT / "code/capacity_readback/planck_posterior/propagate.py"
FLAG_ARTIFACT = (
    ROOT / "code/particles/data/flag_2024_light_quark_ratio_fixture.json"
)
FLAG_LOADER = (
    ROOT
    / "code/particles/scripts/"
    "generate_flag_2024_light_quark_ratio_fixture.py"
)
VUS_ARTIFACT = (
    ROOT / "code/particles/data/pdg_2024_vus_kmu2_fixture.json"
)
VUS_LOADER = (
    ROOT
    / "code/particles/scripts/"
    "generate_pdg_2024_vus_kmu2_fixture.py"
)


@pytest.fixture
def registry() -> dict:
    return json.loads(REGISTRY_PATH.read_text(encoding="utf-8"))


def _entry(registry: dict, dataset_id: str) -> dict:
    return next(
        item
        for item in registry["entries"]
        if item["dataset_id"] == dataset_id
    )


def _load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_registry_pins_every_local_artifact_loader_and_license(registry: dict) -> None:
    summary = provenance.validate_registry(registry)
    assert summary["pass"] is True
    assert summary["entries"] == 12
    assert summary["artifact_pins_checked"] == 12
    assert summary["loader_pins_checked"] == 12
    assert summary["upstream_file_pins_checked"] == 10
    assert summary["license_noassertion_entries"] == 11


def test_mandatory_external_dataset_cannot_disappear(registry: dict) -> None:
    mutant = copy.deepcopy(registry)
    mutant["entries"] = [
        entry
        for entry in mutant["entries"]
        if entry["dataset_id"] != "pdg-2026-wz-running-width-target-fixture"
    ]
    with pytest.raises(
        provenance.ProvenanceError,
        match="does not exactly match the mandatory artifact set",
    ):
        provenance.validate_registry(mutant)


def test_auger_diagnostic_is_registered_as_transcribed_and_exposed(
    registry: dict,
) -> None:
    entry = _entry(registry, "auger-2022-fz12-photon-threshold-diagnostic")
    assert entry["artifact_role"] == (
        "exposed_retrospective_conditional_photon_threshold_diagnostic"
    )
    assert entry["raw_inputs"]["classification"] == (
        "hand_transcribed_published_constants_no_raw_payload"
    )
    assert entry["raw_inputs"]["upstream_files"] == []
    assert entry["loader"]["network_required"] is False
    assert entry["loader"]["deterministic_from_declared_inputs"] is True
    assert "https://arxiv.org/abs/2112.06773" in entry["source"]["urls"]


def test_bd_upstream_pins_are_bound_to_the_source_packet(registry: dict) -> None:
    mutant = copy.deepcopy(registry)
    bd = _entry(
        mutant,
        "bouchard-donagi-threshold-spectrum-literature-packet",
    )
    bd["raw_inputs"]["upstream_files"][0]["sha256"] = "2" * 64
    with pytest.raises(
        provenance.ProvenanceError,
        match="differ between the source packet and provenance registry",
    ):
        provenance.validate_registry(mutant)


def test_artifact_hash_drift_fails_closed(registry: dict) -> None:
    mutant = copy.deepcopy(registry)
    mutant["entries"][0]["artifact_sha256"] = "0" * 64
    with pytest.raises(provenance.ProvenanceError, match="SHA-256 mismatch"):
        provenance.validate_registry(mutant)


def test_license_status_cannot_be_omitted(registry: dict) -> None:
    mutant = copy.deepcopy(registry)
    del mutant["entries"][0]["license"]["expression"]
    with pytest.raises(provenance.ProvenanceError, match="fields must be exactly"):
        provenance.validate_registry(mutant)


def test_transcribed_constants_cannot_claim_an_invented_raw_hash(registry: dict) -> None:
    mutant = copy.deepcopy(registry)
    hadron = _entry(mutant, "knt19-pdg2025-hadronic-spectral-shape")
    hadron["raw_inputs"]["upstream_files"] = [
        {
            "id": "invented-raw-table",
            "url": "https://example.invalid/not-a-retained-payload",
            "sha256": "1" * 64,
            "bytes": 1,
        }
    ]
    with pytest.raises(provenance.ProvenanceError, match="must be empty"):
        provenance.validate_registry(mutant)


def test_hadron_loader_reproduces_scientific_payload_except_declared_timestamp(
    registry: dict,
) -> None:
    entry = _entry(registry, "knt19-pdg2025-hadronic-spectral-shape")
    assert entry["loader"]["declared_variable_metadata_json_pointers"] == [
        "/data_release/retrieved_utc"
    ]
    stored = json.loads(HADRON_ARTIFACT.read_text(encoding="utf-8"))
    rebuilt = _load_module(HADRON_LOADER, "issue_553_hadron_loader").build_payload()
    stored["data_release"].pop("retrieved_utc")
    rebuilt["data_release"].pop("retrieved_utc")
    assert rebuilt == stored


def test_planck_gaussian_approximation_loader_is_byte_exact(
    tmp_path: Path,
    registry: dict,
) -> None:
    entry = _entry(
        registry,
        "planck2018-table2-lambda-to-N-gaussian-approximation",
    )
    assert entry["loader"]["declared_variable_metadata_json_pointers"] == []
    copied_loader = tmp_path / PLANCK_LOADER.name
    shutil.copy2(PLANCK_LOADER, copied_loader)
    environment = os.environ.copy()
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    subprocess.run(
        [sys.executable, str(copied_loader)],
        cwd=tmp_path,
        check=True,
        capture_output=True,
        text=True,
        env=environment,
    )
    rebuilt = tmp_path / PLANCK_ARTIFACT.name
    rebuilt_bytes = rebuilt.read_bytes()
    assert b"\r\n" not in rebuilt_bytes
    assert rebuilt_bytes == PLANCK_ARTIFACT.read_bytes()


def test_flag_loader_is_byte_exact_and_derived_fields_are_stable(
    tmp_path: Path,
    registry: dict,
) -> None:
    entry = _entry(registry, "flag-2024-light-quark-ratio-fixture")
    assert entry["loader"]["declared_variable_metadata_json_pointers"] == []
    copied_loader = (
        tmp_path
        / "code/particles/scripts"
        / FLAG_LOADER.name
    )
    copied_loader.parent.mkdir(parents=True)
    shutil.copy2(FLAG_LOADER, copied_loader)
    environment = os.environ.copy()
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    subprocess.run(
        [sys.executable, str(copied_loader)],
        cwd=tmp_path,
        check=True,
        capture_output=True,
        text=True,
        env=environment,
    )
    rebuilt = (
        tmp_path
        / "code/particles/data"
        / FLAG_ARTIFACT.name
    )
    rebuilt_bytes = rebuilt.read_bytes()
    assert b"\r\n" not in rebuilt_bytes
    assert rebuilt_bytes == FLAG_ARTIFACT.read_bytes()
    payload = json.loads(rebuilt_bytes)
    assert payload["averages"][0]["derived_ms_over_md"]["value"] == "19.9437775"
    assert payload["averages"][1]["derived_ms_over_md"]["value"] == "20.35935"


def test_flag_fixture_cannot_invent_covariance_or_preregistration(
    tmp_path: Path,
    registry: dict,
) -> None:
    flag = _entry(registry, "flag-2024-light-quark-ratio-fixture")
    fixture = json.loads(FLAG_ARTIFACT.read_text(encoding="utf-8"))
    mutant_path = tmp_path / FLAG_ARTIFACT.name
    fixture["derived_quantity"]["input_covariance_available"] = True
    mutant_path.write_text(json.dumps(fixture), encoding="utf-8")
    with pytest.raises(
        provenance.ProvenanceError,
        match="must not invent a covariance",
    ):
        provenance._validate_content_boundary(flag, mutant_path)

    fixture["derived_quantity"]["input_covariance_available"] = False
    fixture["claim_boundary"]["significance_gate_preregistered"] = True
    mutant_path.write_text(json.dumps(fixture), encoding="utf-8")
    with pytest.raises(
        provenance.ProvenanceError,
        match="no-theory-uncertainty/no-fit boundary",
    ):
        provenance._validate_content_boundary(flag, mutant_path)


def test_vus_fixture_loader_is_byte_exact_and_uncertainty_is_pinned(
    tmp_path: Path,
    registry: dict,
) -> None:
    entry = _entry(
        registry,
        "pdg-2024-vus-kmu2-compare-only-fixture",
    )
    assert entry["loader"]["declared_variable_metadata_json_pointers"] == []
    copied_loader = (
        tmp_path
        / "code/particles/scripts"
        / VUS_LOADER.name
    )
    copied_loader.parent.mkdir(parents=True)
    shutil.copy2(VUS_LOADER, copied_loader)
    environment = os.environ.copy()
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    subprocess.run(
        [sys.executable, str(copied_loader)],
        cwd=tmp_path,
        check=True,
        capture_output=True,
        text=True,
        env=environment,
    )
    rebuilt = (
        tmp_path
        / "code/particles/data"
        / VUS_ARTIFACT.name
    )
    rebuilt_bytes = rebuilt.read_bytes()
    assert b"\r\n" not in rebuilt_bytes
    assert rebuilt_bytes == VUS_ARTIFACT.read_bytes()
    payload = json.loads(rebuilt_bytes)
    assert payload["coordinate"]["value"] == "0.2250"
    assert payload["coordinate"]["standard_uncertainty"] == "0.0004"


def test_vus_fixture_cannot_be_retyped_as_a_global_fit(
    tmp_path: Path,
    registry: dict,
) -> None:
    entry = _entry(
        registry,
        "pdg-2024-vus-kmu2-compare-only-fixture",
    )
    payload = json.loads(VUS_ARTIFACT.read_text(encoding="utf-8"))
    payload["claim_boundary"]["global_ckm_fit_value"] = True
    mutant_path = tmp_path / VUS_ARTIFACT.name
    mutant_path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(
        provenance.ProvenanceError,
        match="crossed its compare-only boundary",
    ):
        provenance._validate_content_boundary(entry, mutant_path)
