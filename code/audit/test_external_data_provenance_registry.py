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
    assert summary["entries"] == 9
    assert summary["artifact_pins_checked"] == 9
    assert summary["loader_pins_checked"] == 9
    assert summary["upstream_file_pins_checked"] == 10
    assert summary["license_noassertion_entries"] == 9


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
