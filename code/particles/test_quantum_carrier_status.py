#!/usr/bin/env python3
"""Mutation and replay gates for the issue-552 carrier status packet."""

from __future__ import annotations

import copy
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "code/particles/scripts/build_quantum_carrier_status.py"
INDEPENDENT = REPO_ROOT / "code/particles/verify_quantum_carrier_status_independent.py"
PROJECTION = (
    REPO_ROOT / "code/particles/manifests/quantum_carrier_status_source_projection.json"
)
RECEIPT = REPO_ROOT / "code/particles/runs/status/quantum_carrier_status.json"
MARKDOWN = REPO_ROOT / "code/particles/QUANTUM_CARRIER_STATUS.md"


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _producer():
    return _load(SCRIPT, "build_quantum_carrier_status_test")


def _verifier():
    return _load(INDEPENDENT, "verify_quantum_carrier_status_test")


def _rehash(module, payload: dict) -> dict:
    payload["receipt_sha256"] = module._receipt_hash(payload)
    return payload


def test_typed_classical_vector_is_separate_from_quantum_verdicts() -> None:
    module = _producer()
    payload = module.build_payload()
    assert payload["classical_mode_vector_order"] == ["photon", "gluon", "graviton"]
    assert payload["classical_mode_vector"] == [2, 16, 2]
    assert payload["continuum_spacetime_dimension"] == 4
    assert payload["comparison_values_consumed"] is False
    assert payload["blind_prediction_eligible"] is False
    assert payload["target_named_status_rows"] is True
    assert payload["status"] == "THREE_ROW_EXPLICIT_NOT_EVALUABLE"

    rows = {row["carrier_id"]: row for row in payload["rows"]}
    expected = {
        "photon": (1, "u1_lie_algebra_generator", 2),
        "gluon": (8, "su3_adjoint_generator", 16),
        "graviton": (1, "symmetric_metric_tensor_field", 2),
    }
    for carrier_id, (factor, role, total) in expected.items():
        baseline = rows[carrier_id]["classical_baseline"]
        assert baseline["multiplicity_factor"] == factor
        assert baseline["multiplicity_role"] == role
        assert baseline["exact_total_mode_count"] == total
        assert baseline["continuum_spacetime_dimension"] == 4
        assert baseline["branch_is_additional_input_not_group_output"] is True
        assert "gauge_algebra_dimension" not in baseline
        assert baseline["particle_claim"] is False
        assert rows[carrier_id]["particle_promotion_allowed"] is False
        assert rows[carrier_id]["verdict_class"] == "EXPLICIT_NOT_EVALUABLE"


def test_projection_pins_bounded_source_frontiers_without_exhaustive_claim() -> None:
    module = _producer()
    projection = module.load_and_validate_projection()
    assert projection["comparison_policy"] == {
        "blind_prediction_eligible": False,
        "comparison_data_present": False,
        "comparison_values_consumed": False,
        "laboratory_values_present": False,
        "status_packet_only": True,
        "target_named_rows_present": True,
    }
    boundaries = {
        row["boundary_id"]: row for row in projection["declared_boundary_evidence"]
    }
    assert set(boundaries) == {
        "finite_local_domain",
        "electroweak_source_action",
        "inhabited_einstein_tower",
        "lorentzian_quantum_eft_transfer",
        "qcd_spectral_resource",
    }
    assert boundaries["qcd_spectral_resource"]["evidence_class"] == (
        "bounded_local_resource_receipt"
    )
    assert "limited to the pinned declared corpus" in projection["scope"]


def test_every_row_has_the_full_typed_capability_inventory() -> None:
    module = _producer()
    payload = module.build_payload()
    expected_keys = {
        "state_space",
        "observable_algebra",
        "gauge_quotient",
        "vacuum",
        "spectral_object",
        "physical_current_residue",
        "refinement_control",
        "phase_or_asymptotic_sector",
    }
    for row in payload["rows"]:
        capabilities = row["capabilities"]
        assert set(capabilities) == expected_keys
        assert capabilities["state_space"]["physical_quantum_object_available"] is False
        assert capabilities["observable_algebra"]["classical_object_available"] is True
        assert (
            capabilities["observable_algebra"]["physical_quantum_object_available"]
            is False
        )
        assert (
            capabilities["gauge_quotient"]["physical_quantum_object_available"]
            is False
        )
        assert (
            capabilities["vacuum"]["source_selected_quantum_object_available"]
            is False
        )
        assert (
            capabilities["spectral_object"][
                "positive_physical_quantum_object_available"
            ]
            is False
        )
        assert (
            capabilities["physical_current_residue"][
                "nonzero_positive_residue_available"
            ]
            is False
        )
        assert capabilities["refinement_control"]["available"] is False
        assert (
            capabilities["phase_or_asymptotic_sector"][
                "physical_quantum_sector_available"
            ]
            is False
        )


@pytest.mark.parametrize(
    "mutation",
    [
        "aggregate_status",
        "comparison_value",
        "blind_eligibility",
        "mode_vector",
        "dimension",
        "row_verdict",
        "particle_promotion",
        "positive_residue",
        "capability_status",
        "strongest_statement",
        "multiplicity_factor",
        "multiplicity_role",
        "branch_boundary",
        "blocking_frontier",
        "open_interfaces",
        "resource_boundaries",
        "source_pin",
        "unknown_field",
    ],
)
def test_primary_validator_rejects_semantic_mutations(mutation: str) -> None:
    module = _producer()
    payload = copy.deepcopy(module.build_payload())
    if mutation == "aggregate_status":
        payload["status"] = "SOURCE_POSITIVE"
    elif mutation == "comparison_value":
        payload["comparison_values_consumed"] = True
    elif mutation == "blind_eligibility":
        payload["blind_prediction_eligible"] = True
    elif mutation == "mode_vector":
        payload["classical_mode_vector"] = [2, 15, 2]
    elif mutation == "dimension":
        payload["continuum_spacetime_dimension"] = 5
    elif mutation == "row_verdict":
        payload["rows"][0]["verdict"] = "SOURCE_POSITIVE"
    elif mutation == "particle_promotion":
        payload["rows"][0]["particle_promotion_allowed"] = True
    elif mutation == "positive_residue":
        payload["rows"][0]["capabilities"]["physical_current_residue"][
            "nonzero_positive_residue_available"
        ] = True
    elif mutation == "capability_status":
        payload["rows"][0]["capabilities"]["state_space"]["status"] = (
            "SOURCE_POSITIVE"
        )
    elif mutation == "strongest_statement":
        payload["rows"][0]["strongest_supported_statement"] = (
            "OPH predicts a measured zero photon rest mass."
        )
    elif mutation == "multiplicity_factor":
        payload["rows"][2]["classical_baseline"]["multiplicity_factor"] = 99
    elif mutation == "multiplicity_role":
        payload["rows"][2]["classical_baseline"]["multiplicity_role"] = (
            "gauge_algebra_dimension"
        )
    elif mutation == "branch_boundary":
        payload["rows"][0]["classical_baseline"][
            "branch_is_additional_input_not_group_output"
        ] = False
    elif mutation == "blocking_frontier":
        payload["rows"][0]["blocking_frontier"] = []
    elif mutation == "open_interfaces":
        payload["rows"][0]["open_interfaces"] = []
    elif mutation == "resource_boundaries":
        payload["rows"][1]["resource_boundaries"] = []
    elif mutation == "source_pin":
        payload["source_projection_pin"]["sha256"] = "sha256:" + "0" * 64
    elif mutation == "unknown_field":
        payload["rows"][0]["measured_mass_gev"] = 0
    else:  # pragma: no cover
        raise AssertionError(mutation)
    _rehash(module, payload)
    with pytest.raises(module.CertificateError):
        module.validate_payload(payload)


@pytest.mark.parametrize(
    "mutation",
    [
        "strongest_statement",
        "forbidden_gauge_dimension",
        "multiplicity_factor",
        "capability_status",
        "blocking_frontier",
        "open_interfaces",
        "policy",
        "unknown_field",
    ],
)
def test_independent_reconstruction_rejects_rehashed_semantic_mutations(
    mutation: str,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    verifier = _verifier()
    payload = json.loads(RECEIPT.read_text(encoding="utf-8"))
    if mutation == "strongest_statement":
        payload["rows"][0]["strongest_supported_statement"] = (
            "OPH predicts a measured zero photon rest mass."
        )
    elif mutation == "forbidden_gauge_dimension":
        payload["rows"][2]["classical_baseline"]["gauge_algebra_dimension"] = 1
    elif mutation == "multiplicity_factor":
        payload["rows"][2]["classical_baseline"]["multiplicity_factor"] = 99
    elif mutation == "capability_status":
        payload["rows"][0]["capabilities"]["state_space"]["status"] = (
            "SOURCE_POSITIVE"
        )
    elif mutation == "blocking_frontier":
        payload["rows"][0]["blocking_frontier"] = []
    elif mutation == "open_interfaces":
        payload["rows"][0]["open_interfaces"] = []
    elif mutation == "policy":
        payload["producer_scope"] = "Exhaustive theorem over every OPH completion."
    elif mutation == "unknown_field":
        payload["measured_target"] = 0
    else:  # pragma: no cover
        raise AssertionError(mutation)
    _rehash(verifier, payload)
    altered = tmp_path / "altered_receipt.json"
    altered.write_text(json.dumps(payload), encoding="utf-8")
    monkeypatch.setattr(verifier, "RECEIPT_PATH", altered)
    with pytest.raises(verifier.VerificationError):
        verifier.verify()


def test_projection_drift_is_rejected(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    module = _producer()
    projection = json.loads(module.PROJECTION_PATH.read_text(encoding="utf-8"))
    projection["classical_mode_vector"] = [2, 15, 2]
    altered = tmp_path / "projection.json"
    altered.write_text(json.dumps(projection), encoding="utf-8")
    monkeypatch.setattr(module, "PROJECTION_PATH", altered)
    with pytest.raises(module.CertificateError) as error:
        module.load_and_validate_projection()
    assert error.value.code == "PROJECTION"


def test_bounded_parent_drift_is_rejected(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    module = _producer()
    boundary = json.loads(module.LOCAL_DOMAIN_BOUNDARY_PATH.read_text(encoding="utf-8"))
    boundary["physical_promotion_allowed"] = True
    altered = tmp_path / "local_domain.json"
    altered.write_text(json.dumps(boundary), encoding="utf-8")
    monkeypatch.setattr(module, "LOCAL_DOMAIN_BOUNDARY_PATH", altered)
    with pytest.raises(module.CertificateError) as error:
        module.reconstruct_source_projection()
    assert error.value.code == "LOCAL_DOMAIN_BOUNDARY"


def test_independent_bounded_parent_drift_is_rejected(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    verifier = _verifier()
    boundary = json.loads(verifier.LOCAL_DOMAIN_PATH.read_text(encoding="utf-8"))
    boundary["physical_promotion_allowed"] = True
    altered = tmp_path / "local_domain.json"
    altered.write_text(json.dumps(boundary), encoding="utf-8")
    monkeypatch.setattr(verifier, "LOCAL_DOMAIN_PATH", altered)
    with pytest.raises(verifier.VerificationError):
        verifier.verify()


def test_generated_outputs_are_byte_exact() -> None:
    module = _producer()
    projection = module.reconstruct_source_projection()
    payload = module.build_payload()
    assert PROJECTION.read_bytes() == module._pretty_bytes(projection)
    assert RECEIPT.read_bytes() == module._pretty_bytes(payload)
    assert MARKDOWN.read_bytes() == module.render_markdown(payload).encode("utf-8")
    module.validate_committed(RECEIPT, MARKDOWN)


def test_independent_reconstruction_passes() -> None:
    completed = subprocess.run(
        [sys.executable, str(INDEPENDENT)],
        cwd=REPO_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    assert "QUANTUM_CARRIER_STATUS_INDEPENDENT_VALID" in completed.stdout
