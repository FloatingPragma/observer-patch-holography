from __future__ import annotations

import json
import os
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

import pytest

import fz12_custody_projection as custody_builder
import fz12_observation_map_certificate as certificate


VERIFIER = certificate.HERE / "verify_fz12_observation_map_independent.py"


def resign(value: dict, field: str) -> None:
    body = {key: item for key, item in value.items() if key != field}
    value[field] = certificate.base.tagged_sha256(
        certificate.base.canonical_json_bytes(body)
    )


def write_canonical(path: Path, value: dict) -> None:
    path.write_bytes(certificate.base.canonical_json_bytes(value))


def install_source_mutation(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, mutate
) -> Path:
    source = json.loads(certificate.SOURCE_RECEIPT_PATH.read_text())
    mutate(source)
    resign(source, "receipt_sha256")
    path = tmp_path / "seam_current_edge_prediction_receipt.json"
    write_canonical(path, source)
    monkeypatch.setattr(certificate, "SOURCE_RECEIPT_PATH", path)
    return path


def install_custody_mutation(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, mutate
) -> Path:
    custody = json.loads(certificate.CUSTODY_PROJECTION_PATH.read_text())
    mutate(custody)
    resign(custody, "projection_sha256")
    path = tmp_path / "fz12_custody_projection.json"
    write_canonical(path, custody)
    monkeypatch.setattr(certificate, "CUSTODY_PROJECTION_PATH", path)
    return path


def independent_run(
    *,
    source: Path = certificate.SOURCE_RECEIPT_PATH,
    custody: Path = certificate.CUSTODY_PROJECTION_PATH,
    receipt: Path = certificate.RECEIPT_PATH,
) -> subprocess.CompletedProcess[str]:
    environment = os.environ.copy()
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    return subprocess.run(
        [
            sys.executable,
            str(VERIFIER),
            "--source",
            str(source),
            "--custody-projection",
            str(custody),
            "--receipt",
            str(receipt),
        ],
        cwd=certificate.HERE,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )


def mutated_receipt(tmp_path: Path, mutate) -> Path:
    receipt = json.loads(certificate.RECEIPT_PATH.read_text())
    mutate(receipt)
    resign(receipt, "receipt_sha256")
    path = tmp_path / "fz12_observation_map_receipt.json"
    write_canonical(path, receipt)
    return path


def test_custody_projection_is_canonical_exact_and_data_free() -> None:
    committed = certificate.CUSTODY_PROJECTION_PATH.read_bytes()
    rebuilt = custody_builder.build_projection()
    assert committed == certificate.base.canonical_json_bytes(rebuilt)
    scope = rebuilt["projection_scope"]
    assert scope == {
        "register_row": "FZ-12",
        "source_class": "frozen source theory and custody metadata only",
        "includes_measurement_values": False,
        "includes_comparison_values": False,
        "includes_likelihood_values": False,
        "includes_other_campaign_rows": False,
    }
    assert len(committed) == certificate.EXPECTED_CUSTODY_BYTES
    assert certificate.base.tagged_sha256(committed) == (
        f"sha256:{certificate.EXPECTED_CUSTODY_SHA256}"
    )


def test_producer_has_no_multi_campaign_register_dependency() -> None:
    producer = Path(certificate.__file__).read_text(encoding="utf-8")
    assert "claims/frozen_prediction_register.json" not in producer
    assert "FROZEN_REGISTER_PATH" not in producer
    receipt = certificate.build_receipt()
    assert "frozen_registration" not in receipt["ancestry"]
    assert receipt["ancestry"]["custody_projection"]["row"] == "FZ-12"


def test_exact_fz12_formal_observation_map() -> None:
    result = certificate.build_receipt()["exact_formal_result"]
    assert result["coefficients_in_a_units"] == {
        "omega": {
            "k": "1",
            "a2_k3": "-1/40",
            "a4_k5_isotropic": "19/67200",
            "a4_k5_I6": "-1/25200",
        },
        "radial_group_velocity": {
            "constant": "1",
            "a2_k2": "-3/40",
            "a4_k4_isotropic": "19/13440",
            "a4_k4_I6": "-1/5040",
        },
        "inverse_radial_speed": {
            "constant": "1",
            "a2_k2": "3/40",
            "a4_k4_isotropic": "283/67200",
            "a4_k4_I6": "1/5040",
        },
        "formal_inverse_radial_speed_excess": {
            "a2_k2": "3/40",
            "a4_k4_isotropic": "283/67200",
            "a4_k4_I6": "1/5040",
        },
        "transverse_group_velocity": {"a4_k4_grad_S_I6": "-1/25200"},
    }
    assert result["machine_checks"] == {
        "omega_hat_squared_equals_lambda_hat_mod_q8": True,
        "radial_times_inverse_radial_equals_one_mod_q6": True,
        "formal_inverse_radial_speed_excess_is_inverse_minus_one": True,
        "transverse_coefficient_is_one_over_k_times_spherical_gradient": True,
    }


def test_natural_truncation_orders_are_explicit() -> None:
    receipt = certificate.build_receipt()
    result = receipt["exact_formal_result"]
    assert receipt["formal_input"]["dimensionless_normalization"] == {
        "q": "q = a k",
        "lambda_hat": "lambda_hat = a^2 Lambda",
        "omega_hat": "omega_hat = a omega",
        "normalized_symbol": (
            "lambda_hat = q^2 - q^4/20 + (1/840 - I6/12600) q^6 + O(q^8)"
        ),
    }
    assert result["ring_for_square_root"] == "Q[I6][q] modulo q^8"
    assert result["ring_for_velocity_and_inverse"] == "Q[I6][q] modulo q^6"
    assert "O(a^6 k^7)" in result["frequency"]
    for key in (
        "radial_group_velocity",
        "inverse_radial_speed",
        "formal_inverse_radial_speed_excess",
        "transverse_group_velocity",
    ):
        assert "O(a^6 k^6)" in result[key]


def test_scale_free_coefficients_are_exact() -> None:
    result = certificate.build_receipt()["exact_formal_result"]
    assert result["scale_free_coefficients"] == {
        "omega_k5_isotropic_over_C4_squared": "19/168",
        "omega_k5_I6_over_C4_squared": "-1/63",
        "radial_k4_isotropic_over_C4_squared": "95/168",
        "radial_k4_I6_over_C4_squared": "-5/63",
        "inverse_speed_k4_isotropic_over_C4_squared": "283/168",
        "inverse_speed_k4_I6_over_C4_squared": "5/63",
        "transverse_grad_S_I6_over_C4_squared": "-1/63",
    }


def test_boundary_refuses_time_of_flight_and_dynamics_claims() -> None:
    receipt = certificate.build_receipt()
    boundary = receipt["interpretation_boundary"]
    assert boundary["formal_truncation_only"] is True
    assert boundary["wave_packet_dynamics_proved"] is False
    assert boundary["time_evolution_proved"] is False
    assert boundary["trajectory_or_time_of_flight_map"] is False
    assert "transverse drift" in boundary["statement"]
    assert "inverse_speed_delay_per_unit_coordinate_length" not in json.dumps(receipt)
    assert receipt["formal_input"]["extra_map_premise_discharged_here"] is False


def test_literal_no_data_boundary_has_only_two_target_free_inputs() -> None:
    receipt = certificate.build_receipt()
    exposure = receipt["exposure_boundary"]
    assert exposure["inputs_read"] == [
        "target-free frozen FZ-12 source receipt",
        "FZ-12-only custody projection without measurement values",
    ]
    assert exposure["no_data_inputs_verified"] is True
    assert exposure["comparison_inputs"] == []
    assert not any(
        exposure[key]
        for key in (
            "target_values_read",
            "comparison_data_read",
            "public_measurement_read",
            "comparison_permitted",
        )
    )


def test_committed_receipt_is_canonical_byte_exact_and_self_hashed() -> None:
    committed = certificate.verify_committed_receipt()
    assert (
        certificate.RECEIPT_PATH.read_bytes()
        == certificate.base.canonical_json_bytes(committed)
    )


def test_independent_verifier_is_separate_and_passes_in_subprocess() -> None:
    source = VERIFIER.read_text(encoding="utf-8")
    assert "import fz12_observation_map_certificate" not in source
    assert "from fz12_observation_map_certificate" not in source
    result = independent_run()
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "FZ12_OBSERVATION_MAP_INDEPENDENT_VERIFY_PASS"


@pytest.mark.parametrize("field", ["schema", "status", "receipt_sha256"])
def test_source_type_or_digest_mutation_fails_closed(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, field: str
) -> None:
    source = json.loads(certificate.SOURCE_RECEIPT_PATH.read_text())
    source[field] = "tampered"
    path = tmp_path / "source.json"
    write_canonical(path, source)
    monkeypatch.setattr(certificate, "SOURCE_RECEIPT_PATH", path)
    with pytest.raises(certificate.base.FingerprintError):
        certificate.build_receipt()


@pytest.mark.parametrize("mutation", ["coefficient", "frequency_boundary", "exposure"])
def test_resigned_source_semantic_mutation_fails_fixed_byte_pin(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, mutation: str
) -> None:
    def mutate(source: dict) -> None:
        if mutation == "coefficient":
            source["conditional_physical_candidate"]["coefficients"]["B6_over_a4"] = (
                "1/12600"
            )
        elif mutation == "frequency_boundary":
            source["conditional_physical_candidate"]["symbol_name"] = "omega squared"
        else:
            source["exposure_and_custody_boundary"]["comparison_permitted"] = True

    path = install_source_mutation(monkeypatch, tmp_path, mutate)
    with pytest.raises(certificate.base.FingerprintError, match="fixed byte digest"):
        certificate.build_receipt()
    result = independent_run(source=path)
    assert result.returncode == 1
    assert "source fixed byte digest drift" in result.stderr


@pytest.mark.parametrize(
    "mutation", ["source_commit", "frozen_artifact", "measurement_scope", "source_pin"]
)
def test_resigned_custody_mutations_fail_producer_and_independent_verifier(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, mutation: str
) -> None:
    def mutate(custody: dict) -> None:
        if mutation == "source_commit":
            custody["external_custody_contract"]["source_commit"] = "0" * 40
        elif mutation == "frozen_artifact":
            custody["external_custody_contract"]["artifact_sha256"][
                "seam_current_edge_prediction_frozen_2026-08-02.json"
            ] = "0" * 64
        elif mutation == "measurement_scope":
            custody["projection_scope"]["includes_measurement_values"] = True
        else:
            custody["source_receipt"]["sha256"] = "0" * 64

    path = install_custody_mutation(monkeypatch, tmp_path, mutate)
    with pytest.raises(certificate.base.FingerprintError, match="fixed byte digest"):
        certificate.build_receipt()
    result = independent_run(custody=path)
    assert result.returncode == 1
    assert "custody fixed byte digest drift" in result.stderr


@pytest.mark.parametrize(
    "mutation",
    [
        "coefficient",
        "display_formula",
        "formal_input",
        "machine_check",
        "wave_packet",
        "transverse_warning",
    ],
)
def test_independent_verifier_rejects_resigned_result_and_boundary_mutations(
    tmp_path: Path, mutation: str
) -> None:
    def mutate(receipt: dict) -> None:
        if mutation == "coefficient":
            receipt["exact_formal_result"]["coefficients_in_a_units"]["omega"][
                "a4_k5_I6"
            ] = "1/25200"
        elif mutation == "display_formula":
            receipt["exact_formal_result"]["frequency"] = "omega = k + O(a^6 k^7)"
        elif mutation == "formal_input":
            receipt["formal_input"]["frozen_coefficients"]["B0_over_a4"] = "1/841"
        elif mutation == "machine_check":
            receipt["exact_formal_result"]["machine_checks"][
                "omega_hat_squared_equals_lambda_hat_mod_q8"
            ] = False
        elif mutation == "wave_packet":
            receipt["interpretation_boundary"]["wave_packet_dynamics_proved"] = True
        else:
            receipt["interpretation_boundary"]["statement"] = "formal only"

    result = independent_run(receipt=mutated_receipt(tmp_path, mutate))
    assert result.returncode == 1
    assert "INDEPENDENT_VERIFY_FAIL" in result.stderr


@pytest.mark.parametrize(
    "mutation",
    [
        "extra_physical_claim",
        "ancestry_commit",
        "scope_statement",
        "comparison_state",
        "interpretation_append",
        "extra_nested_key",
    ],
)
def test_independent_verifier_is_schema_closed_for_semantic_mutations(
    tmp_path: Path, mutation: str
) -> None:
    def mutate(receipt: dict) -> None:
        if mutation == "extra_physical_claim":
            receipt["physical_prediction_claim"] = True
        elif mutation == "ancestry_commit":
            receipt["ancestry"]["custody_projection"]["custody_commit"] = "0" * 40
        elif mutation == "scope_statement":
            receipt["certificate_scope"]["statement"] = "this predicts photons"
        elif mutation == "comparison_state":
            receipt["exposure_boundary"]["comparison_state"] = "ELIGIBLE"
        elif mutation == "interpretation_append":
            receipt["interpretation_boundary"]["statement"] += (
                "; this is a measured time delay"
            )
        else:
            receipt["certificate_scope"]["unknown_promotion"] = True

    result = independent_run(receipt=mutated_receipt(tmp_path, mutate))
    assert result.returncode == 1
    assert "INDEPENDENT_VERIFY_FAIL" in result.stderr


def test_wrong_fz12_ray_cannot_emit_the_named_map() -> None:
    with pytest.raises(
        certificate.base.FingerprintError,
        match="FZ-12 formal observation-map coefficient drift",
    ):
        certificate.formal_map(Fraction(-1, 20), Fraction(1, 840), Fraction(1, 12600))


def test_committed_receipt_payload_mutation_fails_closed(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    path = mutated_receipt(
        tmp_path,
        lambda receipt: receipt["exact_formal_result"].update(
            {"formal_inverse_radial_speed_excess": "tampered"}
        ),
    )
    monkeypatch.setattr(certificate, "RECEIPT_PATH", path)
    with pytest.raises(
        certificate.base.FingerprintError, match="ancestry or result drift"
    ):
        certificate.verify_committed_receipt()
