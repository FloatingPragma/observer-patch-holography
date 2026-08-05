from __future__ import annotations

import copy
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

import carrier_frequency_speed_certificate as certificate
import verify_carrier_frequency_speed_independent as independent


VERIFIER = certificate.HERE / "verify_carrier_frequency_speed_independent.py"


def write_canonical(path: Path, value: dict) -> None:
    path.write_bytes(independent.canonical_json_bytes(value))


def resign(value: dict) -> None:
    value["receipt_sha256"] = independent.self_digest(value, "receipt_sha256")


def independent_run(
    *,
    fz11: Path = independent.DEFAULT_FZ11,
    fz12: Path = independent.DEFAULT_FZ12,
    lean: Path = independent.DEFAULT_LEAN,
    receipt: Path = independent.DEFAULT_RECEIPT,
) -> subprocess.CompletedProcess[str]:
    environment = os.environ.copy()
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    return subprocess.run(
        [
            sys.executable,
            str(VERIFIER),
            "--fz11",
            str(fz11),
            "--fz12",
            str(fz12),
            "--lean",
            str(lean),
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
    receipt = copy.deepcopy(certificate.build_receipt())
    mutate(receipt)
    resign(receipt)
    path = tmp_path / "carrier_frequency_speed_receipt.json"
    write_canonical(path, receipt)
    return path


def test_exact_supports_are_reconstructed_independently() -> None:
    producer = certificate.support_certificate()
    verifier = independent.independent_support_checks()
    assert producer == verifier
    assert producer["vertex12"] == {
        "support_count": 12,
        "raw_radius_squared": "5/2+1/2*sqrt5",
        "unit_second_moment": "sum_i (u_i.dot.x)^2 = 4 |x|^2",
        "tight_constant": "4",
        "normalized_symbol_prefactor": "1/(2 a^2)",
    }
    assert producer["edge30"] == {
        "support_count": 30,
        "unit_norm_squared": "1",
        "unit_second_moment": "sum_e (w_e.dot.x)^2 = 10 |x|^2",
        "tight_constant": "10",
        "normalized_symbol_prefactor": "1/(5 a^2)",
    }


def test_global_unit_upper_bound_is_exact() -> None:
    receipt = certificate.build_receipt()
    assert receipt["generic_exact_theorem"] == independent.expected_generic_theorem()
    assert receipt["generic_exact_theorem"]["frequency_bound"] == (
        "|Omega_a(k)-Omega_a(p)| <= |k-p|"
    )
    assert receipt["generic_exact_theorem"]["certified_upper_constant"] == "1"
    assert receipt["contract_mutation_controls"] == (
        independent.expected_contract_mutation_controls()
    )


def test_physical_and_exposure_bridges_remain_open() -> None:
    receipt = certificate.build_receipt()
    assert receipt["physical_boundary"] == independent.expected_physical_boundary()
    assert all(value is False for value in receipt["physical_boundary"].values())
    assert receipt["exposure_boundary"] == independent.expected_exposure_boundary()
    assert receipt["exposure_boundary"]["comparison_inputs"] == []
    assert receipt["branch_bindings"]["new_prediction_payload"] is False


def test_lean_surface_is_fixed_and_placeholder_free() -> None:
    source = independent.verify_lean(independent.DEFAULT_LEAN).decode("utf-8")
    assert "theorem feature_dist_le" in source
    assert "theorem frequency_global_one_lipschitz" in source
    assert "theorem fz11_frequency_global_one_lipschitz" in source
    assert "theorem fz12_frequency_global_one_lipschitz" in source
    assert "sorry" not in source
    assert "admit" not in source
    assert certificate.build_receipt()["lean"] == independent.expected_lean_pin()


def test_committed_receipt_is_canonical_and_byte_exact() -> None:
    expected = certificate.build_receipt()
    assert independent.DEFAULT_RECEIPT.read_bytes() == independent.canonical_json_bytes(
        expected
    )
    parsed = independent.strict_json_loads(
        independent.DEFAULT_RECEIPT.read_bytes(), "committed receipt"
    )
    independent.verify_payload(parsed)


def test_independent_verifier_has_no_producer_import_and_passes() -> None:
    source = VERIFIER.read_text(encoding="utf-8")
    assert "import carrier_frequency_speed_certificate" not in source
    assert "from carrier_frequency_speed_certificate" not in source
    result = independent_run()
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "CARRIER_FREQUENCY_SPEED_INDEPENDENT_PASS"


@pytest.mark.parametrize(
    ("mutation", "expected_error"),
    [
        ("negative_weight", "generic theorem contract drift"),
        ("edge_tight", "support instantiation drift"),
        ("vertex_tight", "support instantiation drift"),
        ("lipschitz", "generic theorem contract drift"),
        ("physical_frequency", "physical boundary drift"),
        ("physical_front", "physical boundary drift"),
        ("comparison", "exposure boundary drift"),
        ("parent_hash", "parent pins drift"),
        ("control_registry", "contract-mutation registry drift"),
        ("frozen_parent", "branch binding drift"),
        ("lean_theorem", "Lean pin drift"),
        ("extra_key", "top-level keys drift"),
    ],
)
def test_resigned_semantic_mutations_fail_closed(
    tmp_path: Path, mutation: str, expected_error: str
) -> None:
    def mutate(receipt: dict) -> None:
        if mutation == "negative_weight":
            receipt["generic_exact_theorem"]["premises"][0] = (
                "a finite support with signed weights"
            )
        elif mutation == "edge_tight":
            receipt["exact_support_instantiations"]["edge30"]["tight_constant"] = "9"
        elif mutation == "vertex_tight":
            receipt["exact_support_instantiations"]["vertex12"]["tight_constant"] = "5"
        elif mutation == "lipschitz":
            receipt["generic_exact_theorem"]["certified_upper_constant"] = "2"
        elif mutation == "physical_frequency":
            receipt["physical_boundary"]["physical_frequency_proved"] = True
        elif mutation == "physical_front":
            receipt["physical_boundary"]["wave_packet_or_signal_front_proved"] = True
        elif mutation == "comparison":
            receipt["exposure_boundary"]["public_measurement_read"] = True
        elif mutation == "parent_hash":
            receipt["parent_pins"][0]["sha256"] = "sha256:" + "0" * 64
        elif mutation == "control_registry":
            receipt["contract_mutation_controls"]["controls"][0] = (
                "allow an unspecified weight mutation"
            )
        elif mutation == "frozen_parent":
            receipt["branch_bindings"]["frozen_bytes_modified"] = True
        elif mutation == "lean_theorem":
            receipt["lean"]["theorems"].append("physical_group_speed")
        elif mutation == "extra_key":
            receipt["audit_mutation"] = True
        else:  # pragma: no cover
            raise AssertionError(mutation)

    path = mutated_receipt(tmp_path, mutate)
    payload = independent.strict_json_loads(path.read_bytes(), "mutated receipt")
    with pytest.raises(independent.VerificationError, match=expected_error):
        independent.verify_payload(payload, check_files=False)


def test_unresigned_receipt_mutation_fails_self_digest() -> None:
    receipt = certificate.build_receipt()
    receipt["issue"] = 705
    with pytest.raises(independent.VerificationError, match="self-digest drift"):
        independent.verify_payload(receipt, check_files=False)


@pytest.mark.parametrize("parent_name", ["fz11", "fz12"])
def test_resigned_parent_mutation_fails_fixed_raw_pin(
    tmp_path: Path, parent_name: str
) -> None:
    source = (
        independent.DEFAULT_FZ11 if parent_name == "fz11" else independent.DEFAULT_FZ12
    )
    parent = json.loads(source.read_text(encoding="ascii"))
    parent["audit_mutation"] = True
    parent["receipt_sha256"] = independent.self_digest(parent, "receipt_sha256")
    path = tmp_path / source.name
    write_canonical(path, parent)
    result = independent_run(**{parent_name: path})
    assert result.returncode == 1
    assert (
        "parent byte-count drift" in result.stderr
        or "parent raw hash drift" in result.stderr
    )


def test_lean_text_mutation_fails_fixed_raw_pin(tmp_path: Path) -> None:
    source = independent.DEFAULT_LEAN.read_bytes()
    path = tmp_path / "CarrierFrequencySpeed.lean"
    path.write_bytes(
        source.replace(b"globally one-Lipschitz", b"globally two-Lipschitz", 1)
    )
    result = independent_run(lean=path)
    assert result.returncode == 1
    assert (
        "Lean byte-count drift" in result.stderr
        or "Lean raw hash drift" in result.stderr
    )


def test_duplicate_keys_and_nonfinite_numbers_are_rejected() -> None:
    with pytest.raises(independent.VerificationError, match="duplicate JSON key"):
        independent.strict_json_loads(b'{"schema":"x","schema":"y"}\n', "duplicate")
    with pytest.raises(independent.VerificationError, match="non-finite JSON constant"):
        independent.strict_json_loads(b'{"value":NaN}\n', "nonfinite")


@pytest.mark.parametrize(
    ("mutation", "expected_error"),
    [
        ("status", "status drift"),
        ("self_digest", "self-digest drift"),
        ("duplicate", "duplicate JSON key"),
        ("nonfinite", "non-finite JSON constant"),
    ],
)
def test_producer_parent_loader_is_typed_and_fail_closed(
    tmp_path: Path, mutation: str, expected_error: str
) -> None:
    contract = dict(certificate.PARENTS[0])
    parent = json.loads(contract["path"].read_text(encoding="ascii"))
    if mutation == "status":
        parent["status"] = "MUTATED_STATUS"
        parent["receipt_sha256"] = independent.self_digest(parent, "receipt_sha256")
        raw = certificate.canonical_json_bytes(parent)
    elif mutation == "self_digest":
        parent["audit_mutation"] = True
        raw = certificate.canonical_json_bytes(parent)
    elif mutation == "duplicate":
        raw = b'{"schema":"x","schema":"y"}\n'
    elif mutation == "nonfinite":
        raw = b'{"value":NaN}\n'
    else:  # pragma: no cover
        raise AssertionError(mutation)
    path = tmp_path / "mutated_parent.json"
    path.write_bytes(raw)
    contract["path"] = path
    contract["bytes"] = len(raw)
    contract["sha256"] = hashlib.sha256(raw).hexdigest()
    with pytest.raises(RuntimeError, match=expected_error):
        certificate.load_parent(contract)
