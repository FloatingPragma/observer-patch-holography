from __future__ import annotations

import copy
import hashlib
import importlib.util
import io
import json
import stat
import sys
from pathlib import Path

import pytest


pytest.importorskip("qiskit")

PROGRAMS = Path(__file__).resolve().parents[1] / "programs"
sys.path.insert(0, str(PROGRAMS))
MODULE_PATH = PROGRAMS / "blind_preregister.py"
SPEC = importlib.util.spec_from_file_location("blind_preregister", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
blind = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = blind
SPEC.loader.exec_module(blind)


def _transport_qpy_sha256(circuit) -> str:
    """Test-only QPY digest demonstrating why it is not a logical binding."""

    from qiskit import qpy

    buffer = io.BytesIO()
    qpy.dump(circuit, buffer)
    return hashlib.sha256(buffer.getvalue()).hexdigest()


@pytest.fixture(scope="module")
def sealed_bundle():
    scope = {
        "cayley_models": ("z5",),
        "cayley_protocols": blind.CAYLEY_PROTOCOLS,
        "mh_spectra": ("z3_cyclic_control",),
        "calibration_codes": (0, 15),
    }
    return blind.build_blinded_preregistration(
        test_seed="blind-preregister-tests",
        _test_scope=scope,
    )


def test_secret_is_32_random_bytes_and_test_mode_is_deterministic() -> None:
    first, first_mode = blind._new_secret(None)
    second, second_mode = blind._new_secret(None)
    assert len(first) == len(second) == 32
    assert first != second
    assert first_mode == second_mode == "production_random"

    deterministic_1, mode_1 = blind._new_secret("fixed")
    deterministic_2, mode_2 = blind._new_secret("fixed")
    assert deterministic_1 == deterministic_2
    assert mode_1 == mode_2 == "deterministic_test_only"


def test_deterministic_bundle_and_acyclic_digest_graph(sealed_bundle) -> None:
    public, reveal = sealed_bundle
    scope = {
        "cayley_models": ("z5",),
        "cayley_protocols": blind.CAYLEY_PROTOCOLS,
        "mh_spectra": ("z3_cyclic_control",),
        "calibration_codes": (0, 15),
    }
    public_again, reveal_again = blind.build_blinded_preregistration(
        test_seed="blind-preregister-tests",
        _test_scope=scope,
    )
    assert public == public_again
    assert reveal == reveal_again
    assert blind.catalog_precommitment(public) == public["catalog_precommitment_sha256"]
    assert (
        blind.manifest_core_hash(public)
        == reveal["sealed_payload"]["public_manifest_core_sha256"]
    )
    secret = bytes.fromhex(reveal["secret_hex"])
    assert (
        blind.secret_commitment(secret, reveal["sealed_payload"])
        == public["secret_commitment_sha256"]
    )
    assert blind.manifest_hash(public) == public["manifest_sha256"]
    assert blind.verify_bundle(public, reveal)


def test_public_manifest_is_opaque_and_private_mapping_is_complete(sealed_bundle) -> None:
    public, reveal = sealed_bundle
    rendered = blind.canonical_json_bytes(public).decode("utf-8")
    for forbidden in (
        "ibm_fez",
        "ibm_kingston",
        "development",
        "heldout",
        "record_gated",
        "open_loop_heat",
        "delayed_record",
        "shuffled_record",
        "inverted_record",
        "z5_cyclic_control",
    ):
        assert forbidden not in rendered
    mappings = reveal["sealed_payload"]["mappings"]
    assert set(mappings) == {
        "model",
        "protocol",
        "state",
        "slots",
        "encoding",
        "layout",
        "backend",
        "backend_role",
    }
    assert set(mappings["protocol"]) == set(blind.CAYLEY_PROTOCOLS)
    assert mappings["state"]["z5"]
    assert mappings["slots"]["z5"]["disturbance"]
    assert mappings["encoding"]
    assert len(public["backend_slots"]) == 2
    assert all(set(slot) == {"role", "backend", "layout"} for slot in public["backend_slots"])


def test_catalog_covers_five_arms_mh_and_diagnostic_calibration(sealed_bundle) -> None:
    public, reveal = sealed_bundle
    descriptors = reveal["sealed_payload"]["circuits"]
    cayley_protocols = {
        descriptor["parameters"]["protocol"]
        for descriptor in descriptors.values()
        if descriptor["family"] == "cayley"
    }
    assert cayley_protocols == set(blind.CAYLEY_PROTOCOLS)
    assert any(descriptor["family"] == "mh" for descriptor in descriptors.values())
    calibrations = [
        descriptor
        for descriptor in descriptors.values()
        if descriptor["family"] == "readout_calibration"
    ]
    assert calibrations
    assert all(
        descriptor["parameters"]["evidentiary_role"] == "diagnostic_only"
        and descriptor["shots"] == 512
        for descriptor in calibrations
    )
    assert all(
        descriptor["shots"] == 192
        for descriptor in descriptors.values()
        if descriptor["family"] in {"cayley", "mh"}
    )
    resources = public["resources"]
    assert resources["transpile_optimization_level"] == 0
    assert resources["account_quota_seconds"] == 600
    assert resources["account_reserve_seconds"] == 180
    assert resources["estimated_qpu_seconds_ceiling"] == 420
    assert resources["total_circuit_instances"] == len(public["circuits"])
    assert len(public["nulls"]) == len(blind.DEFAULT_NULL_SPECS)


def test_two_phase_analysis_binding_preserves_secret_catalog_and_circuits(sealed_bundle) -> None:
    public, reveal = sealed_bundle
    hardened = copy.deepcopy(reveal["sealed_payload"]["analysis_document"])
    hardened.pop("analysis_lock_sha256")
    hardened["lock_revision"] = "synthetic-hardened-test"
    rebound_public, rebound_reveal = blind.bind_analysis_document(public, reveal, hardened)
    assert rebound_public["catalog_precommitment_sha256"] == public[
        "catalog_precommitment_sha256"
    ]
    assert rebound_public["circuits"] == public["circuits"]
    assert rebound_reveal["secret_hex"] == reveal["secret_hex"]
    assert rebound_reveal["sealed_payload"]["circuits"] == reveal["sealed_payload"][
        "circuits"
    ]
    assert rebound_public["analysis_lock_sha256"] != public["analysis_lock_sha256"]
    assert rebound_public["secret_commitment_sha256"] != public[
        "secret_commitment_sha256"
    ]
    assert blind.verify_bundle(rebound_public, rebound_reveal, hardened)


def test_every_arm_mh_and_calibration_has_reproducible_logical_hash(
    sealed_bundle,
) -> None:
    public, reveal = sealed_bundle
    descriptors = reveal["sealed_payload"]["circuits"]
    observed_protocols: set[str] = set()
    observed_families: set[str] = set()
    for opaque_id, descriptor in descriptors.items():
        circuit = blind.rebuild_blinded_circuit(opaque_id, descriptor)
        serialization = blind.canonical_openqasm3_bytes(circuit)
        assert serialization.startswith(b"OPENQASM 3")
        assert serialization.endswith(b"\n")
        assert not serialization.endswith(b"\n\n")
        assert blind.logical_circuit_sha256(circuit) == descriptor[
            "logical_circuit_sha256"
        ]
        observed_families.add(descriptor["family"])
        if descriptor["family"] == "cayley":
            observed_protocols.add(descriptor["parameters"]["protocol"])
    assert observed_protocols == set(blind.CAYLEY_PROTOCOLS)
    assert observed_families == {"cayley", "mh", "readout_calibration"}
    assert all(
        set(row)
        == {"opaque_id", "logical_circuit_sha256", "shots", "backend_role"}
        for row in public["circuits"]
    )
    assert blind.verify_bundle(public, reveal, rebuild_circuits=True)


def test_equal_circuits_with_divergent_qpy_have_one_logical_digest() -> None:
    from qiskit import QuantumCircuit

    first = QuantumCircuit(1, 1, name="same")
    first.metadata = {"alpha": 1, "beta": 2}
    first.h(0)
    first.measure(0, 0)
    second = QuantumCircuit(1, 1, name="same")
    second.metadata = {"beta": 2, "alpha": 1}
    second.h(0)
    second.measure(0, 0)

    assert first == second
    assert _transport_qpy_sha256(first) != _transport_qpy_sha256(second)
    assert blind.canonical_openqasm3_bytes(first) == blind.canonical_openqasm3_bytes(
        second
    )
    assert blind.logical_circuit_sha256(first) == blind.logical_circuit_sha256(second)


def test_logical_serialization_normalization_contract(monkeypatch) -> None:
    from qiskit import qasm3

    monkeypatch.setattr(
        qasm3,
        "dumps",
        lambda _circuit: "OPENQASM 3.0;  \r\nqubit q;\t\r\n\r\n",
    )
    assert blind.canonical_openqasm3_bytes(object()) == b"OPENQASM 3.0;\nqubit q;\n"


def test_catalog_order_does_not_change_logical_bindings() -> None:
    scope = {
        "cayley_models": ("z5",),
        "cayley_protocols": tuple(reversed(blind.CAYLEY_PROTOCOLS)),
        "mh_spectra": ("z3_cyclic_control",),
        "calibration_codes": (15, 0),
    }
    public, _ = blind.build_blinded_preregistration(
        test_seed="blind-preregister-tests",
        _test_scope=scope,
    )
    reversed_hashes = {
        row["opaque_id"]: row["logical_circuit_sha256"] for row in public["circuits"]
    }
    canonical_public, _ = blind.build_blinded_preregistration(
        test_seed="blind-preregister-tests",
        _test_scope={
            "cayley_models": ("z5",),
            "cayley_protocols": blind.CAYLEY_PROTOCOLS,
            "mh_spectra": ("z3_cyclic_control",),
            "calibration_codes": (0, 15),
        },
    )
    canonical_hashes = {
        row["opaque_id"]: row["logical_circuit_sha256"]
        for row in canonical_public["circuits"]
    }
    assert reversed_hashes == canonical_hashes


@pytest.mark.parametrize("target", ["public", "reveal", "analysis"])
def test_tampering_fails_closed(sealed_bundle, target: str) -> None:
    public, reveal = copy.deepcopy(sealed_bundle)
    if target == "public":
        public["resources"]["total_planned_shots"] += 1
    elif target == "reveal":
        descriptor = next(iter(reveal["sealed_payload"]["circuits"].values()))
        descriptor["shots"] += 1
    else:
        reveal["sealed_payload"]["analysis_document"]["comparison_measure"] = "changed"
    with pytest.raises(blind.VerificationError):
        blind.verify_bundle(public, reveal)


def test_reveal_is_outside_checkout_mode_0600_and_exclusive(
    sealed_bundle,
    tmp_path: Path,
) -> None:
    _, reveal = sealed_bundle
    checkout = tmp_path / "checkout"
    checkout.mkdir()
    outside = tmp_path / "private" / "reveal.json"
    written = blind.write_reveal_file(outside, reveal, checkout_root=checkout)
    assert stat.S_IMODE(written.stat().st_mode) == 0o600
    assert json.loads(written.read_text(encoding="utf-8")) == reveal
    with pytest.raises(FileExistsError):
        blind.write_reveal_file(outside, reveal, checkout_root=checkout)
    with pytest.raises(ValueError, match="outside"):
        blind.write_reveal_file(
            checkout / "private-reveal.json",
            reveal,
            checkout_root=checkout,
        )
