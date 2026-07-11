#!/usr/bin/env python3
"""Commit/reveal preregistration for the blinded record-gated benchmark.

This module never connects to IBM and never submits work.  It builds exact
logical circuits, replaces semantic names and metadata with opaque values,
commits the private reveal, and verifies the complete digest graph.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import hmac
import json
import os
import secrets
import stat
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence


SCHEMA_VERSION = "oph.blinded-record-gated-preregistration.v1"
REVEAL_SCHEMA_VERSION = "oph.blinded-record-gated-reveal.v1"
BLINDED_CIRCUIT_SCHEMA_VERSION = "oph.blinded-cayley-logical-circuit.v1"
COMMITMENT_DOMAIN = b"oph-509-blind-v1\0"
OPAQUE_DOMAIN = b"oph-509-opaque-id-v1\0"
SECRET_BYTES = 32
SHA256_HEX_LENGTH = 64
LOGICAL_CIRCUIT_SERIALIZATION = "normalized-openqasm3-utf8-v1"

CAYLEY_PROTOCOLS = (
    "record_gated",
    "open_loop_heat",
    "delayed_record",
    "shuffled_record",
    "inverted_record",
)
MH_SPECTRA = (
    "s3_primary",
    "z3_cyclic_control",
    "z5_cyclic_control",
    "a4_nonabelian_decoy",
    "random_seeded_decoy",
)


@dataclass(frozen=True)
class BackendLayoutSlot:
    role: str
    backend: str
    layout: tuple[int, int, int, int]
    properties_last_update: str

    def validate(self) -> None:
        if not self.role or not self.backend or not self.properties_last_update:
            raise ValueError("backend slot text fields must be nonempty")
        if len(self.layout) != 4 or len(set(self.layout)) != 4:
            raise ValueError("layout must contain four distinct physical qubits")
        if any(not isinstance(qubit, int) or qubit < 0 for qubit in self.layout):
            raise ValueError("layout qubits must be nonnegative integers")


DEFAULT_BACKEND_SLOTS = (
    BackendLayoutSlot(
        role="development",
        backend="ibm_fez",
        layout=(10, 18, 94, 124),
        properties_last_update="2026-07-11T10:15:57+07:00",
    ),
    BackendLayoutSlot(
        role="heldout",
        backend="ibm_kingston",
        layout=(21, 47, 50, 107),
        properties_last_update="2026-07-11T09:49:05+07:00",
    ),
)

DEFAULT_THRESHOLDS: dict[str, Any] = {
    "per_backend_layout_min_conditional_likelihood_ratio_strict": 10.0,
    "pooled_min_conditional_likelihood_ratio_strict": 100.0,
    "null_favor_failure_conditional_likelihood_ratio_inclusive": 100.0,
    "failure_backend_layout_slots_required": 2,
    "simultaneous_prediction_envelope": 0.99,
    "secondary_familywise_alpha": 0.01,
    "mapping_unique_winner_count": 1,
    "required_backend_layout_slots": 2,
    "max_unlisted_exclusions": 0,
    "max_dropped_shots": 0,
    "max_post_reveal_layout_changes": 0,
}

DEFAULT_EXCLUSIONS: dict[str, Any] = {
    "max_unlisted_exclusions": 0,
    "max_dropped_shots": 0,
    "target_dependent_exclusions_allowed": False,
    "post_reveal_layout_selection_allowed": False,
    "invalid_codes": "retain_as_explicit_outcomes",
    "incorrect_readbacks": "retain_as_explicit_outcomes",
    "failed_feedback": "retain_as_explicit_outcomes",
    "failed_submissions": "retain_receipt_and_do_not_replace_from_outcomes",
    "layout_omissions": "invalidate_affected_slot",
    "duplicate_transport_receipts": "deduplicate_only_by_provider_job_id",
}

DEFAULT_NULL_SPECS: dict[str, dict[str, Any]] = {
    "open_loop_heat": {"family": "cayley", "rule": "L_G squared"},
    "delayed_record": {
        "family": "cayley",
        "rule": "gate the second proposal with descent evaluated at the pre-disturbance state",
    },
    "shuffled_record": {
        "family": "cayley",
        "rule": "gate with descent evaluated for (proposal_generator + 1) modulo degree",
    },
    "inverted_record": {
        "family": "cayley",
        "rule": "gate the second proposal with one minus the contemporaneous descent record",
    },
    "state_preparation_only": {
        "family": "likelihood_only",
        "rule": "fit final marginals without a joint heated-record/decision process law",
    },
    "label_layout_only": {
        "family": "global_mapping_component",
        "rule": "frozen shared label and layout mixture, multiplicity corrected",
    },
    "calibration_convolved_noise": {
        "family": "likelihood_only",
        "rule": "frozen robust contamination mixture; basis calibration is diagnostic only",
    },
    "unweighted_sector_boltzmann": {"family": "mh", "kappa": 0.0},
    "full_block_plancherel": {"family": "mh", "kappa": 2.0},
}


class VerificationError(ValueError):
    """Raised when a manifest, reveal, lock, or circuit fails closed."""


def canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")


def sha256_json(value: Any) -> str:
    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()


def _is_sha256(value: Any) -> bool:
    if not isinstance(value, str) or len(value) != SHA256_HEX_LENGTH:
        return False
    return all(character in "0123456789abcdef" for character in value)


def _derive_opaque_key(secret: bytes) -> bytes:
    return hmac.new(secret, b"opaque-id-key", hashlib.sha256).digest()


def _opaque_id(key: bytes, namespace: str, semantic_value: Any) -> str:
    material = (
        OPAQUE_DOMAIN
        + namespace.encode("utf-8")
        + b"\0"
        + canonical_json_bytes(semantic_value)
    )
    return "o_" + hmac.new(key, material, hashlib.sha256).hexdigest()[:32]


def _blind_nonce(key: bytes, opaque_id: str) -> str:
    return hmac.new(key, b"circuit-nonce\0" + opaque_id.encode(), hashlib.sha256).hexdigest()


def _keyed_permutation(key: bytes, namespace: str, values: Sequence[int]) -> tuple[int, ...]:
    return tuple(
        sorted(
            values,
            key=lambda value: hmac.new(
                key,
                namespace.encode() + b"\0" + str(value).encode(),
                hashlib.sha256,
            ).digest(),
        )
    )


def _new_secret(test_seed: str | None) -> tuple[bytes, str]:
    if test_seed is None:
        return secrets.token_bytes(SECRET_BYTES), "production_random"
    if not test_seed:
        raise ValueError("deterministic test seed must be nonempty")
    return (
        hashlib.sha256(b"oph-509-deterministic-test-only\0" + test_seed.encode()).digest(),
        "deterministic_test_only",
    )


def catalog_precommitment(manifest: Mapping[str, Any]) -> str:
    omitted = {
        "catalog_precommitment_sha256",
        "analysis_lock_sha256",
        "secret_commitment_sha256",
        "manifest_sha256",
    }
    return sha256_json({key: value for key, value in manifest.items() if key not in omitted})


def manifest_core_hash(manifest: Mapping[str, Any]) -> str:
    return sha256_json(
        {
            key: value
            for key, value in manifest.items()
            if key not in {"manifest_sha256", "secret_commitment_sha256"}
        }
    )


def manifest_hash(manifest: Mapping[str, Any]) -> str:
    return sha256_json({key: value for key, value in manifest.items() if key != "manifest_sha256"})


def secret_commitment(secret: bytes, sealed_payload: Mapping[str, Any]) -> str:
    if len(secret) != SECRET_BYTES:
        raise ValueError("commitment secret must contain exactly 32 bytes")
    return hashlib.sha256(
        COMMITMENT_DOMAIN + secret + b"\0" + canonical_json_bytes(sealed_payload)
    ).hexdigest()


def _load_circuit_modules() -> tuple[Any, Any, Any]:
    try:
        import generative_repair_circuits as mh_circuits
        import generative_repair_kernel as kernels
        import record_gated_cayley_circuits as cayley_circuits
    except ImportError as exc:  # pragma: no cover - exercised only without optional Qiskit deps
        raise RuntimeError(
            "Qiskit benchmark dependencies are required to build or rebuild sealed circuits"
        ) from exc
    return kernels, cayley_circuits, mh_circuits


def _cayley_semantics(parameters: Mapping[str, Any]) -> dict[str, Any]:
    kernels, _, _ = _load_circuit_modules()
    model_key = str(parameters["model"])
    protocol = str(parameters["protocol"])
    if protocol not in CAYLEY_PROTOCOLS:
        raise ValueError(f"unsupported Cayley protocol {protocol!r}")
    model = kernels.builtin_cayley_models()[model_key]
    initial = int(parameters["initial_state"])
    disturbance_slot = int(parameters["disturbance_slot"])
    second_slot = int(parameters["second_slot"])
    encoding = tuple(int(value) for value in parameters["state_encoding"])
    degree = len(model.generators)
    if initial not in range(model.size):
        raise ValueError("Cayley initial state is outside the model")
    if disturbance_slot not in range(2 * degree) or second_slot not in range(2 * degree):
        raise ValueError("Cayley balanced slot is outside the model")
    if len(encoding) != model.size or len(set(encoding)) != model.size:
        raise ValueError("Cayley encoding must be injective and complete")

    disturbance_generator = (
        None if disturbance_slot < degree else disturbance_slot - degree
    )
    heated = (
        initial
        if disturbance_generator is None
        else model.right_action[disturbance_generator][initial]
    )
    proposal_generator = second_slot % degree
    proposal_target = model.right_action[proposal_generator][heated]
    current_decision = int(model.mismatch[proposal_target] < model.mismatch[heated])

    if protocol == "open_loop_heat":
        open_generator = None if second_slot < degree else second_slot - degree
        final = heated if open_generator is None else model.right_action[open_generator][heated]
        decision = current_decision
    else:
        decision_generator = proposal_generator
        decision_source = heated
        if protocol == "delayed_record":
            decision_source = initial
        elif protocol == "shuffled_record":
            decision_generator = (proposal_generator + 1) % degree
        decision_target = model.right_action[decision_generator][decision_source]
        decision = int(model.mismatch[decision_target] < model.mismatch[decision_source])
        if protocol == "inverted_record":
            decision = 1 - decision
        final = proposal_target if decision else heated

    return {
        "model": model,
        "protocol": protocol,
        "initial_state": initial,
        "disturbance_slot": disturbance_slot,
        "second_slot": second_slot,
        "state_encoding": encoding,
        "heated_state": heated,
        "decision_record": decision,
        "final_state": final,
    }


def _apply_code(circuit: Any, register: Any, code: int, width: int) -> None:
    for bit in range(width):
        if (code >> bit) & 1:
            circuit.x(register[bit])
        else:
            circuit.id(register[bit])


def _apply_xor(circuit: Any, register: Any, source: int, target: int, width: int) -> None:
    mask = source ^ target
    for bit in range(width):
        if (mask >> bit) & 1:
            circuit.x(register[bit])
        else:
            circuit.id(register[bit])


def _build_cayley_circuit(parameters: Mapping[str, Any]) -> Any:
    kernels, cayley_circuits, _ = _load_circuit_modules()
    protocol = str(parameters["protocol"])
    if protocol not in CAYLEY_PROTOCOLS:
        raise ValueError(f"unsupported Cayley protocol {protocol!r}")
    model = kernels.builtin_cayley_models()[str(parameters["model"])]
    recipe = cayley_circuits.build_recipe(
        model,
        protocol,
        int(parameters["initial_state"]),
        int(parameters["disturbance_slot"]),
        int(parameters["second_slot"]),
        tuple(int(value) for value in parameters["state_encoding"]),
    )
    return cayley_circuits.build_circuit(recipe)


def _build_mh_circuit(parameters: Mapping[str, Any]) -> Any:
    kernels, _, mh_circuits = _load_circuit_modules()
    spectrum = kernels.builtin_spectra()[str(parameters["spectrum"])]
    recipe = mh_circuits.edge_recipe(
        spectrum=spectrum,
        beta=float(parameters["beta"]),
        programmed_kappa=float(parameters["programmed_kappa"]),
        semantic_source=int(parameters["semantic_source"]),
        semantic_target=int(parameters["semantic_target"]),
        label_permutation=tuple(int(value) for value in parameters["label_permutation"]),
    )
    return mh_circuits.build_edge_instrument(recipe)


def _build_readout_calibration_circuit(parameters: Mapping[str, Any]) -> Any:
    from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister

    if parameters.get("evidentiary_role") != "diagnostic_only":
        raise ValueError("basis readout calibration must be labeled diagnostic_only")
    width = int(parameters["num_qubits"])
    basis_code = int(parameters["basis_code"])
    if width != 4 or basis_code not in range(2**width):
        raise ValueError("readout calibration requires a four-qubit basis code")
    qubits = QuantumRegister(width, "calibration_q")
    bits = ClassicalRegister(width, "calibration_c")
    circuit = QuantumCircuit(qubits, bits)
    _apply_code(circuit, qubits, basis_code, width)
    circuit.measure(qubits, bits)
    return circuit


def rebuild_blinded_circuit(opaque_id: str, descriptor: Mapping[str, Any]) -> Any:
    """Rebuild an exact logical circuit from a verified private descriptor."""

    required = {
        "family",
        "backend_role",
        "backend_role_opaque_id",
        "shots",
        "logical_circuit_sha256",
        "circuit_name",
        "circuit_metadata",
        "parameters",
    }
    if set(descriptor) != required:
        raise VerificationError(
            f"circuit descriptor keys differ from the frozen schema: {sorted(descriptor)}"
        )
    if descriptor["circuit_name"] != opaque_id:
        raise VerificationError("descriptor circuit name does not match its opaque key")
    family = descriptor["family"]
    if family == "cayley":
        circuit = _build_cayley_circuit(descriptor["parameters"])
    elif family == "mh":
        circuit = _build_mh_circuit(descriptor["parameters"])
    elif family == "readout_calibration":
        circuit = _build_readout_calibration_circuit(descriptor["parameters"])
    else:
        raise VerificationError(f"unsupported blinded circuit family {family!r}")
    circuit.name = opaque_id
    circuit.metadata = dict(descriptor["circuit_metadata"])
    return circuit


def canonical_openqasm3_bytes(circuit: Any) -> bytes:
    """Return the stable logical-circuit serialization used by the seal.

    QPY is a transport format and can encode semantically irrelevant Python
    insertion order (for example, circuit-metadata dictionary order).  It is
    therefore unsuitable as a cross-object or cross-process commitment.  The
    OpenQASM 3 exporter describes the executable logical circuit and omits that
    incidental object state.  Normalization here makes the byte contract
    explicit: LF line endings, no trailing horizontal whitespace, no terminal
    blank lines, exactly one final LF, and UTF-8 encoding.
    """

    from qiskit import qasm3

    exported = qasm3.dumps(circuit).replace("\r\n", "\n").replace("\r", "\n")
    lines = [line.rstrip(" \t") for line in exported.split("\n")]
    while lines and lines[-1] == "":
        lines.pop()
    if not lines:
        raise VerificationError("OpenQASM 3 exporter returned an empty logical circuit")
    return ("\n".join(lines) + "\n").encode("utf-8")


def logical_circuit_sha256(circuit: Any) -> str:
    """Hash the canonical normalized OpenQASM 3 logical circuit."""

    return hashlib.sha256(canonical_openqasm3_bytes(circuit)).hexdigest()


def verify_rebuilt_circuit(opaque_id: str, descriptor: Mapping[str, Any]) -> None:
    observed = logical_circuit_sha256(rebuild_blinded_circuit(opaque_id, descriptor))
    if observed != descriptor["logical_circuit_sha256"]:
        raise VerificationError(
            "logical OpenQASM 3 digest mismatch for "
            f"{opaque_id}: {observed} != {descriptor['logical_circuit_sha256']}"
        )


def _validate_backend_slots(slots: Sequence[BackendLayoutSlot]) -> None:
    if len(slots) != 2:
        raise ValueError("exactly two backend/layout slots are required")
    for slot in slots:
        slot.validate()
    if len({slot.role for slot in slots}) != 2:
        raise ValueError("backend roles must be distinct")
    if len({slot.backend for slot in slots}) != 2:
        raise ValueError("independent slots must use distinct backends")
    if len({slot.layout for slot in slots}) != 2:
        raise ValueError("independent slots must use distinct layouts")


def _descriptor(
    *,
    family: str,
    role: str,
    role_opaque_id: str,
    shots: int,
    parameters: Mapping[str, Any],
    opaque_key: bytes,
) -> tuple[str, dict[str, Any]]:
    semantic_key = {"family": family, "backend_role": role, "parameters": parameters}
    opaque_id = _opaque_id(opaque_key, "circuit", semantic_key)
    metadata = {
        "schema_version": BLINDED_CIRCUIT_SCHEMA_VERSION,
        "opaque_id": opaque_id,
        "blind_nonce": _blind_nonce(opaque_key, opaque_id),
    }
    descriptor: dict[str, Any] = {
        "family": family,
        "backend_role": role,
        "backend_role_opaque_id": role_opaque_id,
        "shots": int(shots),
        "logical_circuit_sha256": "0" * SHA256_HEX_LENGTH,
        "circuit_name": opaque_id,
        "circuit_metadata": metadata,
        "parameters": dict(parameters),
    }
    digest = logical_circuit_sha256(rebuild_blinded_circuit(opaque_id, descriptor))
    descriptor["logical_circuit_sha256"] = digest
    return opaque_id, descriptor


def _public_circuit_entry(opaque_id: str, descriptor: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "opaque_id": opaque_id,
        "logical_circuit_sha256": descriptor["logical_circuit_sha256"],
        "shots": descriptor["shots"],
        "backend_role": descriptor["backend_role_opaque_id"],
    }


def _default_analysis_document(catalog_sha256: str) -> dict[str, Any]:
    return {
        "schema_version": "oph.blinded-record-gated-analysis-lock.v1",
        "catalog_precommitment_sha256": catalog_sha256,
        "primary_endpoint": (
            "pooled heldout joint heated-record, decision, final-state, and leakage likelihood"
        ),
        "comparison_measure": "conditional likelihood ratio with no posterior-odds claim",
        "outcome_policy": "retain every shot and every invalid code as an explicit outcome",
        "calibration_model": {
            "control_rule": (
                "diagnostic basis controls are reported from later count/job/hash-bound data; "
                "they are not a factorized causal channel or a fictional preregistration pass"
            ),
        },
        "multiplicity": "frozen global label/layout components with one shared correction",
        "secondary_correction": "Holm family-wise correction at alpha 0.01",
        "thresholds": dict(DEFAULT_THRESHOLDS),
        "null_specs_sha256": sha256_json(DEFAULT_NULL_SPECS),
    }


def _lock_analysis_document(
    analysis_document: Mapping[str, Any] | None,
    catalog_sha256: str,
) -> tuple[dict[str, Any], str]:
    document = dict(
        _default_analysis_document(catalog_sha256)
        if analysis_document is None
        else analysis_document
    )
    supplied_catalog = document.get("catalog_precommitment_sha256")
    if supplied_catalog is not None and supplied_catalog != catalog_sha256:
        raise ValueError("analysis document binds a different catalog precommitment")
    document["catalog_precommitment_sha256"] = catalog_sha256
    supplied_lock = document.pop("analysis_lock_sha256", None)
    lock_sha256 = sha256_json(document)
    if supplied_lock is not None and supplied_lock != lock_sha256:
        raise ValueError("analysis document self-hash is invalid")
    document["analysis_lock_sha256"] = lock_sha256
    return document, lock_sha256


def _resources(circuits: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    by_role: dict[str, int] = {}
    total_shots = 0
    for descriptor in circuits.values():
        role = str(descriptor["backend_role"])
        by_role[role] = by_role.get(role, 0) + 1
        total_shots += int(descriptor["shots"])
    per_role = set(by_role.values())
    if len(per_role) != 1:
        raise ValueError("the two backend roles must have matched circuit catalogs")
    return {
        "account_quota_seconds": 600,
        "account_reserve_seconds": 180,
        "estimated_qpu_seconds_ceiling": 420,
        "max_execution_time_seconds": 420,
        "max_pubs_per_job": 300,
        "shots_per_circuit": 192,
        "shots_per_balanced_variant": 192,
        "shots_per_mh_edge": 192,
        "shots_per_diagnostic_calibration": 512,
        "backend_slot_count": 2,
        "logical_qubits_per_circuit": 4,
        "classical_bits_per_dynamic_circuit": 7,
        "catalog_entries_per_backend_slot": next(iter(per_role)),
        "total_circuit_instances": len(circuits),
        "total_planned_shots": total_shots,
        "transpile_optimization_level": 0,
        "transpiler_seed": 509,
        "resilience_level": 0,
        "dynamical_decoupling": False,
        "gate_twirling": False,
        "measurement_twirling": False,
        "common_duration_padding": True,
    }


def build_blinded_preregistration(
    analysis_document: Mapping[str, Any] | None = None,
    *,
    backend_slots: Sequence[BackendLayoutSlot] = DEFAULT_BACKEND_SLOTS,
    test_seed: str | None = None,
    _test_scope: Mapping[str, Sequence[Any]] | None = None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Build the public manifest and private reveal without submitting jobs.

    Production mode always builds the full frozen catalog.  ``_test_scope`` is
    deliberately private and is accepted only with deterministic test mode.
    """

    slots = tuple(backend_slots)
    _validate_backend_slots(slots)
    if _test_scope is not None and test_seed is None:
        raise ValueError("a reduced catalog is permitted only in deterministic test mode")
    scope = {
        "cayley_models": ("z5", "s3"),
        "cayley_protocols": CAYLEY_PROTOCOLS,
        "mh_spectra": MH_SPECTRA,
        "calibration_codes": tuple(range(16)),
    }
    if _test_scope is not None:
        scope.update({key: tuple(value) for key, value in _test_scope.items()})

    secret, secret_mode = _new_secret(test_seed)
    opaque_key = _derive_opaque_key(secret)
    kernels, _, _ = _load_circuit_modules()
    cayley_models = kernels.builtin_cayley_models()
    spectra = kernels.builtin_spectra()

    mappings: dict[str, Any] = {
        "model": {},
        "protocol": {},
        "state": {},
        "slots": {},
        "encoding": {},
        "layout": {},
        "backend": {},
        "backend_role": {},
    }
    for model_key in scope["cayley_models"]:
        model = cayley_models[str(model_key)]
        mappings["model"][str(model_key)] = _opaque_id(opaque_key, "model", model_key)
        mappings["state"][str(model_key)] = {
            label: _opaque_id(opaque_key, f"state:{model_key}", label)
            for label in model.states
        }
        degree = len(model.generators)
        mappings["slots"][str(model_key)] = {
            "disturbance": {
                str(index): _opaque_id(
                    opaque_key, f"slot:{model_key}:disturbance", index
                )
                for index in range(2 * degree)
            },
            "second": {
                str(index): _opaque_id(opaque_key, f"slot:{model_key}:second", index)
                for index in range(2 * degree)
            },
        }
    for spectrum_key in scope["mh_spectra"]:
        spectrum = spectra[str(spectrum_key)]
        mappings["model"][str(spectrum_key)] = _opaque_id(
            opaque_key, "model", spectrum_key
        )
        mappings["state"][str(spectrum_key)] = {
            label: _opaque_id(opaque_key, f"state:{spectrum_key}", label)
            for label in spectrum.labels
        }
    for protocol in scope["cayley_protocols"]:
        mappings["protocol"][str(protocol)] = _opaque_id(
            opaque_key, "protocol", protocol
        )

    private_slots: list[dict[str, Any]] = []
    public_slots: list[dict[str, str]] = []
    for slot in slots:
        role_id = _opaque_id(opaque_key, "backend-role", slot.role)
        backend_id = _opaque_id(opaque_key, "backend", slot.backend)
        layout_semantic = {"role": slot.role, "layout": list(slot.layout)}
        layout_id = _opaque_id(opaque_key, "layout", layout_semantic)
        mappings["backend_role"][slot.role] = role_id
        mappings["backend"][slot.backend] = backend_id
        mappings["layout"][json.dumps(layout_semantic, sort_keys=True)] = layout_id
        private_slots.append(
            {
                "role": slot.role,
                "role_opaque_id": role_id,
                "backend": slot.backend,
                "backend_opaque_id": backend_id,
                "layout": list(slot.layout),
                "layout_opaque_id": layout_id,
                "properties_last_update": slot.properties_last_update,
            }
        )
        public_slots.append({"role": role_id, "backend": backend_id, "layout": layout_id})

    encoding_assignments: dict[str, list[int]] = {}
    circuits: dict[str, dict[str, Any]] = {}
    for slot in slots:
        role_id = mappings["backend_role"][slot.role]
        for model_key in scope["cayley_models"]:
            model = cayley_models[str(model_key)]
            encoding = _keyed_permutation(
                opaque_key,
                f"encoding:{slot.role}:{model_key}",
                tuple(range(8)),
            )[: model.size]
            encoding_semantic = f"{slot.role}:{model_key}:{list(encoding)}"
            encoding_id = _opaque_id(opaque_key, "encoding", encoding_semantic)
            mappings["encoding"][encoding_semantic] = encoding_id
            encoding_assignments[encoding_id] = list(encoding)
            degree = len(model.generators)
            for protocol in scope["cayley_protocols"]:
                for initial_state in range(model.size):
                    for disturbance_slot in range(2 * degree):
                        for second_slot in range(2 * degree):
                            parameters = {
                                "model": str(model_key),
                                "protocol": str(protocol),
                                "initial_state": initial_state,
                                "disturbance_slot": disturbance_slot,
                                "second_slot": second_slot,
                                "state_encoding": list(encoding),
                            }
                            opaque_id, descriptor = _descriptor(
                                family="cayley",
                                role=slot.role,
                                role_opaque_id=role_id,
                                shots=192,
                                parameters=parameters,
                                opaque_key=opaque_key,
                            )
                            if opaque_id in circuits:
                                raise RuntimeError("opaque circuit identifier collision")
                            circuits[opaque_id] = descriptor

        for spectrum_key in scope["mh_spectra"]:
            spectrum = spectra[str(spectrum_key)]
            permutation = _keyed_permutation(
                opaque_key,
                f"mh-encoding:{slot.role}:{spectrum_key}",
                tuple(range(spectrum.size)),
            )
            encoding_semantic = f"{slot.role}:{spectrum_key}:{list(permutation)}"
            encoding_id = _opaque_id(opaque_key, "encoding", encoding_semantic)
            mappings["encoding"][encoding_semantic] = encoding_id
            encoding_assignments[encoding_id] = list(permutation)
            for source in range(spectrum.size):
                for target in range(spectrum.size):
                    if source == target:
                        continue
                    parameters = {
                        "spectrum": str(spectrum_key),
                        "beta": float(kernels.DEFAULT_BETA),
                        "programmed_kappa": float(kernels.OPH_KAPPA),
                        "semantic_source": source,
                        "semantic_target": target,
                        "label_permutation": list(permutation),
                    }
                    opaque_id, descriptor = _descriptor(
                        family="mh",
                        role=slot.role,
                        role_opaque_id=role_id,
                        shots=192,
                        parameters=parameters,
                        opaque_key=opaque_key,
                    )
                    circuits[opaque_id] = descriptor

        for basis_code in scope["calibration_codes"]:
            parameters = {
                "basis_code": int(basis_code),
                "num_qubits": 4,
                "evidentiary_role": "diagnostic_only",
            }
            opaque_id, descriptor = _descriptor(
                family="readout_calibration",
                role=slot.role,
                role_opaque_id=role_id,
                shots=512,
                parameters=parameters,
                opaque_key=opaque_key,
            )
            circuits[opaque_id] = descriptor

    public_circuits = [
        _public_circuit_entry(opaque_id, circuits[opaque_id])
        for opaque_id in sorted(circuits)
    ]
    resources = _resources(circuits)
    private_nulls: dict[str, Any] = {}
    public_nulls: list[dict[str, str]] = []
    for name, spec in DEFAULT_NULL_SPECS.items():
        null_id = _opaque_id(opaque_key, "null", name)
        spec_commitment = hashlib.sha256(
            COMMITMENT_DOMAIN
            + secret
            + b"null\0"
            + canonical_json_bytes({"name": name, "spec": spec})
        ).hexdigest()
        private_nulls[null_id] = {"semantic_name": name, "spec": spec}
        public_nulls.append(
            {"opaque_id": null_id, "sealed_spec_sha256": spec_commitment}
        )

    public: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "protocol_id": _opaque_id(opaque_key, "protocol-bundle", "oph-509-cayley"),
        "catalog_precommitment_sha256": "",
        "manifest_sha256": "",
        "analysis_lock_sha256": "",
        "secret_commitment_sha256": "",
        "backend_slots": public_slots,
        "circuits": public_circuits,
        "resources": resources,
        "thresholds": dict(DEFAULT_THRESHOLDS),
        "nulls": sorted(public_nulls, key=lambda row: row["opaque_id"]),
        "exclusions": dict(DEFAULT_EXCLUSIONS),
    }
    catalog_sha256 = catalog_precommitment(public)
    public["catalog_precommitment_sha256"] = catalog_sha256
    locked_analysis, analysis_sha256 = _lock_analysis_document(
        analysis_document, catalog_sha256
    )
    public["analysis_lock_sha256"] = analysis_sha256
    core_sha256 = manifest_core_hash(public)
    sealed_payload = {
        "mode": secret_mode,
        "public_manifest_core_sha256": core_sha256,
        "analysis_document": locked_analysis,
        "mappings": mappings,
        "encoding_assignments": encoding_assignments,
        "backend_slots": private_slots,
        "backend_layout_selection_policy": {
            "backend_rule": "lexicographically first two operational dynamic-circuit backends",
            "layout_rule": "exhaustive minimum role score with no outcome-based selection",
            "role_score": "2*measure_error + x_error + reset_error when exposed",
            "measurement_instruction": "measure",
            "logical_order": ["sector0", "sector1", "sector2", "decision"],
        },
        "circuits": {key: circuits[key] for key in sorted(circuits)},
        "nulls": private_nulls,
        "claim_boundary": (
            "Engineered self-reading record/feedback instrument; not evidence against unrestricted QM."
        ),
    }
    commitment_sha256 = secret_commitment(secret, sealed_payload)
    public["secret_commitment_sha256"] = commitment_sha256
    public["manifest_sha256"] = manifest_hash(public)
    reveal = {
        "schema_version": REVEAL_SCHEMA_VERSION,
        "secret_hex": secret.hex(),
        "opaque_id_key_hex": opaque_key.hex(),
        "sealed_payload": sealed_payload,
    }
    return public, reveal


def bind_analysis_document(
    public_manifest: Mapping[str, Any],
    reveal: Mapping[str, Any],
    analysis_document: Mapping[str, Any],
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Bind a hardened analysis lock without regenerating the blind catalog.

    This is the second phase of sealing: candidate tables may be derived from
    the already-generated opaque catalog, after which this function preserves
    every circuit, mapping, backend assignment, and secret while recomputing
    only the downstream analysis/core/commitment/final-manifest layers.
    """

    verify_bundle(public_manifest, reveal, rebuild_circuits=False)
    updated_public = copy.deepcopy(dict(public_manifest))
    updated_reveal = copy.deepcopy(dict(reveal))
    catalog_sha256 = str(updated_public["catalog_precommitment_sha256"])
    locked_analysis, analysis_sha256 = _lock_analysis_document(
        analysis_document,
        catalog_sha256,
    )
    updated_public["analysis_lock_sha256"] = analysis_sha256
    payload = updated_reveal["sealed_payload"]
    payload["analysis_document"] = locked_analysis
    payload["public_manifest_core_sha256"] = manifest_core_hash(updated_public)
    secret = bytes.fromhex(updated_reveal["secret_hex"])
    updated_public["secret_commitment_sha256"] = secret_commitment(secret, payload)
    updated_public["manifest_sha256"] = manifest_hash(updated_public)
    verify_bundle(
        updated_public,
        updated_reveal,
        locked_analysis,
        rebuild_circuits=False,
    )
    return updated_public, updated_reveal


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def assert_public_manifest_blinded(
    public_manifest: Mapping[str, Any],
    reveal: Mapping[str, Any],
) -> None:
    """Reject lexical disclosure of sealed semantic identifiers."""

    rendered = canonical_json_bytes(public_manifest).decode("utf-8").lower()
    payload = reveal["sealed_payload"]
    mappings = payload["mappings"]
    forbidden: set[str] = set()
    for category in ("model", "protocol", "backend", "backend_role"):
        forbidden.update(str(value).lower() for value in mappings[category].keys())
    for states in mappings["state"].values():
        forbidden.update(
            str(value).lower()
            for value in states
            if len(str(value)) >= 3 or any(character in str(value) for character in "()'")
        )
    leaked = sorted(
        token
        for token in forbidden
        if token
        and json.dumps(token, ensure_ascii=False, separators=(",", ":")) in rendered
    )
    if leaked:
        raise VerificationError(f"public manifest leaks sealed semantic tokens: {leaked}")


def verify_bundle(
    public_manifest: Mapping[str, Any],
    reveal: Mapping[str, Any],
    analysis_document: Mapping[str, Any] | None = None,
    *,
    rebuild_circuits: bool = False,
) -> bool:
    """Verify commitments and, optionally, every logical OpenQASM 3 hash."""

    required_public = {
        "schema_version",
        "protocol_id",
        "catalog_precommitment_sha256",
        "manifest_sha256",
        "analysis_lock_sha256",
        "secret_commitment_sha256",
        "backend_slots",
        "circuits",
        "resources",
        "thresholds",
        "nulls",
        "exclusions",
    }
    _require(set(public_manifest) == required_public, "public manifest keys are not canonical")
    _require(public_manifest["schema_version"] == SCHEMA_VERSION, "public schema mismatch")
    for field in (
        "catalog_precommitment_sha256",
        "manifest_sha256",
        "analysis_lock_sha256",
        "secret_commitment_sha256",
    ):
        _require(_is_sha256(public_manifest[field]), f"{field} is not lowercase SHA-256")
    _require(
        catalog_precommitment(public_manifest)
        == public_manifest["catalog_precommitment_sha256"],
        "catalog precommitment mismatch",
    )
    _require(
        manifest_hash(public_manifest) == public_manifest["manifest_sha256"],
        "final manifest digest mismatch",
    )
    _require(len(public_manifest["backend_slots"]) == 2, "exactly two public slots required")
    for slot in public_manifest["backend_slots"]:
        _require(set(slot) == {"role", "backend", "layout"}, "public slot keys mismatch")
        _require(all(str(value).startswith("o_") for value in slot.values()), "slot is not opaque")

    _require(
        set(reveal) == {
            "schema_version",
            "secret_hex",
            "opaque_id_key_hex",
            "sealed_payload",
        },
        "private reveal keys are not canonical",
    )
    _require(reveal["schema_version"] == REVEAL_SCHEMA_VERSION, "reveal schema mismatch")
    try:
        secret = bytes.fromhex(str(reveal["secret_hex"]))
        supplied_key = bytes.fromhex(str(reveal["opaque_id_key_hex"]))
    except ValueError as exc:
        raise VerificationError("reveal key material is not hexadecimal") from exc
    _require(len(secret) == SECRET_BYTES, "reveal secret is not 32 bytes")
    _require(supplied_key == _derive_opaque_key(secret), "opaque-ID key derivation mismatch")
    payload = reveal["sealed_payload"]
    _require(
        secret_commitment(secret, payload) == public_manifest["secret_commitment_sha256"],
        "secret commitment mismatch",
    )
    _require(
        payload["public_manifest_core_sha256"] == manifest_core_hash(public_manifest),
        "reveal does not bind the final manifest core",
    )

    locked_analysis = dict(payload["analysis_document"])
    supplied_analysis_hash = locked_analysis.pop("analysis_lock_sha256", None)
    computed_analysis_hash = sha256_json(locked_analysis)
    _require(
        supplied_analysis_hash == computed_analysis_hash,
        "analysis document self-hash mismatch",
    )
    _require(
        computed_analysis_hash == public_manifest["analysis_lock_sha256"],
        "public analysis-lock digest mismatch",
    )
    _require(
        locked_analysis.get("catalog_precommitment_sha256")
        == public_manifest["catalog_precommitment_sha256"],
        "analysis lock does not bind the catalog precommitment",
    )
    if analysis_document is not None:
        external = dict(analysis_document)
        external.pop("analysis_lock_sha256", None)
        _require(
            sha256_json(external) == computed_analysis_hash,
            "supplied analysis document differs from the sealed lock",
        )

    private_slots = payload["backend_slots"]
    _require(len(private_slots) == 2, "private reveal must contain two backend slots")
    reconstructed_public_slots = [
        {
            "role": slot["role_opaque_id"],
            "backend": slot["backend_opaque_id"],
            "layout": slot["layout_opaque_id"],
        }
        for slot in private_slots
    ]
    _require(
        reconstructed_public_slots == public_manifest["backend_slots"],
        "private backend/layout assignments differ from public slots",
    )
    _validate_backend_slots(
        tuple(
            BackendLayoutSlot(
                role=str(slot["role"]),
                backend=str(slot["backend"]),
                layout=tuple(int(value) for value in slot["layout"]),
                properties_last_update=str(slot["properties_last_update"]),
            )
            for slot in private_slots
        )
    )

    public_circuits = {row["opaque_id"]: row for row in public_manifest["circuits"]}
    _require(
        len(public_circuits) == len(public_manifest["circuits"]),
        "duplicate opaque circuit identifier",
    )
    private_circuits = payload["circuits"]
    _require(set(public_circuits) == set(private_circuits), "public/private circuit set mismatch")
    role_ids = {slot["role"] for slot in public_manifest["backend_slots"]}
    for opaque_id, descriptor in private_circuits.items():
        public_row = public_circuits[opaque_id]
        _require(
            set(public_row)
            == {"opaque_id", "logical_circuit_sha256", "shots", "backend_role"},
            f"public circuit row keys mismatch for {opaque_id}",
        )
        _require(
            _is_sha256(public_row["logical_circuit_sha256"]),
            "logical circuit hash is malformed",
        )
        _require(int(public_row["shots"]) > 0, "circuit shots must be positive")
        _require(public_row["backend_role"] in role_ids, "unknown opaque backend role")
        _require(
            public_row["logical_circuit_sha256"]
            == descriptor["logical_circuit_sha256"]
            and public_row["shots"] == descriptor["shots"]
            and public_row["backend_role"] == descriptor["backend_role_opaque_id"],
            f"public/private circuit descriptor mismatch for {opaque_id}",
        )
        metadata = descriptor["circuit_metadata"]
        _require(
            metadata
            == {
                "schema_version": BLINDED_CIRCUIT_SCHEMA_VERSION,
                "opaque_id": opaque_id,
                "blind_nonce": metadata.get("blind_nonce"),
            },
            f"blinded circuit metadata mismatch for {opaque_id}",
        )
        _require(_is_sha256(metadata["blind_nonce"]), "blind nonce is malformed")
        if descriptor["family"] == "readout_calibration":
            _require(
                descriptor["parameters"].get("evidentiary_role") == "diagnostic_only",
                "basis calibration may not claim primary causal-channel status",
            )
        if rebuild_circuits:
            verify_rebuilt_circuit(opaque_id, descriptor)

    for public_null in public_manifest["nulls"]:
        null_id = public_null["opaque_id"]
        _require(null_id in payload["nulls"], "public null is absent from reveal")
        private_null = payload["nulls"][null_id]
        expected = hashlib.sha256(
            COMMITMENT_DOMAIN
            + secret
            + b"null\0"
            + canonical_json_bytes(
                {"name": private_null["semantic_name"], "spec": private_null["spec"]}
            )
        ).hexdigest()
        _require(
            public_null["sealed_spec_sha256"] == expected,
            f"null-spec commitment mismatch for {null_id}",
        )
    _require(
        public_manifest["resources"]["backend_slot_count"] == 2,
        "resource contract does not require two backend slots",
    )
    _require(
        public_manifest["resources"]["account_reserve_seconds"] == 180
        and public_manifest["resources"]["estimated_qpu_seconds_ceiling"] == 420,
        "QPU reserve/ceiling contract drifted",
    )
    _require(public_manifest["thresholds"] == DEFAULT_THRESHOLDS, "threshold contract drifted")
    _require(public_manifest["exclusions"] == DEFAULT_EXCLUSIONS, "exclusion contract drifted")
    assert_public_manifest_blinded(public_manifest, reveal)
    return True


def _outside_checkout(path: Path, checkout_root: Path) -> Path:
    checkout = checkout_root.resolve()
    parent = path.parent.resolve()
    target = parent / path.name
    try:
        target.relative_to(checkout)
    except ValueError:
        return target
    raise ValueError("private reveal path must be outside the git checkout")


def write_reveal_file(
    path: Path,
    reveal: Mapping[str, Any],
    *,
    checkout_root: Path,
) -> Path:
    """Exclusively create a private reveal outside the checkout with mode 0600."""

    target = _outside_checkout(Path(path), Path(checkout_root))
    target.parent.mkdir(parents=True, exist_ok=True)
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(target, flags, 0o600)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(canonical_json_bytes(reveal) + b"\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.chmod(target, 0o600)
    except BaseException:
        target.unlink(missing_ok=True)
        raise
    mode = stat.S_IMODE(target.stat().st_mode)
    if mode != 0o600:
        target.unlink(missing_ok=True)
        raise PermissionError(f"private reveal mode is {oct(mode)}, expected 0o600")
    return target


def write_public_manifest(path: Path, manifest: Mapping[str, Any]) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        raise FileExistsError(target)
    target.write_bytes(canonical_json_bytes(manifest) + b"\n")
    return target


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError("JSON document must be an object")
    return value


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Seal or verify the blinded IBM preregistration.")
    subparsers = parser.add_subparsers(dest="command", required=True)
    create = subparsers.add_parser("create")
    create.add_argument("--public-out", type=Path, required=True)
    create.add_argument("--reveal-out", type=Path, required=True)
    create.add_argument("--analysis-document", type=Path)
    create.add_argument("--checkout-root", type=Path, default=Path(__file__).resolve().parents[3])
    create.add_argument("--deterministic-test-seed")
    verify = subparsers.add_parser("verify")
    verify.add_argument("--public", type=Path, required=True)
    verify.add_argument("--reveal", type=Path, required=True)
    verify.add_argument("--analysis-document", type=Path)
    verify.add_argument("--rebuild-circuits", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    if args.command == "create":
        analysis = load_json(args.analysis_document) if args.analysis_document else None
        public, reveal = build_blinded_preregistration(
            analysis,
            test_seed=args.deterministic_test_seed,
        )
        verify_bundle(public, reveal, rebuild_circuits=False)
        write_reveal_file(args.reveal_out, reveal, checkout_root=args.checkout_root)
        write_public_manifest(args.public_out, public)
        print(public["manifest_sha256"])
        return 0
    public = load_json(args.public)
    reveal = load_json(args.reveal)
    analysis = load_json(args.analysis_document) if args.analysis_document else None
    verify_bundle(
        public,
        reveal,
        analysis,
        rebuild_circuits=args.rebuild_circuits,
    )
    print(public["manifest_sha256"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
