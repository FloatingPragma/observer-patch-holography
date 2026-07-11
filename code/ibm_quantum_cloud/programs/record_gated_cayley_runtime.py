#!/usr/bin/env python3
"""Guarded IBM Runtime operator for the blinded record-gated benchmark.

The default command is ``plan``.  It authenticates, snapshots the frozen
backends, compiles the reveal-bound circuits at the preregistered layouts, and
checks the account reserve, but it cannot submit a job.  Submission requires a
separate ``submit`` action, ``--confirm-submit``, and literal confirmations of
both the manifest and analysis-lock SHA-256 digests.

This module deliberately uses SamplerV2 in backend (job) mode.  It never opens
a Session or Batch, and it explicitly disables dynamical decoupling and both
forms of twirling.  Operator artifacts contain hashes, opaque identifiers,
backend snapshots, job identifiers, statuses, metrics, and raw counts; the
private reveal and API token are never copied to them or printed.
"""

from __future__ import annotations

import argparse
import dataclasses
import fcntl
import hashlib
import importlib.metadata
import io
import json
import logging
import math
import os
import re
import sys
from collections import OrderedDict
from dataclasses import dataclass
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Sequence


MANIFEST_SCHEMA = "oph.blinded-record-gated-preregistration.v1"
REVEAL_SCHEMA = "oph.blinded-record-gated-reveal.v1"
RUNTIME_SCHEMA = "oph.record-gated-cayley-runtime.v1"
COMMITMENT_DOMAIN = b"oph-509-blind-v1\0"
HEX64 = re.compile(r"[0-9a-f]{64}\Z")
EXPECTED_RUNTIME_VERSION = "0.47.0"
EXPECTED_QISKIT_MAJOR_MINOR = (2, 5)
LOGICAL_CIRCUIT_SERIALIZATION = "normalized-openqasm3-utf8-v1"
LOGICAL_QUBITS = 4
IBM_SUBJOB_OVERHEAD_SECONDS = 2.0
ESTIMATE_SAFETY_FACTOR = 1.25
# IBM's published baseline describes complete circuits as typically 50--100 us.
# Measurement durations are already present in Target, so reserve 40 us for
# each duration-less classical feed-forward operation and then apply the
# separate 1.25 workload safety factor below.
DYNAMIC_CONTROL_LATENCY_BOUND_SECONDS = 0.00004
UNKNOWN_INSTRUCTION_DURATION_BOUND_SECONDS = 0.0001
FINAL_JOB_STATES = {"DONE", "ERROR", "CANCELLED"}
REDACTED_KEYS = (
    "authorization",
    "credential",
    "password",
    "private_key",
    "secret",
    "token",
)


class RuntimeSafetyError(RuntimeError):
    """A fail-closed operator safety check did not pass."""


@dataclass(frozen=True)
class VerifiedBundle:
    manifest: dict[str, Any]
    manifest_sha256: str
    analysis_lock_sha256: str
    analysis_lock_document: Mapping[str, Any]
    public_manifest_core_sha256: str
    reveal_mode: str
    backend_slots: tuple[dict[str, Any], ...]
    circuit_families: Mapping[str, str]
    circuit_descriptors: Mapping[str, Mapping[str, Any]] = dataclasses.field(
        repr=False, compare=False
    )
    ideal_validation: Mapping[str, Any]
    circuits: tuple[tuple[dict[str, Any], Any], ...]


@dataclass(frozen=True)
class PreparedCircuit:
    opaque_id: str
    family: str
    shots: int
    backend_role: str
    logical_circuit_sha256: str
    compiled_qpy_sha256: str
    compiled_duration_seconds: float
    logical_circuit: Any = dataclasses.field(repr=False, compare=False)
    compiled_circuit: Any = dataclasses.field(repr=False, compare=False)


@dataclass(frozen=True)
class PreparedGroup:
    group_id: str
    backend_role: str
    backend_role_opaque_id: str
    family: str
    backend_name: str
    layout_opaque_id: str
    physical_layout: tuple[int, ...]
    properties_last_update: str
    shots: int
    estimated_qpu_seconds: float
    circuits: tuple[PreparedCircuit, ...]
    compiled_qpy: bytes = dataclasses.field(repr=False, compare=False)


@dataclass(frozen=True)
class PreparedRun:
    plan: dict[str, Any]
    groups: tuple[PreparedGroup, ...]
    backends: Mapping[str, Any] = dataclasses.field(repr=False, compare=False)


@dataclass(frozen=True)
class RuntimeBindings:
    SamplerV2: Any
    SamplerOptions: Any
    transpile: Callable[..., Any]
    qpy_dump: Callable[[Sequence[Any], Any], None]
    qiskit_version: str
    runtime_version: str


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def canonical_json_bytes(value: Any) -> bytes:
    try:
        rendered = json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
            allow_nan=False,
        )
    except (TypeError, ValueError) as exc:
        raise RuntimeSafetyError("document is not canonical finite JSON") from exc
    return rendered.encode("utf-8")


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_json(value: Any) -> str:
    return sha256_bytes(canonical_json_bytes(value))


def operator_source_sha256() -> str:
    """Hash the exact operator source file used to build runtime receipts."""

    try:
        source = Path(__file__).resolve().read_bytes()
    except OSError as exc:
        raise RuntimeSafetyError("could not hash the Runtime operator source") from exc
    return sha256_bytes(source)


def compiled_qpy_artifact(group_id: str, payload: bytes) -> dict[str, str]:
    """Describe a process-local compiled-QPY receipt by content, not job group."""

    group_id = _require_hex64(group_id, "compiled-QPY group_id")
    digest = sha256_bytes(payload)
    return {
        "relative_path": f"runtime_compiled_{group_id}_{digest}.qpy",
        "sha256": digest,
    }


def _without(mapping: Mapping[str, Any], *keys: str) -> dict[str, Any]:
    return {key: value for key, value in mapping.items() if key not in keys}


def _require_hex64(value: Any, label: str) -> str:
    if not isinstance(value, str) or HEX64.fullmatch(value) is None:
        raise RuntimeSafetyError(f"{label} must be a lowercase SHA-256 digest")
    return value


def load_json(path: Path, label: str) -> Any:
    try:
        with path.open("r", encoding="utf-8") as handle:
            return json.load(handle)
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise RuntimeSafetyError(f"could not load {label}") from exc


def _require_mapping(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise RuntimeSafetyError(f"{label} must be a JSON object")
    return value


def _manifest_core(manifest: Mapping[str, Any]) -> dict[str, Any]:
    return _without(manifest, "manifest_sha256", "secret_commitment_sha256")


def _catalog_precommitment_body(manifest: Mapping[str, Any]) -> dict[str, Any]:
    return _without(
        manifest,
        "catalog_precommitment_sha256",
        "analysis_lock_sha256",
        "secret_commitment_sha256",
        "manifest_sha256",
    )


def verify_manifest(manifest: Mapping[str, Any]) -> tuple[str, str]:
    required = {
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
    missing = sorted(required.difference(manifest))
    if missing:
        raise RuntimeSafetyError("public manifest is missing required fields")
    if manifest["schema_version"] != MANIFEST_SCHEMA:
        raise RuntimeSafetyError("unsupported public manifest schema")
    if not isinstance(manifest["protocol_id"], str) or not manifest["protocol_id"]:
        raise RuntimeSafetyError("public manifest protocol_id is invalid")

    claimed_catalog = _require_hex64(
        manifest["catalog_precommitment_sha256"], "catalog_precommitment_sha256"
    )
    if sha256_json(_catalog_precommitment_body(manifest)) != claimed_catalog:
        raise RuntimeSafetyError("public catalog precommitment does not verify")

    claimed_manifest = _require_hex64(manifest["manifest_sha256"], "manifest_sha256")
    computed_manifest = sha256_json(_without(manifest, "manifest_sha256"))
    if claimed_manifest != computed_manifest:
        raise RuntimeSafetyError("public manifest digest does not verify")
    _require_hex64(manifest["analysis_lock_sha256"], "analysis_lock_sha256")
    _require_hex64(manifest["secret_commitment_sha256"], "secret_commitment_sha256")

    slots = manifest["backend_slots"]
    if not isinstance(slots, list) or not slots:
        raise RuntimeSafetyError("backend_slots must be a non-empty list")
    roles: set[str] = set()
    for slot in slots:
        if not isinstance(slot, dict):
            raise RuntimeSafetyError("each backend slot must be an object")
        if set(slot) != {"role", "backend", "layout"}:
            raise RuntimeSafetyError("backend slot lacks role, backend, or layout")
        role, backend, layout = slot["role"], slot["backend"], slot["layout"]
        if not isinstance(role, str) or not role or role in roles:
            raise RuntimeSafetyError("backend slot roles must be unique non-empty strings")
        roles.add(role)
        # All three public values are opaque.  Physical names and qubits exist
        # only in the operator reveal and are validated after commitment check.
        if not isinstance(backend, str) or not backend:
            raise RuntimeSafetyError("backend slot has an invalid opaque backend identifier")
        if not isinstance(layout, str) or not layout:
            raise RuntimeSafetyError("backend slot has an invalid opaque layout identifier")

    circuit_rows = manifest["circuits"]
    if not isinstance(circuit_rows, list) or not circuit_rows:
        raise RuntimeSafetyError("circuits must be a non-empty list")
    opaque_ids: set[str] = set()
    for row in circuit_rows:
        if not isinstance(row, dict) or set(row) != {
            "opaque_id",
            "logical_circuit_sha256",
            "shots",
            "backend_role",
        }:
            raise RuntimeSafetyError("public circuit row lacks required fields")
        opaque_id = row["opaque_id"]
        if not isinstance(opaque_id, str) or not opaque_id or opaque_id in opaque_ids:
            raise RuntimeSafetyError("opaque circuit identifiers must be unique")
        opaque_ids.add(opaque_id)
        _require_hex64(
            row["logical_circuit_sha256"], "circuit logical_circuit_sha256"
        )
        if type(row["shots"]) is not int or row["shots"] <= 0:
            raise RuntimeSafetyError("circuit shots must be a positive integer")
        if row["backend_role"] not in roles:
            raise RuntimeSafetyError("circuit references an unknown backend role")

    resources = _require_mapping(manifest["resources"], "resources")
    _validate_resources(resources)
    if resources["backend_slot_count"] != len(slots):
        raise RuntimeSafetyError("backend_slot_count does not match backend_slots")
    if resources["total_circuit_instances"] != len(circuit_rows):
        raise RuntimeSafetyError("total_circuit_instances does not match circuits")
    if resources["total_planned_shots"] != sum(row["shots"] for row in circuit_rows):
        raise RuntimeSafetyError("total_planned_shots does not match circuits")
    allowed_shots = {
        resources["shots_per_balanced_variant"],
        resources["shots_per_mh_edge"],
        resources["shots_per_diagnostic_calibration"],
    }
    if any(row["shots"] not in allowed_shots for row in circuit_rows):
        raise RuntimeSafetyError("public circuit shots differ from all frozen shot classes")
    counts_by_role = {
        role: sum(row["backend_role"] == role for row in circuit_rows) for role in roles
    }
    catalog_counts = resources["catalog_entries_per_backend_slot"]
    if isinstance(catalog_counts, int):
        if catalog_counts <= 0 or any(count != catalog_counts for count in counts_by_role.values()):
            raise RuntimeSafetyError("catalog_entries_per_backend_slot does not match circuits")
    elif isinstance(catalog_counts, dict):
        if catalog_counts != counts_by_role:
            raise RuntimeSafetyError("catalog_entries_per_backend_slot does not match circuits")
    else:
        raise RuntimeSafetyError("catalog_entries_per_backend_slot has an invalid type")
    return computed_manifest, sha256_json(_manifest_core(manifest))


def _validate_resources(resources: Mapping[str, Any]) -> None:
    required = {
        "account_quota_seconds",
        "account_reserve_seconds",
        "estimated_qpu_seconds_ceiling",
        "max_pubs_per_job",
        "max_execution_time_seconds",
        "shots_per_circuit",
        "shots_per_balanced_variant",
        "shots_per_mh_edge",
        "shots_per_diagnostic_calibration",
        "backend_slot_count",
        "logical_qubits_per_circuit",
        "classical_bits_per_dynamic_circuit",
        "catalog_entries_per_backend_slot",
        "total_circuit_instances",
        "total_planned_shots",
        "transpile_optimization_level",
        "transpiler_seed",
        "resilience_level",
        "dynamical_decoupling",
        "gate_twirling",
        "measurement_twirling",
        "common_duration_padding",
    }
    if not required.issubset(resources):
        raise RuntimeSafetyError("resources lack a required execution or usage guard")
    limit = resources["account_quota_seconds"]
    reserve = resources["account_reserve_seconds"]
    ceiling = resources["estimated_qpu_seconds_ceiling"]
    for value, label in ((limit, "account limit"), (reserve, "account reserve"), (ceiling, "QPU ceiling")):
        if type(value) not in (int, float) or not math.isfinite(value) or value < 0:
            raise RuntimeSafetyError(f"{label} must be a finite non-negative number")
    if limit <= 0 or reserve >= limit or ceiling > limit - reserve:
        raise RuntimeSafetyError("resource ceiling would consume the protected account reserve")
    for key in (
        "max_pubs_per_job",
        "max_execution_time_seconds",
        "shots_per_circuit",
        "shots_per_balanced_variant",
        "shots_per_mh_edge",
        "shots_per_diagnostic_calibration",
        "backend_slot_count",
        "logical_qubits_per_circuit",
        "classical_bits_per_dynamic_circuit",
        "total_circuit_instances",
        "total_planned_shots",
    ):
        if type(resources[key]) is not int or resources[key] <= 0:
            raise RuntimeSafetyError(f"{key} must be a positive integer")
    if resources["logical_qubits_per_circuit"] != LOGICAL_QUBITS:
        raise RuntimeSafetyError("logical_qubits_per_circuit must be four")
    if resources["classical_bits_per_dynamic_circuit"] != 7:
        raise RuntimeSafetyError("classical_bits_per_dynamic_circuit must be seven")
    if resources["transpile_optimization_level"] not in (0, 1, 2, 3):
        raise RuntimeSafetyError("transpile_optimization_level must be 0, 1, 2, or 3")
    if type(resources["transpiler_seed"]) is not int or resources["transpiler_seed"] < 0:
        raise RuntimeSafetyError("transpiler_seed must be a non-negative integer")
    if resources["resilience_level"] != 0:
        raise RuntimeSafetyError("Sampler resilience_level must remain zero")
    for key in ("dynamical_decoupling", "gate_twirling", "measurement_twirling"):
        if resources[key] is not False:
            raise RuntimeSafetyError(f"{key} must be frozen false")
    if resources["common_duration_padding"] is not True:
        raise RuntimeSafetyError("common_duration_padding must be frozen true")


def verify_analysis_lock(
    document: Any, expected_sha256: str, catalog_precommitment_sha256: str | None = None
) -> str:
    expected_sha256 = _require_hex64(expected_sha256, "analysis_lock_sha256")
    document = _require_mapping(document, "analysis-lock document")
    internal = _require_hex64(document.get("analysis_lock_sha256"), "analysis lock self-digest")
    observed = sha256_json(_without(document, "analysis_lock_sha256"))
    if observed != internal or observed != expected_sha256:
        raise RuntimeSafetyError("analysis-lock document digest does not match the manifest")
    if catalog_precommitment_sha256 is not None:
        _require_hex64(catalog_precommitment_sha256, "catalog precommitment")
        bindings = [
            document[key]
            for key in ("catalog_precommitment_sha256", "blind_manifest_commitment")
            if key in document
        ]
        if not bindings or any(value != catalog_precommitment_sha256 for value in bindings):
            raise RuntimeSafetyError("analysis lock is bound to a different catalog precommitment")
    return observed


def verify_reveal(
    reveal: Mapping[str, Any],
    manifest: Mapping[str, Any],
    public_manifest_core_sha256: str,
) -> dict[str, Any]:
    if reveal.get("schema_version") != REVEAL_SCHEMA:
        raise RuntimeSafetyError("unsupported private reveal schema")
    secret_hex = reveal.get("secret_hex")
    if not isinstance(secret_hex, str) or re.fullmatch(r"[0-9a-f]{64}", secret_hex) is None:
        raise RuntimeSafetyError("private reveal secret is malformed")
    payload = _require_mapping(reveal.get("sealed_payload"), "sealed_payload")
    if payload.get("mode") not in {"production_random", "deterministic_test_only"}:
        raise RuntimeSafetyError("private reveal mode is invalid")
    if payload.get("public_manifest_core_sha256") != public_manifest_core_sha256:
        raise RuntimeSafetyError("private reveal is bound to a different manifest core")
    if payload.get("protocol_id") not in (None, manifest["protocol_id"]):
        raise RuntimeSafetyError("private reveal protocol binding does not match")
    computed = sha256_bytes(
        COMMITMENT_DOMAIN
        + bytes.fromhex(secret_hex)
        + b"\0"
        + canonical_json_bytes(payload)
    )
    if computed != manifest["secret_commitment_sha256"]:
        raise RuntimeSafetyError("private reveal commitment does not verify")
    circuits = payload.get("circuits")
    if not isinstance(circuits, dict):
        raise RuntimeSafetyError("private reveal circuit table must be keyed by opaque_id")
    public_ids = {row["opaque_id"] for row in manifest["circuits"]}
    if set(circuits) != public_ids:
        raise RuntimeSafetyError("private reveal circuit identifiers do not match the manifest")
    return payload


def qpy_bytes(circuits: Sequence[Any], qpy_dump: Callable[[Sequence[Any], Any], None]) -> bytes:
    buffer = io.BytesIO()
    qpy_dump(circuits, buffer)
    return buffer.getvalue()


def _descriptor_value(descriptor: Mapping[str, Any], *names: str) -> Any:
    found = [name for name in names if name in descriptor]
    if len(found) != 1:
        raise RuntimeSafetyError(f"reveal circuit descriptor must define exactly one of {', '.join(names)}")
    return descriptor[found[0]]


def _cayley_recipe_from_descriptor(descriptor: Mapping[str, Any]) -> Any:
    try:
        from generative_repair_kernel import builtin_cayley_models
        from record_gated_cayley_circuits import build_recipe
    except ImportError as exc:
        raise RuntimeSafetyError("Qiskit circuit dependencies are unavailable") from exc
    model_key = _descriptor_value(descriptor, "model_key", "model")
    models = builtin_cayley_models()
    if model_key not in models:
        raise RuntimeSafetyError("reveal references an unsupported model")
    encoding_raw = descriptor.get("state_encoding")
    if not isinstance(encoding_raw, list):
        raise RuntimeSafetyError("private reveal contains an invalid state encoding")
    try:
        return build_recipe(
            models[model_key],
            str(descriptor["protocol"]),
            int(descriptor["initial_state"]),
            int(descriptor["disturbance_slot"]),
            int(descriptor["second_slot"]),
            tuple(int(value) for value in encoding_raw),
        )
    except (KeyError, TypeError, ValueError) as exc:
        raise RuntimeSafetyError("private reveal contains an invalid Cayley recipe") from exc


def _mh_recipe_from_descriptor(descriptor: Mapping[str, Any]) -> Any:
    try:
        from generative_repair_circuits import edge_recipe
        from generative_repair_kernel import builtin_spectra
    except ImportError as exc:
        raise RuntimeSafetyError("MH circuit dependencies are unavailable") from exc
    spectrum_key = _descriptor_value(descriptor, "spectrum_key", "spectrum")
    spectra = builtin_spectra()
    if spectrum_key not in spectra:
        raise RuntimeSafetyError("reveal references an unsupported spectrum")
    permutation = descriptor.get("label_permutation")
    if not isinstance(permutation, list):
        raise RuntimeSafetyError("private reveal contains an invalid label permutation")
    try:
        return edge_recipe(
            spectra[spectrum_key],
            float(descriptor["beta"]),
            float(descriptor["programmed_kappa"]),
            int(descriptor["semantic_source"]),
            int(descriptor["semantic_target"]),
            tuple(int(value) for value in permutation),
        )
    except (KeyError, TypeError, ValueError) as exc:
        raise RuntimeSafetyError("private reveal contains an invalid MH recipe") from exc


def rebuild_reveal_circuit(opaque_id: str, descriptor: Mapping[str, Any]) -> Any:
    """Use the preregistration's single canonical blinded-circuit builder."""

    try:
        from blind_preregister import rebuild_blinded_circuit

        return rebuild_blinded_circuit(opaque_id, descriptor)
    except (ImportError, KeyError, TypeError, ValueError) as exc:
        raise RuntimeSafetyError("private reveal contains an invalid circuit descriptor") from exc


def resolve_backend_slots(
    public_slots: Sequence[Mapping[str, Any]], payload: Mapping[str, Any]
) -> tuple[dict[str, Any], ...]:
    private_slots = payload.get("backend_slots")
    if not isinstance(private_slots, list) or len(private_slots) != len(public_slots):
        raise RuntimeSafetyError("private backend slot table does not match the manifest")
    private_by_role_id: dict[str, Mapping[str, Any]] = {}
    for slot in private_slots:
        if not isinstance(slot, dict) or set(slot) != {
            "role",
            "role_opaque_id",
            "backend",
            "backend_opaque_id",
            "layout",
            "layout_opaque_id",
            "properties_last_update",
        }:
            raise RuntimeSafetyError("private backend slot has an unexpected schema")
        role_id = slot["role_opaque_id"]
        if not isinstance(role_id, str) or not role_id or role_id in private_by_role_id:
            raise RuntimeSafetyError("private backend role identifiers are invalid")
        if slot["role"] not in {"development", "heldout"}:
            raise RuntimeSafetyError("private backend semantic role is invalid")
        if not isinstance(slot["backend"], str) or not slot["backend"].startswith("ibm_"):
            raise RuntimeSafetyError("private backend name is invalid")
        layout = slot["layout"]
        if (
            not isinstance(layout, list)
            or len(layout) != LOGICAL_QUBITS
            or any(type(index) is not int or index < 0 for index in layout)
            or len(set(layout)) != LOGICAL_QUBITS
        ):
            raise RuntimeSafetyError("private fixed layout is invalid")
        if not isinstance(slot["properties_last_update"], str):
            raise RuntimeSafetyError("private calibration-selection timestamp is invalid")
        private_by_role_id[role_id] = slot

    resolved = []
    for public in public_slots:
        private = private_by_role_id.get(public["role"])
        if private is None:
            raise RuntimeSafetyError("public backend role has no committed private mapping")
        if (
            private["backend_opaque_id"] != public["backend"]
            or private["layout_opaque_id"] != public["layout"]
        ):
            raise RuntimeSafetyError("private backend or layout mapping does not match public slot")
        resolved.append(
            {
                "public_role": public["role"],
                "public_backend": public["backend"],
                "public_layout": public["layout"],
                "role": private["role"],
                "backend": private["backend"],
                "layout": list(private["layout"]),
                "properties_last_update": private["properties_last_update"],
            }
        )
    if {slot["role"] for slot in resolved} != {"development", "heldout"}:
        raise RuntimeSafetyError("operator reveal must bind development and heldout roles exactly")
    return tuple(resolved)


def validate_ideal_reveal_circuits(
    rows: Sequence[tuple[str, Mapping[str, Any], Any]],
) -> dict[str, Any]:
    """Fail closed when any submitted logical circuit disagrees with its recipe.

    Cayley and calibration circuits are deterministic and use one ideal shot.
    MH circuits use 2048 ideal shots; their record/feedback consistency is
    exact, while acceptance is checked against an eight-sigma binomial band
    (with a one-shot continuity allowance).
    """

    try:
        from qiskit import transpile
        from qiskit_aer import AerSimulator
        from generative_repair_circuits import parse_grouped_counts
        from record_gated_cayley_circuits import parse_counts
    except ImportError as exc:
        raise RuntimeSafetyError("Aer is required for mandatory ideal recipe validation") from exc

    backend = AerSimulator(seed_simulator=509)
    checked = 0
    counts_by_family: dict[str, int] = {}
    for family, shots in (("cayley", 1), ("mh", 2048), ("readout_calibration", 1)):
        selected = [row for row in rows if row[1].get("family") == family]
        if not selected:
            continue
        circuits = [row[2] for row in selected]
        try:
            isa = transpile(circuits, backend=backend, optimization_level=0, seed_transpiler=509)
            result = backend.run(isa, shots=shots).result()
        except Exception as exc:
            raise RuntimeSafetyError("mandatory ideal circuit execution failed") from exc
        for (opaque_id, descriptor, _), compiled in zip(selected, isa):
            parameters = _require_mapping(descriptor.get("parameters"), "circuit parameters")
            try:
                raw_counts = result.get_counts(compiled)
            except Exception as exc:
                raise RuntimeSafetyError("mandatory ideal counts are unavailable") from exc
            if family == "cayley":
                recipe_descriptor = dict(parameters)
                recipe = _cayley_recipe_from_descriptor(recipe_descriptor)
                parsed = parse_counts(raw_counts)
                expected = {
                    "heated": recipe.state_encoding[recipe.heated_state],
                    "decision": recipe.decision_record,
                    "final": recipe.state_encoding[recipe.final_state],
                }
                if sum(item["count"] for item in parsed) != shots or any(
                    any(item[key] != value for key, value in expected.items()) for item in parsed
                ):
                    raise RuntimeSafetyError("ideal Cayley circuit disagrees with its recipe")
            elif family == "mh":
                recipe_descriptor = dict(parameters)
                recipe = _mh_recipe_from_descriptor(recipe_descriptor)
                parsed = parse_grouped_counts(raw_counts)
                accepted = 0
                total = 0
                for item in parsed:
                    total += item["count"]
                    accepted += item["count"] if item["accepted"] == 1 else 0
                    expected_after = (
                        recipe.physical_target if item["accepted"] == 1 else recipe.physical_source
                    )
                    if item["before"] != recipe.physical_source or item["after"] != expected_after:
                        raise RuntimeSafetyError("ideal MH circuit violates record-gated feedback")
                if total != shots:
                    raise RuntimeSafetyError("ideal MH circuit returned the wrong number of shots")
                probability = recipe.expected_acceptance
                tolerance = 8.0 * math.sqrt(max(probability * (1.0 - probability), 1e-12) / shots)
                tolerance += 1.0 / shots
                if abs(accepted / shots - probability) > tolerance:
                    raise RuntimeSafetyError("ideal MH acceptance disagrees with its recipe")
            else:
                basis_code = parameters.get("basis_code")
                num_qubits = parameters.get("num_qubits")
                if type(basis_code) is not int or type(num_qubits) is not int:
                    raise RuntimeSafetyError("readout calibration recipe is invalid")
                if sum(int(value) for value in raw_counts.values()) != shots:
                    raise RuntimeSafetyError("ideal readout calibration returned wrong shot count")
                # The canonical builder names one classical register.  Leading
                # zeros may be retained or omitted by simulators, so compare as
                # integers after removing register separators.
                if any(int(str(key).replace(" ", ""), 2) != basis_code for key in raw_counts):
                    raise RuntimeSafetyError("ideal readout calibration disagrees with its basis code")
            checked += 1
            counts_by_family[family] = counts_by_family.get(family, 0) + 1
            del opaque_id  # never include semantic-adjacent identifiers in this receipt
    if checked != len(rows):
        raise RuntimeSafetyError("one or more circuit families escaped ideal validation")
    return {"circuits_checked": checked, "families": counts_by_family, "passed": True}


def verify_operator_bundle(
    manifest_path: Path,
    reveal_path: Path,
    analysis_lock_path: Path,
    *,
    circuit_rebuilder: Callable[[str, Mapping[str, Any]], Any] = rebuild_reveal_circuit,
    ideal_validator: Callable[
        [Sequence[tuple[str, Mapping[str, Any], Any]]], Mapping[str, Any]
    ] = validate_ideal_reveal_circuits,
    logical_digest: Callable[[Any], str] | None = None,
) -> VerifiedBundle:
    manifest = _require_mapping(load_json(manifest_path, "public manifest"), "public manifest")
    manifest_sha, core_sha = verify_manifest(manifest)
    analysis_document = load_json(analysis_lock_path, "analysis-lock document")
    analysis_sha = verify_analysis_lock(
        analysis_document,
        manifest["analysis_lock_sha256"],
        manifest["catalog_precommitment_sha256"],
    )
    reveal = _require_mapping(load_json(reveal_path, "private reveal"), "private reveal")
    payload = verify_reveal(reveal, manifest, core_sha)
    resolved_slots = resolve_backend_slots(manifest["backend_slots"], payload)

    if logical_digest is None:
        try:
            from blind_preregister import (
                LOGICAL_CIRCUIT_SERIALIZATION as prereg_logical_serialization,
                logical_circuit_sha256,
            )
        except ImportError as exc:
            raise RuntimeSafetyError(
                "canonical logical-circuit hashing is unavailable"
            ) from exc
        if prereg_logical_serialization != LOGICAL_CIRCUIT_SERIALIZATION:
            raise RuntimeSafetyError(
                "preregistration logical-circuit serialization contract drifted"
            )
        logical_digest = logical_circuit_sha256

    public_by_id = {row["opaque_id"]: row for row in manifest["circuits"]}
    resolved_by_public_role = {slot["public_role"]: slot for slot in resolved_slots}
    verified: list[tuple[dict[str, Any], Any]] = []
    ideal_rows: list[tuple[str, Mapping[str, Any], Any]] = []
    for opaque_id in [row["opaque_id"] for row in manifest["circuits"]]:
        descriptor = _require_mapping(payload["circuits"][opaque_id], "reveal circuit descriptor")
        if set(descriptor) != {
            "family",
            "backend_role",
            "backend_role_opaque_id",
            "shots",
            "logical_circuit_sha256",
            "circuit_name",
            "circuit_metadata",
            "parameters",
        }:
            raise RuntimeSafetyError("private circuit descriptor has an unexpected schema")
        family = descriptor["family"]
        if family not in {"cayley", "mh", "readout_calibration"}:
            raise RuntimeSafetyError("private circuit descriptor has an unsupported family")
        expected_family_shots = {
            "cayley": manifest["resources"]["shots_per_balanced_variant"],
            "mh": manifest["resources"]["shots_per_mh_edge"],
            "readout_calibration": manifest["resources"][
                "shots_per_diagnostic_calibration"
            ],
        }[family]
        public = public_by_id[opaque_id]
        slot = resolved_by_public_role.get(public["backend_role"])
        if slot is None:
            raise RuntimeSafetyError("private circuit references an unresolved backend role")
        if (
            descriptor["backend_role_opaque_id"] != public["backend_role"]
            or descriptor["backend_role"] != slot["role"]
            or descriptor["shots"] != public["shots"]
            or descriptor["shots"] != expected_family_shots
            or descriptor["circuit_name"] != opaque_id
        ):
            raise RuntimeSafetyError("private circuit role, shots, or name does not match public row")
        circuit = circuit_rebuilder(opaque_id, descriptor)
        try:
            observed_logical_raw = logical_digest(circuit)
        except Exception as exc:
            raise RuntimeSafetyError(
                "could not serialize the reconstructed logical circuit canonically"
            ) from exc
        observed_logical = _require_hex64(
            observed_logical_raw, "reconstructed logical circuit digest"
        )
        revealed_logical = _require_hex64(
            descriptor.get("logical_circuit_sha256"),
            "reveal logical_circuit_sha256",
        )
        if (
            observed_logical != revealed_logical
            or observed_logical != public["logical_circuit_sha256"]
        ):
            raise RuntimeSafetyError(
                "reconstructed canonical logical-circuit digest does not verify"
            )
        verified.append((public, circuit))
        ideal_rows.append((opaque_id, descriptor, circuit))
    validation = _safe_json(ideal_validator(ideal_rows))
    if validation.get("passed") is not True or validation.get("circuits_checked") != len(ideal_rows):
        raise RuntimeSafetyError("mandatory ideal recipe validation did not cover every circuit")
    return VerifiedBundle(
        manifest=dict(manifest),
        manifest_sha256=manifest_sha,
        analysis_lock_sha256=analysis_sha,
        analysis_lock_document=analysis_document,
        public_manifest_core_sha256=core_sha,
        reveal_mode=payload["mode"],
        backend_slots=resolved_slots,
        circuit_families={row[0]: row[1]["family"] for row in ideal_rows},
        circuit_descriptors={row[0]: row[1] for row in ideal_rows},
        ideal_validation=validation,
        circuits=tuple(verified),
    )


def _parse_version(version: str) -> tuple[int, ...]:
    numbers = []
    for part in version.split("."):
        match = re.match(r"(\d+)", part)
        if match is None:
            break
        numbers.append(int(match.group(1)))
    return tuple(numbers)


def load_runtime_bindings() -> RuntimeBindings:
    try:
        import qiskit
        import qiskit_ibm_runtime
        from qiskit import qpy, transpile
        from qiskit_ibm_runtime import SamplerV2
        from qiskit_ibm_runtime.options import SamplerOptions
    except ImportError as exc:
        raise RuntimeSafetyError("pinned Qiskit Runtime dependencies are unavailable") from exc
    runtime_version = importlib.metadata.version("qiskit-ibm-runtime")
    if runtime_version != EXPECTED_RUNTIME_VERSION:
        raise RuntimeSafetyError(
            f"qiskit-ibm-runtime must be exactly {EXPECTED_RUNTIME_VERSION}"
        )
    qiskit_version = qiskit.__version__
    if _parse_version(qiskit_version)[:2] != EXPECTED_QISKIT_MAJOR_MINOR:
        raise RuntimeSafetyError("Qiskit must be from the frozen 2.5 release line")
    if getattr(qiskit_ibm_runtime, "__version__", runtime_version) != runtime_version:
        raise RuntimeSafetyError("Runtime package version metadata is inconsistent")
    return RuntimeBindings(
        SamplerV2=SamplerV2,
        SamplerOptions=SamplerOptions,
        transpile=transpile,
        qpy_dump=qpy.dump,
        qiskit_version=qiskit_version,
        runtime_version=runtime_version,
    )


def load_api_token(credentials_file: Path) -> str:
    try:
        text = credentials_file.read_text(encoding="utf-8").strip()
    except OSError as exc:
        raise RuntimeSafetyError("could not read the IBM credential file") from exc
    match = re.fullmatch(r"IBM cloud API key:\s*(\S+)", text)
    token = match.group(1) if match else text
    if not token or any(character.isspace() for character in token):
        raise RuntimeSafetyError("IBM credential file does not contain one valid token")
    return token


def create_runtime_service(credentials_file: Path) -> Any:
    try:
        from qiskit_ibm_runtime import QiskitRuntimeService
    except ImportError as exc:
        raise RuntimeSafetyError("Qiskit Runtime is unavailable") from exc
    logging.getLogger("qiskit_ibm_runtime").setLevel(logging.WARNING)
    token = load_api_token(credentials_file)
    try:
        return QiskitRuntimeService(channel="ibm_quantum_platform", token=token)
    except Exception as exc:
        raise RuntimeSafetyError("could not initialize the IBM Quantum Platform service") from exc


def _safe_json(value: Any, *, key: str | None = None) -> Any:
    if key is not None and any(fragment in key.lower() for fragment in REDACTED_KEYS):
        return "<redacted>"
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else str(value)
    if isinstance(value, (datetime, date)):
        return value.isoformat()
    if isinstance(value, Mapping):
        return {
            str(item_key): _safe_json(item_value, key=str(item_key))
            for item_key, item_value in value.items()
        }
    if isinstance(value, (list, tuple, set)):
        return [_safe_json(item) for item in value]
    if hasattr(value, "item"):
        try:
            return _safe_json(value.item())
        except (TypeError, ValueError):
            pass
    if hasattr(value, "to_dict"):
        try:
            return _safe_json(value.to_dict())
        except Exception:  # backend vendor objects are not uniform
            pass
    if hasattr(value, "duration") or hasattr(value, "error"):
        return {
            "duration": _safe_json(getattr(value, "duration", None)),
            "error": _safe_json(getattr(value, "error", None)),
        }
    if hasattr(value, "value") and isinstance(value.value, (str, bool, int, float)):
        return _safe_json(value.value)
    return {"unserialized_type": type(value).__name__}


def _backend_name(backend: Any) -> str:
    value = getattr(backend, "name", None)
    value = value() if callable(value) else value
    if not isinstance(value, str):
        raise RuntimeSafetyError("backend did not expose a stable name")
    return value


def _backend_num_qubits(backend: Any) -> int:
    value = getattr(backend, "num_qubits", None)
    if value is None:
        configuration = backend.configuration()
        value = getattr(configuration, "n_qubits", None)
    if type(value) is not int or value <= 0:
        raise RuntimeSafetyError("backend did not expose a valid qubit count")
    return value


def validate_dynamic_backend(backend: Any, expected_name: str, layout: Sequence[int]) -> None:
    if _backend_name(backend) != expected_name:
        raise RuntimeSafetyError("service resolved a different backend than preregistered")
    status = backend.status()
    if not bool(getattr(status, "operational", False)):
        raise RuntimeSafetyError("preregistered backend is not operational")
    configuration = backend.configuration()
    if bool(getattr(configuration, "simulator", False)):
        raise RuntimeSafetyError("preregistered backend unexpectedly resolved to a simulator")
    if max(layout) >= _backend_num_qubits(backend):
        raise RuntimeSafetyError("preregistered layout is outside the backend")

    supported = set(getattr(backend, "supported_instructions", ()) or ())
    supported.update(getattr(configuration, "supported_instructions", ()) or ())
    target = getattr(backend, "target", None)
    if target is None:
        raise RuntimeSafetyError("backend has no compilation target")
    target_names = set(getattr(target, "operation_names", ()) or ())
    advertised = supported | target_names
    if "if_else" not in advertised or "reset" not in advertised or "measure" not in advertised:
        raise RuntimeSafetyError("backend does not advertise required dynamic-circuit operations")


def backend_snapshot(backend: Any) -> dict[str, Any]:
    target = backend.target
    target_rows: dict[str, Any] = {}
    for operation_name in sorted(set(getattr(target, "operation_names", ()) or ())):
        try:
            properties = target[operation_name]
        except Exception:
            target_rows[operation_name] = {"unavailable": True}
            continue
        rows = []
        if isinstance(properties, Mapping):
            for qargs, instruction_properties in sorted(
                properties.items(), key=lambda item: str(item[0])
            ):
                rows.append(
                    {
                        "qargs": None if qargs is None else [int(value) for value in qargs],
                        "properties": _safe_json(instruction_properties),
                    }
                )
        else:
            rows.append({"qargs": None, "properties": _safe_json(properties)})
        target_rows[operation_name] = rows
    target_snapshot = {
        "num_qubits": _backend_num_qubits(backend),
        "dt": _safe_json(getattr(target, "dt", None)),
        "operations": target_rows,
    }
    try:
        configuration = _safe_json(backend.configuration())
    except Exception:
        configuration = {"unavailable": True}
    try:
        calibration = _safe_json(backend.properties())
    except Exception:
        calibration = {"unavailable": True}
    try:
        status = _safe_json(backend.status())
    except Exception:
        status = {"unavailable": True}
    snapshot = {
        "captured_at_utc": utc_now(),
        "backend_name": _backend_name(backend),
        "target": target_snapshot,
        "target_sha256": sha256_json(target_snapshot),
        "calibration": calibration,
        "calibration_sha256": sha256_json(calibration),
        "configuration": configuration,
        "status": status,
    }
    snapshot["snapshot_sha256"] = sha256_json(snapshot)
    return snapshot


def _used_physical_qubits(circuit: Any) -> set[int]:
    used: set[int] = set()
    for instruction in getattr(circuit, "data", ()):
        qargs = getattr(instruction, "qubits", None)
        if qargs is None and isinstance(instruction, tuple) and len(instruction) >= 2:
            qargs = instruction[1]
        for qubit in qargs or ():
            try:
                used.add(int(circuit.find_bit(qubit).index))
            except Exception as exc:
                raise RuntimeSafetyError("could not audit compiled physical qubits") from exc
    return used


def validate_compiled_circuit(logical: Any, compiled: Any, fixed_layout: Sequence[int]) -> None:
    allowed = set(int(value) for value in fixed_layout)
    used = _used_physical_qubits(compiled)
    if not used or not used.issubset(allowed):
        raise RuntimeSafetyError("compiled circuit escaped the preregistered physical layout")
    logical_counts = dict(logical.count_ops())
    compiled_counts = dict(compiled.count_ops())
    for operation in ("if_else", "reset", "measure"):
        if compiled_counts.get(operation, 0) < logical_counts.get(operation, 0):
            raise RuntimeSafetyError("compilation removed a required dynamic operation")


def _delay_seconds(operation: Any, dt: float | None) -> float:
    duration = getattr(operation, "duration", None)
    if duration is None:
        params = getattr(operation, "params", ())
        duration = params[0] if params else None
    unit = getattr(operation, "unit", "dt")
    if type(duration) not in (int, float):
        raise RuntimeSafetyError("compiled delay has no numeric duration")
    factors = {"s": 1.0, "ms": 1e-3, "us": 1e-6, "ns": 1e-9, "ps": 1e-12}
    if unit == "dt":
        if dt is None:
            raise RuntimeSafetyError("backend target has no dt for compiled delay")
        return float(duration) * dt
    if unit not in factors:
        raise RuntimeSafetyError("compiled delay uses an unsupported time unit")
    return float(duration) * factors[unit]


def _target_instruction_duration(target: Any, name: str, qargs: tuple[int, ...]) -> float | None:
    try:
        table = target[name]
        properties = table.get(qargs) if hasattr(table, "get") else None
        duration = getattr(properties, "duration", None)
        if duration is not None:
            result = float(duration)
            return result if math.isfinite(result) and result >= 0 else None
    except Exception:
        pass
    return None


def _explicit_dynamic_latency(backend: Any) -> float:
    candidates = []
    for owner in (backend, backend.configuration()):
        for name in (
            "dynamic_circuit_latency",
            "conditional_latency",
            "classical_feedforward_latency",
        ):
            value = getattr(owner, name, None)
            if type(value) in (int, float) and 0 < float(value) < 1:
                candidates.append(float(value))
    return max(candidates, default=DYNAMIC_CONTROL_LATENCY_BOUND_SECONDS)


def _recursive_conservative_duration(
    circuit: Any,
    target: Any,
    dynamic_latency: float,
    physical_map: Sequence[int] | None = None,
) -> float:
    dt = getattr(target, "dt", None)
    total = 0.0
    for instruction in getattr(circuit, "data", ()):
        operation = getattr(instruction, "operation", None)
        qbits = getattr(instruction, "qubits", None)
        if operation is None or qbits is None:
            try:
                operation, qbits = instruction[0], instruction[1]
            except Exception as exc:
                raise RuntimeSafetyError("compiled instruction has an unexpected structure") from exc
        local_indices = [int(circuit.find_bit(qubit).index) for qubit in qbits]
        qargs = tuple(
            local_indices[index] if physical_map is None else int(physical_map[local_indices[index]])
            for index in range(len(local_indices))
        )
        name = str(getattr(operation, "name", ""))
        blocks = getattr(operation, "blocks", ())
        if name in {"if_else", "while_loop", "for_loop", "switch_case"}:
            if name != "if_else" or not blocks:
                raise RuntimeSafetyError("compiled circuit contains unsupported control flow")
            block_durations = [
                _recursive_conservative_duration(
                    block,
                    target,
                    dynamic_latency,
                    physical_map=qargs,
                )
                for block in blocks
            ]
            # An absent else block has zero duration.  Static instruction
            # durations are summed, rather than parallelized, so this is an
            # intentional upper bound even when the backend scheduler overlaps
            # operations within a branch.
            if len(block_durations) == 1:
                block_durations.append(0.0)
            total += dynamic_latency + max(block_durations)
        elif name == "barrier":
            continue
        elif name == "delay":
            total += _delay_seconds(operation, dt)
        else:
            duration = _target_instruction_duration(target, name, qargs)
            if duration is None:
                # Virtual phase operations have zero physical duration.  Every
                # other missing target duration receives a deliberately large
                # bound instead of making the plan fail on IBM's duration-less
                # control-flow target entries.
                duration = 0.0 if name in {"rz", "p", "u1"} else UNKNOWN_INSTRUCTION_DURATION_BOUND_SECONDS
            total += duration
    return total


def estimate_compiled_duration(circuit: Any, backend: Any) -> float:
    target = backend.target
    duration = _recursive_conservative_duration(
        circuit, target, _explicit_dynamic_latency(backend)
    )
    if not math.isfinite(duration) or duration <= 0:
        raise RuntimeSafetyError("compiled circuit duration estimate is not positive and finite")
    return duration


def pad_compiled_circuit(
    circuit: Any,
    backend: Any,
    fixed_layout: Sequence[int],
    target_duration_seconds: float,
) -> tuple[Any, float]:
    """Insert a pre-terminal-measure delay to reach a common family envelope."""

    current = estimate_compiled_duration(circuit, backend)
    dt = getattr(backend.target, "dt", None)
    if type(dt) not in (int, float) or float(dt) <= 0:
        raise RuntimeSafetyError("backend target has no positive dt for duration padding")
    dt = float(dt)
    deficit = target_duration_seconds - current
    if deficit <= dt:
        return circuit, current
    cycles = int(math.ceil(deficit / dt))
    try:
        from qiskit.circuit import CircuitInstruction, Delay
    except ImportError as exc:
        raise RuntimeSafetyError("Qiskit delay instructions are unavailable") from exc
    padded = circuit.copy()
    pad_qubit_index = int(fixed_layout[0])
    measurement_indices = []
    for index, instruction in enumerate(padded.data):
        operation = getattr(instruction, "operation", None)
        if getattr(operation, "name", None) != "measure":
            continue
        qargs = getattr(instruction, "qubits", ())
        if any(int(padded.find_bit(qubit).index) == pad_qubit_index for qubit in qargs):
            measurement_indices.append(index)
    if not measurement_indices:
        raise RuntimeSafetyError("could not locate terminal measurement for duration padding")
    boundary = measurement_indices[-1]
    delay = CircuitInstruction(
        operation=Delay(cycles, unit="dt"),
        qubits=(padded.qubits[pad_qubit_index],),
        clbits=(),
    )
    try:
        padded.data.insert(boundary, delay)
    except Exception as exc:
        raise RuntimeSafetyError("could not insert common-duration padding") from exc
    validate_compiled_circuit(circuit, padded, fixed_layout)
    padded_duration = estimate_compiled_duration(padded, backend)
    if padded_duration + dt < target_duration_seconds:
        raise RuntimeSafetyError("common-duration padding did not reach its envelope")
    return padded, padded_duration


def _qpy_blob(circuits: Sequence[Any], bindings: RuntimeBindings) -> bytes:
    return qpy_bytes(circuits, bindings.qpy_dump)


def _usage_snapshot(service: Any) -> dict[str, Any]:
    try:
        usage = _require_mapping(service.usage(), "instance usage")
    except Exception as exc:
        if isinstance(exc, RuntimeSafetyError):
            raise
        raise RuntimeSafetyError("could not query IBM instance usage") from exc
    required = {"usage_consumed_seconds", "usage_remaining_seconds", "usage_limit_seconds"}
    if not required.issubset(usage):
        raise RuntimeSafetyError("IBM instance usage response lacks required limits")
    snapshot = _safe_json(usage)
    snapshot["captured_at_utc"] = utc_now()
    snapshot["snapshot_sha256"] = sha256_json(snapshot)
    return snapshot


def guard_qpu_usage(
    usage: Mapping[str, Any], resources: Mapping[str, Any], estimated_qpu_seconds: float
) -> dict[str, float]:
    try:
        actual_limit = float(usage["usage_limit_seconds"])
        remaining = float(usage["usage_remaining_seconds"])
    except (KeyError, TypeError, ValueError) as exc:
        raise RuntimeSafetyError("instance usage values are invalid") from exc
    configured_limit = float(resources["account_quota_seconds"])
    reserve = float(resources["account_reserve_seconds"])
    ceiling = float(resources["estimated_qpu_seconds_ceiling"])
    if bool(usage.get("usage_limit_reached", False)):
        raise RuntimeSafetyError("IBM instance reports that its usage limit is reached")
    if actual_limit != configured_limit:
        raise RuntimeSafetyError("IBM account limit differs from the preregistered limit")
    if estimated_qpu_seconds > ceiling:
        raise RuntimeSafetyError("compiled workload exceeds the preregistered QPU ceiling")
    if remaining < estimated_qpu_seconds + reserve:
        raise RuntimeSafetyError("compiled workload would consume the protected QPU reserve")
    return {
        "actual_limit_seconds": actual_limit,
        "remaining_seconds": remaining,
        "reserve_seconds": reserve,
        "estimated_qpu_seconds": estimated_qpu_seconds,
        "headroom_after_estimate_seconds": remaining - estimated_qpu_seconds,
    }


def _slot_map(bundle: VerifiedBundle) -> dict[str, dict[str, Any]]:
    return {slot["public_role"]: slot for slot in bundle.backend_slots}


def _chunks(items: Sequence[Any], size: int) -> Iterable[Sequence[Any]]:
    for start in range(0, len(items), size):
        yield items[start : start + size]


def prepare_runtime_run(
    bundle: VerifiedBundle,
    service: Any,
    bindings: RuntimeBindings,
    backend_role: str,
) -> PreparedRun:
    manifest = bundle.manifest
    resources = manifest["resources"]
    slot_map = _slot_map(bundle)
    selected_slots = {
        public_role: slot
        for public_role, slot in slot_map.items()
        if slot["role"] == backend_role
    }
    if len(selected_slots) != 1:
        raise RuntimeSafetyError("requested backend role does not resolve uniquely")
    backends: dict[str, Any] = {}
    snapshots: dict[str, dict[str, Any]] = {}
    for public_role, slot in selected_slots.items():
        try:
            backend = service.backend(slot["backend"])
        except Exception as exc:
            raise RuntimeSafetyError("could not resolve a preregistered backend") from exc
        validate_dynamic_backend(backend, slot["backend"], slot["layout"])
        backends[public_role] = backend
        snapshots[backend_role] = backend_snapshot(backend)

    prepared: list[PreparedCircuit] = []
    for public, logical in bundle.circuits:
        role = public["backend_role"]
        if role not in selected_slots:
            continue
        backend = backends[role]
        layout = slot_map[role]["layout"]
        try:
            compiled = bindings.transpile(
                logical,
                backend=backend,
                optimization_level=resources["transpile_optimization_level"],
                seed_transpiler=resources["transpiler_seed"],
                initial_layout=list(layout),
            )
        except Exception as exc:
            raise RuntimeSafetyError("fixed-layout compilation failed") from exc
        validate_compiled_circuit(logical, compiled, layout)
        compiled_blob = _qpy_blob([compiled], bindings)
        prepared.append(
            PreparedCircuit(
                opaque_id=public["opaque_id"],
                family=bundle.circuit_families[public["opaque_id"]],
                shots=public["shots"],
                backend_role=role,
                logical_circuit_sha256=public["logical_circuit_sha256"],
                compiled_qpy_sha256=sha256_bytes(compiled_blob),
                compiled_duration_seconds=estimate_compiled_duration(compiled, backend),
                logical_circuit=logical,
                compiled_circuit=compiled,
            )
        )

    if not prepared:
        raise RuntimeSafetyError("requested backend role has no preregistered circuits")
    duration_families: OrderedDict[tuple[str, str], list[PreparedCircuit]] = OrderedDict()
    for item in prepared:
        duration_families.setdefault((item.backend_role, item.family), []).append(item)
    padded_by_id: dict[str, PreparedCircuit] = {}
    for (public_role, _family), items in duration_families.items():
        envelope = max(item.compiled_duration_seconds for item in items)
        backend = backends[public_role]
        layout = slot_map[public_role]["layout"]
        durations = []
        for item in items:
            padded, duration = pad_compiled_circuit(
                item.compiled_circuit, backend, layout, envelope
            )
            compiled_blob = _qpy_blob([padded], bindings)
            replacement = dataclasses.replace(
                item,
                compiled_circuit=padded,
                compiled_qpy_sha256=sha256_bytes(compiled_blob),
                compiled_duration_seconds=duration,
            )
            padded_by_id[item.opaque_id] = replacement
            durations.append(duration)
        dt = float(getattr(backend.target, "dt", 0.0) or 0.0)
        if max(durations) - min(durations) > dt + 1e-12:
            raise RuntimeSafetyError("compiled circuits do not share a common duration envelope")
    prepared = [padded_by_id[item.opaque_id] for item in prepared]

    grouped: OrderedDict[tuple[str, str, int], list[PreparedCircuit]] = OrderedDict()
    for item in prepared:
        grouped.setdefault((item.backend_role, item.family, item.shots), []).append(item)
    groups: list[PreparedGroup] = []
    max_pubs = resources["max_pubs_per_job"]
    for (role, family, shots), items in grouped.items():
        for chunk in _chunks(items, max_pubs):
            identifiers = [item.opaque_id for item in chunk]
            group_id = sha256_json(
                {
                    "manifest_sha256": bundle.manifest_sha256,
                    "backend_role": role,
                    "family": family,
                    "shots": shots,
                    "opaque_ids": identifiers,
                }
            )
            backend = backends[role]
            rep_delay = float(getattr(backend, "default_rep_delay", 0.00025) or 0.00025)
            raw_estimate = IBM_SUBJOB_OVERHEAD_SECONDS + sum(
                (rep_delay + item.compiled_duration_seconds) * shots for item in chunk
            )
            guarded_estimate = raw_estimate * ESTIMATE_SAFETY_FACTOR
            compiled_blob = _qpy_blob([item.compiled_circuit for item in chunk], bindings)
            groups.append(
                PreparedGroup(
                    group_id=group_id,
                    backend_role=slot_map[role]["role"],
                    backend_role_opaque_id=role,
                    family=family,
                    backend_name=slot_map[role]["backend"],
                    layout_opaque_id=slot_map[role]["public_layout"],
                    physical_layout=tuple(slot_map[role]["layout"]),
                    properties_last_update=slot_map[role]["properties_last_update"],
                    shots=shots,
                    estimated_qpu_seconds=guarded_estimate,
                    circuits=tuple(chunk),
                    compiled_qpy=compiled_blob,
                )
            )

    estimate = sum(group.estimated_qpu_seconds for group in groups)
    usage = _usage_snapshot(service)
    guard = guard_qpu_usage(usage, resources, estimate)
    plan_groups = []
    for group in groups:
        compiled_artifact = compiled_qpy_artifact(group.group_id, group.compiled_qpy)
        plan_groups.append(
            {
                "group_id": group.group_id,
                "backend_role": group.backend_role,
                "backend_role_opaque_id": group.backend_role_opaque_id,
                "family": group.family,
                "backend_name": group.backend_name,
                "layout_id": group.layout_opaque_id,
                "physical_layout": list(group.physical_layout),
                "properties_last_update": group.properties_last_update,
                "shots": group.shots,
                "opaque_ids": [item.opaque_id for item in group.circuits],
                "logical_circuit_sha256_by_opaque_id": {
                    item.opaque_id: item.logical_circuit_sha256
                    for item in group.circuits
                },
                "compiled_qpy_sha256_by_opaque_id": {
                    item.opaque_id: item.compiled_qpy_sha256
                    for item in group.circuits
                },
                "compiled_qpy_bundle_sha256": compiled_artifact["sha256"],
                "compiled_qpy_artifact": compiled_artifact,
                "compiled_duration_seconds_by_opaque_id": {
                    item.opaque_id: item.compiled_duration_seconds for item in group.circuits
                },
                "estimated_qpu_seconds": group.estimated_qpu_seconds,
            }
        )
    plan = {
        "schema_version": RUNTIME_SCHEMA,
        "artifact_type": "pre_submission_dry_run",
        "created_at_utc": utc_now(),
        "submission_performed": False,
        "operator_source_sha256": operator_source_sha256(),
        "manifest_sha256": bundle.manifest_sha256,
        "analysis_lock_sha256": bundle.analysis_lock_sha256,
        "public_manifest_core_sha256": bundle.public_manifest_core_sha256,
        "private_reveal_commitment_verified": True,
        "logical_circuit_hashes_verified": True,
        "ideal_recipe_validation": bundle.ideal_validation,
        "hardened_analysis_lock_verified": hardened_analysis_lock_status(bundle),
        "selected_backend_role": backend_role,
        "runtime": {
            "qiskit_version": bindings.qiskit_version,
            "qiskit_ibm_runtime_version": bindings.runtime_version,
            "primitive": "SamplerV2",
            "execution_mode": "backend_job",
            "logical_circuit_serialization": LOGICAL_CIRCUIT_SERIALIZATION,
            "sessions": False,
            "batches": False,
            "dynamical_decoupling": False,
            "gate_twirling": False,
            "measurement_twirling": False,
            "measurement_type": "classified",
            "init_qubits_each_shot": True,
        },
        "resource_estimate": {
            "ibm_subjob_overhead_seconds": IBM_SUBJOB_OVERHEAD_SECONDS,
            "safety_factor": ESTIMATE_SAFETY_FACTOR,
            "estimated_qpu_seconds": estimate,
            "guard": guard,
        },
        "instance_usage_before_submission": usage,
        "backend_snapshots": snapshots,
        "groups": plan_groups,
    }
    plan["plan_sha256"] = sha256_json(plan)
    return PreparedRun(plan=plan, groups=tuple(groups), backends=backends)


def validate_submission_confirmation(
    *,
    confirm_submit: bool,
    confirmed_manifest_sha256: str | None,
    confirmed_analysis_lock_sha256: str | None,
    bundle: VerifiedBundle,
    analysis_validator: Callable[..., Any] | None = None,
) -> None:
    if not confirm_submit:
        raise RuntimeSafetyError("submission requires the explicit --confirm-submit flag")
    if confirmed_manifest_sha256 is None or confirmed_analysis_lock_sha256 is None:
        raise RuntimeSafetyError("submission requires literal manifest and analysis-lock digests")
    _require_hex64(confirmed_manifest_sha256, "confirmed manifest digest")
    _require_hex64(confirmed_analysis_lock_sha256, "confirmed analysis-lock digest")
    if confirmed_manifest_sha256 != bundle.manifest_sha256:
        raise RuntimeSafetyError("confirmed manifest digest does not match")
    if confirmed_analysis_lock_sha256 != bundle.analysis_lock_sha256:
        raise RuntimeSafetyError("confirmed analysis-lock digest does not match")
    if bundle.reveal_mode != "production_random":
        raise RuntimeSafetyError("submission refuses deterministic test preregistrations")
    validate_hardened_analysis_lock(bundle, analysis_validator=analysis_validator)


def validate_hardened_analysis_lock(
    bundle: VerifiedBundle,
    *,
    analysis_validator: Callable[..., Any] | None = None,
) -> None:
    if analysis_validator is None:
        try:
            from cayley_blind_likelihood_analysis import validate_analysis_lock
        except ImportError as exc:
            raise RuntimeSafetyError("hardened blinded analysis validator is unavailable") from exc
        analysis_validator = validate_analysis_lock
    try:
        analysis_validator(bundle.analysis_lock_document, verify_code_hash=True)
    except Exception as exc:
        raise RuntimeSafetyError("hardened blinded analysis lock does not verify") from exc
    if (
        bundle.analysis_lock_document.get("catalog_precommitment_sha256")
        != bundle.manifest["catalog_precommitment_sha256"]
    ):
        raise RuntimeSafetyError("hardened analysis lock binds a different catalog")


def hardened_analysis_lock_status(bundle: VerifiedBundle) -> bool:
    try:
        validate_hardened_analysis_lock(bundle)
    except RuntimeSafetyError:
        return False
    return True


def _new_sampler_options(bindings: RuntimeBindings, resources: Mapping[str, Any], group: PreparedGroup) -> Any:
    options = bindings.SamplerOptions()
    options.dynamical_decoupling.enable = False
    options.twirling.enable_gates = False
    options.twirling.enable_measure = False
    options.execution.meas_type = "classified"
    options.execution.init_qubits = True
    # IBM accepts [1, 10800].  Use the preregistered 420-second hard cap,
    # rather than a tiny ceil(estimate) value that can reject a valid job due
    # to control-electronics overhead not represented in local timing.
    maximum = int(resources["max_execution_time_seconds"])
    if maximum not in range(1, 10801):
        raise RuntimeSafetyError("max_execution_time_seconds is outside IBM's accepted range")
    options.max_execution_time = maximum
    tag = f"oph509-{group.group_id[:20]}"
    try:
        options.environment.job_tags = ["oph509", tag]
    except AttributeError as exc:
        raise RuntimeSafetyError("SamplerOptions does not support deterministic job tags") from exc
    if (
        bool(options.dynamical_decoupling.enable)
        or bool(options.twirling.enable_gates)
        or bool(options.twirling.enable_measure)
        or options.execution.meas_type != "classified"
        or options.execution.init_qubits is not True
    ):
        raise RuntimeSafetyError("Sampler options unexpectedly enabled suppression or twirling")
    return options


def _event_hash(event: Mapping[str, Any]) -> str:
    return sha256_json(_without(event, "event_sha256"))


def append_event(path: Path, event: Mapping[str, Any]) -> dict[str, Any]:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(path, os.O_RDWR | os.O_CREAT, 0o600)
    try:
        with os.fdopen(descriptor, "r+", encoding="utf-8", closefd=False) as handle:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
            handle.seek(0)
            lines = [line for line in handle.read().splitlines() if line.strip()]
            previous = None
            for line in lines:
                try:
                    prior = json.loads(line)
                except json.JSONDecodeError as exc:
                    raise RuntimeSafetyError("operator event journal is corrupt") from exc
                if prior.get("previous_event_sha256") != previous:
                    raise RuntimeSafetyError("operator event journal hash chain is corrupt")
                observed = _require_hex64(prior.get("event_sha256"), "prior event digest")
                if observed != _event_hash(prior):
                    raise RuntimeSafetyError("operator event journal hash chain is corrupt")
                previous = observed
            row = dict(event)
            row["previous_event_sha256"] = previous
            row["event_sha256"] = _event_hash(row)
            handle.seek(0, os.SEEK_END)
            handle.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
            handle.flush()
            os.fsync(handle.fileno())
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
            return row
    finally:
        os.close(descriptor)


def write_new_bytes(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    try:
        view = memoryview(payload)
        while view:
            written = os.write(descriptor, view)
            view = view[written:]
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def write_new_json(path: Path, value: Any) -> None:
    write_new_bytes(path, json.dumps(value, indent=2, sort_keys=True).encode("utf-8") + b"\n")


def write_or_verify_bytes(path: Path, payload: bytes) -> None:
    if path.exists():
        try:
            existing = path.read_bytes()
        except OSError as exc:
            raise RuntimeSafetyError("could not verify an existing operator artifact") from exc
        if existing != payload:
            raise RuntimeSafetyError("existing operator artifact differs; refusing overwrite")
        return
    write_new_bytes(path, payload)


def persist_prepared_run(prepared: PreparedRun, output_dir: Path) -> Path:
    claimed_plan_sha = _require_hex64(
        prepared.plan.get("plan_sha256"), "pre-submission plan_sha256"
    )
    if sha256_json(_without(prepared.plan, "plan_sha256")) != claimed_plan_sha:
        raise RuntimeSafetyError("pre-submission plan digest does not verify")
    role = prepared.plan["selected_backend_role"]
    plan_path = output_dir / (
        f"runtime_{role}_pre_submission_plan_{claimed_plan_sha}.json"
    )
    plan_groups = prepared.plan.get("groups")
    if not isinstance(plan_groups, list):
        raise RuntimeSafetyError("pre-submission plan lacks compiled artifact bindings")
    plan_by_group = {
        row.get("group_id"): row for row in plan_groups if isinstance(row, Mapping)
    }
    if len(plan_by_group) != len(plan_groups):
        raise RuntimeSafetyError("pre-submission plan has duplicate or malformed groups")
    prepared_group_ids = {group.group_id for group in prepared.groups}
    if set(plan_by_group) != prepared_group_ids:
        raise RuntimeSafetyError("prepared groups differ from the pre-submission plan")
    for group in prepared.groups:
        planned = plan_by_group.get(group.group_id)
        if planned is None:
            raise RuntimeSafetyError("prepared group is absent from the pre-submission plan")
        expected_artifact = compiled_qpy_artifact(group.group_id, group.compiled_qpy)
        if (
            planned.get("compiled_qpy_artifact") != expected_artifact
            or planned.get("compiled_qpy_bundle_sha256") != expected_artifact["sha256"]
        ):
            raise RuntimeSafetyError(
                "compiled QPY bytes differ from their pre-submission artifact binding"
            )
        write_or_verify_bytes(
            output_dir / expected_artifact["relative_path"], group.compiled_qpy
        )
    plan_payload = (
        json.dumps(prepared.plan, indent=2, sort_keys=True).encode("utf-8") + b"\n"
    )
    write_or_verify_bytes(plan_path, plan_payload)
    return plan_path


def _job_status(job: Any) -> str:
    try:
        status = job.status()
    except Exception:
        return "STATUS_UNAVAILABLE"
    value = getattr(status, "value", status)
    return str(value).upper()


def _job_usage_estimation(job: Any) -> Any:
    try:
        value = job.usage_estimation
        return _safe_json(value() if callable(value) else value)
    except Exception:
        return {"unavailable": True}


def submit_prepared_run(
    prepared: PreparedRun,
    bundle: VerifiedBundle,
    service: Any,
    bindings: RuntimeBindings,
    output_dir: Path,
    submission_journal: Path | None = None,
) -> list[dict[str, Any]]:
    planned_operator_source = _require_hex64(
        prepared.plan.get("operator_source_sha256"), "operator_source_sha256"
    )
    if planned_operator_source != operator_source_sha256():
        raise RuntimeSafetyError(
            "Runtime operator source differs from the pre-submission plan"
        )
    journal = submission_journal or (
        output_dir / f"runtime_{prepared.plan['selected_backend_role']}_submission_events.ndjson"
    )
    resources = bundle.manifest["resources"]
    registered: list[dict[str, Any]] = []
    existing_events = read_event_journal(journal, bundle.manifest_sha256) if journal.exists() else []
    existing_registrations = [
        event for event in existing_events if event.get("event_type") == "job_registered"
    ]
    registered_by_group = {event["group_id"]: event for event in existing_registrations}
    valid_group_ids = {group.group_id for group in prepared.groups}
    if not set(registered_by_group).issubset(valid_group_ids):
        raise RuntimeSafetyError("submission journal contains groups outside the current role plan")
    ambiguous = [
        event
        for event in existing_events
        if event.get("event_type") == "submission_failed_without_registered_job"
    ]
    if ambiguous:
        raise RuntimeSafetyError("journal has an ambiguous submission; recover it by deterministic tag")
    started_groups = {
        event["group_id"]
        for event in existing_events
        if event.get("event_type") == "submission_started"
    }
    if started_groups.difference(registered_by_group):
        raise RuntimeSafetyError("journal has an unregistered started submission; recover by job tag")
    for registration in existing_registrations:
        try:
            prior_job = service.job(registration["job_id"])
        except Exception as exc:
            raise RuntimeSafetyError("could not verify a previously registered job") from exc
        prior_status = _job_status(prior_job)
        append_event(
            journal,
            {
                "schema_version": RUNTIME_SCHEMA,
                "event_type": "pre_next_submission_status_check",
                "timestamp_utc": utc_now(),
                "manifest_sha256": bundle.manifest_sha256,
                "operator_source_sha256": planned_operator_source,
                "group_id": registration["group_id"],
                "job_id": registration["job_id"],
                "status": prior_status,
                "metrics": _job_metrics(prior_job),
            },
        )
        if prior_status != "DONE":
            if prior_status in {"ERROR", "CANCELLED"}:
                raise RuntimeSafetyError("a prior role-scoped job did not complete successfully")
            return []

    pending_groups = [
        group for group in prepared.groups if group.group_id not in registered_by_group
    ]
    if not pending_groups:
        return []
    group = pending_groups[0]
    remaining_estimate = sum(item.estimated_qpu_seconds for item in pending_groups)
    usage = _usage_snapshot(service)
    guard_qpu_usage(usage, resources, remaining_estimate)
    options = _new_sampler_options(bindings, resources, group)
    tag = f"oph509-{group.group_id[:20]}"
    append_event(
        journal,
        {
            "schema_version": RUNTIME_SCHEMA,
            "event_type": "submission_started",
            "timestamp_utc": utc_now(),
            "manifest_sha256": bundle.manifest_sha256,
            "analysis_lock_sha256": bundle.analysis_lock_sha256,
            "operator_source_sha256": planned_operator_source,
            "plan_sha256": prepared.plan["plan_sha256"],
            "group_id": group.group_id,
            "backend_name": group.backend_name,
            "layout_id": group.layout_opaque_id,
            "backend_role": group.backend_role,
            "backend_role_opaque_id": group.backend_role_opaque_id,
            "family": group.family,
            "physical_layout": list(group.physical_layout),
            "properties_last_update": group.properties_last_update,
            "shots": group.shots,
            "opaque_ids": [item.opaque_id for item in group.circuits],
            "circuit_bindings": [
                {
                    "opaque_id": item.opaque_id,
                    "logical_circuit_sha256": item.logical_circuit_sha256,
                    "compiled_qpy_sha256": item.compiled_qpy_sha256,
                }
                for item in group.circuits
            ],
            "job_tag": tag,
            "estimated_qpu_seconds": group.estimated_qpu_seconds,
            "usage_before_submission": usage,
        },
    )
    try:
        sampler = bindings.SamplerV2(
            mode=prepared.backends[group.backend_role_opaque_id], options=options
        )
        job = sampler.run(
            [item.compiled_circuit for item in group.circuits], shots=group.shots
        )
        job_id = job.job_id()
        if not isinstance(job_id, str) or not job_id:
            raise RuntimeSafetyError("Runtime returned an invalid job identifier")
    except Exception as exc:
        append_event(
            journal,
            {
                "schema_version": RUNTIME_SCHEMA,
                "event_type": "submission_failed_without_registered_job",
                "timestamp_utc": utc_now(),
                "manifest_sha256": bundle.manifest_sha256,
                "operator_source_sha256": planned_operator_source,
                "group_id": group.group_id,
                "job_tag": tag,
                "failure_type": type(exc).__name__,
                "operator_action": "recover by deterministic job tag; do not resubmit blindly",
            },
        )
        raise RuntimeSafetyError(
            "submission failed before a job ID was durably registered; recover by job tag"
        ) from exc
    event = append_event(
        journal,
        {
            "schema_version": RUNTIME_SCHEMA,
            "event_type": "job_registered",
            "timestamp_utc": utc_now(),
            "manifest_sha256": bundle.manifest_sha256,
            "analysis_lock_sha256": bundle.analysis_lock_sha256,
            "operator_source_sha256": planned_operator_source,
            "plan_sha256": prepared.plan["plan_sha256"],
            "group_id": group.group_id,
            "job_id": job_id,
            "job_tag": tag,
            "backend_name": group.backend_name,
            "layout_id": group.layout_opaque_id,
            "backend_role": group.backend_role,
            "backend_role_opaque_id": group.backend_role_opaque_id,
            "family": group.family,
            "physical_layout": list(group.physical_layout),
            "properties_last_update": group.properties_last_update,
            "shots": group.shots,
            "opaque_ids": [item.opaque_id for item in group.circuits],
            "circuit_bindings": [
                {
                    "opaque_id": item.opaque_id,
                    "logical_circuit_sha256": item.logical_circuit_sha256,
                    "compiled_qpy_sha256": item.compiled_qpy_sha256,
                }
                for item in group.circuits
            ],
            "status": _job_status(job),
            "usage_estimation": _job_usage_estimation(job),
            "estimated_qpu_seconds": group.estimated_qpu_seconds,
        },
    )
    registered.append(event)
    return registered


def read_event_journal(path: Path, manifest_sha256: str) -> list[dict[str, Any]]:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise RuntimeSafetyError("could not read submission journal") from exc
    previous = None
    events: list[dict[str, Any]] = []
    for line in lines:
        if not line.strip():
            continue
        try:
            event = _require_mapping(json.loads(line), "submission event")
        except json.JSONDecodeError as exc:
            raise RuntimeSafetyError("submission journal contains invalid JSON") from exc
        if event.get("previous_event_sha256") != previous or event.get("event_sha256") != _event_hash(event):
            raise RuntimeSafetyError("submission journal hash chain does not verify")
        previous = event["event_sha256"]
        if event.get("manifest_sha256") != manifest_sha256:
            raise RuntimeSafetyError("submission journal contains a different manifest digest")
        events.append(event)
    return events


def read_registered_jobs(path: Path, manifest_sha256: str) -> list[dict[str, Any]]:
    events = read_event_journal(path, manifest_sha256)
    registered: list[dict[str, Any]] = []
    seen_ids: set[str] = set()
    seen_groups: set[str] = set()
    for event in events:
        if event.get("event_type") != "job_registered":
            continue
        job_id = event.get("job_id")
        group_id = event.get("group_id")
        if (
            not isinstance(job_id, str)
            or not job_id
            or job_id in seen_ids
            or not isinstance(group_id, str)
            or group_id in seen_groups
        ):
            raise RuntimeSafetyError("submission journal has duplicate or invalid jobs")
        seen_ids.add(job_id)
        seen_groups.add(group_id)
        registered.append(event)
    if not registered:
        raise RuntimeSafetyError("submission journal has no registered jobs")
    return registered


def _counts(bit_array: Any) -> dict[str, int]:
    try:
        counts = bit_array.get_counts()
    except Exception as exc:
        raise RuntimeSafetyError("Sampler result does not expose raw counts") from exc
    if not isinstance(counts, Mapping):
        raise RuntimeSafetyError("Sampler raw counts have an unexpected type")
    return {str(key): int(value) for key, value in counts.items()}


def extract_pub_result(
    pub: Any,
    opaque_id: str,
    circuit_binding: Mapping[str, Any] | None = None,
    *,
    family: str | None = None,
    expected_shots: int | None = None,
) -> dict[str, Any]:
    register_counts: dict[str, dict[str, int]] = {}
    data = getattr(pub, "data", None)
    if data is None:
        raise RuntimeSafetyError("Sampler PUB result has no classical data")
    keys_method = getattr(data, "keys", None)
    if not callable(keys_method):
        raise RuntimeSafetyError("Sampler PUB result has no named classical registers")
    register_order = [str(key) for key in keys_method()]
    for key in register_order:
        register_counts[str(key)] = _counts(getattr(data, key))
    try:
        joined = pub.join_data()
    except Exception as exc:
        raise RuntimeSafetyError("Sampler PUB result cannot reconstruct joint counts") from exc
    joint_counts = _counts(joined)
    layout_by_family = {
        "cayley": ({"heated", "decision", "final"}, 7, "final[3]|decision[1]|heated[3]"),
        "mh": ({"before", "accepted", "after"}, 7, "after[3]|accepted[1]|before[3]"),
        "readout_calibration": ({"calibration_c"}, 4, "calibration_c[4]"),
    }
    joined_layout = None
    joined_width = None
    if family is not None:
        if family not in layout_by_family:
            raise RuntimeSafetyError("result references an unsupported circuit family")
        expected_registers, joined_width, joined_layout = layout_by_family[family]
        if set(register_order) != expected_registers:
            raise RuntimeSafetyError("Sampler result classical registers differ from frozen layout")
        if any(
            len(str(key).replace(" ", "")) != joined_width for key in joint_counts
        ):
            raise RuntimeSafetyError("Sampler joined count width differs from frozen layout")
    if expected_shots is not None:
        if sum(joint_counts.values()) != expected_shots or any(
            sum(counts.values()) != expected_shots for counts in register_counts.values()
        ):
            raise RuntimeSafetyError("Sampler result shot totals differ from registration")
    result = {
        "opaque_id": opaque_id,
        "joint_counts": joint_counts,
        "raw_joined_counts_sha256": sha256_json(joint_counts),
        "register_counts": register_counts,
        "classical_register_order": register_order,
        "joined_count_layout": joined_layout,
        "joined_count_width": joined_width,
        "metadata": _safe_json(getattr(pub, "metadata", None)),
    }
    if circuit_binding is not None:
        if circuit_binding.get("opaque_id") != opaque_id:
            raise RuntimeSafetyError("result order differs from registered circuit bindings")
        result["logical_circuit_sha256"] = circuit_binding[
            "logical_circuit_sha256"
        ]
        result["compiled_qpy_sha256"] = circuit_binding["compiled_qpy_sha256"]
    return result


def _job_metrics(job: Any) -> Any:
    try:
        return _safe_json(job.metrics())
    except Exception:
        return {"unavailable": True}


def _job_calibration_snapshot(job: Any) -> dict[str, Any]:
    try:
        properties = _safe_json(job.properties(refresh=False))
    except Exception:
        properties = {"unavailable": True}
    return {
        "properties": properties,
        "properties_sha256": sha256_json(properties),
    }


def harvest_registered_jobs(
    *,
    service: Any,
    bundle: VerifiedBundle,
    submission_journal: Path,
    output_dir: Path,
) -> list[dict[str, Any]]:
    registered = read_registered_jobs(submission_journal, bundle.manifest_sha256)
    harvest_journal = output_dir / "runtime_harvest_events.ndjson"
    harvester_source_sha = operator_source_sha256()
    events: list[dict[str, Any]] = []
    for registration in registered:
        submission_operator_source = _require_hex64(
            registration.get("operator_source_sha256"),
            "registered operator_source_sha256",
        )
        job_id = registration["job_id"]
        try:
            job = service.job(job_id)
        except Exception as exc:
            event = append_event(
                harvest_journal,
                {
                    "schema_version": RUNTIME_SCHEMA,
                    "event_type": "job_lookup_failed",
                    "timestamp_utc": utc_now(),
                    "manifest_sha256": bundle.manifest_sha256,
                    "operator_source_sha256": submission_operator_source,
                    "harvester_source_sha256": harvester_source_sha,
                    "job_id": job_id,
                    "group_id": registration["group_id"],
                    "failure_type": type(exc).__name__,
                },
            )
            events.append(event)
            continue
        status = _job_status(job)
        backend = job.backend()
        if _backend_name(backend) != registration["backend_name"]:
            raise RuntimeSafetyError("registered job resolved to a different backend")
        event_data: dict[str, Any] = {
            "schema_version": RUNTIME_SCHEMA,
            "event_type": "job_harvest",
            "timestamp_utc": utc_now(),
            "manifest_sha256": bundle.manifest_sha256,
            "analysis_lock_sha256": bundle.analysis_lock_sha256,
            "operator_source_sha256": submission_operator_source,
            "harvester_source_sha256": harvester_source_sha,
            "plan_sha256": registration["plan_sha256"],
            "submission_event_sha256": registration["event_sha256"],
            "job_id": job_id,
            "group_id": registration["group_id"],
            "backend_name": registration["backend_name"],
            "layout_id": registration["layout_id"],
            "backend_role": registration["backend_role"],
            "backend_role_opaque_id": registration["backend_role_opaque_id"],
            "family": registration["family"],
            "physical_layout": registration["physical_layout"],
            "properties_last_update": registration["properties_last_update"],
            "shots": registration["shots"],
            "opaque_ids": registration["opaque_ids"],
            "circuit_bindings": registration["circuit_bindings"],
            "status": status,
            "metrics": _job_metrics(job),
            "usage_estimation": _job_usage_estimation(job),
            "calibration_at_execution": _job_calibration_snapshot(job),
            "results": None,
        }
        if status == "DONE":
            try:
                result = job.result()
                if len(result) != len(registration["opaque_ids"]):
                    raise RuntimeSafetyError("Sampler result PUB count differs from registration")
                event_data["results"] = [
                    extract_pub_result(
                        pub,
                        opaque_id,
                        binding,
                        family=registration["family"],
                        expected_shots=registration["shots"],
                    )
                    for pub, opaque_id, binding in zip(
                        result,
                        registration["opaque_ids"],
                        registration["circuit_bindings"],
                    )
                ]
            except RuntimeSafetyError:
                raise
            except Exception as exc:
                raise RuntimeSafetyError("could not harvest a completed Sampler result") from exc
        event = append_event(harvest_journal, event_data)
        events.append(event)
    usage = _usage_snapshot(service)
    append_event(
        harvest_journal,
        {
            "schema_version": RUNTIME_SCHEMA,
            "event_type": "harvest_usage_snapshot",
            "timestamp_utc": utc_now(),
            "manifest_sha256": bundle.manifest_sha256,
            "harvester_source_sha256": harvester_source_sha,
            "instance_usage": usage,
        },
    )
    return events


def require_completed_role_journal(
    service: Any,
    journal: Path,
    manifest_sha256: str,
    expected_group_ids: set[str] | None = None,
) -> None:
    registrations = read_registered_jobs(journal, manifest_sha256)
    if expected_group_ids is not None and {
        registration["group_id"] for registration in registrations
    } != expected_group_ids:
        raise RuntimeSafetyError("prerequisite development role has unsubmitted groups")
    for registration in registrations:
        try:
            job = service.job(registration["job_id"])
        except Exception as exc:
            raise RuntimeSafetyError("could not verify prerequisite role job") from exc
        if _job_status(job) != "DONE":
            raise RuntimeSafetyError("prerequisite development role is not fully complete")


def expected_role_group_ids(bundle: VerifiedBundle, backend_role: str) -> set[str]:
    slot = next(
        (item for item in bundle.backend_slots if item["role"] == backend_role), None
    )
    if slot is None:
        raise RuntimeSafetyError("requested prerequisite role has no committed backend slot")
    rows = [
        row
        for row in bundle.manifest["circuits"]
        if row["backend_role"] == slot["public_role"]
    ]
    grouped: OrderedDict[tuple[str, int], list[dict[str, Any]]] = OrderedDict()
    for row in rows:
        grouped.setdefault(
            (bundle.circuit_families[row["opaque_id"]], row["shots"]), []
        ).append(row)
    result = set()
    for (family, shots), family_rows in grouped.items():
        for chunk in _chunks(
            family_rows, bundle.manifest["resources"]["max_pubs_per_job"]
        ):
            result.add(
                sha256_json(
                    {
                        "manifest_sha256": bundle.manifest_sha256,
                        "backend_role": slot["public_role"],
                        "family": family,
                        "shots": shots,
                        "opaque_ids": [row["opaque_id"] for row in chunk],
                    }
                )
            )
    return result


def _default_credentials_path() -> Path:
    return Path(__file__).resolve().parents[1] / ".api-key"


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plan, explicitly submit, or harvest the blinded record-gated IBM benchmark."
    )
    parser.add_argument("action", nargs="?", choices=("plan", "submit", "harvest"), default="plan")
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--reveal", type=Path, required=True)
    parser.add_argument("--analysis-lock", type=Path, required=True)
    parser.add_argument("--credentials-file", type=Path, default=_default_credentials_path())
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--backend-role", choices=("development", "heldout"), default="development"
    )
    parser.add_argument("--submission-journal", type=Path)
    parser.add_argument("--confirm-submit", action="store_true")
    parser.add_argument("--confirm-manifest-sha256")
    parser.add_argument("--confirm-analysis-lock-sha256")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        bindings = load_runtime_bindings()
        bundle = verify_operator_bundle(
            args.manifest,
            args.reveal,
            args.analysis_lock,
        )
        if args.action == "submit":
            validate_submission_confirmation(
                confirm_submit=args.confirm_submit,
                confirmed_manifest_sha256=args.confirm_manifest_sha256,
                confirmed_analysis_lock_sha256=args.confirm_analysis_lock_sha256,
                bundle=bundle,
            )
        service = create_runtime_service(args.credentials_file)
        if args.action in {"plan", "submit"}:
            prepared = prepare_runtime_run(bundle, service, bindings, args.backend_role)
            plan_path = persist_prepared_run(prepared, args.output_dir)
            if args.action == "plan":
                print(
                    json.dumps(
                        {
                            "action": "plan",
                            "backend_role": args.backend_role,
                            "submission_performed": False,
                            "manifest_sha256": bundle.manifest_sha256,
                            "analysis_lock_sha256": bundle.analysis_lock_sha256,
                            "estimated_qpu_seconds": prepared.plan["resource_estimate"][
                                "estimated_qpu_seconds"
                            ],
                            "plan_path": str(plan_path),
                        },
                        sort_keys=True,
                    )
                )
                return 0
            role_journal = (
                args.output_dir
                / f"runtime_{args.backend_role}_submission_events.ndjson"
            )
            if args.backend_role == "heldout":
                require_completed_role_journal(
                    service,
                    args.output_dir / "runtime_development_submission_events.ndjson",
                    bundle.manifest_sha256,
                    expected_role_group_ids(bundle, "development"),
                )
            jobs = submit_prepared_run(
                prepared,
                bundle,
                service,
                bindings,
                args.output_dir,
                submission_journal=role_journal,
            )
            print(
                json.dumps(
                    {
                        "action": "submit",
                        "backend_role": args.backend_role,
                        "manifest_sha256": bundle.manifest_sha256,
                        "registered_job_ids": [event["job_id"] for event in jobs],
                        "one_group_attempted": bool(jobs),
                        "rerun_required_for_next_group": bool(jobs),
                        "submission_journal": str(role_journal),
                    },
                    sort_keys=True,
                )
            )
            return 0

        journal = args.submission_journal or (
            args.output_dir / f"runtime_{args.backend_role}_submission_events.ndjson"
        )
        events = harvest_registered_jobs(
            service=service,
            bundle=bundle,
            submission_journal=journal,
            output_dir=args.output_dir,
        )
        print(
            json.dumps(
                {
                    "action": "harvest",
                    "backend_role": args.backend_role,
                    "manifest_sha256": bundle.manifest_sha256,
                    "job_statuses": {
                        event.get("job_id", "lookup-failed"): event.get("status", "LOOKUP_FAILED")
                        for event in events
                    },
                    "harvest_journal": str(
                        args.output_dir / "runtime_harvest_events.ndjson"
                    ),
                },
                sort_keys=True,
            )
        )
        return 0
    except RuntimeSafetyError as exc:
        print(f"runtime safety refusal: {exc}", file=sys.stderr)
        return 2
    except Exception as exc:
        # Do not render arbitrary vendor exception text: some HTTP/client
        # exceptions can embed request details.  The class name is sufficient
        # for a safe local diagnostic without risking credential disclosure.
        print(f"runtime internal refusal: {type(exc).__name__}", file=sys.stderr)
        return 3


if __name__ == "__main__":
    raise SystemExit(main())
