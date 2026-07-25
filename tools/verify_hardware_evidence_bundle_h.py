#!/usr/bin/env python3
"""Fail-closed verifier for OPH hardware evidence bundle class H (#325).

The producer is not the verifier.  This module validates the public v1 JSON
Schema, reads and hashes every referenced artifact itself, recomputes the
canonical bundle binding, and evaluates evidentiary predicates.  Values under
``producer_assertions`` are deliberately ignored.  Physical promotion also
requires operator-pinned external evidence: real Ed25519 signatures, signed
pre-run commitments, a signed replay-registry assignment, deterministic
analysis replay, and an independently administered witness.

Verdicts:

* ``INVALID``: malformed schema, unsafe/missing artifacts, or hash/root drift.
* ``INSUFFICIENT``: integrity-valid packet missing an evidentiary predicate.
* ``SUFFICIENT_RELATIVE_TO_DECLARED_THREAT_MODEL``: every internal predicate
  and every applicable operator-pinned external gate passed.  This is an
  evidence-policy result, not proof that nature produced the claimed effect.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable

from jsonschema import Draft202012Validator, FormatChecker

try:
    from tools.hardware_evidence_external import verify_external_evidence
except ModuleNotFoundError:  # Direct execution from tools/.
    from hardware_evidence_external import verify_external_evidence

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SCHEMA = REPO_ROOT / "schemas" / "hardware_evidence_bundle_h_v1.schema.json"

PREDICATES = (
    "raw_capture",
    "calibration_chain",
    "custody",
    "controls",
    "analysis_binding",
    "completeness",
    "replay_protection",
    "device_identity",
    "attestation",
)

ROLE_PREDICATE = {
    "raw_measurement": "raw_capture",
    "calibration": "calibration_chain",
    "calibration_reference": "calibration_chain",
    "custody": "custody",
    "control": "controls",
    "analysis_code": "analysis_binding",
    "claim_text": "analysis_binding",
    "protocol": "analysis_binding",
    "run_schedule": "completeness",
    "device_identity": "device_identity",
    "nonce_reservation": "replay_protection",
    "attestation": "attestation",
}

REQUIRED_ROLES = frozenset(ROLE_PREDICATE)

PROTOCOL_SCHEMA = "oph.hardware_evidence_bundle_h.protocol.v1"
RAW_CAPTURE_SCHEMA = "oph.hardware_evidence_bundle_h.raw_capture.v1"
CONTROL_CAPTURE_SCHEMA = "oph.hardware_evidence_bundle_h.control_capture.v1"
DEVICE_IDENTITY_SCHEMA = "oph.hardware_evidence_bundle_h.device_identity.v1"
CLAIM_TEXT_SCHEMA = "oph.hardware_evidence_bundle_h.structured_claim.v1"
ANALYSIS_REPLAY_COMMAND = [
    "internal:paired_calibrated_mean_max_deviation.v1"
]
CONTROL_KINDS = frozenset({"blank", "detuned_twin", "sham", "synthetic_blank"})
PHYSICAL_CONTROL_KINDS = frozenset({"blank", "detuned_twin", "sham"})


def _canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def recompute_bundle_binding(bundle: dict[str, Any]) -> str:
    """Return the v1 root over the full typed bundle with a blank root field."""
    payload = copy.deepcopy(bundle)
    payload["bundle_binding"]["binding_sha256"] = ""
    return hashlib.sha256(_canonical_bytes(payload)).hexdigest()


def _parse_utc(value: str) -> datetime:
    parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    if parsed.tzinfo is None:
        raise ValueError("timestamp has no timezone")
    return parsed


def _load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _schema_errors(bundle: Any, schema_path: Path) -> list[str]:
    schema = _load_json(schema_path)
    Draft202012Validator.check_schema(schema)
    validator = Draft202012Validator(schema, format_checker=FormatChecker())
    errors = sorted(validator.iter_errors(bundle), key=lambda err: list(err.path))
    rendered: list[str] = []
    for error in errors:
        location = ".".join(str(part) for part in error.absolute_path) or "<root>"
        rendered.append(f"{location}: {error.message}")
    return rendered


def _safe_artifact_path(bundle_dir: Path, relative: str) -> Path | None:
    candidate = Path(relative)
    if candidate.is_absolute() or ".." in candidate.parts:
        return None
    base = bundle_dir.resolve()
    resolved = (base / candidate).resolve()
    if not resolved.is_relative_to(base):
        return None
    return resolved


def _seen_nonces(path: Path | None, registry_id: str) -> tuple[set[str], list[str]]:
    if path is None:
        return set(), []
    try:
        payload = _load_json(path)
    except (OSError, json.JSONDecodeError) as exc:
        return set(), [f"replay registry unreadable: {exc}"]
    if not isinstance(payload, dict):
        return set(), ["replay registry must be an identified JSON object"]
    if payload.get("registry_id") != registry_id:
        return set(), ["replay registry id does not match the bundle"]
    values = payload.get("seen_nonces")
    if not isinstance(values, list):
        return set(), ["replay registry seen_nonces must be a list"]
    return {str(item) for item in values}, []


def verify_bundle(
    bundle_path: Path,
    *,
    schema_path: Path = DEFAULT_SCHEMA,
    replay_registry_path: Path | None = None,
    trust_policy_path: Path | None = None,
    external_evidence_path: Path | None = None,
) -> dict[str, Any]:
    """Verify a class-H bundle and return a deterministic machine report."""
    try:
        bundle = _load_json(bundle_path)
    except (OSError, json.JSONDecodeError) as exc:
        return {
            "schema": "oph.hardware_evidence_bundle_h.verification.v1",
            "verdict": "INVALID",
            "schema_valid": False,
            "integrity_valid": False,
            "predicates": {name: False for name in PREDICATES},
            "rejection_codes": ["BUNDLE_JSON_INVALID"],
            "details": [str(exc)],
            "ignored_producer_assertions": [],
        }

    schema_errors = _schema_errors(bundle, schema_path)
    ignored = sorted((bundle.get("producer_assertions") or {}).keys())
    if schema_errors:
        return {
            "schema": "oph.hardware_evidence_bundle_h.verification.v1",
            "bundle_id": bundle.get("bundle_id"),
            "verdict": "INVALID",
            "schema_valid": False,
            "integrity_valid": False,
            "predicates": {name: False for name in PREDICATES},
            "rejection_codes": ["SCHEMA_INVALID"],
            "details": schema_errors,
            "ignored_producer_assertions": ignored,
        }

    predicates = {name: True for name in PREDICATES}
    invalid_codes: set[str] = set()
    insufficient_codes: set[str] = set()
    details: list[str] = []

    def invalid(code: str, detail: str, predicate: str | None = None) -> None:
        invalid_codes.add(code)
        details.append(f"{code}: {detail}")
        if predicate:
            predicates[predicate] = False

    def insufficient(code: str, detail: str, predicate: str) -> None:
        insufficient_codes.add(code)
        details.append(f"{code}: {detail}")
        predicates[predicate] = False

    bundle_dir = bundle_path.resolve().parent
    artifacts: dict[str, dict[str, Any]] = {}
    artifact_paths: dict[str, Path] = {}
    seen_paths: set[str] = set()
    roles: set[str] = set()
    for row in bundle["artifacts"]:
        artifact_id = row["artifact_id"]
        if artifact_id in artifacts:
            invalid("DUPLICATE_ARTIFACT_ID", artifact_id, ROLE_PREDICATE[row["role"]])
            continue
        if row["path"] in seen_paths:
            invalid("DUPLICATE_ARTIFACT_PATH", row["path"], ROLE_PREDICATE[row["role"]])
        seen_paths.add(row["path"])
        artifacts[artifact_id] = row
        roles.add(row["role"])
        path = _safe_artifact_path(bundle_dir, row["path"])
        if path is None:
            invalid("ARTIFACT_PATH_UNSAFE", row["path"], ROLE_PREDICATE[row["role"]])
            continue
        artifact_paths[artifact_id] = path
        if not path.is_file():
            invalid("ARTIFACT_MISSING", row["path"], ROLE_PREDICATE[row["role"]])
            continue
        actual = sha256_path(path)
        if actual != row["sha256"]:
            invalid(
                "ARTIFACT_HASH_MISMATCH",
                f"{artifact_id}: declared {row['sha256']}, actual {actual}",
                ROLE_PREDICATE[row["role"]],
            )

    for missing_role in sorted(REQUIRED_ROLES - roles):
        predicate = ROLE_PREDICATE[missing_role]
        insufficient(
            f"{missing_role.upper()}_MISSING",
            f"no artifact has role {missing_role}",
            predicate,
        )

    artifact_ids = set(artifacts)
    covered = set(bundle["bundle_binding"]["covered_artifact_ids"])
    if covered != artifact_ids:
        invalid(
            "BINDING_COVERAGE_MISMATCH",
            f"covered={sorted(covered)}, artifacts={sorted(artifact_ids)}",
            "analysis_binding",
        )
    expected_binding = recompute_bundle_binding(bundle)
    if expected_binding != bundle["bundle_binding"]["binding_sha256"]:
        invalid(
            "BUNDLE_BINDING_MISMATCH",
            (
                f"declared {bundle['bundle_binding']['binding_sha256']}, "
                f"actual {expected_binding}"
            ),
            "analysis_binding",
        )

    def require_artifact(
        artifact_id: str,
        role: str,
        *,
        predicate: str,
        code: str,
    ) -> bool:
        row = artifacts.get(artifact_id)
        if row is None or row["role"] != role:
            insufficient(code, f"{artifact_id} is not a {role} artifact", predicate)
            return False
        return True

    def json_artifact(artifact_id: str, predicate: str) -> Any | None:
        path = artifact_paths.get(artifact_id)
        if path is None or not path.is_file():
            return None
        try:
            return _load_json(path)
        except (OSError, json.JSONDecodeError) as exc:
            invalid(
                "ARTIFACT_JSON_INVALID",
                f"{artifact_id}: {exc}",
                predicate,
            )
            return None

    def closed_artifact(
        value: Any,
        *,
        required: set[str],
        label: str,
        code: str,
        predicate: str,
    ) -> bool:
        if not isinstance(value, dict):
            insufficient(code, f"{label} must be a JSON object", predicate)
            return False
        actual = set(value)
        if actual != required:
            insufficient(
                code,
                (
                    f"{label} has missing keys {sorted(required - actual)} "
                    f"and extra keys {sorted(actual - required)}"
                ),
                predicate,
            )
            return False
        return True

    def is_sha256(value: Any) -> bool:
        return (
            isinstance(value, str)
            and len(value) == 64
            and all(character in "0123456789abcdef" for character in value)
        )

    analysis = bundle["analysis_binding"]
    protocol_id = analysis["protocol_artifact_id"]
    protocol_ok = require_artifact(
        protocol_id,
        "protocol",
        predicate="analysis_binding",
        code="PROTOCOL_ARTIFACT_INVALID",
    )
    protocol = (
        json_artifact(protocol_id, "analysis_binding")
        if protocol_ok
        else None
    )
    protocol_keys = {
        "schema",
        "protocol_id",
        "device_id",
        "measurement_channel",
        "raw_unit",
        "reported_unit",
        "firmware_sha256",
        "control_policy",
    }
    if not closed_artifact(
        protocol,
        required=protocol_keys,
        label="protocol artifact",
        code="PROTOCOL_ARTIFACT_FORMAT_INVALID",
        predicate="analysis_binding",
    ):
        protocol = None
    elif (
        protocol["schema"] != PROTOCOL_SCHEMA
        or protocol["protocol_id"] != bundle["claim"]["protocol_id"]
        or protocol["device_id"] != bundle["claim"]["device_id"]
        or not isinstance(protocol["measurement_channel"], str)
        or not protocol["measurement_channel"]
        or not isinstance(protocol["raw_unit"], str)
        or not protocol["raw_unit"]
        or not isinstance(protocol["reported_unit"], str)
        or not protocol["reported_unit"]
        or not is_sha256(protocol["firmware_sha256"])
        or protocol["control_policy"] != "paired_preceding_same_protocol.v1"
    ):
        insufficient(
            "PROTOCOL_ARTIFACT_FORMAT_INVALID",
            "protocol does not define the closed channel, units, firmware, and paired-control policy",
            "analysis_binding",
        )
        protocol = None

    schedule = bundle["run_schedule"]
    schedule_ok = require_artifact(
        schedule["artifact_id"],
        "run_schedule",
        predicate="completeness",
        code="RUN_SCHEDULE_ARTIFACT_INVALID",
    )
    if schedule_ok:
        stored_schedule = json_artifact(schedule["artifact_id"], "completeness")
        if stored_schedule != {"run_ids": schedule["run_ids"]}:
            insufficient(
                "RUN_SCHEDULE_ARTIFACT_MISMATCH",
                "typed run schedule differs from its bound artifact",
                "completeness",
            )

    runs = bundle["runs"]
    run_ids = [row["run_id"] for row in runs]
    if len(set(run_ids)) != len(run_ids):
        insufficient("DUPLICATE_RUN_ID", "run ids are not unique", "completeness")
    if run_ids != schedule["run_ids"]:
        insufficient(
            "SELECTIVE_REPORTING",
            (
                f"scheduled order={schedule['run_ids']}, "
                f"reported order={run_ids}"
            ),
            "completeness",
        )

    claim_device = bundle["claim"]["device_id"]
    claim_protocol = bundle["claim"]["protocol_id"]
    calibrations = {row["calibration_id"]: row for row in bundle["calibrations"]}
    if len(calibrations) != len(bundle["calibrations"]):
        insufficient(
            "DUPLICATE_CALIBRATION_ID",
            "calibration ids are not unique",
            "calibration_chain",
        )
    used_calibration_ids = {row["calibration_id"] for row in runs}
    if set(calibrations) != used_calibration_ids:
        insufficient(
            "CALIBRATION_POPULATION_MISMATCH",
            (
                "declared calibration ids differ from the complete "
                "scheduled-run calibration population"
            ),
            "calibration_chain",
        )

    declared_controls = set(bundle["controls"]["artifact_ids"])
    for control_id in sorted(declared_controls):
        require_artifact(
            control_id,
            "control",
            predicate="controls",
            code="CONTROL_ARTIFACT_INVALID",
        )

    nonces: list[str] = []
    capture_times: list[datetime] = []
    previous_live_capture: datetime | None = None
    used_raw_ids: set[str] = set()
    used_control_ids: set[str] = set()
    for run_index, run in enumerate(runs):
        run_id = run["run_id"]
        control_capture_for_run: datetime | None = None
        nonces.extend(
            [run["control_capture_nonce"], run["capture_nonce"]]
        )
        capture_times.append(_parse_utc(run["captured_utc"]))
        if run["device_id"] != claim_device:
            insufficient(
                "DEVICE_SUBSTITUTION",
                f"run {run_id} names {run['device_id']} instead of {claim_device}",
                "device_identity",
            )
        identity_id = run["device_identity_artifact_id"]
        identity: dict[str, Any] | None = None
        if require_artifact(
            identity_id,
            "device_identity",
            predicate="device_identity",
            code="DEVICE_IDENTITY_ARTIFACT_INVALID",
        ):
            identity_payload = json_artifact(identity_id, "device_identity")
            identity_keys = {
                "schema",
                "device_id",
                "hardware_serial_hash",
                "firmware_sha256",
            }
            if closed_artifact(
                identity_payload,
                required=identity_keys,
                label=f"device identity {identity_id}",
                code="DEVICE_IDENTITY_ARTIFACT_FORMAT_INVALID",
                predicate="device_identity",
            ):
                identity = identity_payload
            if (
                identity is None
                or identity["schema"] != DEVICE_IDENTITY_SCHEMA
                or identity["device_id"] != claim_device
                or not is_sha256(identity["hardware_serial_hash"])
                or not is_sha256(identity["firmware_sha256"])
                or (
                    protocol is not None
                    and identity["firmware_sha256"] != protocol["firmware_sha256"]
                )
            ):
                insufficient(
                    "DEVICE_SUBSTITUTION",
                    (
                        f"identity artifact {identity_id} does not bind the "
                        "declared device, physical mark, and protocol firmware"
                    ),
                    "device_identity",
                )

        if len(run["raw_artifact_ids"]) != 1:
            insufficient(
                "RAW_CHANNEL_COVERAGE_INCOMPLETE",
                f"run {run_id} must carry exactly the one channel supported by v1",
                "raw_capture",
            )
        for raw_id in run["raw_artifact_ids"]:
            used_raw_ids.add(raw_id)
            if not require_artifact(
                raw_id,
                "raw_measurement",
                predicate="raw_capture",
                code="RAW_CAPTURE_ARTIFACT_INVALID",
            ):
                continue
            raw = json_artifact(raw_id, "raw_capture")
            raw_keys = {
                "schema",
                "capture_nonce",
                "channel",
                "raw_unit",
                "device_id",
                "firmware_sha256",
                "hardware_serial_hash",
                "protocol_id",
                "run_id",
                "sequence_index",
                "captured_utc",
                "samples",
            }
            if not closed_artifact(
                raw,
                required=raw_keys,
                label=f"raw capture {raw_id}",
                code="RAW_CAPTURE_FORMAT_INVALID",
                predicate="raw_capture",
            ):
                continue
            raw_matches = (
                raw["schema"] == RAW_CAPTURE_SCHEMA
                and raw["capture_nonce"] == run["capture_nonce"]
                and raw["device_id"] == run["device_id"]
                and raw["protocol_id"] == claim_protocol
                and raw["run_id"] == run_id
                and type(raw["sequence_index"]) is int
                and raw["sequence_index"] == 2 * run_index + 1
                and raw["captured_utc"] == run["captured_utc"]
                and isinstance(raw["samples"], list)
                and bool(raw["samples"])
                and protocol is not None
                and raw["channel"] == protocol["measurement_channel"]
                and raw["raw_unit"] == protocol["raw_unit"]
                and raw["firmware_sha256"] == protocol["firmware_sha256"]
                and identity is not None
                and raw["hardware_serial_hash"] == identity["hardware_serial_hash"]
            )
            if not raw_matches:
                insufficient(
                    "RAW_CAPTURE_BINDING_MISMATCH",
                    (
                        f"{raw_id} does not bind the scheduled run, device mark, "
                        "firmware, protocol channel, units, nonce, and capture time"
                    ),
                    "raw_capture",
                )

        run_controls = set(run["control_artifact_ids"])
        if (
            len(run["control_artifact_ids"]) != 1
            or not run_controls.issubset(declared_controls)
        ):
            insufficient(
                "CONTROL_COVERAGE_INCOMPLETE",
                f"run {run_id} must have exactly one declared paired control in v1",
                "controls",
            )
        for control_id in run_controls:
            used_control_ids.add(control_id)
            if require_artifact(
                control_id,
                "control",
                predicate="controls",
                code="CONTROL_ARTIFACT_INVALID",
            ):
                control = json_artifact(control_id, "controls")
                control_keys = {
                    "schema",
                    "capture_nonce",
                    "control_kind",
                    "channel",
                    "raw_unit",
                    "device_id",
                    "firmware_sha256",
                    "hardware_serial_hash",
                    "protocol_id",
                    "run_id",
                    "sequence_index",
                    "captured_utc",
                    "samples",
                }
                if not closed_artifact(
                    control,
                    required=control_keys,
                    label=f"control capture {control_id}",
                    code="CONTROL_CAPTURE_FORMAT_INVALID",
                    predicate="controls",
                ):
                    continue
                try:
                    control_capture = _parse_utc(control["captured_utc"])
                    live_capture = _parse_utc(run["captured_utc"])
                except (TypeError, ValueError):
                    control_capture = None
                    live_capture = None
                control_matches = (
                    control["schema"] == CONTROL_CAPTURE_SCHEMA
                    and control["capture_nonce"] == run["control_capture_nonce"]
                    and control["control_kind"] in CONTROL_KINDS
                    and control["run_id"] == run_id
                    and control["device_id"] == claim_device
                    and control["protocol_id"] == claim_protocol
                    and type(control["sequence_index"]) is int
                    and control["sequence_index"] == 2 * run_index
                    and isinstance(control["samples"], list)
                    and bool(control["samples"])
                    and control_capture is not None
                    and live_capture is not None
                    and control_capture < live_capture
                    and (
                        previous_live_capture is None
                        or previous_live_capture < control_capture
                    )
                    and protocol is not None
                    and control["channel"] == protocol["measurement_channel"]
                    and control["raw_unit"] == protocol["raw_unit"]
                    and control["firmware_sha256"] == protocol["firmware_sha256"]
                    and identity is not None
                    and control["hardware_serial_hash"]
                    == identity["hardware_serial_hash"]
                    and (
                        bundle["claim_boundary"]["physical_claim"] is not True
                        or control["control_kind"] in PHYSICAL_CONTROL_KINDS
                    )
                )
                if not control_matches:
                    if (
                        bundle["claim_boundary"]["physical_claim"] is True
                        and control["control_kind"] not in PHYSICAL_CONTROL_KINDS
                    ):
                        insufficient(
                            "CONTROL_KIND_NOT_PHYSICAL",
                            (
                                f"{control_id} uses synthetic-only control kind "
                                f"{control['control_kind']}"
                            ),
                            "controls",
                        )
                    insufficient(
                        "CONTROL_BINDING_MISMATCH",
                        (
                            f"{control_id} is not the immediately preceding, "
                            "same-device, same-firmware, same-protocol control for {run_id}"
                        ),
                        "controls",
                    )
                if control_capture is not None:
                    capture_times.append(control_capture)
                    control_capture_for_run = control_capture
        previous_live_capture = _parse_utc(run["captured_utc"])

        calibration = calibrations.get(run["calibration_id"])
        if calibration is None:
            insufficient(
                "CALIBRATION_MISSING_FOR_RUN",
                f"run {run_id} names unknown calibration",
                "calibration_chain",
            )
            continue
        calibration_id = calibration["artifact_id"]
        if require_artifact(
            calibration_id,
            "calibration",
            predicate="calibration_chain",
            code="CALIBRATION_ARTIFACT_INVALID",
        ):
            stored_calibration = json_artifact(calibration_id, "calibration_chain")
            expected_calibration = {
                key: calibration[key]
                for key in (
                    "calibration_id",
                    "device_id",
                    "reference_id",
                    "reference_artifact_id",
                    "reference_uri",
                    "reference_certificate_sha256",
                    "channel_id",
                    "raw_unit",
                    "reported_unit",
                    "transformation",
                    "valid_from_utc",
                    "valid_until_utc",
                )
            }
            if stored_calibration != expected_calibration:
                insufficient(
                    "CALIBRATION_ARTIFACT_MISMATCH",
                    f"{calibration_id} differs from typed calibration",
                    "calibration_chain",
                )
        reference_id = calibration["reference_artifact_id"]
        reference = artifacts.get(reference_id)
        if (
            reference is None
            or reference["role"] != "calibration_reference"
            or reference["sha256"]
            != calibration["reference_certificate_sha256"]
        ):
            insufficient(
                "CALIBRATION_REFERENCE_MISMATCH",
                (
                    f"calibration {calibration['calibration_id']} does not "
                    "bind its declared reference-certificate artifact"
                ),
                "calibration_chain",
            )
        if calibration["device_id"] != run["device_id"]:
            insufficient(
                "CALIBRATION_DEVICE_MISMATCH",
                f"run {run_id} and calibration device differ",
                "calibration_chain",
            )
        if protocol is not None and (
            calibration["channel_id"] != protocol["measurement_channel"]
            or calibration["raw_unit"] != protocol["raw_unit"]
            or calibration["reported_unit"] != protocol["reported_unit"]
        ):
            insufficient(
                "CALIBRATION_PROTOCOL_MISMATCH",
                f"run {run_id} calibration does not transform the protocol channel and units",
                "calibration_chain",
            )
        captured = _parse_utc(run["captured_utc"])
        valid_from = _parse_utc(calibration["valid_from_utc"])
        valid_until = _parse_utc(calibration["valid_until_utc"])
        if (
            control_capture_for_run is None
            or not valid_from <= control_capture_for_run <= valid_until
            or not valid_from <= captured <= valid_until
        ):
            insufficient(
                "CALIBRATION_OUT_OF_WINDOW",
                (
                    f"run {run_id} control or live capture is outside "
                    "calibration validity"
                ),
                "calibration_chain",
            )

    declared_raw_ids = {
        artifact_id
        for artifact_id, row in artifacts.items()
        if row["role"] == "raw_measurement"
    }
    if used_raw_ids != declared_raw_ids:
        insufficient(
            "RAW_CHANNEL_COVERAGE_INCOMPLETE",
            "declared raw captures differ from the scheduled run population",
            "raw_capture",
        )
    if used_control_ids != declared_controls:
        insufficient(
            "CONTROL_COVERAGE_INCOMPLETE",
            "declared controls differ from the interleaved scheduled controls",
            "controls",
        )

    if len(set(nonces)) != len(nonces):
        insufficient(
            "REPLAY_NONCE_REUSED",
            "capture nonce repeats inside the bundle",
            "replay_protection",
        )
    replay_in_scope = bundle["threat_model"]["replay"]["in_scope"]
    if replay_in_scope and replay_registry_path is None:
        insufficient(
            "REPLAY_REGISTRY_REQUIRED",
            "cross-bundle replay cannot be checked without an independent nonce registry",
            "replay_protection",
        )
    seen_nonces, registry_errors = _seen_nonces(
        replay_registry_path,
        bundle["replay_protection"]["registry_id"],
    )
    for error in registry_errors:
        insufficient("REPLAY_REGISTRY_INVALID", error, "replay_protection")
    for nonce in nonces:
        if nonce in seen_nonces:
            insufficient(
                "REPLAY_NONCE_SEEN",
                f"capture nonce {nonce} is present in the independent registry",
                "replay_protection",
            )
    reservation_id = bundle["replay_protection"]["reservation_artifact_id"]
    if require_artifact(
        reservation_id,
        "nonce_reservation",
        predicate="replay_protection",
        code="REPLAY_RESERVATION_ARTIFACT_INVALID",
    ):
        reservation = json_artifact(reservation_id, "replay_protection")
        expected_reservation = {
            "schema": "oph.hardware_evidence_bundle_h.nonce_reservation.v1",
            "registry_id": bundle["replay_protection"]["registry_id"],
            "bundle_id": bundle["bundle_id"],
            "device_id": claim_device,
            "protocol_id": claim_protocol,
            "capture_nonces": nonces,
        }
        if reservation != expected_reservation:
            insufficient(
                "REPLAY_RESERVATION_MISMATCH",
                (
                    "bound nonce reservation must name every ordered raw and "
                    "control capture nonce"
                ),
                "replay_protection",
            )

    custody = bundle["custody"]
    if require_artifact(
        custody["artifact_id"],
        "custody",
        predicate="custody",
        code="CUSTODY_ARTIFACT_INVALID",
    ):
        stored_custody = json_artifact(custody["artifact_id"], "custody")
        expected_custody = {
            "continuous_declared": custody["continuous_declared"],
            "data_artifact_ids": custody["data_artifact_ids"],
            "segments": custody["segments"],
        }
        if stored_custody != expected_custody:
            insufficient(
                "CUSTODY_ARTIFACT_MISMATCH",
                "typed custody chain differs from its bound artifact",
                "custody",
            )
    if not custody["continuous_declared"]:
        insufficient("CUSTODY_BREAK", "custody is not declared continuous", "custody")
    custody_data = set(custody["data_artifact_ids"])
    required_custody_data = artifact_ids - {custody["artifact_id"]}
    if custody_data != required_custody_data:
        insufficient(
            "CUSTODY_DATA_COVERAGE_INCOMPLETE",
            (
                f"custody data={sorted(custody_data)}, "
                f"required={sorted(required_custody_data)}"
            ),
            "custody",
        )
    segments: list[tuple[datetime, datetime]] = []
    for segment in custody["segments"]:
        if segment["device_id"] != claim_device:
            insufficient(
                "DEVICE_SUBSTITUTION",
                "custody segment names a different device",
                "device_identity",
            )
            continue
        start = _parse_utc(segment["start_utc"])
        end = _parse_utc(segment["end_utc"])
        if end < start:
            insufficient("CUSTODY_BREAK", "custody segment ends before it starts", "custody")
            continue
        segments.append((start, end))
    if capture_times and segments:
        chain_start = min(capture_times)
        chain_end = _parse_utc(bundle["created_utc"])
        cursor = chain_start
        for start, end in sorted(segments):
            if end < cursor:
                continue
            if start > cursor:
                break
            cursor = max(cursor, end)
            if cursor >= chain_end:
                break
        if cursor < chain_end:
            insufficient(
                "CUSTODY_BREAK",
                "custody segments do not continuously cover capture through publication",
                "custody",
            )

    for field, role in (
        ("analysis_artifact_id", "analysis_code"),
        ("claim_text_artifact_id", "claim_text"),
        ("protocol_artifact_id", "protocol"),
    ):
        require_artifact(
            analysis[field],
            role,
            predicate="analysis_binding",
            code=f"{role.upper()}_ARTIFACT_INVALID",
        )
    if not analysis["frozen_before_unblinding"]:
        insufficient(
            "ANALYSIS_NOT_FROZEN",
            "analysis was not frozen before unblinding",
            "analysis_binding",
        )
    if analysis["replay_command"] != ANALYSIS_REPLAY_COMMAND:
        insufficient(
            "ANALYSIS_REPLAY_COMMAND_UNSUPPORTED",
            "analysis must name the closed verifier-internal replay operation",
            "analysis_binding",
        )
    required_inputs = {
        artifact_id
        for artifact_id, row in artifacts.items()
        if row["role"]
        in {
            "raw_measurement",
            "calibration",
            "calibration_reference",
            "custody",
            "control",
            "protocol",
            "run_schedule",
            "device_identity",
        }
    }
    if required_inputs != set(analysis["input_artifact_ids"]):
        insufficient(
            "ANALYSIS_INPUT_POPULATION_MISMATCH",
            (
                "analysis inputs must equal the raw, calibration/reference, "
                "custody, control, protocol, schedule, and device evidence"
            ),
            "analysis_binding",
        )
    claim_payload = json_artifact(
        analysis["claim_text_artifact_id"],
        "analysis_binding",
    )
    if claim_payload != {
        "schema": CLAIM_TEXT_SCHEMA,
        "claim": bundle["claim"],
    }:
        insufficient(
            "CLAIM_TEXT_MISMATCH",
            (
                "bound claim text must be the closed structured claim and "
                "match E, M, U, device, protocol, conditions, and boundary flags exactly"
            ),
            "analysis_binding",
        )

    attestation = bundle["attestation"]
    for artifact_id in attestation["artifact_ids"]:
        require_artifact(
            artifact_id,
            "attestation",
            predicate="attestation",
            code="ATTESTATION_ARTIFACT_INVALID",
        )
    signer_index = {row["signer_id"]: row for row in bundle["signers"]}
    if len(signer_index) != len(bundle["signers"]):
        insufficient("DUPLICATE_SIGNER_ID", "signer ids are not unique", "attestation")
    used_signers: list[dict[str, Any]] = []
    for signer_id in attestation["signer_ids"]:
        signer = signer_index.get(signer_id)
        if signer is None:
            insufficient(
                "ATTESTATION_SIGNER_MISSING",
                f"unknown signer {signer_id}",
                "attestation",
            )
            continue
        used_signers.append(signer)
        if not set(attestation["artifact_ids"]).issubset(
            set(signer["signed_artifact_ids"])
        ):
            insufficient(
                "ATTESTATION_SIGNATURE_COVERAGE_INCOMPLETE",
                f"signer {signer_id} does not cover every attestation artifact",
                "attestation",
            )
        if signer["compromised"]:
            insufficient(
                "SIGNER_COMPROMISED",
                f"attestation signer {signer_id} is marked compromised",
                "attestation",
            )
    attestation_mode = attestation["mode"]
    independent_parties = set(attestation["party_ids"]) - {bundle["claimant_id"]}
    signers_from_independent_parties = [
        row
        for row in used_signers
        if row["party_id"] in independent_parties and not row["compromised"]
    ]
    mode_passes = False
    if attestation_mode == "independent_reproduction":
        mode_passes = (
            attestation["independent_of_claimant"]
            and attestation["fresh_device_or_run_series"]
            and bool(attestation["artifact_ids"])
            and bool(independent_parties)
            and bool(signers_from_independent_parties)
        )
    elif attestation_mode == "independent_end_to_end_witness":
        mode_passes = (
            attestation["independent_of_claimant"]
            and attestation["witnessed_end_to_end"]
            and bool(attestation["artifact_ids"])
            and bool(independent_parties)
            and bool(signers_from_independent_parties)
        )
    elif attestation_mode == "threat_model_argument":
        mode_passes = (
            not bundle["claim"]["extraordinary_effect"]
            and bool(attestation["artifact_ids"])
        )
    elif attestation_mode == "same_capture_multi_signature":
        insufficient(
            "SIGNATURES_NOT_INDEPENDENT_ATTESTATION",
            "multiple signatures over one capture do not constitute reproduction or witnessing",
            "attestation",
        )
    if not mode_passes:
        insufficient(
            "ATTESTATION_REQUIRED",
            f"attestation mode {attestation_mode} does not satisfy the claim's rule",
            "attestation",
        )

    attestation_artifacts = [
        json_artifact(artifact_id, "attestation")
        for artifact_id in attestation["artifact_ids"]
    ]
    for payload in attestation_artifacts:
        if payload is None:
            continue
        if (
            not isinstance(payload, dict)
            or payload.get("bundle_id") != bundle["bundle_id"]
            or payload.get("mode") != attestation_mode
        ):
            insufficient(
                "ATTESTATION_ARTIFACT_MISMATCH",
                "attestation artifact does not bind bundle id and mode",
                "attestation",
            )

    if not bundle["claim_boundary"]["physical_claim"]:
        insufficient(
            "NON_PHYSICAL_FIXTURE",
            "bundle is explicitly a non-physical contract fixture",
            "attestation",
        )

    try:
        external = verify_external_evidence(
            bundle=bundle,
            artifacts=artifacts,
            artifact_paths=artifact_paths,
            replay_registry_path=replay_registry_path,
            trust_policy_path=trust_policy_path,
            external_evidence_path=external_evidence_path,
        )
    except Exception as exc:
        # Trust/evidence files are untrusted inputs. Any parser or type failure
        # is an INVALID report, never a traceback or an INSUFFICIENT exit code.
        external = {
            "invalid_codes": {"EXTERNAL_EVIDENCE_VERIFICATION_EXCEPTION"},
            "insufficient_codes": set(),
            "details": [
                (
                    "EXTERNAL_EVIDENCE_VERIFICATION_EXCEPTION: "
                    f"{type(exc).__name__}: {exc}"
                )
            ],
            "predicate_failures": {"attestation"},
            "open_gates": [
                "TRUST_ROOT_SIGNATURE_VERIFICATION_OPEN",
                "FRESH_REPRODUCTION_BUNDLE_VERIFICATION_OPEN",
                "ANALYSIS_TO_CLAIM_EXECUTION_OPEN",
                "EXTERNAL_RUN_SCHEDULE_COMMITMENT_OPEN",
                "DEVICE_CUSTODY_PROVENANCE_VERIFICATION_OPEN",
                "REPLAY_REGISTRY_AUTHORITY_OPEN",
            ],
            "verified_anchor_ids": [],
        }
    invalid_codes.update(external["invalid_codes"])
    insufficient_codes.update(external["insufficient_codes"])
    details.extend(external["details"])
    for predicate in external["predicate_failures"]:
        predicates[predicate] = False
    open_gates = external["open_gates"]

    if invalid_codes:
        verdict = "INVALID"
    elif insufficient_codes:
        verdict = "INSUFFICIENT"
    else:
        verdict = "SUFFICIENT_RELATIVE_TO_DECLARED_THREAT_MODEL"

    return {
        "schema": "oph.hardware_evidence_bundle_h.verification.v1",
        "bundle_id": bundle["bundle_id"],
        "verdict": verdict,
        "schema_valid": True,
        "integrity_valid": not invalid_codes,
        "predicates": predicates,
        "rejection_codes": sorted(invalid_codes | insufficient_codes),
        "details": sorted(details),
        "ignored_producer_assertions": ignored,
        "open_gates": open_gates,
        "claim_boundary": bundle["claim_boundary"],
        "verified_external_anchor_ids": external["verified_anchor_ids"],
        "verifier_boundary": (
            "The verdict evaluates evidence sufficiency relative to the "
            "operator-pinned trust policy and declared threat model. It does "
            "not prove the physical truth of the claim or protect against "
            "collusion by every independent authority."
        ),
    }


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle", type=Path)
    parser.add_argument("--schema", type=Path, default=DEFAULT_SCHEMA)
    parser.add_argument(
        "--replay-registry",
        type=Path,
        help=(
            "independently maintained registry snapshot with registry_id, "
            "seen_nonces, and nonce assignments"
        ),
    )
    parser.add_argument(
        "--trust-policy",
        type=Path,
        help="operator-pinned Ed25519 trust roots; never read from the producer bundle",
    )
    parser.add_argument(
        "--external-evidence",
        type=Path,
        help=(
            "artifact signatures, schedule/protocol/nonce/analysis "
            "preregistration, witness, and replay-authority evidence"
        ),
    )
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args(argv)

    report = verify_bundle(
        args.bundle,
        schema_path=args.schema,
        replay_registry_path=args.replay_registry,
        trust_policy_path=args.trust_policy,
        external_evidence_path=args.external_evidence,
    )
    encoded = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.json_out:
        args.json_out.write_text(encoded, encoding="utf-8")
    print(encoded, end="")
    return {
        "SUFFICIENT_RELATIVE_TO_DECLARED_THREAT_MODEL": 0,
        "INSUFFICIENT": 1,
        "INVALID": 2,
    }[report["verdict"]]


if __name__ == "__main__":
    raise SystemExit(main())
