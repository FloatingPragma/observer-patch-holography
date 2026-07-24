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
    "custody": "custody",
    "control": "controls",
    "analysis_code": "analysis_binding",
    "claim_text": "analysis_binding",
    "protocol": "analysis_binding",
    "run_schedule": "completeness",
    "device_identity": "device_identity",
    "attestation": "attestation",
}

REQUIRED_ROLES = frozenset(ROLE_PREDICATE)


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
    if set(run_ids) != set(schedule["run_ids"]):
        insufficient(
            "SELECTIVE_REPORTING",
            (
                f"scheduled={sorted(schedule['run_ids'])}, "
                f"reported={sorted(set(run_ids))}"
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

    declared_controls = set(bundle["controls"]["artifact_ids"])
    for control_id in sorted(declared_controls):
        require_artifact(
            control_id,
            "control",
            predicate="controls",
            code="CONTROL_ARTIFACT_INVALID",
        )

    nonces: list[str] = []
    for run in runs:
        run_id = run["run_id"]
        nonces.append(run["capture_nonce"])
        if run["device_id"] != claim_device:
            insufficient(
                "DEVICE_SUBSTITUTION",
                f"run {run_id} names {run['device_id']} instead of {claim_device}",
                "device_identity",
            )
        identity_id = run["device_identity_artifact_id"]
        if require_artifact(
            identity_id,
            "device_identity",
            predicate="device_identity",
            code="DEVICE_IDENTITY_ARTIFACT_INVALID",
        ):
            identity = json_artifact(identity_id, "device_identity")
            if not isinstance(identity, dict) or identity.get("device_id") != claim_device:
                insufficient(
                    "DEVICE_SUBSTITUTION",
                    f"identity artifact {identity_id} does not bind {claim_device}",
                    "device_identity",
                )

        for raw_id in run["raw_artifact_ids"]:
            if not require_artifact(
                raw_id,
                "raw_measurement",
                predicate="raw_capture",
                code="RAW_CAPTURE_ARTIFACT_INVALID",
            ):
                continue
            raw = json_artifact(raw_id, "raw_capture")
            expected_raw = {
                "capture_nonce": run["capture_nonce"],
                "device_id": run["device_id"],
                "run_id": run_id,
            }
            if not isinstance(raw, dict) or any(
                raw.get(key) != value for key, value in expected_raw.items()
            ):
                insufficient(
                    "RAW_CAPTURE_BINDING_MISMATCH",
                    f"{raw_id} does not bind run/device/nonce",
                    "raw_capture",
                )

        run_controls = set(run["control_artifact_ids"])
        if not run_controls or not run_controls.issubset(declared_controls):
            insufficient(
                "CONTROL_COVERAGE_INCOMPLETE",
                f"run {run_id} lacks declared controls",
                "controls",
            )
        for control_id in run_controls:
            if require_artifact(
                control_id,
                "control",
                predicate="controls",
                code="CONTROL_ARTIFACT_INVALID",
            ):
                control = json_artifact(control_id, "controls")
                if not isinstance(control, dict) or control.get("run_id") != run_id:
                    insufficient(
                        "CONTROL_BINDING_MISMATCH",
                        f"{control_id} does not bind run {run_id}",
                        "controls",
                    )

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
                "calibration_id": calibration["calibration_id"],
                "device_id": calibration["device_id"],
                "reference_id": calibration["reference_id"],
                "valid_from_utc": calibration["valid_from_utc"],
                "valid_until_utc": calibration["valid_until_utc"],
            }
            if stored_calibration != expected_calibration:
                insufficient(
                    "CALIBRATION_ARTIFACT_MISMATCH",
                    f"{calibration_id} differs from typed calibration",
                    "calibration_chain",
                )
        if calibration["device_id"] != run["device_id"]:
            insufficient(
                "CALIBRATION_DEVICE_MISMATCH",
                f"run {run_id} and calibration device differ",
                "calibration_chain",
            )
        captured = _parse_utc(run["captured_utc"])
        valid_from = _parse_utc(calibration["valid_from_utc"])
        valid_until = _parse_utc(calibration["valid_until_utc"])
        if not valid_from <= captured <= valid_until:
            insufficient(
                "CALIBRATION_OUT_OF_WINDOW",
                f"run {run_id} capture is outside calibration validity",
                "calibration_chain",
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
    if runs and segments:
        chain_start = min(_parse_utc(run["captured_utc"]) for run in runs)
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

    analysis = bundle["analysis_binding"]
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
    required_inputs = {
        artifact_id
        for artifact_id, row in artifacts.items()
        if row["role"]
        in {
            "raw_measurement",
            "calibration",
            "custody",
            "control",
            "protocol",
            "run_schedule",
            "device_identity",
        }
    }
    if not required_inputs.issubset(set(analysis["input_artifact_ids"])):
        insufficient(
            "ANALYSIS_INPUT_COVERAGE_INCOMPLETE",
            "analysis inputs omit raw, calibration, custody, control, protocol, schedule, or device evidence",
            "analysis_binding",
        )
    claim_path = artifact_paths.get(analysis["claim_text_artifact_id"])
    if claim_path and claim_path.is_file():
        claim_text = claim_path.read_text(encoding="utf-8")
        if bundle["claim"]["effect_statement"] not in claim_text:
            insufficient(
                "CLAIM_TEXT_MISMATCH",
                "structured effect statement is absent from bound claim text",
                "analysis_binding",
            )
    protocol_payload = json_artifact(analysis["protocol_artifact_id"], "analysis_binding")
    if not isinstance(protocol_payload, dict) or protocol_payload.get("protocol_id") != claim_protocol:
        insufficient(
            "PROTOCOL_BINDING_MISMATCH",
            "protocol artifact does not bind the structured protocol id",
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
        help="identified {registry_id, seen_nonces} JSON object maintained independently",
    )
    parser.add_argument(
        "--trust-policy",
        type=Path,
        help="operator-pinned Ed25519 trust roots; never read from the producer bundle",
    )
    parser.add_argument(
        "--external-evidence",
        type=Path,
        help="signature, preregistration, witness, and replay-authority evidence",
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
