"""External-evidence verification for hardware claim class H.

The bundle producer does not choose the trust roots.  The verifier operator
supplies a closed trust-policy file separately from the bundle and its external
evidence.  Ed25519 signatures then bind the canonical bundle root, artifacts,
pre-run commitments, and replay-registry snapshot to those pinned roots.
"""

from __future__ import annotations

import base64
import hashlib
import json
from datetime import datetime
from fractions import Fraction
from pathlib import Path
from typing import Any

from cryptography.exceptions import InvalidSignature
from cryptography.hazmat.primitives.asymmetric.ed25519 import Ed25519PublicKey

TRUST_SCHEMA = "oph.hardware_evidence_bundle_h.trust_policy.v1"
EVIDENCE_SCHEMA = "oph.hardware_evidence_bundle_h.external_evidence.v1"

OPEN_GATES = (
    "TRUST_ROOT_SIGNATURE_VERIFICATION_OPEN",
    "FRESH_REPRODUCTION_BUNDLE_VERIFICATION_OPEN",
    "ANALYSIS_TO_CLAIM_EXECUTION_OPEN",
    "EXTERNAL_RUN_SCHEDULE_COMMITMENT_OPEN",
    "DEVICE_CUSTODY_PROVENANCE_VERIFICATION_OPEN",
    "REPLAY_REGISTRY_AUTHORITY_OPEN",
)

ANCHOR_ROLES = {
    "claimant",
    "measurement_authority",
    "calibration_authority",
    "device_authority",
    "custody_authority",
    "preregistration_authority",
    "replay_authority",
    "independent_attestor",
}

ARTIFACT_AUTHORITY = {
    "raw_measurement": "measurement_authority",
    "calibration": "calibration_authority",
    "custody": "custody_authority",
    "control": "measurement_authority",
    "analysis_code": "claimant",
    "claim_text": "claimant",
    "protocol": "claimant",
    "device_identity": "device_authority",
    "attestation": "independent_attestor",
}


def canonical_bytes(value: Any) -> bytes:
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


def _parse_utc(value: str) -> datetime:
    parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    if parsed.tzinfo is None:
        raise ValueError("timestamp has no timezone")
    return parsed


def _read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _closed_object(
    value: Any,
    *,
    required: set[str],
    label: str,
) -> str | None:
    if not isinstance(value, dict):
        return f"{label} must be an object"
    actual = set(value)
    missing = required - actual
    extra = actual - required
    if missing or extra:
        return (
            f"{label} has missing keys {sorted(missing)} "
            f"and extra keys {sorted(extra)}"
        )
    return None


def _decode_public_key(value: str) -> Ed25519PublicKey:
    raw = base64.b64decode(value, validate=True)
    if len(raw) != 32:
        raise ValueError("Ed25519 public key must contain 32 bytes")
    return Ed25519PublicKey.from_public_bytes(raw)


def _verify_signature(
    key: Ed25519PublicKey,
    signature_b64: str,
    payload: dict[str, Any],
) -> bool:
    try:
        signature = base64.b64decode(signature_b64, validate=True)
        key.verify(signature, canonical_bytes(payload))
    except (InvalidSignature, ValueError):
        return False
    return True


def signature_payload(
    *,
    subject_kind: str,
    subject_id: str,
    sha256: str,
    key_id: str,
) -> dict[str, str]:
    """Return the exact signed message for a root or artifact signature."""
    return {
        "key_id": key_id,
        "sha256": sha256,
        "subject_id": subject_id,
        "subject_kind": subject_kind,
    }


def commitment_payload(
    *,
    kind: str,
    artifact_id: str,
    sha256: str,
    committed_utc: str,
    key_id: str,
    policy_id: str,
    bundle_id: str,
    claimant_id: str,
    device_id: str,
    protocol_id: str,
) -> dict[str, str]:
    """Return the exact signed message for a preregistration commitment."""
    return {
        "artifact_id": artifact_id,
        "bundle_id": bundle_id,
        "claimant_id": claimant_id,
        "committed_utc": committed_utc,
        "device_id": device_id,
        "key_id": key_id,
        "kind": kind,
        "policy_id": policy_id,
        "protocol_id": protocol_id,
        "sha256": sha256,
    }


def registry_payload(
    *,
    registry_id: str,
    snapshot_utc: str,
    seen_nonces: list[str],
    assignments: list[dict[str, str]],
    key_id: str,
) -> dict[str, Any]:
    """Return the exact signed replay-registry snapshot message."""
    return {
        "assignments": assignments,
        "key_id": key_id,
        "registry_id": registry_id,
        "seen_nonces": seen_nonces,
        "snapshot_utc": snapshot_utc,
    }


def _analysis_result(
    *,
    bundle: dict[str, Any],
    artifacts: dict[str, dict[str, Any]],
    artifact_paths: dict[str, Path],
) -> tuple[bool, str]:
    analysis_id = bundle["analysis_binding"]["analysis_artifact_id"]
    path = artifact_paths.get(analysis_id)
    if path is None:
        return False, "analysis artifact is unavailable"
    try:
        recipe = _read_json(path)
    except (OSError, json.JSONDecodeError) as exc:
        return False, f"analysis artifact is not a declarative JSON recipe: {exc}"
    error = _closed_object(
        recipe,
        required={
            "schema",
            "raw_artifact_ids",
            "sample_field",
            "effect_statement",
            "unit",
            "uncertainty_kind",
        },
        label="analysis recipe",
    )
    if error:
        return False, error
    if recipe["schema"] != "oph.hardware_evidence_bundle_h.analysis.mean_max_deviation.v1":
        return False, "analysis recipe schema is unsupported"
    if recipe["uncertainty_kind"] != "max_absolute_deviation":
        return False, "analysis uncertainty rule is unsupported"
    raw_ids = recipe["raw_artifact_ids"]
    if (
        not isinstance(raw_ids, list)
        or not raw_ids
        or len(set(raw_ids)) != len(raw_ids)
        or not all(isinstance(item, str) for item in raw_ids)
    ):
        return False, "analysis raw_artifact_ids must be a nonempty unique string list"
    declared_raw_ids = {
        artifact_id
        for artifact_id, row in artifacts.items()
        if row["role"] == "raw_measurement"
    }
    scheduled_raw_ids = {
        artifact_id
        for run in bundle["runs"]
        for artifact_id in run["raw_artifact_ids"]
    }
    if set(raw_ids) != declared_raw_ids or set(raw_ids) != scheduled_raw_ids:
        return False, (
            "analysis recipe must consume every and only the raw artifacts "
            "bound to the reported run population"
        )
    if not set(raw_ids).issubset(set(bundle["analysis_binding"]["input_artifact_ids"])):
        return False, "analysis recipe consumes an undeclared input artifact"
    samples: list[Fraction] = []
    sample_field = recipe["sample_field"]
    if not isinstance(sample_field, str) or not sample_field:
        return False, "analysis sample_field must be a nonempty string"
    for artifact_id in raw_ids:
        raw_path = artifact_paths.get(artifact_id)
        if raw_path is None:
            return False, f"analysis raw artifact {artifact_id} is unavailable"
        try:
            raw = _read_json(raw_path)
        except (OSError, json.JSONDecodeError) as exc:
            return False, f"analysis raw artifact {artifact_id} is invalid: {exc}"
        values = raw.get(sample_field) if isinstance(raw, dict) else None
        if not isinstance(values, list) or not values:
            return False, f"analysis raw artifact {artifact_id} has no samples"
        try:
            samples.extend(Fraction(str(value)) for value in values)
        except (ValueError, ZeroDivisionError):
            return False, f"analysis raw artifact {artifact_id} has a nonnumeric sample"
    mean = sum(samples, Fraction(0)) / len(samples)
    uncertainty = max(abs(value - mean) for value in samples)
    claim = bundle["claim"]
    try:
        claimed_mean = Fraction(str(claim["magnitude"]["value"]))
        claimed_uncertainty = Fraction(str(claim["uncertainty"]["value"]))
    except (ValueError, ZeroDivisionError):
        return False, "structured claim has a nonnumeric magnitude or uncertainty"
    expected = (
        recipe["effect_statement"],
        recipe["unit"],
        recipe["uncertainty_kind"],
        mean,
        uncertainty,
    )
    actual = (
        claim["effect_statement"],
        claim["magnitude"]["unit"],
        claim["uncertainty"]["kind"],
        claimed_mean,
        claimed_uncertainty,
    )
    if actual != expected or claim["uncertainty"]["unit"] != recipe["unit"]:
        return False, (
            "executed analysis result does not equal the structured "
            "effect, magnitude, unit, and uncertainty"
        )
    return True, (
        f"declarative analysis reproduced mean={mean} and "
        f"max_absolute_deviation={uncertainty}"
    )


def verify_external_evidence(
    *,
    bundle: dict[str, Any],
    artifacts: dict[str, dict[str, Any]],
    artifact_paths: dict[str, Path],
    replay_registry_path: Path | None,
    trust_policy_path: Path | None,
    external_evidence_path: Path | None,
) -> dict[str, Any]:
    """Verify operator-pinned external evidence for one internally valid bundle."""
    result: dict[str, Any] = {
        "invalid_codes": set(),
        "insufficient_codes": set(),
        "details": [],
        "predicate_failures": set(),
        "open_gates": [],
        "verified_anchor_ids": [],
    }

    def invalid(code: str, detail: str, predicate: str | None = None) -> None:
        result["invalid_codes"].add(code)
        result["details"].append(f"{code}: {detail}")
        if predicate:
            result["predicate_failures"].add(predicate)

    def insufficient(code: str, detail: str, predicate: str) -> None:
        result["insufficient_codes"].add(code)
        result["details"].append(f"{code}: {detail}")
        result["predicate_failures"].add(predicate)

    if trust_policy_path is None or external_evidence_path is None:
        missing = []
        if trust_policy_path is None:
            missing.append("operator trust policy")
        if external_evidence_path is None:
            missing.append("external evidence")
        for code, predicate in (
            (OPEN_GATES[0], "attestation"),
            (OPEN_GATES[1], "attestation"),
            (OPEN_GATES[2], "analysis_binding"),
            (OPEN_GATES[3], "completeness"),
            (OPEN_GATES[4], "device_identity"),
            (OPEN_GATES[5], "replay_protection"),
        ):
            insufficient(code, f"missing {' and '.join(missing)}", predicate)
        result["open_gates"] = list(OPEN_GATES)
        return result

    try:
        policy = _read_json(trust_policy_path)
        evidence = _read_json(external_evidence_path)
    except (OSError, json.JSONDecodeError) as exc:
        invalid("EXTERNAL_EVIDENCE_JSON_INVALID", str(exc))
        result["open_gates"] = list(OPEN_GATES)
        return result

    policy_error = _closed_object(
        policy,
        required={"schema", "policy_id", "anchors"},
        label="trust policy",
    )
    evidence_error = _closed_object(
        evidence,
        required={"schema", "policy_id", "signatures", "commitments", "nonce_registry"},
        label="external evidence",
    )
    if policy_error or evidence_error:
        invalid(
            "EXTERNAL_EVIDENCE_SHAPE_INVALID",
            "; ".join(item for item in (policy_error, evidence_error) if item),
        )
        result["open_gates"] = list(OPEN_GATES)
        return result
    if policy["schema"] != TRUST_SCHEMA or evidence["schema"] != EVIDENCE_SCHEMA:
        invalid("EXTERNAL_EVIDENCE_SCHEMA_INVALID", "unsupported trust/evidence schema")
        result["open_gates"] = list(OPEN_GATES)
        return result
    if evidence["policy_id"] != policy["policy_id"]:
        invalid("TRUST_POLICY_ID_MISMATCH", "external evidence names another policy")
        result["open_gates"] = list(OPEN_GATES)
        return result

    anchors: dict[str, dict[str, Any]] = {}
    public_keys: dict[str, Ed25519PublicKey] = {}
    anchor_required = {
        "key_id",
        "party_id",
        "organization_id",
        "role",
        "public_key_ed25519",
        "valid_from_utc",
        "valid_until_utc",
        "revoked",
    }
    for index, anchor in enumerate(policy["anchors"] if isinstance(policy["anchors"], list) else []):
        error = _closed_object(
            anchor,
            required=anchor_required,
            label=f"anchor[{index}]",
        )
        if error:
            invalid("TRUST_ANCHOR_INVALID", error, "attestation")
            continue
        key_id = anchor["key_id"]
        if not isinstance(key_id, str) or not key_id or key_id in anchors:
            invalid("TRUST_ANCHOR_INVALID", "anchor key ids must be unique strings", "attestation")
            continue
        if anchor["role"] not in ANCHOR_ROLES or not isinstance(anchor["revoked"], bool):
            invalid("TRUST_ANCHOR_INVALID", f"anchor {key_id} has an invalid role/status", "attestation")
            continue
        try:
            valid_from = _parse_utc(anchor["valid_from_utc"])
            valid_until = _parse_utc(anchor["valid_until_utc"])
            key = _decode_public_key(anchor["public_key_ed25519"])
        except (TypeError, ValueError) as exc:
            invalid("TRUST_ANCHOR_INVALID", f"anchor {key_id}: {exc}", "attestation")
            continue
        evaluation_time = _parse_utc(bundle["created_utc"])
        if anchor["revoked"] or not valid_from <= evaluation_time <= valid_until:
            insufficient(
                "TRUST_ANCHOR_NOT_VALID",
                f"anchor {key_id} is revoked or outside its validity interval",
                "attestation",
            )
            continue
        anchors[key_id] = anchor
        public_keys[key_id] = key

    if not anchors:
        insufficient(OPEN_GATES[0], "no valid operator-pinned trust anchor", "attestation")
        result["open_gates"] = list(OPEN_GATES)
        return result

    signatures: list[dict[str, Any]] = (
        evidence["signatures"] if isinstance(evidence["signatures"], list) else []
    )
    verified_signatures: list[tuple[dict[str, Any], dict[str, Any]]] = []
    signature_required = {
        "subject_kind",
        "subject_id",
        "sha256",
        "key_id",
        "signature_base64",
    }
    valid_subjects: dict[tuple[str, str], str] = {
        ("bundle_binding", bundle["bundle_id"]): bundle["bundle_binding"]["binding_sha256"],
    }
    valid_subjects.update(
        {
            ("artifact", artifact_id): row["sha256"]
            for artifact_id, row in artifacts.items()
        }
    )
    registry_digest = (
        sha256_path(replay_registry_path)
        if replay_registry_path is not None and replay_registry_path.is_file()
        else None
    )
    if registry_digest is not None:
        valid_subjects[("replay_registry", bundle["replay_protection"]["registry_id"])] = registry_digest

    for index, row in enumerate(signatures):
        error = _closed_object(
            row,
            required=signature_required,
            label=f"signature[{index}]",
        )
        if error:
            invalid("SIGNATURE_ENVELOPE_INVALID", error, "attestation")
            continue
        subject = (row["subject_kind"], row["subject_id"])
        expected_digest = valid_subjects.get(subject)
        if expected_digest is None or row["sha256"] != expected_digest:
            invalid(
                "SIGNATURE_SUBJECT_MISMATCH",
                f"signature[{index}] does not bind a current subject digest",
                "attestation",
            )
            continue
        anchor = anchors.get(row["key_id"])
        if anchor is None:
            insufficient(
                "SIGNATURE_KEY_UNTRUSTED",
                f"signature[{index}] key is absent from the operator policy",
                "attestation",
            )
            continue
        payload = signature_payload(
            subject_kind=row["subject_kind"],
            subject_id=row["subject_id"],
            sha256=row["sha256"],
            key_id=row["key_id"],
        )
        if not _verify_signature(public_keys[row["key_id"]], row["signature_base64"], payload):
            invalid(
                "SIGNATURE_INVALID",
                f"signature[{index}] failed Ed25519 verification",
                "attestation",
            )
            continue
        verified_signatures.append((row, anchor))

    def matching_signatures(
        subject_kind: str,
        subject_id: str,
        role: str,
    ) -> list[tuple[dict[str, Any], dict[str, Any]]]:
        return [
            (row, anchor)
            for row, anchor in verified_signatures
            if row["subject_kind"] == subject_kind
            and row["subject_id"] == subject_id
            and anchor["role"] == role
        ]

    claimant_roots = [
        anchor
        for row, anchor in matching_signatures(
            "bundle_binding", bundle["bundle_id"], "claimant"
        )
        if anchor["party_id"] == bundle["claimant_id"]
    ]
    if not claimant_roots:
        insufficient(
            OPEN_GATES[0],
            "no valid claimant signature binds the canonical bundle root",
            "attestation",
        )

    provenance_missing = False
    attestor_matches: list[tuple[dict[str, Any], dict[str, Any]]] = []
    for artifact_id, row in artifacts.items():
        role = ARTIFACT_AUTHORITY.get(row["role"])
        if role is None:
            continue
        matches = matching_signatures("artifact", artifact_id, role)
        if not matches:
            provenance_missing = True
            insufficient(
                "ARTIFACT_PROVENANCE_SIGNATURE_MISSING",
                f"{artifact_id} lacks a valid {role} signature",
                {
                    "calibration": "calibration_chain",
                    "custody": "custody",
                    "device_identity": "device_identity",
                    "raw_measurement": "raw_capture",
                    "control": "controls",
                    "analysis_code": "analysis_binding",
                    "claim_text": "analysis_binding",
                    "protocol": "analysis_binding",
                    "attestation": "attestation",
                }[row["role"]],
            )
        if row["role"] == "attestation":
            attestor_matches.extend(matches)
    if provenance_missing:
        insufficient(
            OPEN_GATES[4],
            "device, measurement, calibration, custody, and bundle provenance is incomplete",
            "device_identity",
        )

    claimant_organizations = {
        anchor["organization_id"] for anchor in claimant_roots
    }
    signer_index = {
        row["signer_id"]: row
        for row in bundle["signers"]
    }
    declared_attestor_keys = {
        row["key_id"]
        for signer_id in bundle["attestation"]["signer_ids"]
        if (row := signer_index.get(signer_id)) is not None
        and row["party_id"] in bundle["attestation"]["party_ids"]
        and not row["compromised"]
    }
    independent_attestors = [
        anchor
        for _, anchor in attestor_matches
        if anchor["key_id"] in declared_attestor_keys
        and anchor["party_id"] != bundle["claimant_id"]
        and anchor["organization_id"] not in claimant_organizations
    ]
    witness_root_signatures = [
        anchor
        for _, anchor in matching_signatures(
            "bundle_binding",
            bundle["bundle_id"],
            "independent_attestor",
        )
        if anchor["key_id"] in declared_attestor_keys
        and anchor["party_id"] != bundle["claimant_id"]
        and anchor["organization_id"] not in claimant_organizations
    ]
    attestation_mode = bundle["attestation"]["mode"]
    if attestation_mode == "independent_end_to_end_witness":
        if (
            not independent_attestors
            or not witness_root_signatures
            or not bundle["attestation"]["witnessed_end_to_end"]
        ):
            insufficient(
                OPEN_GATES[1],
                (
                    "no independently administered witness signs both the "
                    "attestation artifact and the canonical bundle root"
                ),
                "attestation",
            )
    elif attestation_mode == "independent_reproduction":
        insufficient(
            OPEN_GATES[1],
            "a separate fresh-run bundle is named but not resolved by this evidence packet",
            "attestation",
        )
    else:
        insufficient(
            OPEN_GATES[1],
            "physical promotion requires a verified independent witness or reproduction",
            "attestation",
        )

    commitments = (
        evidence["commitments"] if isinstance(evidence["commitments"], list) else []
    )
    commitment_required = {
        "kind",
        "artifact_id",
        "sha256",
        "committed_utc",
        "key_id",
        "policy_id",
        "bundle_id",
        "claimant_id",
        "device_id",
        "protocol_id",
        "signature_base64",
    }
    expected_commitments = {
        "run_schedule": bundle["run_schedule"]["artifact_id"],
        "analysis": bundle["analysis_binding"]["analysis_artifact_id"],
    }
    verified_commitment_kinds: set[str] = set()
    verified_commitment_anchor_ids: set[str] = set()
    first_capture = min(_parse_utc(run["captured_utc"]) for run in bundle["runs"])
    for index, row in enumerate(commitments):
        error = _closed_object(
            row,
            required=commitment_required,
            label=f"commitment[{index}]",
        )
        if error:
            invalid("COMMITMENT_INVALID", error, "completeness")
            continue
        expected_id = expected_commitments.get(row["kind"])
        artifact = artifacts.get(row["artifact_id"])
        if (
            expected_id != row["artifact_id"]
            or artifact is None
            or row["sha256"] != artifact["sha256"]
            or row["policy_id"] != policy["policy_id"]
            or row["bundle_id"] != bundle["bundle_id"]
            or row["claimant_id"] != bundle["claimant_id"]
            or row["device_id"] != bundle["claim"]["device_id"]
            or row["protocol_id"] != bundle["claim"]["protocol_id"]
        ):
            invalid(
                "COMMITMENT_SUBJECT_MISMATCH",
                f"commitment[{index}] does not bind the declared artifact",
                "completeness",
            )
            continue
        anchor = anchors.get(row["key_id"])
        if anchor is None or anchor["role"] != "preregistration_authority":
            insufficient(
                "COMMITMENT_KEY_UNTRUSTED",
                f"commitment[{index}] lacks a preregistration authority",
                "completeness",
            )
            continue
        try:
            committed = _parse_utc(row["committed_utc"])
        except (TypeError, ValueError) as exc:
            invalid("COMMITMENT_INVALID", str(exc), "completeness")
            continue
        payload = commitment_payload(
            kind=row["kind"],
            artifact_id=row["artifact_id"],
            sha256=row["sha256"],
            committed_utc=row["committed_utc"],
            key_id=row["key_id"],
            policy_id=row["policy_id"],
            bundle_id=row["bundle_id"],
            claimant_id=row["claimant_id"],
            device_id=row["device_id"],
            protocol_id=row["protocol_id"],
        )
        if not _verify_signature(public_keys[row["key_id"]], row["signature_base64"], payload):
            invalid("COMMITMENT_SIGNATURE_INVALID", f"commitment[{index}]", "completeness")
            continue
        if committed >= first_capture:
            insufficient(
                "COMMITMENT_NOT_PRE_RUN",
                f"{row['kind']} was not committed before capture",
                "completeness" if row["kind"] == "run_schedule" else "analysis_binding",
            )
            continue
        anchor_valid_from = _parse_utc(anchor["valid_from_utc"])
        anchor_valid_until = _parse_utc(anchor["valid_until_utc"])
        if not anchor_valid_from <= committed <= anchor_valid_until:
            insufficient(
                "COMMITMENT_KEY_NOT_VALID_AT_EVENT",
                f"{row['kind']} signer was not valid at commitment time",
                "completeness" if row["kind"] == "run_schedule" else "analysis_binding",
            )
            continue
        verified_commitment_kinds.add(row["kind"])
        verified_commitment_anchor_ids.add(row["key_id"])
    if verified_commitment_kinds != set(expected_commitments):
        insufficient(
            OPEN_GATES[3],
            "run schedule and analysis require signed pre-run commitments",
            "completeness",
        )

    analysis_ok, analysis_detail = _analysis_result(
        bundle=bundle,
        artifacts=artifacts,
        artifact_paths=artifact_paths,
    )
    if not analysis_ok:
        insufficient(OPEN_GATES[2], analysis_detail, "analysis_binding")
    else:
        result["details"].append(f"ANALYSIS_REPLAY_OK: {analysis_detail}")

    registry = evidence["nonce_registry"]
    registry_error = _closed_object(
        registry,
        required={
            "registry_id",
            "snapshot_utc",
            "seen_nonces",
            "assignments",
            "key_id",
            "signature_base64",
        },
        label="nonce_registry",
    )
    registry_ok = False
    replay_anchor_verified: str | None = None
    if registry_error:
        invalid("REPLAY_REGISTRY_EVIDENCE_INVALID", registry_error, "replay_protection")
    elif registry["registry_id"] != bundle["replay_protection"]["registry_id"]:
        invalid(
            "REPLAY_REGISTRY_EVIDENCE_INVALID",
            "registry id differs from the bundle",
            "replay_protection",
        )
    else:
        anchor = anchors.get(registry["key_id"])
        assignments = registry["assignments"]
        payload = registry_payload(
            registry_id=registry["registry_id"],
            snapshot_utc=registry["snapshot_utc"],
            seen_nonces=registry["seen_nonces"],
            assignments=assignments,
            key_id=registry["key_id"],
        )
        if anchor is None or anchor["role"] != "replay_authority":
            insufficient(
                "REPLAY_REGISTRY_KEY_UNTRUSTED",
                "registry signer is not the pinned replay authority",
                "replay_protection",
            )
        elif not _verify_signature(
            public_keys[registry["key_id"]],
            registry["signature_base64"],
            payload,
        ):
            invalid(
                "REPLAY_REGISTRY_SIGNATURE_INVALID",
                "registry snapshot signature failed",
                "replay_protection",
            )
        elif replay_registry_path is None or not replay_registry_path.is_file():
            insufficient(
                "REPLAY_REGISTRY_REQUIRED",
                "the signed registry snapshot file is unavailable",
                "replay_protection",
            )
        else:
            replay_anchor_verified = registry["key_id"]
            try:
                stored_registry = _read_json(replay_registry_path)
            except (OSError, json.JSONDecodeError) as exc:
                invalid("REPLAY_REGISTRY_INVALID", str(exc), "replay_protection")
                stored_registry = None
            if isinstance(stored_registry, dict) and (
                stored_registry.get("registry_id") != registry["registry_id"]
                or stored_registry.get("seen_nonces") != registry["seen_nonces"]
                or stored_registry.get("assignments") != assignments
            ):
                invalid(
                    "REPLAY_REGISTRY_SNAPSHOT_MISMATCH",
                    "signed assignments differ from the supplied registry",
                    "replay_protection",
                )
            elif isinstance(stored_registry, dict):
                try:
                    snapshot = _parse_utc(registry["snapshot_utc"])
                    last_capture = max(
                        _parse_utc(run["captured_utc"]) for run in bundle["runs"]
                    )
                    publication = _parse_utc(bundle["created_utc"])
                except (TypeError, ValueError) as exc:
                    invalid(
                        "REPLAY_REGISTRY_EVIDENCE_INVALID",
                        f"registry snapshot time is invalid: {exc}",
                        "replay_protection",
                    )
                    snapshot = None
                if snapshot is not None and not last_capture < snapshot <= publication:
                    insufficient(
                        "REPLAY_REGISTRY_SNAPSHOT_STALE",
                        "registry snapshot must follow capture and precede publication",
                        "replay_protection",
                    )
                if snapshot is not None:
                    anchor_valid_from = _parse_utc(anchor["valid_from_utc"])
                    anchor_valid_until = _parse_utc(anchor["valid_until_utc"])
                    if not anchor_valid_from <= snapshot <= anchor_valid_until:
                        insufficient(
                            "REPLAY_KEY_NOT_VALID_AT_EVENT",
                            "replay authority was not valid at snapshot time",
                            "replay_protection",
                        )
                        snapshot = None
                seen_nonces = registry["seen_nonces"]
                if (
                    not isinstance(seen_nonces, list)
                    or not all(isinstance(nonce, str) for nonce in seen_nonces)
                    or len(seen_nonces) != len(set(seen_nonces))
                ):
                    invalid(
                        "REPLAY_REGISTRY_EVIDENCE_INVALID",
                        "seen_nonces must be a unique string list",
                        "replay_protection",
                    )
                    snapshot = None
                assignments_by_nonce: dict[str, list[dict[str, str]]] = {}
                for assignment in assignments if isinstance(assignments, list) else []:
                    error = _closed_object(
                        assignment,
                        required={
                            "nonce",
                            "bundle_id",
                            "binding_sha256",
                            "registered_utc",
                            "state",
                            "consumed_utc",
                        },
                        label="nonce assignment",
                    )
                    if error:
                        invalid(
                            "REPLAY_REGISTRY_ASSIGNMENT_INVALID",
                            error,
                            "replay_protection",
                        )
                        continue
                    assignments_by_nonce.setdefault(assignment["nonce"], []).append(assignment)
                registry_ok = snapshot is not None and last_capture < snapshot <= publication
                for run in bundle["runs"]:
                    matches = assignments_by_nonce.get(run["capture_nonce"], [])
                    assignment_ok = False
                    if len(matches) == 1:
                        assignment = matches[0]
                        try:
                            registered = _parse_utc(assignment["registered_utc"])
                            consumed = _parse_utc(assignment["consumed_utc"])
                            captured = _parse_utc(run["captured_utc"])
                            assignment_ok = (
                                assignment["bundle_id"] == bundle["bundle_id"]
                                and assignment["binding_sha256"]
                                == bundle["bundle_binding"]["binding_sha256"]
                                and assignment["state"] == "consumed"
                                and registered < captured <= consumed
                                and snapshot is not None
                                and consumed <= snapshot
                            )
                        except (TypeError, ValueError):
                            assignment_ok = False
                    if not assignment_ok:
                        registry_ok = False
                        insufficient(
                            "REPLAY_NONCE_NOT_ATOMICALLY_CONSUMED",
                            (
                                f"nonce {run['capture_nonce']} lacks one "
                                "pre-run reservation and post-capture consume receipt"
                            ),
                            "replay_protection",
                        )
    if not registry_ok:
        insufficient(
            OPEN_GATES[5],
            "authoritative replay-registry verification is incomplete",
            "replay_protection",
        )

    result["verified_anchor_ids"] = sorted(
        {anchor["key_id"] for _, anchor in verified_signatures}
        | verified_commitment_anchor_ids
        | ({replay_anchor_verified} if replay_anchor_verified else set())
    )
    result["open_gates"] = sorted(
        code for code in OPEN_GATES if code in result["insufficient_codes"]
    )
    return result
