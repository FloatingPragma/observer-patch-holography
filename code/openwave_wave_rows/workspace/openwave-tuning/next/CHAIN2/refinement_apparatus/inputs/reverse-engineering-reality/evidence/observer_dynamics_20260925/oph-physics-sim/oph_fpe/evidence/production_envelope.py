"""Immutable, fail-closed production-envelope inventory replay.

This module implements the P0 evidence boundary shared by the physical A5/SM
and W/Z source-to-pole programs.  It verifies bytes, identities, provenance,
runtime subject/output bindings, and envelope ancestry.  It intentionally does
*not* verify scientific claims: the scientific-checker registry is empty, a
producer ``PASS`` is untrusted, and promotion is always false.

Only an on-disk manifest is accepted.  No manifest-provided code is executed.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import stat
from types import MappingProxyType
from typing import Any, Mapping, Sequence


PRODUCTION_BUNDLE_MANIFEST_SCHEMA = "oph.production-evidence.bundle-manifest/1.0.0"
PRODUCTION_ENVELOPE_SCHEMA = "oph.production-evidence.envelope/1.0.0"
FREEZE_ANCHOR_SCHEMA = "oph.production-evidence.freeze-anchor/1.0.0"
PRODUCTION_BUNDLE_REPORT_SCHEMA = "oph.production-evidence.bundle-report/1.0.0"
PRODUCTION_BUNDLE_REPORT_ARTIFACT_TYPE = "OPH_PRODUCTION_BUNDLE_INVENTORY_REPLAY"

COMMON_STAGE = "COMMON_STAGE"
WZ_SOURCE_TO_POLE = "WZ_SOURCE_TO_POLE"
PROFILES = frozenset({COMMON_STAGE, WZ_SOURCE_TO_POLE})
CANONICAL_JSON_POLICY = "OPH_CANONICAL_JSON_V1"

ALLOWED_STATUSES = frozenset(
    {"PASS", "OPEN", "UNRESOLVED", "FAIL", "NOT_APPLICABLE"}
)
WZ_SHARED_HASH_FIELDS = (
    "action_ast_hash",
    "field_census_hash",
    "scheme_hash",
    "fj_convention_hash",
    "term_mask_hash",
    "analytic_sheet_hash",
    "units_basis_hash",
)
TARGET_FIREWALL_FALSE_FIELDS = (
    "source_ancestry_contains_target",
    "target_used_to_select_candidate",
    "target_used_to_tune_producer",
    "target_used_to_tune_checker",
    "target_bytes_available_to_producer",
    "target_bytes_available_to_checker",
    "comparison_process_can_mutate_bundle",
)

P0_RECEIPT_KEYS = (
    "P0_MANIFEST_ON_DISK_RECEIPT",
    "P0_MANIFEST_SCHEMA_RECEIPT",
    "P0_MANIFEST_SELF_HASH_RECEIPT",
    "P0_FILE_CONTAINMENT_RECEIPT",
    "P0_FILE_HASH_REPLAY_RECEIPT",
    "P0_FILE_BYTE_COUNT_REPLAY_RECEIPT",
    "P0_FILE_MEDIA_TYPE_RECEIPT",
    "P0_EXACT_FILE_CENSUS_RECEIPT",
    "P0_DUPLICATE_KEY_AND_FINITE_JSON_RECEIPT",
    "P0_ENVELOPE_SCHEMA_RECEIPT",
    "P0_EXACT_SCHEMA_BINDING_RECEIPT",
    "P0_RUNTIME_SUBJECT_BINDING_RECEIPT",
    "P0_OUTPUT_BINDING_RECEIPT",
    "P0_PROVENANCE_BINDING_RECEIPT",
    "P0_FREEZE_PREOUTCOME_ANCHOR_RECEIPT",
    "P0_PARENT_ENVELOPE_RESOLUTION_RECEIPT",
    "P0_ANCESTRY_DAG_RECEIPT",
    "P0_SINGLE_SOURCE_BRANCH_FREEZE_FAMILY_RECEIPT",
    "P0_PRODUCER_CHECKER_SEPARATION_RECEIPT",
    "P0_TARGET_FIREWALL_DECLARATION_RECEIPT",
    "P0_WZ_SHARED_HASH_FAMILY_RECEIPT",
    "P0_IMMUTABLE_INVENTORY_REPLAY_RECEIPT",
    "SCIENTIFIC_REPLAY_RECEIPT",
    "PROMOTION_RECEIPT",
)

REGISTERED_SCIENTIFIC_CHECKERS: Mapping[str, object] = MappingProxyType({})

MAX_MANIFEST_BYTES = 8 * 1024 * 1024
MAX_ARTIFACT_BYTES = 512 * 1024 * 1024
MAX_FILES = 8192
MAX_ENVELOPES = 4096

_SHA256_RE = re.compile(r"^sha256:[0-9a-f]{64}$")
_IDENTIFIER_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.:+/\-]{0,191}$")
_MEDIA_TYPE_RE = re.compile(r"^[a-z0-9][a-z0-9.+-]{0,63}/[a-z0-9][a-z0-9.+-]{0,95}$")
_UTC_RE = re.compile(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z$")
_JSON_MEDIA_TYPES = frozenset({"application/json", "application/schema+json"})

_MANIFEST_FIELDS = frozenset(
    {"schema", "bundle_id", "manifest_payload_sha256", "files", "envelope_paths"}
)
_FILE_FIELDS = frozenset({"path", "sha256", "byte_count", "media_type"})
_REF_FIELDS = _FILE_FIELDS
_ENVELOPE_FIELDS = frozenset(
    {
        "envelope_schema",
        "profile",
        "schema_id",
        "schema_version",
        "schema_sha256",
        "schema_ref",
        "artifact_id",
        "stage_id",
        "receipt_type",
        "claim_lane",
        "claim_scope",
        "branch_id",
        "freeze_id",
        "source_root_hash",
        "subject_type",
        "subject_canonicalization",
        "subject_ref",
        "subject_digest",
        "output_ref",
        "output_digest",
        "producer",
        "checker",
        "freeze_anchor_ref",
        "parent_receipts",
        "shared_contract_hashes",
        "numeric_precision_and_rounding",
        "target_firewall",
        "status",
        "blockers",
        "generated_utc",
    }
)
_PRODUCER_FIELDS = frozenset(
    {"producer_id", "source_tree_ref", "executable_ref", "environment_lock_ref"}
)
_CHECKER_FIELDS = frozenset(
    {"checker_id", "source_tree_ref", "executable_ref", "checker_independence_class"}
)
_PARENT_FIELDS = frozenset(
    {"artifact_id", "envelope_ref", "subject_digest", "output_digest"}
)
_NUMERIC_FIELDS = frozenset(
    {"arithmetic", "precision_bits", "rounding_mode", "numeric_backend"}
)
_TARGET_FIELDS = frozenset(
    set(TARGET_FIREWALL_FALSE_FIELDS) | {"comparison_boundary", "exposure_status"}
)
_ANCHOR_FIELDS = frozenset(
    {"schema", "source_root_hash", "branch_id", "freeze_id", "commitment_phase", "frozen_utc"}
)
_ENVELOPE_RECEIPT_KEYS = (
    "P0_ENVELOPE_SCHEMA_RECEIPT",
    "P0_EXACT_SCHEMA_BINDING_RECEIPT",
    "P0_RUNTIME_SUBJECT_BINDING_RECEIPT",
    "P0_OUTPUT_BINDING_RECEIPT",
    "P0_PROVENANCE_BINDING_RECEIPT",
    "P0_FREEZE_PREOUTCOME_ANCHOR_RECEIPT",
    "P0_PARENT_ENVELOPE_RESOLUTION_RECEIPT",
    "P0_ANCESTRY_DAG_RECEIPT",
    "P0_SINGLE_SOURCE_BRANCH_FREEZE_FAMILY_RECEIPT",
    "P0_PRODUCER_CHECKER_SEPARATION_RECEIPT",
    "P0_TARGET_FIREWALL_DECLARATION_RECEIPT",
    "P0_WZ_SHARED_HASH_FAMILY_RECEIPT",
)


class ProductionEnvelopeError(ValueError):
    """Malformed or unreplayable production evidence."""


@dataclass(frozen=True)
class _FileSpec:
    path: str
    sha256: str
    byte_count: int
    media_type: str


@dataclass
class _EnvelopeState:
    path: str
    file_sha256: str | None = None
    payload: Mapping[str, Any] | None = None
    artifact_id: str | None = None
    stage_id: str | None = None
    receipt_type: str | None = None
    profile: str | None = None
    claim_lane: str | None = None
    claim_scope: str | None = None
    source_root_hash: str | None = None
    branch_id: str | None = None
    freeze_id: str | None = None
    freeze_anchor_sha256: str | None = None
    frozen_utc: str | None = None
    generated_utc: str | None = None
    producer_status: str | None = None
    subject_digest: str | None = None
    output_digest: str | None = None
    parent_rows: list[Mapping[str, Any]] = field(default_factory=list)
    parent_artifact_ids: list[str] = field(default_factory=list)
    shared_contract_hashes: dict[str, str] = field(default_factory=dict)
    checks: dict[str, bool] = field(
        default_factory=lambda: {key: False for key in _ENVELOPE_RECEIPT_KEYS}
    )
    blockers: list[str] = field(default_factory=list)


class _Store:
    def __init__(self, manifest_path: Path, specs: Mapping[str, _FileSpec]):
        self.manifest_path = manifest_path
        self.root = manifest_path.parent.resolve(strict=True)
        self.specs = dict(specs)
        self.rows: dict[str, dict[str, Any]] = {}
        self.raw: dict[str, bytes] = {}
        self.json_values: dict[str, Any] = {}
        self.blockers: list[str] = []
        self.containment = True
        self.hashes = True
        self.byte_counts = True
        self.media_types = True
        self.exact_census = True
        self.json_safety = True

    def verify(self) -> None:
        self._verify_census()
        for relative in sorted(self.specs):
            self._verify_file(self.specs[relative])

    def _verify_census(self) -> None:
        seen: set[str] = set()
        symlinks: list[str] = []
        manifest_resolved = self.manifest_path.resolve(strict=True)
        for current, dirs, files in os.walk(self.root, followlinks=False):
            current_path = Path(current)
            for name in list(dirs):
                candidate = current_path / name
                if candidate.is_symlink():
                    symlinks.append(candidate.relative_to(self.root).as_posix())
            for name in files:
                candidate = current_path / name
                relative = candidate.relative_to(self.root).as_posix()
                if candidate.is_symlink():
                    symlinks.append(relative)
                    seen.add(relative)
                    continue
                try:
                    candidate_resolved = candidate.resolve(strict=True)
                except OSError as exc:
                    self.containment = False
                    self.blockers.append(
                        f"bundle_entry_resolution_failed:{relative}:{type(exc).__name__}"
                    )
                    seen.add(relative)
                    continue
                if candidate_resolved != manifest_resolved:
                    seen.add(relative)
        expected = set(self.specs)
        extras = sorted(seen - expected)
        missing = sorted(expected - seen)
        if extras:
            self.exact_census = False
            self.blockers.extend(f"unlisted_extra_file:{path}" for path in extras)
        if missing:
            self.exact_census = False
            self.blockers.extend(f"listed_file_missing_from_bundle:{path}" for path in missing)
        if symlinks:
            self.containment = False
            self.blockers.extend(f"bundle_symlink_forbidden:{path}" for path in sorted(symlinks))

    def _verify_file(self, spec: _FileSpec) -> None:
        row: dict[str, Any] = {
            "path": spec.path,
            "declared_sha256": spec.sha256,
            "declared_byte_count": spec.byte_count,
            "media_type": spec.media_type,
            "passed": False,
            "blockers": [],
        }
        try:
            path = _contained_regular_file(self.root, spec.path)
            raw = path.read_bytes()
            if len(raw) > MAX_ARTIFACT_BYTES:
                raise ProductionEnvelopeError("artifact_exceeds_size_limit")
            actual_hash = _sha256_bytes(raw)
            row["actual_sha256"] = actual_hash
            row["actual_byte_count"] = len(raw)
            mismatches: list[str] = []
            if actual_hash != spec.sha256:
                self.hashes = False
                mismatches.append("declared_sha256_mismatch")
            if len(raw) != spec.byte_count:
                self.byte_counts = False
                mismatches.append("declared_byte_count_mismatch")
            if mismatches:
                raise ProductionEnvelopeError("+".join(mismatches))
            self.raw[spec.path] = raw
            if spec.media_type in _JSON_MEDIA_TYPES or spec.media_type.endswith("+json"):
                try:
                    self.json_values[spec.path] = _strict_json_loads(raw)
                except (UnicodeDecodeError, json.JSONDecodeError, ProductionEnvelopeError) as exc:
                    self.json_safety = False
                    raise ProductionEnvelopeError(f"unsafe_json:{exc}") from exc
            row["passed"] = True
        except (OSError, ProductionEnvelopeError) as exc:
            blocker = f"file:{spec.path}:{exc}"
            row["blockers"].append(str(exc))
            self.blockers.append(blocker)
            if "symlink" in str(exc) or "contain" in str(exc) or "regular" in str(exc):
                self.containment = False
            if "sha256" in str(exc):
                self.hashes = False
            if "byte_count" in str(exc):
                self.byte_counts = False
        self.rows[spec.path] = row

    def resolve_ref(self, value: Any, label: str) -> tuple[_FileSpec, bytes]:
        ref = _mapping(value, label)
        _exact_fields(ref, _REF_FIELDS, label)
        path = _relative_path(ref.get("path"), f"{label}.path")
        declared = _FileSpec(
            path=path,
            sha256=_hash(ref.get("sha256"), f"{label}.sha256"),
            byte_count=_byte_count(ref.get("byte_count"), f"{label}.byte_count"),
            media_type=_media_type(ref.get("media_type"), f"{label}.media_type"),
        )
        spec = self.specs.get(path)
        if spec is None:
            raise ProductionEnvelopeError(f"{label}_path_not_in_manifest")
        if declared != spec:
            raise ProductionEnvelopeError(f"{label}_does_not_match_manifest_row")
        if path not in self.raw or self.rows.get(path, {}).get("passed") is not True:
            raise ProductionEnvelopeError(f"{label}_bytes_not_verified")
        return spec, self.raw[path]

    def resolve_json_ref(self, value: Any, label: str) -> tuple[_FileSpec, Any]:
        spec, _ = self.resolve_ref(value, label)
        if spec.media_type not in _JSON_MEDIA_TYPES and not spec.media_type.endswith("+json"):
            raise ProductionEnvelopeError(f"{label}_must_have_json_media_type")
        if spec.path not in self.json_values:
            raise ProductionEnvelopeError(f"{label}_json_not_available")
        return spec, self.json_values[spec.path]


def canonical_json_bytes(value: Any) -> bytes:
    """Return the frozen OPH canonical JSON encoding.

    This is a frozen equivalent rather than a claim of full RFC 8785 support:
    UTF-8, sorted keys, no insignificant whitespace, no ASCII escaping, and no
    non-finite numbers.  The policy name is hash-bound in every envelope.
    """

    _validate_finite(value)
    try:
        text = json.dumps(
            value,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        )
    except (TypeError, ValueError) as exc:
        raise ProductionEnvelopeError(f"value_not_canonical_json:{exc}") from exc
    return text.encode("utf-8")


def canonical_json_sha256(value: Any) -> str:
    """Hash :func:`canonical_json_bytes` with an explicit algorithm prefix."""

    return _sha256_bytes(canonical_json_bytes(value))


def verify_production_bundle_manifest(
    manifest_path: str | Path | Mapping[str, Any],
) -> dict[str, Any]:
    """Replay one immutable production bundle and return a stable report.

    A mapping is deliberately rejected even if it has the correct shape.  The
    manifest and all referenced bytes must already exist in one closed bundle.
    """

    if isinstance(manifest_path, Mapping):
        return _failure_report(
            ["manifest_must_be_on_disk_file:not_a_mapping"], manifest_path=None
        )
    try:
        path = _regular_nonsymlink_file(Path(manifest_path))
        if path.stat().st_size > MAX_MANIFEST_BYTES:
            raise ProductionEnvelopeError("manifest_exceeds_size_limit")
        raw = path.read_bytes()
        payload = _strict_json_loads(raw)
        manifest = _mapping(payload, "manifest")
        _exact_fields(manifest, _MANIFEST_FIELDS, "manifest")
        if manifest.get("schema") != PRODUCTION_BUNDLE_MANIFEST_SCHEMA:
            raise ProductionEnvelopeError("manifest_schema_mismatch")
        bundle_id = _identifier(manifest.get("bundle_id"), "bundle_id")
        expected_payload_hash = _hash(
            manifest.get("manifest_payload_sha256"), "manifest_payload_sha256"
        )
        unhashed = dict(manifest)
        unhashed.pop("manifest_payload_sha256", None)
        manifest_self_hash_valid = expected_payload_hash == canonical_json_sha256(unhashed)
        specs = _parse_file_specs(manifest.get("files"), path.name)
        envelope_paths = _parse_envelope_paths(manifest.get("envelope_paths"), specs)
    except (
        OSError,
        UnicodeDecodeError,
        json.JSONDecodeError,
        ProductionEnvelopeError,
        TypeError,
    ) as exc:
        return _failure_report(
            [f"manifest_parse_failed:{type(exc).__name__}:{exc}"],
            manifest_path=str(manifest_path),
        )

    store = _Store(path, specs)
    store.verify()
    states = [_parse_envelope(envelope_path, store) for envelope_path in envelope_paths]
    relationship_blockers, envelope_order, relationship_checks = _verify_relationships(
        states, store
    )

    receipts = {key: False for key in P0_RECEIPT_KEYS}
    receipts.update(
        {
            "P0_MANIFEST_ON_DISK_RECEIPT": True,
            "P0_MANIFEST_SCHEMA_RECEIPT": True,
            "P0_MANIFEST_SELF_HASH_RECEIPT": manifest_self_hash_valid,
            "P0_FILE_CONTAINMENT_RECEIPT": store.containment,
            "P0_FILE_HASH_REPLAY_RECEIPT": store.hashes
            and len(store.raw) == len(specs),
            "P0_FILE_BYTE_COUNT_REPLAY_RECEIPT": store.byte_counts
            and len(store.raw) == len(specs),
            "P0_FILE_MEDIA_TYPE_RECEIPT": store.media_types,
            "P0_EXACT_FILE_CENSUS_RECEIPT": store.exact_census,
            "P0_DUPLICATE_KEY_AND_FINITE_JSON_RECEIPT": store.json_safety,
        }
    )
    per_envelope_keys = (
        "P0_ENVELOPE_SCHEMA_RECEIPT",
        "P0_EXACT_SCHEMA_BINDING_RECEIPT",
        "P0_RUNTIME_SUBJECT_BINDING_RECEIPT",
        "P0_OUTPUT_BINDING_RECEIPT",
        "P0_PROVENANCE_BINDING_RECEIPT",
        "P0_FREEZE_PREOUTCOME_ANCHOR_RECEIPT",
        "P0_PRODUCER_CHECKER_SEPARATION_RECEIPT",
        "P0_TARGET_FIREWALL_DECLARATION_RECEIPT",
    )
    for key in per_envelope_keys:
        receipts[key] = bool(states) and all(state.checks.get(key) is True for state in states)
    receipts.update(relationship_checks)
    inventory_keys = tuple(
        key
        for key in P0_RECEIPT_KEYS
        if key
        not in {
            "P0_IMMUTABLE_INVENTORY_REPLAY_RECEIPT",
            "SCIENTIFIC_REPLAY_RECEIPT",
            "PROMOTION_RECEIPT",
        }
    )
    inventory_passed = all(receipts[key] is True for key in inventory_keys)
    receipts["P0_IMMUTABLE_INVENTORY_REPLAY_RECEIPT"] = inventory_passed
    receipts["SCIENTIFIC_REPLAY_RECEIPT"] = False
    receipts["PROMOTION_RECEIPT"] = False

    inventory_blockers = list(store.blockers) + relationship_blockers
    if not manifest_self_hash_valid:
        inventory_blockers.append("manifest_payload_sha256_mismatch")
    envelope_reports: dict[str, dict[str, Any]] = {}
    seen_keys: set[str] = set()
    for index, state in enumerate(states):
        state_inventory = inventory_passed and all(state.checks.values())
        key = state.artifact_id or f"INVALID_ENVELOPE_{index}"
        if key in seen_keys:
            key = f"{key}__DUPLICATE_{index}"
        seen_keys.add(key)
        inventory_blockers.extend(state.blockers)
        envelope_reports[key] = _envelope_report(state, state_inventory)

    inventory_blockers = sorted(set(inventory_blockers))
    open_requirements = ["no_registered_independent_scientific_checker"]
    return {
        "schema": PRODUCTION_BUNDLE_REPORT_SCHEMA,
        "artifact_type": PRODUCTION_BUNDLE_REPORT_ARTIFACT_TYPE,
        "manifest_path": path.name,
        "manifest_sha256": _sha256_bytes(raw),
        "manifest_payload_sha256": expected_payload_hash,
        "bundle_id": bundle_id,
        "canonical_json_policy": CANONICAL_JSON_POLICY,
        "envelope_order": envelope_order,
        "envelopes": envelope_reports,
        "file_rows": {key: store.rows[key] for key in sorted(store.rows)},
        "receipts": receipts,
        "inventory_replay_passed": inventory_passed,
        "scientific_replay_passed": False,
        "promotion_allowed": False,
        "passed": False,
        "status": "OPEN" if inventory_passed else "FAIL",
        "inventory_blockers": inventory_blockers,
        "open_requirements": open_requirements,
        "blockers": sorted(set(inventory_blockers + open_requirements)),
        "registered_scientific_checker_count": len(REGISTERED_SCIENTIFIC_CHECKERS),
        "ignored_producer_statuses": {
            key: row["producer_status"] for key, row in envelope_reports.items()
        },
        "claim_boundary": (
            "This report proves only immutable P0 inventory and binding consistency. "
            "Producer statuses and output booleans are untrusted. With no registered "
            "independent scientific checker, SCIENTIFIC_REPLAY and PROMOTION remain false."
        ),
    }


def _parse_envelope(path: str, store: _Store) -> _EnvelopeState:
    state = _EnvelopeState(path=path)
    try:
        spec = store.specs[path]
        state.file_sha256 = spec.sha256
        if spec.media_type != "application/json":
            raise ProductionEnvelopeError("envelope_media_type_must_be_application_json")
        payload = store.json_values.get(path)
        envelope = _mapping(payload, f"envelope:{path}")
        state.payload = envelope
        _exact_fields(envelope, _ENVELOPE_FIELDS, f"envelope:{path}")
        if envelope.get("envelope_schema") != PRODUCTION_ENVELOPE_SCHEMA:
            raise ProductionEnvelopeError("envelope_schema_mismatch")
        state.profile = _enum(envelope.get("profile"), PROFILES, "profile")
        state.artifact_id = _identifier(envelope.get("artifact_id"), "artifact_id")
        state.stage_id = _identifier(envelope.get("stage_id"), "stage_id")
        state.receipt_type = _identifier(envelope.get("receipt_type"), "receipt_type")
        state.claim_lane = _identifier(envelope.get("claim_lane"), "claim_lane")
        state.claim_scope = _identifier(envelope.get("claim_scope"), "claim_scope")
        state.branch_id = _identifier(envelope.get("branch_id"), "branch_id")
        state.freeze_id = _identifier(envelope.get("freeze_id"), "freeze_id")
        state.source_root_hash = _hash(envelope.get("source_root_hash"), "source_root_hash")
        _identifier(envelope.get("schema_id"), "schema_id")
        _identifier(envelope.get("schema_version"), "schema_version")
        _identifier(envelope.get("subject_type"), "subject_type")
        if envelope.get("subject_canonicalization") != CANONICAL_JSON_POLICY:
            raise ProductionEnvelopeError("subject_canonicalization_policy_mismatch")
        state.subject_digest = _hash(envelope.get("subject_digest"), "subject_digest")
        state.output_digest = _hash(envelope.get("output_digest"), "output_digest")
        state.producer_status = _enum(envelope.get("status"), ALLOWED_STATUSES, "status")
        _string_list(envelope.get("blockers"), "blockers")
        state.generated_utc = _utc(envelope.get("generated_utc"), "generated_utc")
        _validate_numeric_policy(envelope.get("numeric_precision_and_rounding"))
        parents = envelope.get("parent_receipts")
        if not isinstance(parents, list):
            raise ProductionEnvelopeError("parent_receipts_must_be_array")
        state.parent_rows = [_mapping(row, "parent_receipt") for row in parents]
        hashes = _mapping(envelope.get("shared_contract_hashes"), "shared_contract_hashes")
        state.shared_contract_hashes = {
            _identifier(key, "shared_contract_hash_key"): _hash(
                value, f"shared_contract_hashes.{key}"
            )
            for key, value in hashes.items()
        }
        state.checks["P0_ENVELOPE_SCHEMA_RECEIPT"] = True
    except (ProductionEnvelopeError, TypeError) as exc:
        state.blockers.append(f"envelope:{path}:schema:{exc}")
        return state

    envelope = state.payload
    assert envelope is not None
    try:
        schema_spec, schema_value = store.resolve_json_ref(
            envelope["schema_ref"], f"{state.artifact_id}.schema_ref"
        )
        if schema_spec.media_type != "application/schema+json":
            raise ProductionEnvelopeError("schema_ref_media_type_mismatch")
        if envelope["schema_sha256"] != schema_spec.sha256:
            raise ProductionEnvelopeError("schema_sha256_mismatch")
        schema_object = _mapping(schema_value, "resolved_schema")
        if schema_object.get("$id") != envelope["schema_id"]:
            raise ProductionEnvelopeError("resolved_schema_id_mismatch")
        if schema_object.get("schema_version") != envelope["schema_version"]:
            raise ProductionEnvelopeError("resolved_schema_version_mismatch")
        state.checks["P0_EXACT_SCHEMA_BINDING_RECEIPT"] = True
    except (ProductionEnvelopeError, TypeError) as exc:
        state.blockers.append(f"envelope:{state.artifact_id}:schema_binding:{exc}")

    try:
        _, subject = store.resolve_json_ref(
            envelope["subject_ref"], f"{state.artifact_id}.subject_ref"
        )
        if canonical_json_sha256(subject) != state.subject_digest:
            raise ProductionEnvelopeError("runtime_subject_digest_mismatch")
        state.checks["P0_RUNTIME_SUBJECT_BINDING_RECEIPT"] = True
    except (ProductionEnvelopeError, TypeError) as exc:
        state.blockers.append(f"envelope:{state.artifact_id}:subject_binding:{exc}")

    try:
        _, output = store.resolve_json_ref(
            envelope["output_ref"], f"{state.artifact_id}.output_ref"
        )
        if canonical_json_sha256(output) != state.output_digest:
            raise ProductionEnvelopeError("output_digest_mismatch")
        state.checks["P0_OUTPUT_BINDING_RECEIPT"] = True
    except (ProductionEnvelopeError, TypeError) as exc:
        state.blockers.append(f"envelope:{state.artifact_id}:output_binding:{exc}")

    producer_refs: list[_FileSpec] = []
    checker_refs: list[_FileSpec] = []
    try:
        producer = _mapping(envelope["producer"], "producer")
        checker = _mapping(envelope["checker"], "checker")
        _exact_fields(producer, _PRODUCER_FIELDS, "producer")
        _exact_fields(checker, _CHECKER_FIELDS, "checker")
        producer_id = _identifier(producer.get("producer_id"), "producer_id")
        checker_id = _identifier(checker.get("checker_id"), "checker_id")
        _identifier(checker.get("checker_independence_class"), "checker_independence_class")
        for key in ("source_tree_ref", "executable_ref", "environment_lock_ref"):
            spec, _ = store.resolve_ref(producer[key], f"{state.artifact_id}.producer.{key}")
            producer_refs.append(spec)
        for key in ("source_tree_ref", "executable_ref"):
            spec, _ = store.resolve_ref(checker[key], f"{state.artifact_id}.checker.{key}")
            checker_refs.append(spec)
        state.checks["P0_PROVENANCE_BINDING_RECEIPT"] = True
        producer_hashes = {item.sha256 for item in producer_refs[:2]}
        checker_hashes = {item.sha256 for item in checker_refs}
        if producer_id == checker_id:
            raise ProductionEnvelopeError("producer_and_checker_id_identical")
        if producer_hashes & checker_hashes:
            raise ProductionEnvelopeError("producer_checker_source_or_executable_hash_shared")
        if state.file_sha256 in checker_hashes:
            raise ProductionEnvelopeError("envelope_cannot_be_its_own_checker")
        state.checks["P0_PRODUCER_CHECKER_SEPARATION_RECEIPT"] = True
    except (ProductionEnvelopeError, TypeError) as exc:
        state.blockers.append(f"envelope:{state.artifact_id}:provenance:{exc}")

    try:
        anchor_spec, anchor_value = store.resolve_json_ref(
            envelope["freeze_anchor_ref"], f"{state.artifact_id}.freeze_anchor_ref"
        )
        state.freeze_anchor_sha256 = anchor_spec.sha256
        anchor = _mapping(anchor_value, "freeze_anchor")
        _exact_fields(anchor, _ANCHOR_FIELDS, "freeze_anchor")
        if anchor.get("schema") != FREEZE_ANCHOR_SCHEMA:
            raise ProductionEnvelopeError("freeze_anchor_schema_mismatch")
        if anchor.get("source_root_hash") != state.source_root_hash:
            raise ProductionEnvelopeError("freeze_anchor_source_root_mismatch")
        if anchor.get("branch_id") != state.branch_id:
            raise ProductionEnvelopeError("freeze_anchor_branch_mismatch")
        if anchor.get("freeze_id") != state.freeze_id:
            raise ProductionEnvelopeError("freeze_anchor_id_mismatch")
        if anchor.get("commitment_phase") != "pre_outcome":
            raise ProductionEnvelopeError("freeze_anchor_not_pre_outcome")
        state.frozen_utc = _utc(anchor.get("frozen_utc"), "freeze_anchor.frozen_utc")
        if state.generated_utc is None or not (
            _utc_instant(state.frozen_utc) < _utc_instant(state.generated_utc)
        ):
            raise ProductionEnvelopeError("freeze_anchor_not_before_envelope_generation")
        state.checks["P0_FREEZE_PREOUTCOME_ANCHOR_RECEIPT"] = True
    except (ProductionEnvelopeError, TypeError) as exc:
        state.blockers.append(f"envelope:{state.artifact_id}:freeze_anchor:{exc}")

    try:
        firewall = _mapping(envelope["target_firewall"], "target_firewall")
        _exact_fields(firewall, _TARGET_FIELDS, "target_firewall")
        for key in TARGET_FIREWALL_FALSE_FIELDS:
            if firewall.get(key) is not False:
                raise ProductionEnvelopeError(f"target_firewall_field_not_false:{key}")
        if firewall.get("comparison_boundary") != "read_only_separate_process":
            raise ProductionEnvelopeError("target_comparison_boundary_not_read_only")
        exposure = firewall.get("exposure_status")
        if exposure not in {"not_targeted", "prospective_blind", "post_exposure_validation"}:
            raise ProductionEnvelopeError("target_exposure_status_invalid")
        if state.profile == WZ_SOURCE_TO_POLE and exposure != "post_exposure_validation":
            raise ProductionEnvelopeError("wz_requires_post_exposure_validation")
        state.checks["P0_TARGET_FIREWALL_DECLARATION_RECEIPT"] = True
    except (ProductionEnvelopeError, TypeError) as exc:
        state.blockers.append(f"envelope:{state.artifact_id}:target_firewall:{exc}")

    try:
        if state.profile == WZ_SOURCE_TO_POLE:
            missing = [
                key for key in WZ_SHARED_HASH_FIELDS if key not in state.shared_contract_hashes
            ]
            if missing:
                raise ProductionEnvelopeError(f"wz_shared_hash_fields_missing:{missing}")
        verified_file_hashes = {
            spec.sha256 for path, spec in store.specs.items() if path in store.raw
        }
        unresolved = sorted(
            key
            for key, digest in state.shared_contract_hashes.items()
            if digest not in verified_file_hashes
        )
        if unresolved:
            raise ProductionEnvelopeError(
                f"shared_contract_hashes_do_not_resolve_to_bundle_bytes:{unresolved}"
            )
        state.checks["P0_WZ_SHARED_HASH_FAMILY_RECEIPT"] = True
    except ProductionEnvelopeError as exc:
        state.blockers.append(f"envelope:{state.artifact_id}:wz_hashes:{exc}")
    return state


def _verify_relationships(
    states: Sequence[_EnvelopeState], store: _Store
) -> tuple[list[str], list[str], dict[str, bool]]:
    blockers: list[str] = []
    checks = {
        "P0_PARENT_ENVELOPE_RESOLUTION_RECEIPT": False,
        "P0_ANCESTRY_DAG_RECEIPT": False,
        "P0_SINGLE_SOURCE_BRANCH_FREEZE_FAMILY_RECEIPT": False,
        "P0_WZ_SHARED_HASH_FAMILY_RECEIPT": False,
    }
    by_id: dict[str, _EnvelopeState] = {}
    duplicate_ids: set[str] = set()
    for state in states:
        if state.artifact_id is None:
            continue
        if state.artifact_id in by_id:
            duplicate_ids.add(state.artifact_id)
        by_id[state.artifact_id] = state
    for artifact_id in sorted(duplicate_ids):
        blockers.append(f"duplicate_envelope_artifact_id:{artifact_id}")

    graph: dict[str, set[str]] = {artifact_id: set() for artifact_id in by_id}
    parents_valid = not duplicate_ids and len(by_id) == len(states)
    for state in states:
        state_parent_valid = state.artifact_id is not None
        for index, raw_parent in enumerate(state.parent_rows):
            try:
                _exact_fields(raw_parent, _PARENT_FIELDS, "parent_receipt")
                parent_id = _identifier(raw_parent.get("artifact_id"), "parent.artifact_id")
                parent_subject = _hash(raw_parent.get("subject_digest"), "parent.subject_digest")
                parent_output = _hash(raw_parent.get("output_digest"), "parent.output_digest")
                parent = by_id.get(parent_id)
                if parent is None:
                    raise ProductionEnvelopeError("parent_envelope_not_in_manifest_envelope_paths")
                if parent_id == state.artifact_id:
                    raise ProductionEnvelopeError("envelope_cannot_parent_itself")
                # Build the declared identity graph before checking its byte
                # bindings.  A cycle must remain visible even when changing an
                # envelope also makes a content-hash reference stale.
                assert state.artifact_id is not None
                graph[state.artifact_id].add(parent_id)
                ref_spec, _ = store.resolve_ref(
                    raw_parent.get("envelope_ref"),
                    f"{state.artifact_id}.parent_receipts[{index}].envelope_ref",
                )
                if ref_spec.path != parent.path or ref_spec.sha256 != parent.file_sha256:
                    raise ProductionEnvelopeError("parent_envelope_file_binding_mismatch")
                if parent_subject != parent.subject_digest:
                    raise ProductionEnvelopeError("parent_subject_digest_mismatch")
                if parent_output != parent.output_digest:
                    raise ProductionEnvelopeError("parent_output_digest_mismatch")
                if (
                    parent.source_root_hash != state.source_root_hash
                    or parent.branch_id != state.branch_id
                    or parent.freeze_id != state.freeze_id
                ):
                    raise ProductionEnvelopeError("parent_mixed_source_branch_or_freeze_family")
                state.parent_artifact_ids.append(parent_id)
            except (ProductionEnvelopeError, TypeError) as exc:
                state_parent_valid = False
                parents_valid = False
                blocker = f"envelope:{state.artifact_id}:parent:{exc}"
                state.blockers.append(blocker)
                blockers.append(blocker)
        state.checks["P0_PARENT_ENVELOPE_RESOLUTION_RECEIPT"] = state_parent_valid

    acyclic, order = _topological_parent_first(graph)
    if not acyclic:
        blockers.append("parent_envelope_ancestry_cycle")
    checks["P0_PARENT_ENVELOPE_RESOLUTION_RECEIPT"] = parents_valid
    checks["P0_ANCESTRY_DAG_RECEIPT"] = acyclic and not duplicate_ids and len(by_id) == len(states)
    for state in states:
        state.checks["P0_ANCESTRY_DAG_RECEIPT"] = checks["P0_ANCESTRY_DAG_RECEIPT"]

    roots = {state.source_root_hash for state in states if state.source_root_hash is not None}
    branches = {state.branch_id for state in states if state.branch_id is not None}
    freezes = {state.freeze_id for state in states if state.freeze_id is not None}
    freeze_anchors = {
        state.freeze_anchor_sha256
        for state in states
        if state.freeze_anchor_sha256 is not None
    }
    one_family = (
        bool(states)
        and len(by_id) == len(states)
        and len(roots) == 1
        and len(branches) == 1
        and len(freezes) == 1
        and len(freeze_anchors) == 1
    )
    if not one_family:
        blockers.append("bundle_must_have_one_source_root_branch_freeze_family")
    checks["P0_SINGLE_SOURCE_BRANCH_FREEZE_FAMILY_RECEIPT"] = one_family
    for state in states:
        state.checks["P0_SINGLE_SOURCE_BRANCH_FREEZE_FAMILY_RECEIPT"] = one_family

    hash_values: dict[str, set[str]] = {}
    for state in states:
        for key, value in state.shared_contract_hashes.items():
            hash_values.setdefault(key, set()).add(value)
    shared = all(len(values) == 1 for values in hash_values.values())
    if not shared:
        blockers.append("shared_contract_hash_family_mismatch")
    wz_rows_valid = all(
        state.checks.get("P0_WZ_SHARED_HASH_FAMILY_RECEIPT") is True
        for state in states
        if state.profile == WZ_SOURCE_TO_POLE
    )
    checks["P0_WZ_SHARED_HASH_FAMILY_RECEIPT"] = shared and wz_rows_valid
    for state in states:
        state.checks["P0_WZ_SHARED_HASH_FAMILY_RECEIPT"] = shared and (
            state.profile != WZ_SOURCE_TO_POLE
            or state.checks.get("P0_WZ_SHARED_HASH_FAMILY_RECEIPT") is True
        )
    return blockers, order, checks


def _topological_parent_first(graph: Mapping[str, set[str]]) -> tuple[bool, list[str]]:
    indegree = {node: len(parents) for node, parents in graph.items()}
    children: dict[str, set[str]] = {node: set() for node in graph}
    for child, parents in graph.items():
        for parent in parents:
            children.setdefault(parent, set()).add(child)
    ready = sorted(node for node, degree in indegree.items() if degree == 0)
    order: list[str] = []
    while ready:
        node = ready.pop(0)
        order.append(node)
        for child in sorted(children.get(node, ())):
            indegree[child] -= 1
            if indegree[child] == 0:
                ready.append(child)
                ready.sort()
    return len(order) == len(graph), order if len(order) == len(graph) else sorted(graph)


def _parse_file_specs(value: Any, manifest_name: str) -> dict[str, _FileSpec]:
    if not isinstance(value, list) or not value or len(value) > MAX_FILES:
        raise ProductionEnvelopeError("files_must_be_nonempty_bounded_array")
    specs: dict[str, _FileSpec] = {}
    for index, row_value in enumerate(value):
        row = _mapping(row_value, f"files[{index}]")
        _exact_fields(row, _FILE_FIELDS, f"files[{index}]")
        spec = _FileSpec(
            path=_relative_path(row.get("path"), f"files[{index}].path"),
            sha256=_hash(row.get("sha256"), f"files[{index}].sha256"),
            byte_count=_byte_count(row.get("byte_count"), f"files[{index}].byte_count"),
            media_type=_media_type(row.get("media_type"), f"files[{index}].media_type"),
        )
        if spec.path == manifest_name:
            raise ProductionEnvelopeError("manifest_must_not_list_itself")
        if spec.path in specs:
            raise ProductionEnvelopeError(f"duplicate_manifest_file_path:{spec.path}")
        specs[spec.path] = spec
    return specs


def _parse_envelope_paths(value: Any, specs: Mapping[str, _FileSpec]) -> list[str]:
    if not isinstance(value, list) or not value or len(value) > MAX_ENVELOPES:
        raise ProductionEnvelopeError("envelope_paths_must_be_nonempty_bounded_array")
    paths = [_relative_path(item, "envelope_path") for item in value]
    if len(set(paths)) != len(paths):
        raise ProductionEnvelopeError("duplicate_envelope_path")
    for path in paths:
        if path not in specs:
            raise ProductionEnvelopeError(f"envelope_path_not_listed:{path}")
    return paths


def _validate_numeric_policy(value: Any) -> None:
    policy = _mapping(value, "numeric_precision_and_rounding")
    _exact_fields(policy, _NUMERIC_FIELDS, "numeric_precision_and_rounding")
    arithmetic = _enum(
        policy.get("arithmetic"), {"exact", "floating", "interval", "mixed"}, "arithmetic"
    )
    precision = policy.get("precision_bits")
    if arithmetic == "exact":
        if precision is not None:
            raise ProductionEnvelopeError("exact_arithmetic_precision_bits_must_be_null")
    elif isinstance(precision, bool) or not isinstance(precision, int) or not 1 <= precision <= 1_000_000:
        raise ProductionEnvelopeError("precision_bits_must_be_positive_integer")
    _identifier(policy.get("rounding_mode"), "rounding_mode")
    _identifier(policy.get("numeric_backend"), "numeric_backend")


def _envelope_report(state: _EnvelopeState, inventory_passed: bool) -> dict[str, Any]:
    receipts = dict(state.checks)
    receipts["P0_IMMUTABLE_INVENTORY_REPLAY_RECEIPT"] = inventory_passed
    receipts["SCIENTIFIC_REPLAY_RECEIPT"] = False
    receipts["PROMOTION_RECEIPT"] = False
    blockers = sorted(set(state.blockers))
    return {
        "artifact_id": state.artifact_id,
        "stage_id": state.stage_id,
        "receipt_type": state.receipt_type,
        "profile": state.profile,
        "claim_lane": state.claim_lane,
        "claim_scope": state.claim_scope,
        "source_root_hash": state.source_root_hash,
        "branch_id": state.branch_id,
        "freeze_id": state.freeze_id,
        "freeze_anchor_sha256": state.freeze_anchor_sha256,
        "frozen_utc": state.frozen_utc,
        "generated_utc": state.generated_utc,
        "producer_status": state.producer_status,
        "producer_status_trusted": False,
        "subject_digest": state.subject_digest,
        "output_digest": state.output_digest,
        "envelope_path": state.path,
        "envelope_sha256": state.file_sha256,
        "parent_artifact_ids": sorted(set(state.parent_artifact_ids)),
        "shared_contract_hashes": dict(sorted(state.shared_contract_hashes.items())),
        "receipts": receipts,
        "inventory_replay_passed": inventory_passed,
        "scientific_replay_passed": False,
        "promotion_allowed": False,
        "evidence_class": "IMMUTABLE_INVENTORY_ONLY",
        "blockers": blockers,
    }


def _failure_report(blockers: Sequence[str], *, manifest_path: str | None) -> dict[str, Any]:
    receipts = {key: False for key in P0_RECEIPT_KEYS}
    inventory = sorted(set(blockers))
    open_requirements = ["no_registered_independent_scientific_checker"]
    return {
        "schema": PRODUCTION_BUNDLE_REPORT_SCHEMA,
        "artifact_type": PRODUCTION_BUNDLE_REPORT_ARTIFACT_TYPE,
        "manifest_path": manifest_path,
        "manifest_sha256": None,
        "manifest_payload_sha256": None,
        "bundle_id": None,
        "canonical_json_policy": CANONICAL_JSON_POLICY,
        "envelope_order": [],
        "envelopes": {},
        "file_rows": {},
        "receipts": receipts,
        "inventory_replay_passed": False,
        "scientific_replay_passed": False,
        "promotion_allowed": False,
        "passed": False,
        "status": "FAIL",
        "inventory_blockers": inventory,
        "open_requirements": open_requirements,
        "blockers": sorted(set(inventory + open_requirements)),
        "registered_scientific_checker_count": 0,
        "ignored_producer_statuses": {},
        "claim_boundary": (
            "The on-disk production bundle did not pass immutable inventory replay. "
            "Scientific replay and promotion are false."
        ),
    }


def _strict_json_loads(raw: bytes) -> Any:
    def reject_constant(value: str) -> None:
        raise ProductionEnvelopeError(f"nonfinite_json_number:{value}")

    def unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                raise ProductionEnvelopeError(f"duplicate_json_key:{key}")
            result[key] = value
        return result

    value = json.loads(
        raw.decode("utf-8"),
        object_pairs_hook=unique_object,
        parse_constant=reject_constant,
    )
    _validate_finite(value)
    return value


def _validate_finite(value: Any) -> None:
    if isinstance(value, float) and not math.isfinite(value):
        raise ProductionEnvelopeError("nonfinite_json_number")
    if isinstance(value, Mapping):
        for key, child in value.items():
            if not isinstance(key, str):
                raise ProductionEnvelopeError("json_object_key_not_string")
            _validate_finite(child)
    elif isinstance(value, (list, tuple)):
        for child in value:
            _validate_finite(child)


def _regular_nonsymlink_file(path: Path) -> Path:
    if path.is_symlink():
        raise ProductionEnvelopeError("manifest_path_is_symlink")
    resolved = path.resolve(strict=True)
    if not stat.S_ISREG(resolved.stat().st_mode):
        raise ProductionEnvelopeError("manifest_path_is_not_regular_file")
    return resolved


def _contained_regular_file(root: Path, relative: str) -> Path:
    normalized = _relative_path(relative, "artifact_path")
    current = root
    for part in Path(normalized).parts:
        current = current / part
        if current.is_symlink():
            raise ProductionEnvelopeError("artifact_path_contains_symlink")
    resolved = current.resolve(strict=True)
    try:
        resolved.relative_to(root)
    except ValueError as exc:
        raise ProductionEnvelopeError("artifact_path_escapes_bundle") from exc
    if not stat.S_ISREG(resolved.stat().st_mode):
        raise ProductionEnvelopeError("artifact_path_is_not_regular_file")
    return resolved


def _mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ProductionEnvelopeError(f"{label}_must_be_object")
    return value


def _exact_fields(value: Mapping[str, Any], expected: frozenset[str], label: str) -> None:
    actual = set(value)
    missing = sorted(expected - actual)
    unexpected = sorted(actual - expected)
    if missing:
        raise ProductionEnvelopeError(f"{label}_missing_fields:{missing}")
    if unexpected:
        raise ProductionEnvelopeError(f"{label}_unexpected_fields:{unexpected}")


def _identifier(value: Any, label: str) -> str:
    if not isinstance(value, str) or _IDENTIFIER_RE.fullmatch(value) is None:
        raise ProductionEnvelopeError(f"{label}_must_be_bounded_identifier")
    return value


def _hash(value: Any, label: str) -> str:
    if not isinstance(value, str) or _SHA256_RE.fullmatch(value) is None:
        raise ProductionEnvelopeError(f"{label}_must_be_sha256")
    return value


def _relative_path(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value or len(value) > 512 or "\\" in value:
        raise ProductionEnvelopeError(f"{label}_must_be_bounded_relative_posix_path")
    path = Path(value)
    if path.is_absolute() or ".." in path.parts or "." in path.parts or path.as_posix() != value:
        raise ProductionEnvelopeError(f"{label}_must_be_normalized_relative_posix_path")
    return value


def _byte_count(value: Any, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not 0 <= value <= MAX_ARTIFACT_BYTES:
        raise ProductionEnvelopeError(f"{label}_must_be_bounded_nonnegative_integer")
    return value


def _media_type(value: Any, label: str) -> str:
    if not isinstance(value, str) or _MEDIA_TYPE_RE.fullmatch(value) is None:
        raise ProductionEnvelopeError(f"{label}_must_be_media_type")
    return value


def _enum(value: Any, allowed: set[str] | frozenset[str], label: str) -> str:
    if not isinstance(value, str) or value not in allowed:
        raise ProductionEnvelopeError(f"{label}_invalid")
    return value


def _utc(value: Any, label: str) -> str:
    if not isinstance(value, str) or _UTC_RE.fullmatch(value) is None:
        raise ProductionEnvelopeError(f"{label}_must_be_second_precision_utc")
    try:
        _utc_instant(value)
    except ValueError as exc:
        raise ProductionEnvelopeError(f"{label}_must_be_valid_utc") from exc
    return value


def _utc_instant(value: str) -> datetime:
    return datetime.strptime(value, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=timezone.utc)


def _string_list(value: Any, label: str) -> list[str]:
    if not isinstance(value, list) or any(
        not isinstance(item, str) or not item or len(item) > 512 for item in value
    ):
        raise ProductionEnvelopeError(f"{label}_must_be_bounded_string_array")
    if len(set(value)) != len(value):
        raise ProductionEnvelopeError(f"{label}_contains_duplicates")
    return value


def _sha256_bytes(raw: bytes) -> str:
    return "sha256:" + hashlib.sha256(raw).hexdigest()


__all__ = [
    "ALLOWED_STATUSES",
    "CANONICAL_JSON_POLICY",
    "COMMON_STAGE",
    "FREEZE_ANCHOR_SCHEMA",
    "P0_RECEIPT_KEYS",
    "PRODUCTION_BUNDLE_MANIFEST_SCHEMA",
    "PRODUCTION_BUNDLE_REPORT_ARTIFACT_TYPE",
    "PRODUCTION_BUNDLE_REPORT_SCHEMA",
    "PRODUCTION_ENVELOPE_SCHEMA",
    "PROFILES",
    "REGISTERED_SCIENTIFIC_CHECKERS",
    "TARGET_FIREWALL_FALSE_FIELDS",
    "WZ_SHARED_HASH_FIELDS",
    "WZ_SOURCE_TO_POLE",
    "ProductionEnvelopeError",
    "canonical_json_bytes",
    "canonical_json_sha256",
    "verify_production_bundle_manifest",
]
