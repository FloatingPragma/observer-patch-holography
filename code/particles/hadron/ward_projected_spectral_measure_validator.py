#!/usr/bin/env python3
"""Strict fail-closed validator for the Ward-projected spectral-measure export.

This module is the single typed validator/construction boundary for the
production export contract (`ward_projected_spectral_measure.schema.json`).
It positively enforces every semantic requirement of the contract; anything
unknown, missing, malformed, contradictory, or non-finite fails closed.
Schema conformance alone is never acceptance: the semantic layer runs on
top of the JSON Schema and each check appends a typed reason, and a payload
is accepted only when the reason list is empty.

Positively enforced requirements (ids in SEMANTIC_REQUIREMENTS):

- finite_numeric_values: every number anywhere in the payload is finite
  (NaN/Inf rejected); booleans are never accepted where numbers are typed.
- nonempty_level_support: at least one ensemble block with at least one level.
- unique_level_identifiers: level_id values are unique across the payload.
- residue_level_reference_integrity: residue rows and levels are in exact
  one-to-one correspondence (every residue references an existing level,
  no dangling or duplicate references, every level carries a residue).
- positive_energies: E_n > 0 for every level.
- s_equals_energy_squared: s_n = E_n^2 within the declared relative
  tolerance S_EQUALS_ENERGY_SQUARED_REL_TOL.
- nonnegative_weights_and_residues: w_n >= 0 and residues >= 0.
- weight_residue_consistency: the declared relation between level weights
  and Ward-projected residues (residue_n = weight_n within
  WEIGHT_RESIDUE_REL_TOL; both are the coefficient Z_n of
  delta(s - E_n^2) in the zero-momentum transverse channel).
- complete_typed_provenance: a closed provenance manifest with at least one
  typed source input (kind from the fixed allowlist) and the two
  no-leak booleans exactly false.
- no_measured_hvp_or_target_inputs: external_targets_used must be an empty
  list, and a recursive scan over every key and string value rejects any
  occurrence of a forbidden target token or a downstream-forbidden key
  (the downstream FORBIDDEN_SOURCE_KEYS set is imported, not copied).
- covariance_dimension_and_symmetry: square matrix matching the declared
  dimension and row basis, symmetric within COVARIANCE_SYMMETRY_REL_TOL.
- covariance_positive_semidefinite: min eigenvalue >= -COVARIANCE_PSD_REL_TOL
  relative to the matrix scale.
- covariance_rows_cover_all_levels: the covariance row basis must consist of
  exactly the rows {level_id}_energy and {level_id}_weight for every declared
  level, each exactly once (joint covariance completeness: a payload cannot
  add a level without extending the covariance).
- complete_budget_rows: exactly the seven required budget rows.
- finite_ordered_bound_intervals: every budget carries bound_interval with
  finite decimals 0 <= lo <= hi.
- transport_moment_certificate_complete: kernel, Delta_had_image interval,
  quadrature_error_bound and tail_bound, finite and correctly ordered
  (this is the block the downstream source-transport validator requires).
- strict_branch_typing: flavors exactly "2+1", qed from the closed enum;
  quenched, surrogate, target-calibrated and compare-only payloads are
  excluded from acceptance.
- format_version_allowlisted: the format version must be a known version
  (currently exactly 2); unknown versions fail closed.
- channel_allowlisted: every level block channel must come from the closed
  channel allowlist (currently exactly "U(1)_Q_vector").
- typed_normalization_convention: current_normalization must be a typed
  object with a convention from the closed allowlist and a certificate
  pointer to premise_certificates.conserved_current_or_matching; arbitrary
  normalization strings fail closed.
- typed_premise_certificate_references: the physical premises (gauge-quotient
  ensemble, conserved current / normalization matching, Ward identity,
  reflection positivity / transfer) must be typed references
  {artifact, path, sha256} with the correct const artifact type, a safe
  repo-relative path that resolves to an existing JSON file whose
  LF-normalized sha256 matches the pin, whose declared artifact type matches,
  and whose external_targets_used list is empty.
- reflection_positivity_certificate_required: positivity_status
  "certified_positive" is only accepted together with the required
  positivity_certificate pointer to
  premise_certificates.reflection_positivity_transfer and a resolvable
  reflection-positivity/transfer premise certificate.
- provenance_ensemble_cross_reference: the set of ensemble_ids used by the
  level blocks must equal the set of source_ensemble identifiers declared in
  provenance.source_inputs (no undeclared or unused ensembles).
- referenced_artifacts_source_only: every oph_source_artifact provenance
  identifier must be a safe repo-relative path inside the declared
  source-lane roots (ALLOWED_SOURCE_ARTIFACT_ROOTS) that resolves to an
  existing file; references outside the source lane (comparison artifacts,
  observational data, documentation) are rejected by location, and JSON
  files that declare guards.source_only false, an empirical row class, or
  compare-only guards are additionally rejected by content.
- source_artifact_positive_certification: a source-artifact reference is
  accepted only when the referenced file POSITIVELY certifies itself as a
  typed source-only artifact: it must be JSON, parse to an object, declare
  a nonempty typed "artifact" identifier, and explicitly declare
  "external_targets_used": []. Mere absence of negative markers is never
  sufficient; non-JSON files and unattested files fail closed.
- guard_flags_exactly_false: the three guard booleans are literally False
  (strings, 0, None, or missing fail).
- ward_projection_required: projection.ward_projected is literally True and
  the lane is U(1)_Q.
- closed_object_boundaries: unknown keys at any object boundary are
  rejected (additionalProperties: false everywhere in the schema).
- positivity_status_allowlisted: rho_had_or_measure.positivity_status must
  be exactly "certified_positive"; unknown status strings fail closed.

Certificate references are resolved against base_dir (default: the repo
root); hashes are computed over LF-normalized bytes so the pin is invariant
under platform line-ending conversion. The committed gate specimens under
code/particles/runs/hadron/gate_specimens/ declare
specimen_for_gate_testing: true; the gate accepts them (contract testing),
but the physical availability verdict in the verifier requires every premise
certificate to declare both specimen_for_gate_testing: false and
promotion_allowed: true.

The gate is downstream-compatible by construction: an accepted payload
embeds as `source_measure` into the source-transport payload contract of
code/P_derivation/thomson_spectral_transport.py and passes its
_validate_source_measure checks; verify_issue_317_spectral_measure_packet.py
executes that implication as an integration witness.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
import sys
from dataclasses import dataclass, field
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
P_DERIVATION = ROOT / "P_derivation"
if str(P_DERIVATION) not in sys.path:
    sys.path.insert(0, str(P_DERIVATION))

from thomson_spectral_transport import (  # noqa: E402
    FORBIDDEN_SOURCE_KEYS as DOWNSTREAM_FORBIDDEN_KEYS,
    SOURCE_SYSTEMATICS_BUDGETS,
)

SCHEMA_PATH = HERE / "ward_projected_spectral_measure.schema.json"
REPO_ROOT = ROOT.parent

FORBIDDEN_TARGETS = (
    "CODATA_ALPHA",
    "MUON_G_MINUS_2",
    "EE_TO_HADRONS",
    "RARE_DECAY_DATA",
    "HADRON_MASS_TARGETS",
    "PDG_QCD_FITS",
)

# Broader token set for the recursive scan: full forbidden targets plus the
# bare compilation names, so target provenance cannot hide behind shortened
# labels.
SCAN_TOKENS = FORBIDDEN_TARGETS + ("CODATA", "NIST", "PDG")

# Declared numerical tolerances of the contract.
S_EQUALS_ENERGY_SQUARED_REL_TOL = 1e-12
WEIGHT_RESIDUE_REL_TOL = 1e-12
COVARIANCE_SYMMETRY_REL_TOL = 1e-12
COVARIANCE_PSD_REL_TOL = 1e-12

ALLOWED_SOURCE_INPUT_KINDS = (
    "bare_lattice_parameter",
    "source_ensemble",
    "oph_source_artifact",
    "declared_convention",
    "backend_provenance",
)

ALLOWED_POSITIVITY_STATUS = ("certified_positive",)

ALLOWED_FLAVORS = ("2+1",)
ALLOWED_QED = ("off", "explicitly_superseded")
ALLOWED_FORMAT_VERSIONS = (2,)
ALLOWED_CHANNELS = ("U(1)_Q_vector",)
ALLOWED_NORMALIZATION_CONVENTIONS = (
    "conserved_current_ZV_equals_1",
    "conserved_local_matching_certificate",
)
NORMALIZATION_CERTIFICATE_POINTER = "premise_certificates.conserved_current_or_matching"
POSITIVITY_CERTIFICATE_POINTER = "premise_certificates.reflection_positivity_transfer"

# Typed premise-certificate references: the physical premises P1-P4 (and the
# P6 normalization matching) must be carried as hash-pinned references to
# certificate artifacts of the correct type, not as prose assertions.
PREMISE_CERTIFICATE_TYPES = {
    "gauge_quotient_ensemble": "oph_qcd_source_ensemble_certificate",
    "conserved_current_or_matching": "oph_ward_current_source_certificate",
    "ward_identity": "oph_ward_identity_certificate",
    "reflection_positivity_transfer": "oph_reflection_positivity_transfer_certificate",
}

# Committed typed specimens for exercising the gate's reference checks. They
# declare specimen_for_gate_testing: true, so the physical availability
# verdict rejects any payload backed by them.
SPECIMEN_CERTIFICATE_PATHS = {
    "gauge_quotient_ensemble": (
        "code/particles/runs/hadron/gate_specimens/"
        "specimen_gauge_quotient_ensemble_certificate.json"
    ),
    "conserved_current_or_matching": (
        "code/particles/runs/hadron/gate_specimens/"
        "specimen_ward_current_source_certificate.json"
    ),
    "ward_identity": (
        "code/particles/runs/hadron/gate_specimens/"
        "specimen_ward_identity_certificate.json"
    ),
    "reflection_positivity_transfer": (
        "code/particles/runs/hadron/gate_specimens/"
        "specimen_reflection_positivity_transfer_certificate.json"
    ),
}

_HEX64 = re.compile(r"^[0-9a-f]{64}$")

# Closed set of repo locations from which oph_source_artifact provenance
# references may be drawn: the production source-backend lane and the
# committed gate specimens. Everything else fails closed by location,
# regardless of content.
ALLOWED_SOURCE_ARTIFACT_ROOTS = (
    "code/particles/runs/qcd/hadron_source_backend/",
    "code/particles/runs/hadron/gate_specimens/",
)

SPECIMEN_GATE_ROOT = "code/particles/runs/hadron/gate_specimens/"
SPECIMEN_SOURCE_ARTIFACT_PATH = (
    "code/particles/runs/hadron/gate_specimens/specimen_source_artifact.json"
)

SEMANTIC_REQUIREMENTS = (
    "finite_numeric_values",
    "nonempty_level_support",
    "unique_level_identifiers",
    "residue_level_reference_integrity",
    "positive_energies",
    "s_equals_energy_squared",
    "nonnegative_weights_and_residues",
    "weight_residue_consistency",
    "complete_typed_provenance",
    "no_measured_hvp_or_target_inputs",
    "covariance_dimension_and_symmetry",
    "covariance_positive_semidefinite",
    "covariance_rows_cover_all_levels",
    "complete_budget_rows",
    "finite_ordered_bound_intervals",
    "transport_moment_certificate_complete",
    "strict_branch_typing",
    "format_version_allowlisted",
    "channel_allowlisted",
    "typed_normalization_convention",
    "typed_premise_certificate_references",
    "reflection_positivity_certificate_required",
    "provenance_ensemble_cross_reference",
    "referenced_artifacts_source_only",
    "source_artifact_positive_certification",
    "guard_flags_exactly_false",
    "ward_projection_required",
    "closed_object_boundaries",
    "positivity_status_allowlisted",
)


@dataclass(frozen=True)
class ValidationResult:
    accepted: bool
    reasons: tuple[str, ...] = field(default_factory=tuple)

    def to_json(self) -> dict[str, Any]:
        return {"accepted": self.accepted, "reasons": list(self.reasons)}


def load_schema() -> dict[str, Any]:
    return json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))


def lf_sha256(path: Path) -> str:
    """sha256 over LF-normalized bytes (invariant under CRLF checkouts)."""
    data = path.read_bytes().replace(b"\r\n", b"\n")
    return hashlib.sha256(data).hexdigest()


def specimen_certificate_reference(key: str, base_dir: Path | None = None) -> dict[str, Any]:
    """Typed reference to the committed gate specimen for a premise key."""
    base = base_dir if base_dir is not None else REPO_ROOT
    path_text = SPECIMEN_CERTIFICATE_PATHS[key]
    return {
        "artifact": PREMISE_CERTIFICATE_TYPES[key],
        "path": path_text,
        "sha256": lf_sha256(base / path_text),
    }


def _safe_relative_path(text: Any) -> bool:
    """Accept only repo-relative forward-slash paths without traversal."""
    if not isinstance(text, str) or not text:
        return False
    if "\\" in text or text.startswith("/") or ":" in text:
        return False
    return ".." not in text.split("/")


def _is_number(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def _finite_number(value: Any, where: str, reasons: list[str]) -> float | None:
    if not _is_number(value):
        reasons.append(f"not_a_number:{where}")
        return None
    if not math.isfinite(value):
        reasons.append(f"nonfinite_number:{where}")
        return None
    return float(value)


def _scan_finiteness(payload: Any, where: str, reasons: list[str]) -> None:
    """Reject NaN/Inf anywhere in the payload, at any nesting depth."""
    if isinstance(payload, dict):
        for key, value in payload.items():
            _scan_finiteness(value, f"{where}.{key}" if where else str(key), reasons)
    elif isinstance(payload, list):
        for index, value in enumerate(payload):
            _scan_finiteness(value, f"{where}[{index}]", reasons)
    elif _is_number(payload) and not math.isfinite(payload):
        reasons.append(f"nonfinite_number:{where}")


def _scan_forbidden_tokens(payload: Any, where: str, reasons: list[str]) -> None:
    """Recursive forbidden-token scan over every key and string value.

    Rejects (a) any key or string value containing a scan token (the full
    forbidden targets plus the bare compilation names, case-insensitive),
    and (b) any key equal to a downstream-forbidden source key. This closes
    the nested-provenance leak class, e.g.
    {"measured_hvp_input": {"source": "EE_TO_HADRONS"}}, independently of
    whether the schema boundary already rejects the unknown key.
    """
    if isinstance(payload, dict):
        for key, value in payload.items():
            key_text = str(key)
            key_where = f"{where}.{key_text}" if where else key_text
            upper = key_text.upper()
            for token in SCAN_TOKENS:
                if token in upper:
                    reasons.append(f"TARGET_LEAK_DETECTED:key:{key_where}")
            if key_text in DOWNSTREAM_FORBIDDEN_KEYS:
                reasons.append(f"downstream_forbidden_key:{key_where}")
            _scan_forbidden_tokens(value, key_where, reasons)
    elif isinstance(payload, list):
        for index, value in enumerate(payload):
            _scan_forbidden_tokens(value, f"{where}[{index}]", reasons)
    elif isinstance(payload, str):
        upper = payload.upper()
        for token in SCAN_TOKENS:
            if token in upper:
                reasons.append(f"TARGET_LEAK_DETECTED:value:{where}")


def _decimal_bound(value: Any, where: str, reasons: list[str]) -> Decimal | None:
    """Parse a nonnegative finite decimal from a number or string."""
    if isinstance(value, bool) or value is None:
        reasons.append(f"invalid_decimal:{where}")
        return None
    if _is_number(value) and not math.isfinite(value):
        reasons.append(f"nonfinite_number:{where}")
        return None
    try:
        parsed = Decimal(str(value))
    except (InvalidOperation, ValueError):
        reasons.append(f"invalid_decimal:{where}")
        return None
    if not parsed.is_finite():
        reasons.append(f"nonfinite_number:{where}")
        return None
    if parsed < 0:
        reasons.append(f"negative_decimal:{where}")
        return None
    return parsed


def _check_interval(value: Any, where: str, reasons: list[str]) -> None:
    if not isinstance(value, dict):
        reasons.append(f"required_interval_missing:{where}")
        return
    lo = _decimal_bound(value.get("lo"), f"{where}.lo", reasons)
    hi = _decimal_bound(value.get("hi"), f"{where}.hi", reasons)
    if lo is not None and hi is not None and lo > hi:
        reasons.append(f"invalid_interval_order:{where}")


def _check_schema(payload: dict[str, Any], schema: dict[str, Any], reasons: list[str]) -> None:
    import jsonschema

    validator = jsonschema.Draft202012Validator(schema)
    for error in sorted(validator.iter_errors(payload), key=lambda e: e.json_path):
        reasons.append(f"schema_violation:{error.json_path}:{error.message}")


def _check_provenance(payload: dict[str, Any], base_dir: Path, reasons: list[str]) -> None:
    targets = payload.get("external_targets_used")
    if not isinstance(targets, list):
        reasons.append("external_targets_used_not_a_list")
    elif targets:
        reasons.append(
            "TARGET_LEAK_DETECTED:external_targets_used:" + ",".join(str(t) for t in targets)
        )
    provenance = payload.get("provenance")
    if not isinstance(provenance, dict):
        reasons.append("provenance_manifest_missing")
        return
    source_inputs = provenance.get("source_inputs")
    if not isinstance(source_inputs, list) or not source_inputs:
        reasons.append("provenance_source_inputs_missing_or_empty")
    else:
        for index, item in enumerate(source_inputs):
            if not isinstance(item, dict):
                reasons.append(f"provenance_source_input_not_typed:[{index}]")
                continue
            kind = item.get("kind")
            if kind not in ALLOWED_SOURCE_INPUT_KINDS:
                reasons.append(
                    f"provenance_source_input_kind_not_allowlisted:[{index}]:{kind!r}"
                )
            identifier = item.get("identifier")
            if not isinstance(identifier, str) or not identifier:
                reasons.append(f"provenance_source_input_identifier_missing:[{index}]")
            elif kind == "oph_source_artifact":
                _check_source_artifact_reference(identifier, index, base_dir, reasons)
    if provenance.get("measured_hvp_input_present") is not False:
        reasons.append("provenance_measured_hvp_input_present_not_false")
    if provenance.get("target_calibration_present") is not False:
        reasons.append("provenance_target_calibration_present_not_false")


def _check_source_artifact_reference(
    identifier: str, index: int, base_dir: Path, reasons: list[str]
) -> None:
    """Resolve an oph_source_artifact provenance reference fail-closed.

    Acceptance is a POSITIVE certification, never the mere absence of
    negative markers: the identifier must be a safe repo-relative path
    inside the declared source-lane roots, resolve to an existing JSON
    file, parse to an object, declare a nonempty typed "artifact"
    identifier, and explicitly declare "external_targets_used": []. On top
    of that, negative markers (source_only false, empirical row class,
    compare-only guard) reject regardless. Non-JSON files, files outside
    the source lane, and unattested files all fail closed, so comparison
    or target artifacts cannot be disguised as OPH source artifacts.
    """
    where = f"provenance.source_inputs[{index}]"
    if not _safe_relative_path(identifier):
        reasons.append(f"source_artifact_reference_invalid:{where}:{identifier!r}")
        return
    if not identifier.startswith(ALLOWED_SOURCE_ARTIFACT_ROOTS):
        reasons.append(f"source_artifact_location_not_allowlisted:{where}:{identifier}")
        return
    file_path = base_dir / identifier
    if not file_path.is_file():
        reasons.append(f"source_artifact_reference_absent:{where}:{identifier}")
        return
    if file_path.suffix != ".json":
        reasons.append(f"source_artifact_reference_not_json:{where}:{identifier}")
        return
    try:
        content = json.loads(file_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        reasons.append(f"source_artifact_reference_unreadable:{where}:{identifier}")
        return
    if not isinstance(content, dict):
        reasons.append(f"source_artifact_not_positively_certified:{where}:{identifier}")
        return
    artifact = content.get("artifact")
    if not isinstance(artifact, str) or not artifact:
        reasons.append(f"source_artifact_reference_untyped:{where}:{identifier}")
    if content.get("external_targets_used") != []:
        reasons.append(
            f"source_artifact_external_targets_not_declared_empty:{where}:{identifier}"
        )
    guards = content.get("guards") if isinstance(content.get("guards"), dict) else {}
    row_class = content.get("row_class")
    if guards.get("source_only") is False or row_class == "oph_plus_empirical_hadron_closure":
        reasons.append(f"referenced_artifact_not_source_only:{where}:{identifier}")
    if guards.get("compare_only_external_endpoint") is True:
        reasons.append(f"referenced_artifact_compare_only:{where}:{identifier}")


def _check_provenance_cross_reference(payload: dict[str, Any], reasons: list[str]) -> None:
    """Ensemble ids used by level blocks must match declared source ensembles."""
    provenance = payload.get("provenance")
    source_inputs = provenance.get("source_inputs") if isinstance(provenance, dict) else None
    declared: set[str] = set()
    if isinstance(source_inputs, list):
        for item in source_inputs:
            if isinstance(item, dict) and item.get("kind") == "source_ensemble":
                identifier = item.get("identifier")
                if isinstance(identifier, str) and identifier:
                    declared.add(identifier)
    blocks = payload.get("finite_volume_levels")
    used: set[str] = set()
    if isinstance(blocks, list):
        for block in blocks:
            if isinstance(block, dict):
                ensemble_id = block.get("ensemble_id")
                if isinstance(ensemble_id, str) and ensemble_id:
                    used.add(ensemble_id)
    for missing in sorted(used - declared):
        reasons.append(f"ensemble_not_declared_in_provenance:{missing}")
    for unused in sorted(declared - used):
        reasons.append(f"declared_source_ensemble_unused:{unused}")


def _check_premise_certificates(
    payload: dict[str, Any], base_dir: Path, reasons: list[str]
) -> None:
    """Typed, hash-pinned premise-certificate references (P1-P4, P6)."""
    block = payload.get("premise_certificates")
    if not isinstance(block, dict):
        reasons.append("premise_certificates_missing")
        return
    for key in block:
        if key not in PREMISE_CERTIFICATE_TYPES:
            reasons.append(f"unknown_premise_certificate:{key}")
    for key, artifact in PREMISE_CERTIFICATE_TYPES.items():
        entry = block.get(key)
        where = f"premise_certificates.{key}"
        if not isinstance(entry, dict):
            reasons.append(f"premise_certificate_missing:{key}")
            continue
        if entry.get("artifact") != artifact:
            reasons.append(f"premise_certificate_artifact_mismatch:{key}")
        path_text = entry.get("path")
        if not _safe_relative_path(path_text):
            reasons.append(f"premise_certificate_path_invalid:{key}")
            continue
        sha = entry.get("sha256")
        if not isinstance(sha, str) or not _HEX64.match(sha):
            reasons.append(f"premise_certificate_sha256_invalid:{key}")
            continue
        file_path = base_dir / path_text
        if not file_path.is_file():
            reasons.append(f"premise_certificate_file_absent:{key}:{path_text}")
            continue
        if lf_sha256(file_path) != sha:
            reasons.append(f"premise_certificate_hash_mismatch:{key}")
            continue
        try:
            content = json.loads(file_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            reasons.append(f"premise_certificate_unreadable:{key}")
            continue
        if not isinstance(content, dict) or content.get("artifact") != artifact:
            reasons.append(f"premise_certificate_content_artifact_mismatch:{key}")
        if isinstance(content, dict) and content.get("external_targets_used") != []:
            reasons.append(f"premise_certificate_external_targets_not_empty:{key}")


def source_artifact_specimen_flags(
    payload: dict[str, Any], base_dir: Path | None = None
) -> dict[str, bool]:
    """Which oph_source_artifact provenance references are gate specimens.

    Fails closed: references under the gate-specimens root, unreadable
    references, and references whose content declares
    specimen_for_gate_testing true all count as specimen (True).
    """
    base = base_dir if base_dir is not None else REPO_ROOT
    provenance = payload.get("provenance") if isinstance(payload, dict) else None
    inputs = provenance.get("source_inputs") if isinstance(provenance, dict) else None
    flags: dict[str, bool] = {}
    if not isinstance(inputs, list):
        return flags
    for item in inputs:
        if not (isinstance(item, dict) and item.get("kind") == "oph_source_artifact"):
            continue
        identifier = item.get("identifier")
        if not _safe_relative_path(identifier):
            flags[str(identifier)] = True
            continue
        if identifier.startswith(SPECIMEN_GATE_ROOT):
            flags[identifier] = True
            continue
        try:
            content = json.loads((base / identifier).read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            flags[identifier] = True
            continue
        flags[identifier] = (
            not isinstance(content, dict)
            or content.get("specimen_for_gate_testing") is True
        )
    return flags


def premise_certificate_physical_readiness_flags(
    payload: dict[str, Any], base_dir: Path | None = None
) -> dict[str, bool]:
    """Which premise certificates are not ready for physical availability.

    Fails closed: unresolvable references, specimens, and certificates whose
    own promotion gate is not literally true all return True.  Merely changing
    ``specimen_for_gate_testing`` to false must not promote a specimen that
    still declares ``promotion_allowed: false`` or carries no physical
    certificate content.
    """
    base = base_dir if base_dir is not None else REPO_ROOT
    block = payload.get("premise_certificates")
    flags: dict[str, bool] = {}
    for key in PREMISE_CERTIFICATE_TYPES:
        entry = block.get(key) if isinstance(block, dict) else None
        path_text = entry.get("path") if isinstance(entry, dict) else None
        if not _safe_relative_path(path_text):
            flags[key] = True
            continue
        file_path = base / path_text
        try:
            content = json.loads(file_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            flags[key] = True
            continue
        flags[key] = (
            not isinstance(content, dict)
            or content.get("specimen_for_gate_testing") is not False
            or content.get("promotion_allowed") is not True
        )
    return flags


def _check_levels_and_residues(payload: dict[str, Any], reasons: list[str]) -> set[str]:
    blocks = payload.get("finite_volume_levels")
    if not isinstance(blocks, list) or not blocks:
        reasons.append("finite_volume_levels_missing_or_empty")
        blocks = []
    level_weights: dict[str, float] = {}
    seen_ids: set[str] = set()
    for b_index, block in enumerate(blocks):
        if not isinstance(block, dict):
            reasons.append(f"level_block_not_typed:[{b_index}]")
            continue
        channel = block.get("channel")
        if channel not in ALLOWED_CHANNELS:
            reasons.append(
                f"channel_not_allowlisted:finite_volume_levels[{b_index}]:{channel!r}"
            )
        levels = block.get("levels")
        if not isinstance(levels, list) or not levels:
            reasons.append(f"levels_missing_or_empty:[{b_index}]")
            continue
        for l_index, level in enumerate(levels):
            where = f"finite_volume_levels[{b_index}].levels[{l_index}]"
            if not isinstance(level, dict):
                reasons.append(f"level_not_typed:{where}")
                continue
            level_id = level.get("level_id")
            if not isinstance(level_id, str) or not level_id:
                reasons.append(f"level_id_missing:{where}")
                level_id = None
            elif level_id in seen_ids:
                reasons.append(f"duplicate_level_id:{level_id}")
            else:
                seen_ids.add(level_id)
            energy = _finite_number(level.get("energy"), f"{where}.energy", reasons)
            s_value = _finite_number(level.get("s"), f"{where}.s", reasons)
            weight = _finite_number(level.get("weight"), f"{where}.weight", reasons)
            if energy is not None and energy <= 0.0:
                reasons.append(f"nonpositive_energy:{where}")
            if weight is not None and weight < 0.0:
                reasons.append(f"negative_level_weight:{where}")
            if energy is not None and s_value is not None and energy > 0.0:
                expected = energy * energy
                if abs(s_value - expected) > S_EQUALS_ENERGY_SQUARED_REL_TOL * max(1.0, expected):
                    reasons.append(f"s_not_equal_energy_squared:{where}")
            for error_key in ("energy_jackknife_error", "weight_jackknife_error"):
                if error_key in level:
                    err = _finite_number(level[error_key], f"{where}.{error_key}", reasons)
                    if err is not None and err < 0.0:
                        reasons.append(f"negative_error:{where}.{error_key}")
            if level_id is not None and weight is not None:
                level_weights[level_id] = weight

    residues = payload.get("ward_projected_residues")
    if not isinstance(residues, list) or not residues:
        reasons.append("ward_projected_residues_missing_or_empty")
        residues = []
    residue_ids: set[str] = set()
    for r_index, row in enumerate(residues):
        where = f"ward_projected_residues[{r_index}]"
        if not isinstance(row, dict):
            reasons.append(f"residue_row_not_typed:{where}")
            continue
        level_id = row.get("level_id")
        if not isinstance(level_id, str) or not level_id:
            reasons.append(f"residue_level_id_missing:{where}")
            level_id = None
        elif level_id in residue_ids:
            reasons.append(f"duplicate_residue_reference:{level_id}")
        else:
            residue_ids.add(level_id)
        if level_id is not None and level_id not in seen_ids:
            reasons.append(f"residue_references_unknown_level:{level_id}")
        residue = _finite_number(row.get("residue"), f"{where}.residue", reasons)
        if residue is not None and residue < 0.0:
            reasons.append(f"negative_residue:{where}")
        normalization = row.get("current_normalization")
        if not isinstance(normalization, dict):
            reasons.append(f"normalization_not_typed:{where}")
        else:
            convention = normalization.get("convention")
            if convention not in ALLOWED_NORMALIZATION_CONVENTIONS:
                reasons.append(
                    f"normalization_convention_not_allowlisted:{where}:{convention!r}"
                )
            if normalization.get("certificate") != NORMALIZATION_CERTIFICATE_POINTER:
                reasons.append(f"normalization_certificate_pointer_missing:{where}")
        if (
            residue is not None
            and level_id is not None
            and level_id in level_weights
        ):
            weight = level_weights[level_id]
            if abs(residue - weight) > WEIGHT_RESIDUE_REL_TOL * max(1.0, abs(weight)):
                reasons.append(f"weight_residue_inconsistent:{level_id}")
    for level_id in sorted(seen_ids - residue_ids):
        reasons.append(f"level_without_residue:{level_id}")
    return seen_ids


def _check_measure_block(payload: dict[str, Any], reasons: list[str]) -> None:
    rho = payload.get("rho_had_or_measure")
    if not isinstance(rho, dict):
        reasons.append("rho_had_or_measure_missing")
        return
    if rho.get("support_variable") != "s":
        reasons.append("support_variable_not_s")
    if not rho.get("pushforward_rule"):
        reasons.append("pushforward_rule_missing")
    if rho.get("positivity_status") not in ALLOWED_POSITIVITY_STATUS:
        reasons.append(
            f"positivity_status_not_allowlisted:{rho.get('positivity_status')!r}"
        )
    if rho.get("positivity_certificate") != POSITIVITY_CERTIFICATE_POINTER:
        reasons.append("positivity_certificate_pointer_missing")


def _check_covariance(
    payload: dict[str, Any], level_ids: set[str], reasons: list[str]
) -> None:
    covariance = payload.get("covariance")
    if not isinstance(covariance, dict):
        reasons.append("covariance_missing")
        return
    dimension = covariance.get("dimension")
    matrix = covariance.get("matrix")
    row_basis = covariance.get("row_basis")
    if not isinstance(dimension, int) or isinstance(dimension, bool) or dimension < 1:
        reasons.append("covariance_dimension_invalid")
        return
    if not isinstance(row_basis, list) or len(row_basis) != dimension:
        reasons.append("covariance_row_basis_dimension_mismatch")
    if isinstance(row_basis, list) and level_ids:
        rows = [str(row) for row in row_basis]
        if len(set(rows)) != len(rows):
            reasons.append("covariance_row_basis_duplicates")
        expected = {f"{level_id}_energy" for level_id in level_ids} | {
            f"{level_id}_weight" for level_id in level_ids
        }
        if set(rows) != expected:
            reasons.append("covariance_rows_do_not_cover_levels")
    if not isinstance(matrix, list) or len(matrix) != dimension:
        reasons.append("covariance_matrix_dimension_mismatch")
        return
    rows: list[list[float]] = []
    for r_index, row in enumerate(matrix):
        if not isinstance(row, list) or len(row) != dimension:
            reasons.append(f"covariance_matrix_dimension_mismatch:row[{r_index}]")
            return
        parsed_row: list[float] = []
        for c_index, value in enumerate(row):
            entry = _finite_number(value, f"covariance.matrix[{r_index}][{c_index}]", reasons)
            if entry is None:
                return
            parsed_row.append(entry)
        rows.append(parsed_row)
    arr = np.array(rows, dtype=float)
    scale = float(np.max(np.abs(arr))) if arr.size else 0.0
    tol_scale = max(1.0, scale)
    if float(np.max(np.abs(arr - arr.T))) > COVARIANCE_SYMMETRY_REL_TOL * tol_scale:
        reasons.append("covariance_not_symmetric")
        return
    min_eig = float(np.linalg.eigvalsh(0.5 * (arr + arr.T)).min())
    if min_eig < -COVARIANCE_PSD_REL_TOL * tol_scale:
        reasons.append("covariance_not_positive_semidefinite")


def _check_transport_certificate(payload: dict[str, Any], reasons: list[str]) -> None:
    certificate = payload.get("transport_moment_certificate")
    if not isinstance(certificate, dict):
        reasons.append("transport_moment_certificate_missing")
        return
    kernel = certificate.get("kernel")
    if not isinstance(kernel, str) or not kernel:
        reasons.append("transport_moment_certificate_kernel_missing")
    _check_interval(
        certificate.get("Delta_had_image"),
        "transport_moment_certificate.Delta_had_image",
        reasons,
    )
    _decimal_bound(
        certificate.get("quadrature_error_bound"),
        "transport_moment_certificate.quadrature_error_bound",
        reasons,
    )
    _decimal_bound(
        certificate.get("tail_bound"),
        "transport_moment_certificate.tail_bound",
        reasons,
    )


def _check_budgets(payload: dict[str, Any], reasons: list[str]) -> None:
    systematics = payload.get("systematics")
    if not isinstance(systematics, dict):
        reasons.append("systematics_missing")
        return
    for name in SOURCE_SYSTEMATICS_BUDGETS:
        budget = systematics.get(name)
        if not isinstance(budget, dict):
            reasons.append(f"budget_missing:{name}")
            continue
        if "bound_interval" not in budget:
            reasons.append(f"budget_bound_interval_missing:{name}")
            continue
        _check_interval(budget["bound_interval"], f"systematics.{name}.bound_interval", reasons)
    for name in systematics:
        if name not in SOURCE_SYSTEMATICS_BUDGETS:
            reasons.append(f"unknown_budget_row:{name}")


def _check_typing_and_guards(payload: dict[str, Any], reasons: list[str]) -> None:
    version = payload.get("format_version")
    if isinstance(version, bool) or version not in ALLOWED_FORMAT_VERSIONS:
        reasons.append(f"format_version_not_allowlisted:{version!r}")
    branch = payload.get("branch")
    if not isinstance(branch, dict):
        reasons.append("branch_missing")
    else:
        if branch.get("flavors") not in ALLOWED_FLAVORS:
            reasons.append(f"branch_flavors_not_allowlisted:{branch.get('flavors')!r}")
        if branch.get("qed") not in ALLOWED_QED:
            reasons.append(f"branch_qed_not_allowlisted:{branch.get('qed')!r}")
    projection = payload.get("projection")
    if not isinstance(projection, dict):
        reasons.append("projection_missing")
    else:
        if projection.get("lane") != "U(1)_Q":
            reasons.append("projection_lane_not_u1q")
        if projection.get("ward_projected") is not True:
            reasons.append("ward_projection_dropped")
    guards = payload.get("guards")
    if not isinstance(guards, dict):
        reasons.append("guards_missing")
        return
    for name in ("stable_channel_only", "surrogate_hadron_artifact", "compare_only_external_endpoint"):
        if guards.get(name) is not False:
            reasons.append(f"guard_not_literal_false:{name}")


def validate_production_payload(
    payload: Any,
    schema: dict[str, Any] | None = None,
    base_dir: Path | None = None,
) -> ValidationResult:
    """Strict fail-closed acceptance decision for a production export payload.

    base_dir is the directory against which certificate and source-artifact
    references are resolved (default: the repository root).
    """
    if not isinstance(payload, dict):
        return ValidationResult(False, ("payload_not_an_object",))
    schema = schema or load_schema()
    base = base_dir if base_dir is not None else REPO_ROOT
    reasons: list[str] = []
    _check_schema(payload, schema, reasons)
    _scan_finiteness(payload, "", reasons)
    _scan_forbidden_tokens(payload, "", reasons)
    _check_provenance(payload, base, reasons)
    _check_provenance_cross_reference(payload, reasons)
    _check_premise_certificates(payload, base, reasons)
    level_ids = _check_levels_and_residues(payload, reasons)
    _check_measure_block(payload, reasons)
    _check_covariance(payload, level_ids, reasons)
    _check_transport_certificate(payload, reasons)
    _check_budgets(payload, reasons)
    _check_typing_and_guards(payload, reasons)
    deduped = tuple(dict.fromkeys(reasons))
    return ValidationResult(accepted=not deduped, reasons=deduped)
