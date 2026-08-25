"""Executable replay verifier for INS-03 v2 primitive transcripts.

Semantics under verification:
``oph.ins03.primitive_generation_semantics.v2_proposed``, the PROPOSED and
unfrozen grammar of ``PRIMITIVE_GENERATION_SEMANTICS.md`` in this directory.
Transcript schema: ``oph.sim.ins03_primitive_transcript.v2_proposed``.

Given a transcript JSON file, this module re-executes every primitive with
exact arithmetic from the committed starting objects (the committed
eight-context effect table, count literals, and run-state diagonal of the v1
validator's reference constants), checks each step's declared intermediate
state digest, and checks the final replayed export byte-for-byte against the
transcript's claimed export under canonical serialization.  All checks are
fail-closed with named error codes; no floating-point value enters any
computation.

Verdict grammar, extending v1 exactly:

* ``TRANSCRIPT_REPLAY_VERIFIED_UNAUTHENTICATED``: full replay match on a
  ``synthetic``-class transcript.  The verdict states a derivation fact
  only.  It is explicitly NOT production, NOT provenance, and NOT custody:
  replay proves that the claimed export is computable from the committed
  starting objects and the transcribed arguments under the proposed
  grammar, and proves nothing about who or what produced the transcript.
* ``NONCONFORMANT``: any check fails; the report lists every named error.
* ``PRODUCER_AUTHENTICATION_UNIMPLEMENTED`` is preserved from v1 as the
  fail-closed error for every producer or custody claim: the reserved
  ``producer`` transcript class, a ``producer`` provenance class in the
  claimed export, an ``evidential`` label, or an ``authenticated_binding``
  block (whose named verification hook is unimplemented; see
  ``AUTHENTICATED_BINDING_SPEC.md``).
* The class marker ``SYNTHETIC`` is carried in the report for a
  ``synthetic``-class transcript; the shipped sample is so marked and can
  never present as production.

Exact arithmetic, encodings, canonical serialization, digests, and the
committed reference constants are reused by import from the v1 module
``code/phase_instrument_export/ins03_export_validator.py``; they are
single-sourced there with their Lean sources cited.

Self-digests are not custody: the transcript's ``transcript_digest_sha256``
and the export's ``custody_digest_sha256`` are integrity conveniences that
whoever edits the content recomputes.  They bind the content to no agent,
key, time, or execution.  The audit ruling of
``plan/audits/RER_POST_R2020_ACTUAL_COMMITS_DEEP_AUDIT_2026-08-24.md``
(finding F3) governs: production, provenance, and custody claims require
authenticated generation semantics, transcripts, and independent replay,
and until those exist every such claim fails closed here.

What is not proved here: no producer exists, no owner key exists, no run
exists, and no freeze event exists; the semantics document is a proposal
and accepts only the ``PROPOSED_UNFROZEN`` status string.  A verified
replay discharges no register row: PR-03, PR-64, and PR-65 are open, PR-04
stays a declared row, and OL-C5 stays ``partial``.  Nothing here arms an
instrument, draws a seed, freezes a rule, scores a comparison, or supports
any physical claim of the OPH program.

Pure standard-library Python plus the v1 module.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# v1 import: exact arithmetic, encodings, committed reference constants
# ---------------------------------------------------------------------------

_V1_PATH = Path(__file__).resolve().parent.parent / "ins03_export_validator.py"


def _load_v1():
    existing = sys.modules.get("ins03_export_validator")
    if existing is not None and getattr(existing, "__file__", None):
        if Path(existing.__file__).resolve() == _V1_PATH:
            return existing
    spec = importlib.util.spec_from_file_location(
        "ins03_export_validator", _V1_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules["ins03_export_validator"] = module
    spec.loader.exec_module(module)
    return module


v1 = _load_v1()

# Re-exported committed reference constants (single-sourced in v1, where the
# Lean sources are cited; the tests check identity and the literals).
COMMITTED_CONTEXTS = v1.COMMITTED_CONTEXTS
COMMITTED_EFFECT0 = v1.COMMITTED_EFFECT0
COMMITTED_COUNTS = v1.COMMITTED_COUNTS
COMMITTED_RUN_STATE_DIAGONAL = v1.COMMITTED_RUN_STATE_DIAGONAL
COMMITTED_OFFDIAG_NORMSQ_BOUND = v1.COMMITTED_OFFDIAG_NORMSQ_BOUND

VERIFIER_ID = "ins03_transcript_replay_verifier.v2_proposed"
TRANSCRIPT_SCHEMA_ID = "oph.sim.ins03_primitive_transcript.v2_proposed"
SEMANTICS_ID = "oph.ins03.primitive_generation_semantics.v2_proposed"
SEMANTICS_STATUS_PROPOSED = "PROPOSED_UNFROZEN"
STATE_SCHEMA_ID = "oph.ins03.replay_state.v2_proposed"

REPLAY_VERDICT = "TRANSCRIPT_REPLAY_VERIFIED_UNAUTHENTICATED"
NONCONFORMANT_VERDICT = "NONCONFORMANT"
SYNTHETIC_CLASS_MARKER = "SYNTHETIC"

TRANSCRIPT_CLASSES = ("synthetic", "producer")

DECLARED_FIELDS = ("carrier_id", "source_record_ids", "operations",
                   "provenance", "labels", "derived_for_display")
EMIT_REQUIRED_DECLARATIONS = ("carrier_id", "source_record_ids",
                              "operations", "provenance", "labels")

PRIMITIVES = ("load_committed_effect_table", "prepare_state",
              "assemble_lueders_kraus", "readout_table",
              "bind_declaration", "emit_export")

TRANSCRIPT_REQUIRED = ("schema", "semantics_id", "semantics_status",
                       "transcript_class", "steps", "claimed_export",
                       "claimed_export_digest_sha256",
                       "transcript_digest_sha256")
TRANSCRIPT_OPTIONAL = ("authenticated_binding",)

STEP_REQUIRED = ("index", "primitive", "args", "state_digest_sha256")

BOUNDARY_STATEMENT = (
    "A verified replay proves derivation only: the claimed export is "
    "computable from the committed starting objects and the transcribed "
    "arguments under the PROPOSED, unfrozen grammar. It is not production, "
    "not provenance, and not custody; origin requires the authenticated "
    "binding of AUTHENTICATED_BINDING_SPEC.md, which is unimplemented. "
    "in_class = false in any operations list is not provenance. Register "
    "rows PR-03, PR-64, and PR-65 are open."
)

V1_STATIC_VERDICT_NOTE = (
    "Informational static-fixture reading of the claimed export under the "
    "v1 checker; not a replay gate, and not production, provenance, or "
    "custody. The v1 checker pins the committed fixture literals; v2 "
    "replay deliberately does not."
)


# ---------------------------------------------------------------------------
# Replay machine
# ---------------------------------------------------------------------------


class ReplayError(Exception):
    """A named, fail-closed replay failure at one step."""

    def __init__(self, code: str, where: str, detail: str) -> None:
        super().__init__(f"{code} at {where}: {detail}")
        self.code = code
        self.where = where
        self.detail = detail


class ReplayState:
    """The section D machine state of PRIMITIVE_GENERATION_SEMANTICS.md."""

    def __init__(self) -> None:
        self.effects: dict[str, dict[str, Any]] | None = None
        self.preparation: dict[str, Any] | None = None
        self.outcome_maps: dict[str, Any] = {}
        self.summed_channel: dict[str, Any] = {}
        self.readback: dict[str, Any] = {}
        self.declarations: dict[str, Any] = {}
        self.export: dict[str, Any] | None = None


def encode_state(state: ReplayState) -> dict[str, Any]:
    """The canonical JSON encoding of the machine state (section D)."""
    effects = None
    if state.effects is not None:
        effects = {
            context: {"0": v1.encode_matrix(pair["0"]),
                      "1": v1.encode_matrix(pair["1"])}
            for context, pair in state.effects.items()
        }
    preparation = None
    if state.preparation is not None:
        prep = state.preparation
        preparation = {
            "rho_00": [prep["rho_00"].numerator, prep["rho_00"].denominator],
            "rho_11": [prep["rho_11"].numerator, prep["rho_11"].denominator],
            "rho_01": prep["rho_01"].encode(),
        }
    return {
        "schema": STATE_SCHEMA_ID,
        "effects": effects,
        "preparation": preparation,
        "outcome_maps": state.outcome_maps,
        "summed_channel": state.summed_channel,
        "readback": state.readback,
        "declarations": state.declarations,
        "export": state.export,
    }


def state_digest(state: ReplayState) -> str:
    text = v1.canonical_text(encode_state(state))
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def compute_transcript_digest(transcript: dict[str, Any]) -> str:
    """SHA-256 over the canonical text of the transcript object with the
    ``transcript_digest_sha256`` field removed, no trailing line break.
    A self-digest: an integrity convenience, never custody."""
    body = {key: value for key, value in transcript.items()
            if key != "transcript_digest_sha256"}
    return hashlib.sha256(
        v1.canonical_text(body).encode("utf-8")).hexdigest()


def compute_claimed_export_digest(claimed_export: Any) -> str:
    return hashlib.sha256(
        v1.canonical_text(claimed_export).encode("utf-8")).hexdigest()


# ---------------------------------------------------------------------------
# The named authenticated-binding verification hook (UNIMPLEMENTED)
# ---------------------------------------------------------------------------


def verify_authenticated_binding(node: Any, where: str,
                                 errors: "v1.Errors") -> None:
    """Named verification hook of ``AUTHENTICATED_BINDING_SPEC.md``.

    UNIMPLEMENTED and failing closed: no owner key exists, no signature
    scheme is registered, and no verification is performed here.  Every
    transcript carrying an ``authenticated_binding`` block fails with
    ``PRODUCER_AUTHENTICATION_UNIMPLEMENTED`` until the owner implements
    and freezes the binding of that specification.  A self-digest cannot
    substitute for this hook: it binds content to no agent, key, time, or
    execution.
    """
    del node
    errors.add(
        "PRODUCER_AUTHENTICATION_UNIMPLEMENTED", where,
        "the authenticated-binding verification hook is unimplemented; no "
        "owner key exists and no signature can be verified, so the claim "
        "fails closed")


# ---------------------------------------------------------------------------
# Argument helpers
# ---------------------------------------------------------------------------


def _require_arg_keys(args: Any, keys: set[str], where: str) -> None:
    if not isinstance(args, dict) or set(args.keys()) != keys:
        expected = sorted(keys) if keys else "{}"
        raise ReplayError("ARGUMENT_OUTSIDE_GRAMMAR", where,
                          f"arguments must be exactly {expected}")


def _decode_or_raise(decoder, node: Any, where: str, what: str):
    collector = v1.Errors()
    value = decoder(node, where, collector)
    if value is None or collector.items:
        details = "; ".join(item["detail"] for item in collector.items)
        raise ReplayError("ARGUMENT_OUTSIDE_GRAMMAR", where,
                          f"{what}: {details or 'malformed encoding'}")
    return value


def _contains_float(node: Any) -> bool:
    if isinstance(node, float):
        return True
    if isinstance(node, dict):
        return any(_contains_float(value) for value in node.values())
    if isinstance(node, list):
        return any(_contains_float(value) for value in node)
    return False


def _check_no_export(state: ReplayState, where: str) -> None:
    if state.export is not None:
        raise ReplayError("STEP_PRECONDITION_UNMET", where,
                          "a step after emit_export is outside the grammar")


# ---------------------------------------------------------------------------
# Primitive executors (section E of the semantics document)
# ---------------------------------------------------------------------------


def _exec_load_committed_effect_table(state: ReplayState, args: Any,
                                      where: str) -> None:
    _require_arg_keys(args, set(), where)
    _check_no_export(state, where)
    if state.effects is not None:
        raise ReplayError("STEP_PRECONDITION_UNMET", where,
                          "the committed effect table is loaded once")
    identity = v1.mat_identity(2)
    state.effects = {
        context: {"0": COMMITTED_EFFECT0[context],
                  "1": v1.mat_sub(identity, COMMITTED_EFFECT0[context])}
        for context in COMMITTED_CONTEXTS
    }


def _exec_prepare_state(state: ReplayState, args: Any, where: str) -> None:
    _require_arg_keys(args, {"rho_00", "rho_11", "rho_01"}, where)
    _check_no_export(state, where)
    if state.preparation is not None:
        raise ReplayError("STEP_PRECONDITION_UNMET", where,
                          "the preparation is set once")
    rho_00 = _decode_or_raise(v1.decode_rational, args["rho_00"],
                              f"{where}.rho_00", "rho_00")
    rho_11 = _decode_or_raise(v1.decode_rational, args["rho_11"],
                              f"{where}.rho_11", "rho_11")
    rho_01 = _decode_or_raise(v1.decode_c3, args["rho_01"],
                              f"{where}.rho_01", "rho_01")
    if rho_00 + rho_11 != 1:
        raise ReplayError("ARGUMENT_OUTSIDE_GRAMMAR", where,
                          "diagonal entries must sum to exactly one")
    if rho_00 < 0 or rho_11 < 0:
        raise ReplayError("ARGUMENT_OUTSIDE_GRAMMAR", where,
                          "diagonal entries must be nonnegative")
    norm_sq = rho_01.norm_sq()
    product = v1.Q3(rho_00 * rho_11, Fraction(0))
    if not v1.q3_le(norm_sq, product):
        raise ReplayError("ARGUMENT_OUTSIDE_GRAMMAR", where,
                          "|rho_01|^2 <= rho_00 * rho_11 fails exactly")
    state.preparation = {"rho_00": rho_00, "rho_11": rho_11,
                         "rho_01": rho_01}


def _exec_assemble_lueders_kraus(state: ReplayState, args: Any,
                                 where: str) -> None:
    _require_arg_keys(args, {"context"}, where)
    _check_no_export(state, where)
    if state.effects is None:
        raise ReplayError("STEP_PRECONDITION_UNMET", where,
                          "the committed effect table is not loaded")
    context = args["context"]
    if not isinstance(context, str) or context not in state.effects:
        raise ReplayError("ARGUMENT_OUTSIDE_GRAMMAR", where,
                          f"unknown context {context!r}")
    if context in state.outcome_maps:
        raise ReplayError("STEP_PRECONDITION_UNMET", where,
                          f"context {context!r} is assembled once")
    effect0 = state.effects[context]["0"]
    effect1 = state.effects[context]["1"]
    if not v1.mat_eq(effect0, v1.mat_dagger(effect0)) \
            or not v1.mat_eq(v1.mat_mul(effect0, effect0), effect0):
        raise ReplayError(
            "ARGUMENT_OUTSIDE_GRAMMAR", where,
            "the Lueders singleton assembly is defined for projector "
            "effects only")
    dimension = 2
    families = {"0": [effect0], "1": [effect1]}
    context_maps: dict[str, Any] = {}
    for outcome, family in families.items():
        induced = v1.mat_zero(dimension)
        for matrix in family:
            induced = v1.mat_add(
                induced, v1.mat_mul(v1.mat_dagger(matrix), matrix))
        context_maps[outcome] = {
            "kraus": [v1.encode_matrix(matrix) for matrix in family],
            "effect_from_kraus": v1.encode_matrix(induced),
            "declared_effect": v1.encode_matrix(induced),
            "effect_residual": v1.encode_matrix(v1.mat_zero(dimension)),
            "trace_nonincreasing": v1._trace_entries([family], dimension),
        }
    state.outcome_maps[context] = context_maps

    both_families = [families["0"], families["1"]]
    normalization = v1.mat_zero(dimension)
    for family in both_families:
        for matrix in family:
            normalization = v1.mat_add(
                normalization, v1.mat_mul(v1.mat_dagger(matrix), matrix))
    state.summed_channel[context] = {
        "kraus_normalization": v1.encode_matrix(normalization),
        "kraus_normalization_residual": v1.encode_matrix(
            v1.mat_sub(normalization, v1.mat_identity(dimension))),
        "trace_checks": v1._trace_entries(both_families, dimension),
    }


def _exec_readout_table(state: ReplayState, args: Any, where: str) -> None:
    _require_arg_keys(args, {"context", "count_0", "count_1", "mass"}, where)
    _check_no_export(state, where)
    context = args["context"]
    if not isinstance(context, str):
        raise ReplayError("ARGUMENT_OUTSIDE_GRAMMAR", where,
                          "context must be a string")
    if state.preparation is None:
        raise ReplayError("STEP_PRECONDITION_UNMET", where,
                          "no preparation exists")
    if context not in state.outcome_maps or state.effects is None:
        raise ReplayError("STEP_PRECONDITION_UNMET", where,
                          f"context {context!r} is not assembled")
    if context in state.readback:
        raise ReplayError("STEP_PRECONDITION_UNMET", where,
                          f"context {context!r} is read out once")
    count_0, count_1, mass = args["count_0"], args["count_1"], args["mass"]
    for name, value in (("count_0", count_0), ("count_1", count_1),
                        ("mass", mass)):
        if not (isinstance(value, int) and not isinstance(value, bool)):
            raise ReplayError("ARGUMENT_OUTSIDE_GRAMMAR", where,
                              f"{name} must be an integer")
    if mass <= 0 or count_0 < 0 or count_1 < 0 \
            or count_0 > mass or count_1 > mass:
        raise ReplayError("ARGUMENT_OUTSIDE_GRAMMAR", where,
                          "counts must lie in [0, mass] with positive mass")
    if count_0 + count_1 != mass:
        raise ReplayError("ARGUMENT_OUTSIDE_GRAMMAR", where,
                          "outcome counts must sum to the context mass")
    prep = state.preparation
    rho: v1.Matrix = [
        [v1.c3_rational(prep["rho_00"]), prep["rho_01"]],
        [prep["rho_01"].conjugate(), v1.c3_rational(prep["rho_11"])],
    ]
    entries: dict[str, Any] = {}
    for outcome, count in (("0", count_0), ("1", count_1)):
        effect = state.effects[context][outcome]
        trace = v1.mat_trace(v1.mat_mul(rho, effect))
        if not trace.im.is_zero() or trace.re.b != 0:
            raise ReplayError(
                "ARGUMENT_OUTSIDE_GRAMMAR", where,
                "the primitive is defined only where Tr(rho E) is a plain "
                "rational")
        residual = trace.re.a - Fraction(count, mass)
        entries[outcome] = {
            "outcome_symbol": outcome,
            "count": count,
            "mass": mass,
            "compatibility_residual": [residual.numerator,
                                       residual.denominator],
        }
    state.readback[context] = entries


def _exec_bind_declaration(state: ReplayState, args: Any, where: str) -> None:
    _require_arg_keys(args, {"field", "value"}, where)
    _check_no_export(state, where)
    field = args["field"]
    if field not in DECLARED_FIELDS:
        raise ReplayError("ARGUMENT_OUTSIDE_GRAMMAR", where,
                          f"unknown declared field {field!r}")
    if field in state.declarations:
        raise ReplayError("STEP_PRECONDITION_UNMET", where,
                          f"declared field {field!r} is bound once")
    value = args["value"]
    if field == "carrier_id":
        if not isinstance(value, str) or not value:
            raise ReplayError("ARGUMENT_OUTSIDE_GRAMMAR", where,
                              "carrier_id must be a nonempty string")
    elif field == "source_record_ids":
        if not isinstance(value, list) or not all(
                isinstance(item, str) and item for item in value):
            raise ReplayError("ARGUMENT_OUTSIDE_GRAMMAR", where,
                              "source_record_ids must be a list of nonempty "
                              "strings")
    elif field == "operations":
        if not isinstance(value, list) or not all(
                isinstance(entry, dict)
                and isinstance(entry.get("operation"), str)
                and isinstance(entry.get("carrier_before"), str)
                and isinstance(entry.get("carrier_after"), str)
                and isinstance(entry.get("in_class"), bool)
                for entry in value):
            raise ReplayError(
                "ARGUMENT_OUTSIDE_GRAMMAR", where,
                "operations must be a list of objects with string "
                "operation, carrier_before, carrier_after and boolean "
                "in_class")
    elif field == "provenance":
        if not isinstance(value, dict) \
                or not isinstance(value.get("simulator_commit"), str) \
                or not v1.HEX40_RE.match(value["simulator_commit"]):
            raise ReplayError(
                "ARGUMENT_OUTSIDE_GRAMMAR", where,
                "provenance must be an object whose simulator_commit is 40 "
                "lowercase hex digits")
    elif field == "labels":
        if not isinstance(value, dict) \
                or not isinstance(value.get("exploratory"), bool) \
                or not isinstance(value.get("evidential"), bool) \
                or not isinstance(value.get("claim_boundary"), str) \
                or not value["claim_boundary"]:
            raise ReplayError(
                "ARGUMENT_OUTSIDE_GRAMMAR", where,
                "labels must carry boolean exploratory, boolean evidential, "
                "and a nonempty claim_boundary")
        if value["evidential"]:
            raise ReplayError(
                "PRODUCER_AUTHENTICATION_UNIMPLEMENTED", where,
                "an evidential label is a custody claim; no authenticator "
                "exists, so the claim fails closed")
    elif field == "derived_for_display":
        if not isinstance(value, dict):
            raise ReplayError("ARGUMENT_OUTSIDE_GRAMMAR", where,
                              "derived_for_display must be an object")
    state.declarations[field] = value


def _exec_emit_export(state: ReplayState, args: Any, where: str) -> None:
    _require_arg_keys(args, set(), where)
    if state.export is not None:
        raise ReplayError("STEP_PRECONDITION_UNMET", where,
                          "the export is emitted once")
    if state.effects is None:
        raise ReplayError("EMIT_PRECONDITION_UNMET", where,
                          "the committed effect table is not loaded")
    if state.preparation is None:
        raise ReplayError("EMIT_PRECONDITION_UNMET", where,
                          "no preparation exists")
    for context in state.effects:
        if context not in state.outcome_maps:
            raise ReplayError("EMIT_PRECONDITION_UNMET", where,
                              f"context {context!r} is not assembled")
        if context not in state.readback:
            raise ReplayError("EMIT_PRECONDITION_UNMET", where,
                              f"context {context!r} has no readback")
    for field in EMIT_REQUIRED_DECLARATIONS:
        if field not in state.declarations:
            raise ReplayError("EMIT_PRECONDITION_UNMET", where,
                              f"declared field {field!r} is not bound")

    prep = state.preparation
    rho_00, rho_11, rho_01 = prep["rho_00"], prep["rho_11"], prep["rho_01"]
    rho_00_enc = [rho_00.numerator, rho_00.denominator]
    rho_11_enc = [rho_11.numerator, rho_11.denominator]
    rho_01_enc = rho_01.encode()
    diagonal_product = rho_00 * rho_11
    carrier_id = state.declarations["carrier_id"]
    source_record_ids = state.declarations["source_record_ids"]
    simulator_commit = state.declarations["provenance"]["simulator_commit"]

    preparation = {
        "rho_00": rho_00_enc,
        "rho_11": rho_11_enc,
        "rho_01": rho_01_enc,
        "positivity_certificate": {
            "offdiag_norm_sq": rho_01.norm_sq().encode(),
            "diagonal_product": [diagonal_product.numerator,
                                 diagonal_product.denominator],
            "holds": True,
        },
        "record_diagonal": rho_01.is_zero(),
        "record_diagonal_offdiag": rho_01.encode(),
        "carrier_id": carrier_id,
        "source_record_ids": source_record_ids,
        "operations": state.declarations["operations"],
        "preparation_content_sha256": v1.compute_preparation_content_sha256(
            carrier_id, simulator_commit, source_record_ids,
            rho_00_enc, rho_01_enc, rho_11_enc),
    }

    export: dict[str, Any] = {
        "schema": v1.SCHEMA_ID,
        "provenance_class": "synthetic",
        "contexts": list(state.effects.keys()),
        "carrier_dimension": 2,
        "outcome_maps": state.outcome_maps,
        "summed_channel": state.summed_channel,
        "readback": state.readback,
        "preparation": preparation,
        "provenance": state.declarations["provenance"],
        "labels": state.declarations["labels"],
    }
    if "derived_for_display" in state.declarations:
        export["derived_for_display"] = \
            state.declarations["derived_for_display"]
    export["custody_digest_sha256"] = v1.compute_custody_digest(export)
    state.export = export


PRIMITIVE_EXECUTORS = {
    "load_committed_effect_table": _exec_load_committed_effect_table,
    "prepare_state": _exec_prepare_state,
    "assemble_lueders_kraus": _exec_assemble_lueders_kraus,
    "readout_table": _exec_readout_table,
    "bind_declaration": _exec_bind_declaration,
    "emit_export": _exec_emit_export,
}


def execute_primitive(state: ReplayState, primitive: str, args: Any,
                      where: str) -> None:
    """Execute one primitive on the state; raises ReplayError fail-closed."""
    executor = PRIMITIVE_EXECUTORS.get(primitive)
    if executor is None:
        raise ReplayError("UNKNOWN_PRIMITIVE", where,
                          f"primitive {primitive!r} is outside the grammar")
    scan_target: Any = args
    if primitive == "bind_declaration" and isinstance(args, dict) \
            and args.get("field") == "derived_for_display":
        scan_target = {key: value for key, value in args.items()
                       if key != "value"}
    if _contains_float(scan_target):
        raise ReplayError(
            "ARGUMENT_OUTSIDE_GRAMMAR", where,
            "floats are outside the grammar in every argument except the "
            "derived_for_display value")
    executor(state, args, where)


# ---------------------------------------------------------------------------
# The verifier
# ---------------------------------------------------------------------------


def verify_transcript_text(text: str,
                           transcript_path: str = "<text>") -> dict[str, Any]:
    """Verify one transcript document given as its file text.  Returns the
    report object; never raises on malformed content."""
    errors = v1.Errors()
    report: dict[str, Any] = {
        "verifier": VERIFIER_ID,
        "transcript_path": transcript_path,
        "verdict": NONCONFORMANT_VERDICT,
        "transcript_class": None,
        "class_marker": None,
        "semantics_status": None,
        "boundary": BOUNDARY_STATEMENT,
        "errors": errors.items,
    }

    def _reject_constant(name: str) -> None:
        raise ValueError(f"non-finite number {name!r}")

    try:
        parsed = json.loads(text, parse_constant=_reject_constant)
    except ValueError as error:
        errors.add("JSON_PARSE", "<document>", str(error))
        return report
    if not isinstance(parsed, dict):
        errors.add("JSON_PARSE", "<document>",
                   "top-level value must be an object")
        return report

    # Top-level field inventory.
    for field in TRANSCRIPT_REQUIRED:
        if field not in parsed:
            errors.add("MISSING_FIELD", "<top-level>",
                       f"required field {field!r} is absent")
    for field in parsed:
        if field not in TRANSCRIPT_REQUIRED + TRANSCRIPT_OPTIONAL:
            errors.add("UNEXPECTED_FIELD", "<top-level>",
                       f"field {field!r} is outside the transcript schema")

    # Canonical serialization of the file text.
    stripped = text[:-1] if text.endswith("\n") else text
    try:
        if stripped != v1.canonical_text(parsed):
            errors.add("NONCANONICAL_SERIALIZATION", "<document>",
                       "file text differs from json.dumps(obj, "
                       "sort_keys=True, indent=2, ensure_ascii=False)")
    except (TypeError, ValueError) as error:
        errors.add("NONCANONICAL_SERIALIZATION", "<document>", str(error))

    # Schema, semantics id, status, class.
    if parsed.get("schema") != TRANSCRIPT_SCHEMA_ID:
        errors.add("TRANSCRIPT_SCHEMA_MISMATCH", "schema",
                   f"expected {TRANSCRIPT_SCHEMA_ID!r}, found "
                   f"{parsed.get('schema')!r}")
    if parsed.get("semantics_id") != SEMANTICS_ID:
        errors.add("SEMANTICS_ID_MISMATCH", "semantics_id",
                   f"expected {SEMANTICS_ID!r}, found "
                   f"{parsed.get('semantics_id')!r}")
    semantics_status = parsed.get("semantics_status")
    report["semantics_status"] = semantics_status \
        if isinstance(semantics_status, str) else None
    if semantics_status != SEMANTICS_STATUS_PROPOSED:
        errors.add(
            "SEMANTICS_STATUS_INVALID", "semantics_status",
            f"only {SEMANTICS_STATUS_PROPOSED!r} is accepted; a claimed "
            "freeze is unverifiable because no freeze registration exists, "
            "so it fails closed")

    transcript_class = parsed.get("transcript_class")
    if transcript_class not in TRANSCRIPT_CLASSES:
        errors.add("TRANSCRIPT_CLASS_INVALID", "transcript_class",
                   f"expected one of {TRANSCRIPT_CLASSES}, found "
                   f"{transcript_class!r}")
        transcript_class = None
    report["transcript_class"] = transcript_class
    if transcript_class == "synthetic":
        report["class_marker"] = SYNTHETIC_CLASS_MARKER
    if transcript_class == "producer":
        errors.add(
            "PRODUCER_AUTHENTICATION_UNIMPLEMENTED", "transcript_class",
            "the producer transcript class claims production; no "
            "authenticator exists, so the claim fails closed")

    # The named authenticated-binding hook (unimplemented, fail-closed).
    if "authenticated_binding" in parsed:
        verify_authenticated_binding(parsed["authenticated_binding"],
                                     "authenticated_binding", errors)

    # Producer or custody claims inside the claimed export.
    claimed_export = parsed.get("claimed_export")
    if isinstance(claimed_export, dict):
        if claimed_export.get("provenance_class") == "producer":
            errors.add(
                "PRODUCER_AUTHENTICATION_UNIMPLEMENTED",
                "claimed_export.provenance_class",
                "a producer provenance class claims production; no "
                "authenticator exists, so the claim fails closed")
        labels = claimed_export.get("labels")
        if isinstance(labels, dict) and labels.get("evidential") is True:
            errors.add(
                "PRODUCER_AUTHENTICATION_UNIMPLEMENTED",
                "claimed_export.labels",
                "an evidential label is a custody claim; no authenticator "
                "exists, so the claim fails closed")

    # Transcript self-digest (integrity convenience, never custody).
    declared_transcript_digest = parsed.get("transcript_digest_sha256")
    if not isinstance(declared_transcript_digest, str) \
            or not v1.HEX64_RE.match(declared_transcript_digest):
        errors.add("TRANSCRIPT_DIGEST_FORMAT", "transcript_digest_sha256",
                   "expected 64 lowercase hex digits")
    else:
        try:
            recomputed = compute_transcript_digest(parsed)
        except (TypeError, ValueError):
            recomputed = None
        if recomputed is not None \
                and declared_transcript_digest != recomputed:
            errors.add("TRANSCRIPT_DIGEST_MISMATCH",
                       "transcript_digest_sha256",
                       "declared transcript digest differs from the "
                       "canonical recomputation")

    # Claimed-export digests.
    if isinstance(claimed_export, dict):
        embedded = claimed_export.get("custody_digest_sha256")
        if not isinstance(embedded, str) or not v1.HEX64_RE.match(embedded):
            errors.add("CLAIMED_EXPORT_DIGEST_FORMAT",
                       "claimed_export.custody_digest_sha256",
                       "expected 64 lowercase hex digits")
        else:
            try:
                recomputed = v1.compute_custody_digest(claimed_export)
            except (TypeError, ValueError):
                recomputed = None
            if recomputed is not None and embedded != recomputed:
                errors.add("CLAIMED_EXPORT_DIGEST_MISMATCH",
                           "claimed_export.custody_digest_sha256",
                           "embedded custody digest differs from the "
                           "canonical recomputation")
        declared_export_digest = parsed.get("claimed_export_digest_sha256")
        if not isinstance(declared_export_digest, str) \
                or not v1.HEX64_RE.match(declared_export_digest):
            errors.add("CLAIMED_EXPORT_DIGEST_FORMAT",
                       "claimed_export_digest_sha256",
                       "expected 64 lowercase hex digits")
        else:
            try:
                recomputed = compute_claimed_export_digest(claimed_export)
            except (TypeError, ValueError):
                recomputed = None
            if recomputed is not None \
                    and declared_export_digest != recomputed:
                errors.add("CLAIMED_EXPORT_DIGEST_MISMATCH",
                           "claimed_export_digest_sha256",
                           "declared claimed-export digest differs from "
                           "the canonical recomputation")
    elif "claimed_export" in parsed:
        errors.add("MISSING_FIELD", "claimed_export",
                   "claimed_export must be an object")
        claimed_export = None

    # Steps: structure, then replay.
    steps = parsed.get("steps")
    replay_completed = False
    state = ReplayState()
    if not isinstance(steps, list) or not steps:
        errors.add("STEP_FORMAT", "steps",
                   "steps must be a nonempty ordered list")
    else:
        replay_completed = True
        for position, step in enumerate(steps):
            where = f"steps[{position}]"
            if not isinstance(step, dict) \
                    or set(step.keys()) != set(STEP_REQUIRED):
                errors.add("STEP_FORMAT", where,
                           f"each step carries exactly {STEP_REQUIRED}")
                replay_completed = False
                break
            if step["index"] != position:
                errors.add("STEP_INDEX_NONCONSECUTIVE", where,
                           f"declared index {step['index']!r} differs from "
                           f"list position {position}")
                replay_completed = False
                break
            declared_digest = step["state_digest_sha256"]
            digest_ok = isinstance(declared_digest, str) \
                and bool(v1.HEX64_RE.match(declared_digest))
            if not digest_ok:
                errors.add("STEP_DIGEST_FORMAT", where,
                           "expected 64 lowercase hex digits")
            try:
                execute_primitive(state, step["primitive"], step["args"],
                                  where)
            except ReplayError as failure:
                errors.add(failure.code, failure.where, failure.detail)
                replay_completed = False
                break
            if digest_ok:
                recomputed = state_digest(state)
                if recomputed != declared_digest:
                    errors.add(
                        "STEP_DIGEST_MISMATCH", where,
                        f"declared intermediate digest at step {position} "
                        "differs from the exact replay recomputation")

    # Final export comparison.
    if replay_completed:
        if state.export is None:
            errors.add("EMIT_ABSENT", "steps",
                       "the replayed transcript emits no export")
        elif isinstance(claimed_export, dict):
            if v1.canonical_text(state.export) \
                    != v1.canonical_text(claimed_export):
                errors.add("FINAL_EXPORT_MISMATCH", "claimed_export",
                           "the replayed export differs from the claimed "
                           "export under canonical serialization")

    # Informational, non-gating v1 reading of the claimed export.
    if isinstance(claimed_export, dict):
        try:
            v1_report = v1.validate_export_text(
                v1.canonical_text(claimed_export))
            report["v1_static_verdict"] = v1_report["verdict"]
        except (TypeError, ValueError):
            report["v1_static_verdict"] = None
        report["v1_static_verdict_note"] = V1_STATIC_VERDICT_NOTE

    if not errors.items and transcript_class == "synthetic":
        report["verdict"] = REPLAY_VERDICT
    return report


def verify_file(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as error:
        return {
            "verifier": VERIFIER_ID,
            "transcript_path": str(path),
            "verdict": NONCONFORMANT_VERDICT,
            "transcript_class": None,
            "class_marker": None,
            "semantics_status": None,
            "boundary": BOUNDARY_STATEMENT,
            "errors": [{"code": "FILE_UNREADABLE", "where": str(path),
                        "detail": str(error)}],
        }
    return verify_transcript_text(text, str(path))


# ---------------------------------------------------------------------------
# Synthetic sample transcript builder
# ---------------------------------------------------------------------------
#
# The builder derives the committed static fixture of the v1 module
# (byte-identical to sample_conforming_export.json) through the grammar:
# load the committed effect table, prepare the committed run-state diagonal,
# assemble and read out every committed context, bind the marked synthetic
# declarations of the v1 sample verbatim, and emit.  The transcript is
# marked synthetic in its class, in the emitted provenance class, and in
# every placeholder value; it demonstrates the machinery and can never
# present as production.  The builder certifies nothing: it exists so that
# the sample is regenerable and so that the tests can mutate a
# known-verified document.


def build_sample_synthetic_transcript() -> dict[str, Any]:
    fixture = v1.build_synthetic_sample()
    diagonal = COMMITTED_RUN_STATE_DIAGONAL
    steps_spec: list[tuple[str, dict[str, Any]]] = [
        ("load_committed_effect_table", {}),
        ("prepare_state", {
            "rho_00": [diagonal[0].numerator, diagonal[0].denominator],
            "rho_11": [diagonal[1].numerator, diagonal[1].denominator],
            "rho_01": v1.C3_ZERO.encode(),
        }),
    ]
    for context in COMMITTED_CONTEXTS:
        steps_spec.append(("assemble_lueders_kraus", {"context": context}))
    for context in COMMITTED_CONTEXTS:
        count_0, count_1, mass = COMMITTED_COUNTS[context]
        steps_spec.append(("readout_table", {
            "context": context, "count_0": count_0, "count_1": count_1,
            "mass": mass,
        }))
    steps_spec.extend([
        ("bind_declaration", {"field": "carrier_id",
                              "value": fixture["preparation"]["carrier_id"]}),
        ("bind_declaration", {
            "field": "source_record_ids",
            "value": fixture["preparation"]["source_record_ids"]}),
        ("bind_declaration", {"field": "operations",
                              "value": fixture["preparation"]["operations"]}),
        ("bind_declaration", {"field": "provenance",
                              "value": fixture["provenance"]}),
        ("bind_declaration", {"field": "labels",
                              "value": fixture["labels"]}),
        ("bind_declaration", {"field": "derived_for_display",
                              "value": fixture["derived_for_display"]}),
        ("emit_export", {}),
    ])

    state = ReplayState()
    steps: list[dict[str, Any]] = []
    for position, (primitive, args) in enumerate(steps_spec):
        execute_primitive(state, primitive, args, f"steps[{position}]")
        steps.append({
            "index": position,
            "primitive": primitive,
            "args": args,
            "state_digest_sha256": state_digest(state),
        })
    assert state.export is not None
    transcript: dict[str, Any] = {
        "schema": TRANSCRIPT_SCHEMA_ID,
        "semantics_id": SEMANTICS_ID,
        "semantics_status": SEMANTICS_STATUS_PROPOSED,
        "transcript_class": "synthetic",
        "steps": steps,
        "claimed_export": state.export,
        "claimed_export_digest_sha256":
            compute_claimed_export_digest(state.export),
    }
    transcript["transcript_digest_sha256"] = \
        compute_transcript_digest(transcript)
    return transcript


def sample_transcript_file_text() -> str:
    """The canonical file text of the shipped sample transcript, with one
    trailing newline for the file form; digests are computed without it."""
    return v1.canonical_text(build_sample_synthetic_transcript()) + "\n"


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Fail-closed replay verifier for the "
            f"{TRANSCRIPT_SCHEMA_ID} transcript interface; a verified "
            "replay is a derivation fact only, and every producer or "
            "custody claim is unsupported."))
    parser.add_argument("transcript_path",
                        help="path to the transcript JSON file")
    args = parser.parse_args(argv)
    report = verify_file(args.transcript_path)
    print(json.dumps(report, indent=2))
    if report["errors"] and report["errors"][0]["code"] in ("FILE_UNREADABLE",
                                                            "JSON_PARSE"):
        return 2
    return 0 if report["verdict"] == REPLAY_VERDICT else 1


if __name__ == "__main__":
    sys.exit(main())
