"""Static committed-fixture checker for the INS-03 v1 export interface.

Schema under validation: ``oph.sim.ins03_phase_instrument_export.v1``, the
simulator-side export specified in section A of
``plan/INS03_SOURCE_BOUND_PHASE_INSTRUMENT_DESIGN.md`` (OPH meta planning
surface).  In the typed Lean binding interface
``Lean/EventAlgebra/SourceBoundInstrumentInterface.lean``, the current
placeholder fields are not determined by the committed objects and do not
authenticate custody.  A richer in-corpus or external data-bearing
construction remains viable.  This v1 checker receives or authenticates no
such data.  It checks a
static synthetic transcription of the committed fixture, fail-closed with
named error codes:

* schema id and version (the ``.v1`` suffix of the schema string);
* canonical serialization: ``json.dumps(obj, sort_keys=True, indent=2,
  ensure_ascii=False)``, no absolute paths, no timestamps, and no floats
  outside a ``derived_for_display`` block;
* well-formed exact encodings: every exact rational is a two-integer array
  ``[numerator, denominator]`` in lowest terms with positive denominator
  (``oph.sim.qm_observer_viz.v1`` convention, ``oph_fpe/qm_observer/DESIGN.md``
  section 7); every scalar of ``Q(sqrt(3), i)`` is the four-string list
  ``[re_rational, re_sqrt3, im_rational, im_sqrt3]`` written by ``C3.encode``
  of ``oph_fpe/quantum/phase_operation.py``; matrices follow ``encode_matrix``
  (nested two-by-two lists) for dimension two and the ``encode_matrix_n``
  extension (row-major entries beside the integer ``n``) otherwise;
* required fields per context and outcome, per design section A: the Kraus
  families (field 1), ``effect_from_kraus`` against ``declared_effect`` with
  the exact residual matrix (field 2), the legacy-named
  ``trace_nonincreasing`` diagnostic trace values on the frozen spanning set
  of matrix units (field 3), the summed-channel Kraus
  normalization and trace checks (field 4), the readback counts with the
  exact compatibility residual (field 5), the preparation coordinatization
  with its positivity certificate, operation list, and content hash
  (field 6), provenance (field 7), and labels (field 8);
* dimensional consistency of every matrix against the declared carrier
  dimension;
* exact algebraic conformance in rational and symbolic ``Q(sqrt(3))``
  arithmetic, with no floating-point value in any check: each Kraus family's
  induced effect ``sum_k K_k^dagger K_k`` equals the declared effect and is
  positive semidefinite; each outcome's complement equals the other induced
  effect and is positive semidefinite; per-context effects sum to the
  identity; declared fixture-state outcome traces ``Tr(rho E)`` equal the
  committed fixture frequencies ``count/mass``.  The matrix-unit trace data
  are recomputation diagnostics only: their differences need not be zero and
  do not certify an order inequality;
* the committed-table cross-check: where the export declares the committed
  eight contexts, its effect table, run state, and count table must equal
  the committed literals, hard-coded below as reference constants with their
  sources cited;
* the custody digest: the canonical SHA-256 over the content, recomputed and
  compared with the export's declared digest.

Verdict grammar: ``STATIC_COMMITTED_FIXTURE_CONFORMANT`` for a conformant
export whose ``provenance_class`` is ``synthetic``, and ``NONCONFORMANT``
otherwise.  The reserved ``producer`` class always fails with
``PRODUCER_AUTHENTICATION_UNIMPLEMENTED``.  Marker-free strings and
self-declared hashes cannot authenticate a run, so v1 intentionally has no
producer success verdict.

Digest convention: the digest is the lowercase hexadecimal SHA-256 of the
UTF-8 bytes of the canonical serialization of the export object with the
``custody_digest_sha256`` field removed, with no trailing line break.  The
design document names the committed ``canonical_sha256`` of
``oph_fpe/core/charged_response.py`` together with the indent-2 canonical
text and the convention of ``plan/SIM_ALIGNMENT_2026-08-20.md`` (the JSON
text without a trailing newline).  The committed function re-serializes with
compact separators, so the two conventions differ on which byte stream is
hashed; this validator binds the digest to the indent-2 canonical text, and
``VALIDATOR_CONTRACT.md`` records that binding.

Committed reference constants and their sources:

* the eight context names, the effect table, the run-state diagonal
  ``diag(111/179, 68/179)``, and the count literals ``(111, 68, 179)``,
  ``(315, 401, 716)``, ``(179, 179, 358)`` transcribe
  ``Lean/EventAlgebra/LuedersPhaseInstrument.lean``
  (``luedersPhaseInstrument_run_table``, ``committedRunState_eq_literal``,
  ``committedContextEffect_diagonal_eq``, ``sourcePhaseLift_entries``,
  ``binaryFrequency_rotated_run``) and
  ``Lean/EventAlgebra/SourceBoundInstrumentInterface.lean``
  (``SourceBoundDeterminedData.publicTable_literal``,
  ``modelFrequency_phase_zero`` with the phase count literals
  ``(179, 179)``), read against the simulator transcription in
  ``oph_fpe/quantum/phase_operation.py`` (``named_effects``,
  ``record_diagonal_state``, ``outcome_table``);
* the eight outcome-0 frequencies ``111/179`` (three entries), ``315/716``
  (four entries), ``1/2`` (one entry);
* the off-diagonal modulus bound ``7548/32041`` of the committed diagonal,
  the specialization recorded in design section A field 6
  (``EventAlgebra.prep_offdiag_normSq_le``);
* the constructor list of ``EventAlgebra.SourceReachability.Reachable``
  from ``Lean/EventAlgebra/SourceReachabilityDelimitation.lean``, as
  transcribed in design section A field 6.

What is not proved here: the checker certifies static-fixture schema and exact
algebraic conformance of one JSON document and the integrity of its declared
digest, and nothing else.  It cannot certify source production, provenance,
custody, or run reality.  The current placeholder fields in
``SourceBoundInstrumentInterface.lean`` are freely stipulable and
non-authenticating (``committed_corpus_does_not_determine_binding``,
``binding_digest_free_parameter``), so a passing validation is attainable by
freely stipulated data and discharges no register row.  No run exists; the
shipped sample is synthetic and marked; register rows PR-03, PR-64, and
PR-65 are open, and nothing here discharges the source-production row.  The
frozen decision rule of design section C, including the ``TOL_READBACK``
band for produced sampling residuals, binds at freeze time and is outside
this checker; the readback checks here are the exact-equality reading that a
deterministic committed-table fixture satisfies.  A future finite-run v2
validator must be separately frozen and implemented.

Pure standard-library Python.  All checks run in exact arithmetic over
``fractions.Fraction`` and a minimal exact ``Q(sqrt(3))`` implementation;
float values are detected as data and never enter a computation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any

SCHEMA_ID = "oph.sim.ins03_phase_instrument_export.v1"
VALIDATOR_ID = "ins03_export_validator.v1"
STATIC_FIXTURE_VERDICT = "STATIC_COMMITTED_FIXTURE_CONFORMANT"

PROVENANCE_CLASSES = ("synthetic", "producer")
SYNTHETIC_MARKER = "SYNTHETIC_PLACEHOLDER"
NULL_COMMIT = "0" * 40
NULL_SHA256 = "0" * 64

HEX40_RE = re.compile(r"^[0-9a-f]{40}$")
HEX64_RE = re.compile(r"^[0-9a-f]{64}$")
TIMESTAMP_RE = re.compile(r"\d{4}-\d{2}-\d{2}[T ]\d{2}:\d{2}")

# Constructor list of EventAlgebra.SourceReachability.Reachable
# (Lean/EventAlgebra/SourceReachabilityDelimitation.lean), transcribed from
# design section A field 6.
REACHABLE_CONSTRUCTORS = (
    "base247", "base88", "baseTriple",
    "step247", "step88", "stepTriple",
    "swap", "anchor247", "anchor88",
    "margSnd247", "margFst247", "margSnd88", "margFst88",
    "margSndTriple", "margMiddle",
    "expectL247", "expectR247", "expectL88", "expectR88",
    "scal", "prodMarg", "mix",
)

# Excluded-input set floor of design section C.2 (c), carried in the sample
# provenance block.
EXCLUDED_INPUT_FLOOR = (
    "oph_fpe/quantum/phase_operation.py",
    "oph_fpe/qm_observer/tables.py",
    "code/phase_operation_producer/PHASE_OPERATION_RECEIPT.v1.json",
)


# ---------------------------------------------------------------------------
# Exact arithmetic: Q, Q(sqrt(3)), and Q(sqrt(3), i)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Q3:
    """Exact element ``a + b*sqrt(3)`` of the real field Q(sqrt(3))."""

    a: Fraction = Fraction(0)
    b: Fraction = Fraction(0)

    def __add__(self, other: "Q3") -> "Q3":
        return Q3(self.a + other.a, self.b + other.b)

    def __sub__(self, other: "Q3") -> "Q3":
        return Q3(self.a - other.a, self.b - other.b)

    def __neg__(self) -> "Q3":
        return Q3(-self.a, -self.b)

    def __mul__(self, other: "Q3") -> "Q3":
        return Q3(
            self.a * other.a + 3 * self.b * other.b,
            self.a * other.b + self.b * other.a,
        )

    def is_zero(self) -> bool:
        return self.a == 0 and self.b == 0

    def is_nonneg(self) -> bool:
        """Exact sign of ``a + b*sqrt(3)`` by integer comparison only."""
        if self.b == 0:
            return self.a >= 0
        if self.a >= 0 and self.b >= 0:
            return True
        if self.a <= 0 and self.b <= 0:
            return False
        if self.a >= 0:  # b < 0: nonneg iff a >= |b|*sqrt(3)
            return self.a * self.a >= 3 * self.b * self.b
        # a < 0, b > 0: nonneg iff b*sqrt(3) >= |a|
        return 3 * self.b * self.b >= self.a * self.a

    def encode(self) -> list[str]:
        return [str(self.a), str(self.b)]


Q3_ZERO = Q3()


def q3_le(x: Q3, y: Q3) -> bool:
    return (y - x).is_nonneg()


@dataclass(frozen=True)
class C3:
    """Exact element ``re + i*im`` of Q(sqrt(3), i) with ``re, im`` in Q3."""

    re: Q3 = Q3_ZERO
    im: Q3 = Q3_ZERO

    def __add__(self, other: "C3") -> "C3":
        return C3(self.re + other.re, self.im + other.im)

    def __sub__(self, other: "C3") -> "C3":
        return C3(self.re - other.re, self.im - other.im)

    def __neg__(self) -> "C3":
        return C3(-self.re, -self.im)

    def __mul__(self, other: "C3") -> "C3":
        return C3(
            self.re * other.re - self.im * other.im,
            self.re * other.im + self.im * other.re,
        )

    def conjugate(self) -> "C3":
        return C3(self.re, -self.im)

    def is_zero(self) -> bool:
        return self.re.is_zero() and self.im.is_zero()

    def norm_sq(self) -> Q3:
        return self.re * self.re + self.im * self.im

    def encode(self) -> list[str]:
        """The reference scalar encoding
        ``[re_rational, re_sqrt3, im_rational, im_sqrt3]`` over Q, matching
        ``C3.encode`` of ``oph_fpe/quantum/phase_operation.py``."""
        return [str(self.re.a), str(self.re.b), str(self.im.a), str(self.im.b)]


C3_ZERO = C3()


def c3_rational(value: Fraction | int) -> C3:
    return C3(Q3(Fraction(value), Fraction(0)), Q3_ZERO)


def c3_sqrt3(value: Fraction | int) -> C3:
    return C3(Q3(Fraction(0), Fraction(value)), Q3_ZERO)


def c3_imag(value: Fraction | int) -> C3:
    return C3(Q3_ZERO, Q3(Fraction(value), Fraction(0)))


Matrix = list[list[C3]]


def mat_identity(n: int) -> Matrix:
    return [[c3_rational(1) if i == j else C3_ZERO for j in range(n)]
            for i in range(n)]


def mat_zero(n: int) -> Matrix:
    return [[C3_ZERO for _ in range(n)] for _ in range(n)]


def mat_add(x: Matrix, y: Matrix) -> Matrix:
    n = len(x)
    return [[x[i][j] + y[i][j] for j in range(n)] for i in range(n)]


def mat_sub(x: Matrix, y: Matrix) -> Matrix:
    n = len(x)
    return [[x[i][j] - y[i][j] for j in range(n)] for i in range(n)]


def mat_mul(x: Matrix, y: Matrix) -> Matrix:
    n = len(x)
    return [
        [sum((x[i][k] * y[k][j] for k in range(n)), C3_ZERO) for j in range(n)]
        for i in range(n)
    ]


def mat_dagger(x: Matrix) -> Matrix:
    n = len(x)
    return [[x[j][i].conjugate() for j in range(n)] for i in range(n)]


def mat_trace(x: Matrix) -> C3:
    return sum((x[i][i] for i in range(len(x))), C3_ZERO)


def mat_eq(x: Matrix, y: Matrix) -> bool:
    return x == y


def mat_is_zero(x: Matrix) -> bool:
    return all(entry.is_zero() for row in x for entry in row)


def mat_is_pos_semidefinite_2x2(x: Matrix) -> bool:
    """Exact Hermitian PSD test for the two-dimensional v1 carrier.

    A Hermitian two-by-two matrix is positive semidefinite exactly when its
    two diagonal entries and determinant are nonnegative.  All comparisons
    stay in ``Q(sqrt(3))``.  The v1 preparation field already rejects every
    carrier dimension other than two, so this is the complete PSD test needed
    by a conformant v1 document.
    """
    if len(x) != 2 or any(len(row) != 2 for row in x):
        return False
    if not mat_eq(x, mat_dagger(x)):
        return False
    diagonal = (x[0][0], x[1][1])
    if any(not entry.im.is_zero() or not entry.re.is_nonneg()
           for entry in diagonal):
        return False
    determinant = x[0][0] * x[1][1] - x[0][1] * x[1][0]
    return determinant.im.is_zero() and determinant.re.is_nonneg()


def mat_unit(n: int, i: int, j: int) -> Matrix:
    unit = mat_zero(n)
    unit[i][j] = c3_rational(1)
    return unit


def encode_matrix(x: Matrix) -> list[list[list[str]]]:
    """The ``encode_matrix`` convention of ``phase_operation.py`` for the
    two-by-two shape: nested lists of four-string scalar entries."""
    n = len(x)
    return [[x[i][j].encode() for j in range(n)] for i in range(n)]


def encode_matrix_n(x: Matrix) -> dict[str, Any]:
    """The producer-side ``encode_matrix_n`` extension of design section A:
    the per-entry four-string encoding in row-major order with the integer
    dimension carried beside the matrix."""
    return {"n": len(x), "entries": encode_matrix(x)}


# ---------------------------------------------------------------------------
# Committed reference constants (sources cited in the module docstring)
# ---------------------------------------------------------------------------

COMMITTED_CONTEXTS = (
    "web_diagonal",
    "web_conjugated_0",
    "web_conjugated_1",
    "web_conjugated_2",
    "web_conjugated_3",
    "web_conjugated_4",
    "web_conjugated_5",
    "phase",
)


def _record_projector() -> Matrix:
    return [[c3_rational(1), C3_ZERO], [C3_ZERO, C3_ZERO]]


def _rotated_projector(sign: int) -> Matrix:
    off = c3_sqrt3(Fraction(sign, 4))
    return [
        [c3_rational(Fraction(1, 4)), off],
        [off, c3_rational(Fraction(3, 4))],
    ]


def _phase_projector() -> Matrix:
    return [
        [c3_rational(Fraction(1, 2)), c3_imag(Fraction(-1, 2))],
        [c3_imag(Fraction(1, 2)), c3_rational(Fraction(1, 2))],
    ]


COMMITTED_EFFECT0: dict[str, Matrix] = {
    "web_diagonal": _record_projector(),
    "web_conjugated_0": _record_projector(),
    "web_conjugated_1": _record_projector(),
    "web_conjugated_2": _rotated_projector(-1),
    "web_conjugated_3": _rotated_projector(-1),
    "web_conjugated_4": _rotated_projector(+1),
    "web_conjugated_5": _rotated_projector(+1),
    "phase": _phase_projector(),
}

# (count_0, count_1, mass) per committed context.
COMMITTED_COUNTS: dict[str, tuple[int, int, int]] = {
    "web_diagonal": (111, 68, 179),
    "web_conjugated_0": (111, 68, 179),
    "web_conjugated_1": (111, 68, 179),
    "web_conjugated_2": (315, 401, 716),
    "web_conjugated_3": (315, 401, 716),
    "web_conjugated_4": (315, 401, 716),
    "web_conjugated_5": (315, 401, 716),
    "phase": (179, 179, 358),
}

# The eight outcome-0 entries of the committed table, in context order.
COMMITTED_OUTCOME0_FREQUENCIES = (
    Fraction(111, 179), Fraction(111, 179), Fraction(111, 179),
    Fraction(315, 716), Fraction(315, 716), Fraction(315, 716),
    Fraction(315, 716), Fraction(1, 2),
)

COMMITTED_RUN_STATE_DIAGONAL = (Fraction(111, 179), Fraction(68, 179))

# |rho_01|^2 <= 111*68/179^2 on the committed diagonal
# (EventAlgebra.prep_offdiag_normSq_le specialization).
COMMITTED_OFFDIAG_NORMSQ_BOUND = Fraction(7548, 32041)


# ---------------------------------------------------------------------------
# Canonical serialization and digest
# ---------------------------------------------------------------------------


def canonical_text(obj: Any) -> str:
    return json.dumps(obj, sort_keys=True, indent=2, ensure_ascii=False)


def compute_custody_digest(obj: dict[str, Any]) -> str:
    """SHA-256 over the canonical text of the export object with the
    ``custody_digest_sha256`` field removed, no trailing line break."""
    body = {k: v for k, v in obj.items() if k != "custody_digest_sha256"}
    return hashlib.sha256(canonical_text(body).encode("utf-8")).hexdigest()


def compute_preparation_content_sha256(
    carrier_id: str,
    producer_commit: str,
    source_record_ids: list[str],
    rho_00: list[int],
    rho_01: list[str],
    rho_11: list[int],
) -> str:
    """The canonical preparation body of design section C.2 (a): source
    record ids, carrier id, state coordinates, and producer commit, with no
    context field inside the hashed body."""
    body = {
        "carrier_id": carrier_id,
        "producer_commit": producer_commit,
        "source_record_ids": source_record_ids,
        "state": {"rho_00": rho_00, "rho_01": rho_01, "rho_11": rho_11},
    }
    return hashlib.sha256(canonical_text(body).encode("utf-8")).hexdigest()


# ---------------------------------------------------------------------------
# Error collection
# ---------------------------------------------------------------------------


class Errors:
    def __init__(self) -> None:
        self.items: list[dict[str, str]] = []

    def add(self, code: str, where: str, detail: str) -> None:
        self.items.append({"code": code, "where": where, "detail": detail})

    def codes(self) -> set[str]:
        return {item["code"] for item in self.items}


# ---------------------------------------------------------------------------
# Decoders (fail-closed; return None and record an error on malformation)
# ---------------------------------------------------------------------------


def _is_int(value: Any) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


def decode_rational(node: Any, where: str, errors: Errors) -> Fraction | None:
    if (not isinstance(node, list) or len(node) != 2
            or not _is_int(node[0]) or not _is_int(node[1])):
        errors.add("RATIONAL_ENCODING", where,
                   "expected a two-integer array [numerator, denominator]")
        return None
    numerator, denominator = node
    if denominator <= 0:
        errors.add("RATIONAL_ENCODING", where, "denominator must be positive")
        return None
    value = Fraction(numerator, denominator)
    if value.numerator != numerator or value.denominator != denominator:
        errors.add("RATIONAL_ENCODING", where, "array is not in lowest terms")
        return None
    return value


def _decode_fraction_string(text: Any, where: str, errors: Errors) -> Fraction | None:
    if not isinstance(text, str):
        errors.add("SCALAR_ENCODING", where, "scalar component is not a string")
        return None
    try:
        value = Fraction(text)
    except (ValueError, ZeroDivisionError):
        errors.add("SCALAR_ENCODING", where,
                   f"component {text!r} is not a rational literal")
        return None
    if str(value) != text:
        errors.add("SCALAR_ENCODING", where,
                   f"component {text!r} is not in canonical form")
        return None
    return value


def decode_q3(node: Any, where: str, errors: Errors) -> Q3 | None:
    if not isinstance(node, list) or len(node) != 2:
        errors.add("SCALAR_ENCODING", where,
                   "expected a two-string list [rational, sqrt3]")
        return None
    parts = [_decode_fraction_string(part, where, errors) for part in node]
    if any(part is None for part in parts):
        return None
    return Q3(parts[0], parts[1])


def decode_c3(node: Any, where: str, errors: Errors) -> C3 | None:
    if not isinstance(node, list) or len(node) != 4:
        errors.add("SCALAR_ENCODING", where,
                   "expected the four-string list "
                   "[re_rational, re_sqrt3, im_rational, im_sqrt3]")
        return None
    parts = [_decode_fraction_string(part, where, errors) for part in node]
    if any(part is None for part in parts):
        return None
    return C3(Q3(parts[0], parts[1]), Q3(parts[2], parts[3]))


def decode_matrix(node: Any, dimension: int, where: str,
                  errors: Errors) -> Matrix | None:
    """Decode ``encode_matrix`` (plain nested lists, dimension two only) or
    the ``encode_matrix_n`` extension (object with ``n`` and row-major
    ``entries``, dimensions other than two)."""
    if isinstance(node, dict):
        if dimension == 2:
            errors.add("MATRIX_ENCODING", where,
                       "dimension two uses the plain encode_matrix form")
            return None
        if set(node.keys()) != {"n", "entries"}:
            errors.add("MATRIX_ENCODING", where,
                       "encode_matrix_n object must carry exactly n and entries")
            return None
        if not _is_int(node["n"]) or node["n"] != dimension:
            errors.add("DIMENSION_MISMATCH", where,
                       f"declared n {node.get('n')!r} does not match carrier "
                       f"dimension {dimension}")
            return None
        rows = node["entries"]
    else:
        if dimension != 2:
            errors.add("MATRIX_ENCODING", where,
                       "dimensions other than two use the encode_matrix_n form")
            return None
        rows = node
    if (not isinstance(rows, list) or len(rows) != dimension
            or not all(isinstance(row, list) and len(row) == dimension
                       for row in rows)):
        errors.add("DIMENSION_MISMATCH", where,
                   f"expected a {dimension}x{dimension} matrix")
        return None
    matrix: Matrix = []
    for i, row in enumerate(rows):
        decoded_row: list[C3] = []
        for j, entry in enumerate(row):
            scalar = decode_c3(entry, f"{where}[{i}][{j}]", errors)
            if scalar is None:
                return None
            decoded_row.append(scalar)
        matrix.append(decoded_row)
    return matrix


# ---------------------------------------------------------------------------
# Content scans (floats, absolute paths, timestamps)
# ---------------------------------------------------------------------------


def _scan_content(node: Any, path: str, errors: Errors) -> None:
    if isinstance(node, dict):
        for key, value in node.items():
            if key == "derived_for_display":
                continue
            _scan_content(value, f"{path}.{key}" if path else key, errors)
    elif isinstance(node, list):
        for index, value in enumerate(node):
            _scan_content(value, f"{path}[{index}]", errors)
    elif isinstance(node, float):
        errors.add("FLOAT_OUTSIDE_DISPLAY", path,
                   "float value outside a derived_for_display block")
    elif isinstance(node, str):
        if node.startswith("/") or node.startswith("~") \
                or re.match(r"^[A-Za-z]:\\", node):
            errors.add("NONCANONICAL_PATH", path,
                       "absolute path outside a derived_for_display block")
        if TIMESTAMP_RE.search(node):
            errors.add("NONCANONICAL_TIMESTAMP", path,
                       "timestamp outside a derived_for_display block")


def _scan_synthetic_markers(node: Any, path: str, found: list[str]) -> None:
    if isinstance(node, dict):
        for key, value in node.items():
            if key == "derived_for_display":
                continue
            _scan_synthetic_markers(value, f"{path}.{key}" if path else key, found)
    elif isinstance(node, list):
        for index, value in enumerate(node):
            _scan_synthetic_markers(value, f"{path}[{index}]", found)
    elif isinstance(node, str):
        if SYNTHETIC_MARKER in node or node in (NULL_COMMIT, NULL_SHA256) \
                or "example.invalid" in node:
            found.append(path)


# ---------------------------------------------------------------------------
# Section validators
# ---------------------------------------------------------------------------

TOP_LEVEL_REQUIRED = (
    "schema", "provenance_class", "contexts", "carrier_dimension",
    "outcome_maps", "summed_channel", "readback", "preparation",
    "provenance", "labels", "custody_digest_sha256",
)
TOP_LEVEL_OPTIONAL = ("derived_for_display",)


def _require(mapping: Any, field: str, where: str, errors: Errors) -> Any:
    if not isinstance(mapping, dict) or field not in mapping:
        errors.add("MISSING_FIELD", where, f"required field {field!r} is absent")
        return None
    return mapping[field]


def _validate_outcome_map(node: Any, dimension: int, where: str,
                          errors: Errors) -> dict[str, Any] | None:
    """Fields 1 through 3 of design section A for one context and outcome.
    Returns decoded pieces for the cross-checks, or None."""
    if not isinstance(node, dict):
        errors.add("MISSING_FIELD", where, "outcome map is not an object")
        return None
    for field in ("kraus", "effect_from_kraus", "declared_effect",
                  "effect_residual", "trace_nonincreasing"):
        if field not in node:
            errors.add("MISSING_FIELD", where,
                       f"required field {field!r} is absent")
            return None

    kraus_nodes = node["kraus"]
    if not isinstance(kraus_nodes, list) or not kraus_nodes:
        errors.add("MISSING_FIELD", f"{where}.kraus",
                   "Kraus family must be a nonempty ordered list of matrices")
        return None
    kraus: list[Matrix] = []
    for index, entry in enumerate(kraus_nodes):
        matrix = decode_matrix(entry, dimension, f"{where}.kraus[{index}]", errors)
        if matrix is None:
            return None
        kraus.append(matrix)

    effect_from_kraus = decode_matrix(
        node["effect_from_kraus"], dimension, f"{where}.effect_from_kraus", errors)
    declared_effect = decode_matrix(
        node["declared_effect"], dimension, f"{where}.declared_effect", errors)
    effect_residual = decode_matrix(
        node["effect_residual"], dimension, f"{where}.effect_residual", errors)
    if effect_from_kraus is None or declared_effect is None \
            or effect_residual is None:
        return None

    # Field 2: induced effect from the Kraus family, exact.
    induced = mat_zero(dimension)
    for matrix in kraus:
        induced = mat_add(induced, mat_mul(mat_dagger(matrix), matrix))
    if not mat_eq(induced, effect_from_kraus):
        errors.add("EFFECT_FROM_KRAUS_MISMATCH", where,
                   "sum_k K_k^dagger K_k differs from the declared "
                   "effect_from_kraus")
    if not mat_eq(effect_from_kraus, declared_effect):
        errors.add("EFFECT_DECLARED_MISMATCH", where,
                   "effect_from_kraus differs from declared_effect")
    if dimension == 2 and not mat_is_pos_semidefinite_2x2(effect_from_kraus):
        errors.add("EFFECT_NOT_PSD", where,
                   "effect_from_kraus is not positive semidefinite under the "
                   "exact Hermitian two-by-two test")
    recomputed_residual = mat_sub(effect_from_kraus, declared_effect)
    if not mat_eq(effect_residual, recomputed_residual) \
            or not mat_is_zero(effect_residual):
        errors.add("EFFECT_RESIDUAL_MISMATCH", where,
                   "effect_residual must equal the recomputed difference and "
                   "be exactly zero")

    # Field 3: diagnostic trace values on the frozen spanning set of matrix
    # units, row-major order.  The schema key is retained for v1 compatibility
    # but is legacy-named: matrix units span the linear map, not the PSD cone,
    # and their trace differences are not required to vanish.  Outcome trace
    # nonincrease is certified below at the context level by E_i >= 0 and
    # I - E_i = E_(1-i) >= 0.
    trace_nodes = node["trace_nonincreasing"]
    expected_units = [(i, j) for i in range(dimension) for j in range(dimension)]
    ok_shape = (isinstance(trace_nodes, list)
                and len(trace_nodes) == len(expected_units))
    if ok_shape:
        for entry, (i, j) in zip(trace_nodes, expected_units):
            if not (isinstance(entry, dict)
                    and entry.get("basis_unit") == [i, j]):
                ok_shape = False
                break
    if not ok_shape:
        errors.add("SPANNING_SET_INCOMPLETE", f"{where}.trace_nonincreasing",
                   "expected one entry per matrix unit in row-major order")
    else:
        for entry, (i, j) in zip(trace_nodes, expected_units):
            entry_where = f"{where}.trace_nonincreasing[{i},{j}]"
            declared_output = decode_c3(entry.get("trace_output"),
                                        f"{entry_where}.trace_output", errors)
            declared_input = decode_c3(entry.get("trace_input"),
                                       f"{entry_where}.trace_input", errors)
            declared_diff = decode_c3(entry.get("difference"),
                                      f"{entry_where}.difference", errors)
            if None in (declared_output, declared_input, declared_diff):
                continue
            unit = mat_unit(dimension, i, j)
            output_matrix = mat_zero(dimension)
            for matrix in kraus:
                output_matrix = mat_add(
                    output_matrix,
                    mat_mul(mat_mul(matrix, unit), mat_dagger(matrix)))
            recomputed_output = mat_trace(output_matrix)
            recomputed_input = mat_trace(unit)
            if (declared_output != recomputed_output
                    or declared_input != recomputed_input
                    or declared_diff != recomputed_output - recomputed_input):
                errors.add("TRACE_CHECK_MISMATCH", entry_where,
                           "diagnostic trace values differ from the exact "
                           "recomputation over the Kraus family; no zero-"
                           "difference condition is imposed")

    return {"kraus": kraus, "effect_from_kraus": effect_from_kraus,
            "declared_effect": declared_effect}


def _validate_summed_channel(node: Any, outcome_data: dict[str, dict[str, Any]],
                             dimension: int, where: str, errors: Errors) -> None:
    """Field 4 of design section A for one context."""
    if not isinstance(node, dict):
        errors.add("MISSING_FIELD", where, "summed_channel entry is not an object")
        return
    for field in ("kraus_normalization", "kraus_normalization_residual",
                  "trace_checks"):
        if field not in node:
            errors.add("MISSING_FIELD", where,
                       f"required field {field!r} is absent")
            return
    normalization = decode_matrix(node["kraus_normalization"], dimension,
                                  f"{where}.kraus_normalization", errors)
    residual = decode_matrix(node["kraus_normalization_residual"], dimension,
                             f"{where}.kraus_normalization_residual", errors)
    if normalization is None or residual is None:
        return
    identity = mat_identity(dimension)
    recomputed = mat_zero(dimension)
    kraus_families = [data["kraus"] for data in outcome_data.values()]
    for family in kraus_families:
        for matrix in family:
            recomputed = mat_add(recomputed, mat_mul(mat_dagger(matrix), matrix))
    if not mat_eq(recomputed, normalization):
        errors.add("KRAUS_NORMALIZATION_MISMATCH", where,
                   "sum_i sum_k K_ik^dagger K_ik differs from the declared "
                   "kraus_normalization")
    if not mat_eq(normalization, identity):
        errors.add("KRAUS_NORMALIZATION_NOT_IDENTITY", where,
                   "declared kraus_normalization differs from the identity")
    if not mat_eq(residual, mat_sub(normalization, identity)) \
            or not mat_is_zero(residual):
        errors.add("KRAUS_NORMALIZATION_RESIDUAL", where,
                   "kraus_normalization_residual must equal the recomputed "
                   "difference and be exactly zero")

    trace_nodes = node["trace_checks"]
    expected_units = [(i, j) for i in range(dimension) for j in range(dimension)]
    ok_shape = (isinstance(trace_nodes, list)
                and len(trace_nodes) == len(expected_units))
    if ok_shape:
        for entry, (i, j) in zip(trace_nodes, expected_units):
            if not (isinstance(entry, dict)
                    and entry.get("basis_unit") == [i, j]):
                ok_shape = False
                break
    if not ok_shape:
        errors.add("SPANNING_SET_INCOMPLETE", f"{where}.trace_checks",
                   "expected one entry per matrix unit in row-major order")
        return
    for entry, (i, j) in zip(trace_nodes, expected_units):
        entry_where = f"{where}.trace_checks[{i},{j}]"
        declared_output = decode_c3(entry.get("trace_output"),
                                    f"{entry_where}.trace_output", errors)
        declared_input = decode_c3(entry.get("trace_input"),
                                   f"{entry_where}.trace_input", errors)
        declared_diff = decode_c3(entry.get("difference"),
                                  f"{entry_where}.difference", errors)
        if None in (declared_output, declared_input, declared_diff):
            continue
        unit = mat_unit(dimension, i, j)
        output_matrix = mat_zero(dimension)
        for family in kraus_families:
            for matrix in family:
                output_matrix = mat_add(
                    output_matrix,
                    mat_mul(mat_mul(matrix, unit), mat_dagger(matrix)))
        recomputed_output = mat_trace(output_matrix)
        recomputed_input = mat_trace(unit)
        if (declared_output != recomputed_output
                or declared_input != recomputed_input
                or declared_diff != recomputed_output - recomputed_input):
            errors.add("SUMMED_TRACE_MISMATCH", entry_where,
                       "declared summed-channel trace values differ from the "
                       "exact recomputation")


def _validate_preparation(node: Any, simulator_commit: str | None,
                          errors: Errors) -> dict[str, Any] | None:
    """Field 6 of design section A.  Returns the decoded state or None."""
    where = "preparation"
    if not isinstance(node, dict):
        errors.add("MISSING_FIELD", where, "preparation is not an object")
        return None
    for field in ("rho_00", "rho_11", "rho_01", "positivity_certificate",
                  "record_diagonal", "record_diagonal_offdiag", "carrier_id",
                  "source_record_ids", "operations",
                  "preparation_content_sha256"):
        if field not in node:
            errors.add("MISSING_FIELD", where,
                       f"required field {field!r} is absent")
            return None

    rho_00 = decode_rational(node["rho_00"], f"{where}.rho_00", errors)
    rho_11 = decode_rational(node["rho_11"], f"{where}.rho_11", errors)
    rho_01 = decode_c3(node["rho_01"], f"{where}.rho_01", errors)
    if rho_00 is None or rho_11 is None or rho_01 is None:
        return None
    if rho_00 + rho_11 != 1:
        errors.add("PREP_TRACE_NOT_ONE", where,
                   "diagonal entries must sum to exactly one")
    if rho_00 < 0 or rho_11 < 0:
        errors.add("PREP_DIAGONAL_NEGATIVE", where,
                   "diagonal entries must be nonnegative")

    certificate = node["positivity_certificate"]
    norm_sq = rho_01.norm_sq()
    product = Q3(rho_00 * rho_11, Fraction(0))
    if not (isinstance(certificate, dict)
            and set(certificate.keys()) == {"offdiag_norm_sq",
                                            "diagonal_product", "holds"}):
        errors.add("MISSING_FIELD", f"{where}.positivity_certificate",
                   "expected offdiag_norm_sq, diagonal_product, holds")
    else:
        declared_norm_sq = decode_q3(
            certificate["offdiag_norm_sq"],
            f"{where}.positivity_certificate.offdiag_norm_sq", errors)
        declared_product = decode_rational(
            certificate["diagonal_product"],
            f"{where}.positivity_certificate.diagonal_product", errors)
        if declared_norm_sq is not None and declared_norm_sq != norm_sq:
            errors.add("PREP_CERTIFICATE_MISMATCH", where,
                       "offdiag_norm_sq differs from the recomputed |rho_01|^2")
        if declared_product is not None and declared_product != rho_00 * rho_11:
            errors.add("PREP_CERTIFICATE_MISMATCH", where,
                       "diagonal_product differs from rho_00 * rho_11")
        if certificate["holds"] is not True:
            errors.add("PREP_POSITIVITY_VIOLATION", where,
                       "the positivity certificate must declare holds = true")
    if not q3_le(norm_sq, product):
        errors.add("PREP_POSITIVITY_VIOLATION", where,
                   "|rho_01|^2 <= rho_00 * rho_11 fails exactly")

    record_diagonal = node["record_diagonal"]
    declared_offdiag = decode_c3(node["record_diagonal_offdiag"],
                                 f"{where}.record_diagonal_offdiag", errors)
    if not isinstance(record_diagonal, bool) or declared_offdiag is None \
            or record_diagonal != rho_01.is_zero() \
            or declared_offdiag != rho_01:
        errors.add("PREP_RECORD_DIAGONAL_INCONSISTENT", where,
                   "record_diagonal and record_diagonal_offdiag must certify "
                   "the exported off-diagonal coordinate")

    if not isinstance(node["carrier_id"], str) or not node["carrier_id"]:
        errors.add("MISSING_FIELD", f"{where}.carrier_id",
                   "carrier_id must be a nonempty string")
    record_ids = node["source_record_ids"]
    if not isinstance(record_ids, list) \
            or not all(isinstance(r, str) and r for r in record_ids):
        errors.add("MISSING_FIELD", f"{where}.source_record_ids",
                   "source_record_ids must be a list of nonempty strings")

    operations = node["operations"]
    if not isinstance(operations, list):
        errors.add("MISSING_FIELD", f"{where}.operations",
                   "operations must be a list")
    else:
        for index, entry in enumerate(operations):
            op_where = f"{where}.operations[{index}]"
            if not (isinstance(entry, dict)
                    and isinstance(entry.get("operation"), str)
                    and isinstance(entry.get("carrier_before"), str)
                    and isinstance(entry.get("carrier_after"), str)
                    and isinstance(entry.get("in_class"), bool)):
                errors.add("MISSING_FIELD", op_where,
                           "expected operation, carrier_before, carrier_after, "
                           "in_class")
                continue
            in_constructor_list = entry["operation"] in REACHABLE_CONSTRUCTORS
            if entry["in_class"] != in_constructor_list:
                errors.add("PREP_OPERATION_CLASS_MISMATCH", op_where,
                           "in_class must equal membership in the constructor "
                           "list of EventAlgebra.SourceReachability.Reachable")
            if not in_constructor_list and not (
                    isinstance(entry.get("description"), str)
                    and entry["description"]):
                errors.add("PREP_OPERATION_DESCRIPTION_MISSING", op_where,
                           "an operation matching no constructor must carry "
                           "its own description")

    declared_hash = node["preparation_content_sha256"]
    if not isinstance(declared_hash, str) or not HEX64_RE.match(declared_hash):
        errors.add("PREP_HASH_MISMATCH", where,
                   "preparation_content_sha256 must be 64 lowercase hex digits")
    elif simulator_commit is not None and isinstance(node["carrier_id"], str) \
            and isinstance(record_ids, list):
        recomputed = compute_preparation_content_sha256(
            node["carrier_id"], simulator_commit, record_ids,
            node["rho_00"], node["rho_01"], node["rho_11"])
        if declared_hash != recomputed:
            errors.add("PREP_HASH_MISMATCH", where,
                       "preparation_content_sha256 differs from the "
                       "recomputed canonical-body hash")

    return {"rho_00": rho_00, "rho_11": rho_11, "rho_01": rho_01}


def _validate_provenance(node: Any, errors: Errors) -> str | None:
    """Field 7 of design section A.  Returns the simulator commit or None."""
    where = "provenance"
    if not isinstance(node, dict):
        errors.add("MISSING_FIELD", where, "provenance is not an object")
        return None
    for field in ("producer_modules", "simulator_commit", "repository_url",
                  "rer_commit", "run_id", "input_inventory",
                  "runtime_read_log", "import_graph_independence"):
        if field not in node:
            errors.add("MISSING_FIELD", where,
                       f"required field {field!r} is absent")
            return None
    modules = node["producer_modules"]
    if not isinstance(modules, list) or not all(
            isinstance(m, dict) and isinstance(m.get("path"), str)
            and isinstance(m.get("sha256"), str)
            and HEX64_RE.match(m["sha256"]) for m in modules):
        errors.add("PROVENANCE_FORMAT", f"{where}.producer_modules",
                   "expected a list of {path, sha256} objects with 64-hex "
                   "digests")
    for field in ("simulator_commit", "rer_commit"):
        if not isinstance(node[field], str) or not HEX40_RE.match(node[field]):
            errors.add("PROVENANCE_FORMAT", f"{where}.{field}",
                       "expected 40 lowercase hex digits")
    for field in ("repository_url", "run_id"):
        if not isinstance(node[field], str) or not node[field]:
            errors.add("PROVENANCE_FORMAT", f"{where}.{field}",
                       "expected a nonempty string")
    for field in ("input_inventory", "runtime_read_log"):
        if not isinstance(node[field], list) or not all(
                isinstance(entry, str) for entry in node[field]):
            errors.add("PROVENANCE_FORMAT", f"{where}.{field}",
                       "expected a list of strings")
    report = node["import_graph_independence"]
    if not (isinstance(report, dict)
            and isinstance(report.get("excluded_modules"), list)
            and _is_int(report.get("edge_count"))
            and _is_int(report.get("dynamic_import_count"))):
        errors.add("PROVENANCE_FORMAT", f"{where}.import_graph_independence",
                   "expected excluded_modules, integer edge_count, integer "
                   "dynamic_import_count")
    commit = node["simulator_commit"]
    return commit if isinstance(commit, str) else None


def _validate_labels(node: Any, provenance_class: str | None,
                     errors: Errors) -> None:
    """Field 8 of design section A."""
    where = "labels"
    if not isinstance(node, dict):
        errors.add("MISSING_FIELD", where, "labels is not an object")
        return
    for field in ("exploratory", "evidential", "claim_boundary"):
        if field not in node:
            errors.add("MISSING_FIELD", where,
                       f"required field {field!r} is absent")
            return
    if not isinstance(node["exploratory"], bool) \
            or not isinstance(node["evidential"], bool) \
            or not isinstance(node["claim_boundary"], str) \
            or not node["claim_boundary"]:
        errors.add("LABELS_INVALID", where,
                   "expected boolean exploratory, boolean evidential, and a "
                   "nonempty verbatim claim_boundary")
        return
    if provenance_class == "synthetic" and node["evidential"]:
        errors.add("SYNTHETIC_EVIDENTIAL_CONFLICT", where,
                   "a synthetic export must not be labeled evidential")


# ---------------------------------------------------------------------------
# The validator
# ---------------------------------------------------------------------------


def validate_export_text(text: str, export_path: str = "<text>") -> dict[str, Any]:
    """Validate one export document given as its file text.  Returns the
    report object; never raises on malformed content."""
    errors = Errors()
    report: dict[str, Any] = {
        "validator": VALIDATOR_ID,
        "export_path": export_path,
        "verdict": "NONCONFORMANT",
        "provenance_class": None,
        "errors": errors.items,
    }

    def _fail(code: str, where: str, detail: str) -> dict[str, Any]:
        errors.add(code, where, detail)
        return report

    def _reject_constant(name: str) -> None:
        raise ValueError(f"non-finite number {name!r}")

    try:
        parsed = json.loads(text, parse_constant=_reject_constant)
    except ValueError as error:
        return _fail("JSON_PARSE", "<document>", str(error))
    if not isinstance(parsed, dict):
        return _fail("JSON_PARSE", "<document>",
                     "top-level value must be an object")

    # Top-level field inventory.
    for field in TOP_LEVEL_REQUIRED:
        if field not in parsed:
            errors.add("MISSING_FIELD", "<top-level>",
                       f"required field {field!r} is absent")
    for field in parsed:
        if field not in TOP_LEVEL_REQUIRED + TOP_LEVEL_OPTIONAL:
            errors.add("UNEXPECTED_FIELD", "<top-level>",
                       f"field {field!r} is outside the schema")

    # Schema id and version.
    if parsed.get("schema") != SCHEMA_ID:
        errors.add("SCHEMA_ID_MISMATCH", "schema",
                   f"expected {SCHEMA_ID!r}, found {parsed.get('schema')!r}")

    provenance_class = parsed.get("provenance_class")
    if provenance_class not in PROVENANCE_CLASSES:
        errors.add("PROVENANCE_CLASS_INVALID", "provenance_class",
                   f"expected one of {PROVENANCE_CLASSES}, found "
                   f"{provenance_class!r}")
        provenance_class = None
    report["provenance_class"] = provenance_class

    # Canonical serialization of the file text.
    stripped = text[:-1] if text.endswith("\n") else text
    try:
        if stripped != canonical_text(parsed):
            errors.add("NONCANONICAL_SERIALIZATION", "<document>",
                       "file text differs from json.dumps(obj, sort_keys=True, "
                       "indent=2, ensure_ascii=False)")
    except (TypeError, ValueError) as error:
        errors.add("NONCANONICAL_SERIALIZATION", "<document>", str(error))

    # Custody digest.
    declared_digest = parsed.get("custody_digest_sha256")
    if not isinstance(declared_digest, str) or not HEX64_RE.match(declared_digest):
        errors.add("DIGEST_FORMAT", "custody_digest_sha256",
                   "expected 64 lowercase hex digits")
    else:
        try:
            recomputed_digest = compute_custody_digest(parsed)
        except (TypeError, ValueError):
            recomputed_digest = None
        if recomputed_digest is not None and declared_digest != recomputed_digest:
            errors.add("DIGEST_MISMATCH", "custody_digest_sha256",
                       "declared digest differs from the canonical SHA-256 "
                       "over the content")

    # Floats, absolute paths, timestamps outside derived_for_display.
    _scan_content(parsed, "", errors)

    # Contexts and carrier dimension.
    contexts = parsed.get("contexts")
    if not (isinstance(contexts, list) and contexts
            and all(isinstance(c, str) for c in contexts)
            and len(set(contexts)) == len(contexts)):
        errors.add("CONTEXT_COVERAGE", "contexts",
                   "expected a nonempty list of distinct context names")
        contexts = None
    dimension = parsed.get("carrier_dimension")
    if not _is_int(dimension) or dimension < 1:
        errors.add("DIMENSION_MISMATCH", "carrier_dimension",
                   "expected a positive integer")
        dimension = None

    # Provenance first: the preparation hash body consumes the commit.
    simulator_commit = None
    if "provenance" in parsed:
        simulator_commit = _validate_provenance(parsed["provenance"], errors)

    # Preparation (field 6).  The design fixes the two-by-two
    # coordinatization; a carrier dimension other than two has no committed
    # preparation surface here.
    prep = None
    if "preparation" in parsed:
        if dimension is not None and dimension != 2:
            errors.add("CARRIER_DIMENSION_UNSUPPORTED", "carrier_dimension",
                       "field 6 states the two-by-two coordinatization; "
                       "dimension is not 2")
        else:
            prep = _validate_preparation(parsed["preparation"],
                                         simulator_commit, errors)

    # Outcome maps (fields 1 through 3) and per-context effect sums.
    outcome_data: dict[str, dict[str, dict[str, Any]]] = {}
    if "outcome_maps" in parsed and contexts is not None and dimension is not None:
        outcome_maps = parsed["outcome_maps"]
        if not isinstance(outcome_maps, dict) \
                or set(outcome_maps.keys()) != set(contexts):
            errors.add("CONTEXT_COVERAGE", "outcome_maps",
                       "outcome_maps keys must equal the declared context list")
        else:
            for context in contexts:
                context_node = outcome_maps[context]
                if not isinstance(context_node, dict) \
                        or set(context_node.keys()) != {"0", "1"}:
                    errors.add("OUTCOME_COVERAGE", f"outcome_maps.{context}",
                               "expected exactly the outcome keys '0' and '1'")
                    continue
                decoded: dict[str, dict[str, Any]] = {}
                for outcome in ("0", "1"):
                    result = _validate_outcome_map(
                        context_node[outcome], dimension,
                        f"outcome_maps.{context}.{outcome}", errors)
                    if result is not None:
                        decoded[outcome] = result
                if set(decoded.keys()) == {"0", "1"}:
                    outcome_data[context] = decoded
                    identity = mat_identity(dimension)
                    effect0 = decoded["0"]["effect_from_kraus"]
                    effect1 = decoded["1"]["effect_from_kraus"]
                    total = mat_add(effect0, effect1)
                    if not mat_eq(total, identity):
                        errors.add("CONTEXT_SUM_NOT_IDENTITY",
                                   f"outcome_maps.{context}",
                                   "the two induced effects do not sum to the "
                                   "identity")
                    for outcome, effect, other_effect in (
                            ("0", effect0, effect1),
                            ("1", effect1, effect0)):
                        complement = mat_sub(identity, effect)
                        complement_where = \
                            f"outcome_maps.{context}.{outcome}"
                        if not mat_eq(complement, other_effect):
                            errors.add(
                                "EFFECT_COMPLEMENT_MISMATCH",
                                complement_where,
                                "identity minus this induced effect differs "
                                "from the other outcome's induced effect")
                        if dimension == 2 and not \
                                mat_is_pos_semidefinite_2x2(complement):
                            errors.add(
                                "EFFECT_COMPLEMENT_NOT_PSD",
                                complement_where,
                                "identity minus this induced effect is not "
                                "positive semidefinite under the exact "
                                "Hermitian two-by-two test; outcome trace "
                                "nonincrease is therefore uncertified")

    # Summed channels (field 4).
    if "summed_channel" in parsed and contexts is not None \
            and dimension is not None:
        summed = parsed["summed_channel"]
        if not isinstance(summed, dict) or set(summed.keys()) != set(contexts):
            errors.add("CONTEXT_COVERAGE", "summed_channel",
                       "summed_channel keys must equal the declared context "
                       "list")
        else:
            for context in contexts:
                if context in outcome_data:
                    _validate_summed_channel(
                        summed[context], outcome_data[context], dimension,
                        f"summed_channel.{context}", errors)

    # Readback (field 5).
    readback_counts: dict[str, tuple[int, int, int]] = {}
    if "readback" in parsed and contexts is not None:
        readback = parsed["readback"]
        if not isinstance(readback, dict) \
                or set(readback.keys()) != set(contexts):
            errors.add("CONTEXT_COVERAGE", "readback",
                       "readback keys must equal the declared context list")
        else:
            for context in contexts:
                context_node = readback[context]
                where = f"readback.{context}"
                if not isinstance(context_node, dict) \
                        or set(context_node.keys()) != {"0", "1"}:
                    errors.add("OUTCOME_COVERAGE", where,
                               "expected exactly the outcome keys '0' and '1'")
                    continue
                masses: list[int] = []
                counts: list[int] = []
                for outcome in ("0", "1"):
                    entry = context_node[outcome]
                    entry_where = f"{where}.{outcome}"
                    if not (isinstance(entry, dict)
                            and isinstance(entry.get("outcome_symbol"), str)
                            and _is_int(entry.get("count"))
                            and _is_int(entry.get("mass"))
                            and "compatibility_residual" in entry):
                        errors.add("MISSING_FIELD", entry_where,
                                   "expected outcome_symbol, integer count, "
                                   "integer mass, compatibility_residual")
                        continue
                    count, mass = entry["count"], entry["mass"]
                    if mass <= 0 or count < 0 or count > mass:
                        errors.add("READBACK_COUNT_INVALID", entry_where,
                                   "count must lie in [0, mass] with positive "
                                   "mass")
                        continue
                    masses.append(mass)
                    counts.append(count)
                    declared_residual = decode_rational(
                        entry["compatibility_residual"],
                        f"{entry_where}.compatibility_residual", errors)
                    if prep is not None and context in outcome_data:
                        effect = outcome_data[context][outcome]["effect_from_kraus"]
                        rho: Matrix = [
                            [c3_rational(prep["rho_00"]), prep["rho_01"]],
                            [prep["rho_01"].conjugate(),
                             c3_rational(prep["rho_11"])],
                        ]
                        trace = mat_trace(mat_mul(rho, effect))
                        if not trace.im.is_zero() or trace.re.b != 0:
                            errors.add("BORN_WEIGHT_NOT_RATIONAL", entry_where,
                                       "Tr(rho E) is not a plain rational")
                        else:
                            weight = trace.re.a
                            frequency = Fraction(count, mass)
                            if weight != frequency:
                                errors.add(
                                    "READBACK_TRACE_MISMATCH", entry_where,
                                    "declared run-state outcome trace differs "
                                    "from the declared frequency")
                            if declared_residual is not None \
                                    and declared_residual != weight - frequency:
                                errors.add(
                                    "READBACK_RESIDUAL_MISMATCH", entry_where,
                                    "compatibility_residual differs from the "
                                    "recomputed Tr(rho E) - count/mass")
                if len(masses) == 2:
                    if masses[0] != masses[1]:
                        errors.add("MASS_INCONSISTENT", where,
                                   "the two outcomes declare different masses")
                    elif counts[0] + counts[1] != masses[0]:
                        errors.add("COUNT_MASS_MISMATCH", where,
                                   "outcome counts do not sum to the context "
                                   "mass")
                    else:
                        readback_counts[context] = (counts[0], counts[1],
                                                    masses[0])

    # Labels (field 8).
    if "labels" in parsed:
        _validate_labels(parsed["labels"], provenance_class, errors)

    # The v1 checker has no authenticator and therefore no producer success
    # path.  Marker scanning remains useful diagnostic output, but even a
    # marker-free self-declaration fails closed: syntax, self-hashes, and a
    # canonical document digest do not establish that a producer ran.
    if provenance_class == "producer":
        found: list[str] = []
        _scan_synthetic_markers(parsed, "", found)
        for path in found:
            errors.add("SYNTHETIC_MARKER_IN_PRODUCER", path,
                       "a producer export must not carry synthetic "
                       "placeholder values")
        errors.add(
            "PRODUCER_AUTHENTICATION_UNIMPLEMENTED", "provenance_class",
            "the v1 checker is restricted to the static committed fixture; "
            "producer provenance is self-asserted and cannot receive a "
            "conformant or evidential verdict")

    # Committed-table cross-check: triggered when the export declares the
    # committed contexts.  A context set that overlaps the committed names
    # without equalling them is refused outright, so a renamed context can
    # never route a committed-shaped export around the cross-check.
    if contexts is not None and sorted(contexts) != sorted(COMMITTED_CONTEXTS) \
            and set(contexts) & set(COMMITTED_CONTEXTS):
        errors.add("PARTIAL_COMMITTED_CONTEXT_SET", "contexts",
                   "the context set overlaps the committed context names "
                   "without equalling them; a committed-shaped export must "
                   "declare exactly the committed context set")
    if contexts is not None and sorted(contexts) == sorted(COMMITTED_CONTEXTS):
        for context in COMMITTED_CONTEXTS:
            committed_effect = COMMITTED_EFFECT0[context]
            if context in outcome_data:
                declared = outcome_data[context]["0"]["declared_effect"]
                if not mat_eq(declared, committed_effect):
                    errors.add("COMMITTED_EFFECT_MISMATCH",
                               f"outcome_maps.{context}.0",
                               "declared_effect differs from the committed "
                               "table literal")
                declared_compl = outcome_data[context]["1"]["declared_effect"]
                if not mat_eq(declared_compl,
                              mat_sub(mat_identity(2), committed_effect)):
                    errors.add("COMMITTED_EFFECT_MISMATCH",
                               f"outcome_maps.{context}.1",
                               "declared_effect differs from the committed "
                               "complement literal")
            if context in readback_counts:
                if readback_counts[context] != COMMITTED_COUNTS[context]:
                    errors.add("COMMITTED_FREQUENCY_MISMATCH",
                               f"readback.{context}",
                               "counts and mass differ from the committed "
                               "count literals")
        if prep is not None:
            expected_00, expected_11 = COMMITTED_RUN_STATE_DIAGONAL
            if (prep["rho_00"] != expected_00 or prep["rho_11"] != expected_11
                    or not prep["rho_01"].is_zero()):
                errors.add("COMMITTED_RUN_STATE_MISMATCH", "preparation",
                           "preparation differs from the committed run-state "
                           "literal diag(111/179, 68/179)")

    if not errors.items and provenance_class == "synthetic":
        report["verdict"] = STATIC_FIXTURE_VERDICT
    return report


def validate_file(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as error:
        return {
            "validator": VALIDATOR_ID,
            "export_path": str(path),
            "verdict": "NONCONFORMANT",
            "provenance_class": None,
            "errors": [{"code": "FILE_UNREADABLE", "where": str(path),
                        "detail": str(error)}],
        }
    return validate_export_text(text, str(path))


# ---------------------------------------------------------------------------
# Synthetic sample builder
# ---------------------------------------------------------------------------
#
# The builder assembles the shipped sample_conforming_export.json from the
# committed objects only: the Lueders Kraus family of each committed context
# is the singleton {E} for outcome 0 and {1 - E} for outcome 1 (the effects
# are projectors, so each effect is its own square root), the preparation is
# the committed run-state diagonal, and the counts are the committed count
# literals.  Every producer-side field the corpus does not fix carries a
# clearly marked SYNTHETIC_PLACEHOLDER value, and provenance_class is
# "synthetic".  The builder certifies nothing: it exists so that the sample
# is regenerable and so that the tests can mutate a known-conformant
# document.


def _trace_entries(families: list[list[Matrix]], dimension: int) -> list[dict[str, Any]]:
    entries = []
    for i in range(dimension):
        for j in range(dimension):
            unit = mat_unit(dimension, i, j)
            output_matrix = mat_zero(dimension)
            for family in families:
                for matrix in family:
                    output_matrix = mat_add(
                        output_matrix,
                        mat_mul(mat_mul(matrix, unit), mat_dagger(matrix)))
            output = mat_trace(output_matrix)
            unit_trace = mat_trace(unit)
            entries.append({
                "basis_unit": [i, j],
                "trace_output": output.encode(),
                "trace_input": unit_trace.encode(),
                "difference": (output - unit_trace).encode(),
            })
    return entries


def build_synthetic_sample() -> dict[str, Any]:
    dimension = 2
    identity = mat_identity(dimension)
    simulator_commit = NULL_COMMIT

    outcome_maps: dict[str, Any] = {}
    summed_channel: dict[str, Any] = {}
    readback: dict[str, Any] = {}
    prep_00, prep_11 = COMMITTED_RUN_STATE_DIAGONAL
    rho: Matrix = [
        [c3_rational(prep_00), C3_ZERO],
        [C3_ZERO, c3_rational(prep_11)],
    ]

    for context in COMMITTED_CONTEXTS:
        effect0 = COMMITTED_EFFECT0[context]
        effect1 = mat_sub(identity, effect0)
        families = {"0": [effect0], "1": [effect1]}
        context_maps: dict[str, Any] = {}
        for outcome, family in families.items():
            induced = mat_zero(dimension)
            for matrix in family:
                induced = mat_add(induced, mat_mul(mat_dagger(matrix), matrix))
            context_maps[outcome] = {
                "kraus": [encode_matrix(matrix) for matrix in family],
                "effect_from_kraus": encode_matrix(induced),
                "declared_effect": encode_matrix(induced),
                "effect_residual": encode_matrix(mat_zero(dimension)),
                "trace_nonincreasing": _trace_entries([family], dimension),
            }
        outcome_maps[context] = context_maps

        both_families = [families["0"], families["1"]]
        normalization = mat_zero(dimension)
        for family in both_families:
            for matrix in family:
                normalization = mat_add(normalization,
                                        mat_mul(mat_dagger(matrix), matrix))
        summed_channel[context] = {
            "kraus_normalization": encode_matrix(normalization),
            "kraus_normalization_residual":
                encode_matrix(mat_sub(normalization, identity)),
            "trace_checks": _trace_entries(both_families, dimension),
        }

        count0, count1, mass = COMMITTED_COUNTS[context]
        entries: dict[str, Any] = {}
        for outcome, count, effect in (("0", count0, effect0),
                                       ("1", count1, effect1)):
            trace = mat_trace(mat_mul(rho, effect))
            residual = trace.re.a - Fraction(count, mass)
            entries[outcome] = {
                "outcome_symbol": outcome,
                "count": count,
                "mass": mass,
                "compatibility_residual": [residual.numerator,
                                           residual.denominator],
            }
        readback[context] = entries

    rho_00_enc = [prep_00.numerator, prep_00.denominator]
    rho_11_enc = [prep_11.numerator, prep_11.denominator]
    rho_01_enc = C3_ZERO.encode()
    carrier_id = f"{SYNTHETIC_MARKER}_CARRIER"
    source_record_ids = [f"{SYNTHETIC_MARKER}_RECORD_0"]
    diagonal_product = prep_00 * prep_11
    preparation = {
        "rho_00": rho_00_enc,
        "rho_11": rho_11_enc,
        "rho_01": rho_01_enc,
        "positivity_certificate": {
            "offdiag_norm_sq": Q3_ZERO.encode(),
            "diagonal_product": [diagonal_product.numerator,
                                 diagonal_product.denominator],
            "holds": True,
        },
        "record_diagonal": True,
        "record_diagonal_offdiag": C3_ZERO.encode(),
        "carrier_id": carrier_id,
        "source_record_ids": source_record_ids,
        "operations": [
            {
                "operation": f"{SYNTHETIC_MARKER}_OPERATION",
                "carrier_before": carrier_id,
                "carrier_after": carrier_id,
                "in_class": False,
                "description": (
                    "Placeholder entry: the name matches no constructor of "
                    "EventAlgebra.SourceReachability.Reachable and no source "
                    "operation was executed."
                ),
            }
        ],
        "preparation_content_sha256": compute_preparation_content_sha256(
            carrier_id, simulator_commit, source_record_ids,
            rho_00_enc, rho_01_enc, rho_11_enc),
    }

    provenance = {
        "producer_modules": [
            {"path": f"{SYNTHETIC_MARKER}/no_producer_module.py",
             "sha256": NULL_SHA256}
        ],
        "simulator_commit": simulator_commit,
        "repository_url": f"https://example.invalid/{SYNTHETIC_MARKER}",
        "rer_commit": NULL_COMMIT,
        "run_id": f"{SYNTHETIC_MARKER}_RUN",
        "input_inventory": [],
        "runtime_read_log": [],
        "import_graph_independence": {
            "excluded_modules": list(EXCLUDED_INPUT_FLOOR),
            "edge_count": 0,
            "dynamic_import_count": 0,
            "statement": (f"{SYNTHETIC_MARKER}: no producer executed and no "
                          "import graph exists."),
        },
    }

    labels = {
        "exploratory": True,
        "evidential": False,
        "claim_boundary": (
            "Synthetic conformance sample assembled from committed objects. "
            "No run exists, no seed is drawn, no freeze event exists. The "
            "Kraus families, effects, counts, and preparation transcribe the "
            "committed table of Lean/EventAlgebra/LuedersPhaseInstrument.lean; "
            "the custody fields are placeholders. Nothing here certifies "
            "source production, provenance, or custody, and nothing "
            "discharges register rows PR-03, PR-64, or PR-65."
        ),
    }

    sample: dict[str, Any] = {
        "schema": SCHEMA_ID,
        "provenance_class": "synthetic",
        "contexts": list(COMMITTED_CONTEXTS),
        "carrier_dimension": dimension,
        "outcome_maps": outcome_maps,
        "summed_channel": summed_channel,
        "readback": readback,
        "preparation": preparation,
        "provenance": provenance,
        "labels": labels,
        "derived_for_display": {
            "note": "Float renderings sit only in this block.",
            "phase_outcome0_frequency_float": 0.5,
        },
    }
    sample["custody_digest_sha256"] = compute_custody_digest(sample)
    return sample


def sample_file_text() -> str:
    """The canonical file text of the shipped sample, with one trailing
    newline for the file form; the digest is computed without it."""
    return canonical_text(build_synthetic_sample()) + "\n"


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=("Fail-closed static committed-fixture checker for the "
                     f"{SCHEMA_ID} export interface; producer claims are "
                     "unsupported."))
    parser.add_argument("export_path", help="path to the export JSON file")
    args = parser.parse_args(argv)
    report = validate_file(args.export_path)
    print(json.dumps(report, indent=2))
    if report["errors"] and report["errors"][0]["code"] in ("FILE_UNREADABLE",
                                                            "JSON_PARSE"):
        return 2
    return 0 if report["verdict"] == STATIC_FIXTURE_VERDICT else 1


if __name__ == "__main__":
    sys.exit(main())
