#!/usr/bin/env python3
"""Deterministic producer for the PR-04 recorded-decision phase-operation receipt.

Register row PR-04 carries a recorded decision (dated 2026-08-18, lane issue
number 730): the exact phase lift ``I/2 - (2*sqrt(3)/3)*i*(Q*P - P*Q)``
declared in ``code/born_context_phase_lift`` is a declared architecture
operation under the row's axiomatize disposition. This producer computes the
complete outcome table of that declared operation and of the committed web
contexts on the committed record-diagonal run state, in exact arithmetic over
``Q(sqrt(3), i)``, and writes one receipt JSON.

Count semantics, stated exactly: this is an explicitly exhaustive
deterministic semantics, one of the two validation routes the committed
boundary text of ``Lean/EventAlgebra/OperationalPhaseInstrument.lean`` names
as acceptable. Each context count pair is the exact Born weight of the
context effect under the committed state, scaled to the least positive
integer multiple of the committed run mass 179 that clears the weight's
denominator. No sampling, no randomness, and no floating-point arithmetic
occurs anywhere; rerunning the producer reproduces the receipt byte for byte.

Nonclaims: the operation is declared, never source-produced; the counts are
produced by the declared semantics, never measured; no laboratory or
emergent-instrument claim follows, and the physical-attachment premises of
the register stay open.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any, Sequence

HERE = Path(__file__).resolve().parent
PAYLOAD_PATH = HERE.parent / "born_context_phase_lift" / "BORN_CONTEXT_WEB_PAYLOAD.v1.json"
RECEIPT_PATH = HERE / "PHASE_OPERATION_RECEIPT.v1.json"

EXPECTED_PAYLOAD_SHA256 = (
    "71a06f1c15192123cd09feb2386da702b572c8ac57c9b7633f5aa60c5d404e22"
)
EXPECTED_RUN_COMMIT = "b39b78faf894894ebe573571e0902ccfaaeac32a"
RUN_MASS = 179
RUN_COUNTS = (111, 68)
DECISION_DATE = "2026-08-18"

RECEIPT_SCHEMA = "oph.receipt.pr04_phase_operation.v1"

SEMANTICS_DECLARATION = (
    "Explicitly exhaustive deterministic semantics: every count pair is the "
    "exact Born weight of the context effect under the committed "
    "record-diagonal run state diag(111/179, 68/179), computed in exact "
    "arithmetic over Q(sqrt(3), i) and scaled to the least positive integer "
    "multiple of the committed run mass 179 that makes both counts integers. "
    "No sampling, no randomness, no statistical estimate, and no "
    "floating-point arithmetic enters at any step."
)

DECISION_STATEMENT = (
    "Recorded decision on register row PR-04 (disposition axiomatize): the "
    "committed exact lift I/2 - (2*sqrt(3)/3)*i*(Q*P - P*Q) of "
    "code/born_context_phase_lift is a declared architecture operation. The "
    "information cost is one declared non-diagonal operation with the exact "
    "matrix [[1/2, -i/2], [i/2, 1/2]], equal to the Pauli +Y projector."
)

NONCLAIMS = (
    "The operation is declared, never source-produced. The counts are "
    "produced by the declared deterministic semantics, never measured. No "
    "laboratory or emergent-instrument claim follows; the physical-attachment "
    "premises of the register stay open."
)


class ProductionError(ValueError):
    """Raised when a pinned input or an exact invariant fails."""


def need(condition: bool, message: str) -> None:
    if not condition:
        raise ProductionError(message)


@dataclass(frozen=True)
class R3:
    """An exact element ``a + b*sqrt(3)`` of the real field Q(sqrt(3))."""

    a: Fraction = Fraction(0)
    b: Fraction = Fraction(0)

    @staticmethod
    def of(value: Any) -> "R3":
        if isinstance(value, R3):
            return value
        return R3(Fraction(value), Fraction(0))

    def __add__(self, other: "R3") -> "R3":
        return R3(self.a + other.a, self.b + other.b)

    def __sub__(self, other: "R3") -> "R3":
        return R3(self.a - other.a, self.b - other.b)

    def __neg__(self) -> "R3":
        return R3(-self.a, -self.b)

    def __mul__(self, other: "R3") -> "R3":
        return R3(
            self.a * other.a + 3 * self.b * other.b,
            self.a * other.b + self.b * other.a,
        )

    def is_zero(self) -> bool:
        return self.a == 0 and self.b == 0


@dataclass(frozen=True)
class X3:
    """An exact element ``re + i*im`` of Q(sqrt(3), i)."""

    re: R3 = R3()
    im: R3 = R3()

    @staticmethod
    def of(value: Any) -> "X3":
        if isinstance(value, X3):
            return value
        return X3(R3.of(value), R3())

    def __add__(self, other: "X3") -> "X3":
        return X3(self.re + other.re, self.im + other.im)

    def __sub__(self, other: "X3") -> "X3":
        return X3(self.re - other.re, self.im - other.im)

    def __mul__(self, other: "X3") -> "X3":
        return X3(
            self.re * other.re - self.im * other.im,
            self.re * other.im + self.im * other.re,
        )

    def conj(self) -> "X3":
        return X3(self.re, -self.im)

    def encode(self) -> list[str]:
        return [
            str(self.re.a),
            str(self.re.b),
            str(self.im.a),
            str(self.im.b),
        ]


Mat = tuple[tuple[X3, X3], tuple[X3, X3]]


def mat(rows: Sequence[Sequence[Any]]) -> Mat:
    need(len(rows) == 2 and all(len(r) == 2 for r in rows), "expected 2x2 matrix")
    return tuple(tuple(X3.of(x) for x in r) for r in rows)  # type: ignore[return-value]


IDENTITY = mat(((1, 0), (0, 1)))
RECORD = mat(((1, 0), (0, 0)))


def madd(x: Mat, y: Mat) -> Mat:
    return mat([[x[i][j] + y[i][j] for j in range(2)] for i in range(2)])


def msub(x: Mat, y: Mat) -> Mat:
    return mat([[x[i][j] - y[i][j] for j in range(2)] for i in range(2)])


def mmul(x: Mat, y: Mat) -> Mat:
    return mat(
        [
            [x[i][0] * y[0][j] + x[i][1] * y[1][j] for j in range(2)]
            for i in range(2)
        ]
    )


def mscale(s: X3, x: Mat) -> Mat:
    return mat([[s * x[i][j] for j in range(2)] for i in range(2)])


def mdagger(x: Mat) -> Mat:
    return mat([[x[j][i].conj() for j in range(2)] for i in range(2)])


def mtranspose(x: Mat) -> Mat:
    return mat([[x[j][i] for j in range(2)] for i in range(2)])


def mtrace(x: Mat) -> X3:
    return x[0][0] + x[1][1]


def encode_matrix(x: Mat) -> list[list[list[str]]]:
    return [[x[i][j].encode() for j in range(2)] for i in range(2)]


def decode_payload_scalar(value: Any) -> X3:
    need(isinstance(value, list) and len(value) == 2, "bad payload scalar")
    return X3(R3(Fraction(value[0]), Fraction(value[1])), R3())


def decode_payload_matrix(value: Any) -> Mat:
    need(isinstance(value, list) and len(value) == 2, "bad payload matrix rows")
    need(all(isinstance(r, list) and len(r) == 2 for r in value), "bad payload matrix")
    return mat([[decode_payload_scalar(value[i][j]) for j in range(2)] for i in range(2)])


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def phase_lift(commutator: Mat, coefficient: R3 | None = None) -> Mat:
    """The declared operation ``I/2 - c*i*[Q, P]`` with committed
    ``c = 2*sqrt(3)/3``; the coefficient knob exists for mutation tests."""

    if coefficient is None:
        coefficient = R3(Fraction(0), Fraction(2, 3))
    return msub(
        mscale(X3.of(Fraction(1, 2)), IDENTITY),
        mscale(X3(R3(), coefficient), commutator),
    )


def exact_rational_weight(state: Mat, effect: Mat) -> Fraction:
    """The Born weight ``Tr(state * effect)``, demanded to be a plain
    rational in [0, 1]."""

    weight = mtrace(mmul(state, effect))
    need(weight.im.is_zero(), "Born weight has nonzero imaginary part")
    need(weight.re.b == 0, "Born weight has nonzero sqrt(3) part")
    value = weight.re.a
    need(0 <= value <= 1, "Born weight outside [0, 1]")
    return value


def deterministic_counts(weight: Fraction) -> tuple[int, int, int, int]:
    """Counts of the exhaustive deterministic semantics: the least positive
    multiple ``k`` of the run mass with ``weight * k * RUN_MASS`` integral,
    the mass ``k * RUN_MASS``, and the two integer counts."""

    k = 1
    while (weight * (k * RUN_MASS)).denominator != 1:
        k += 1
        need(k <= weight.denominator, "no admissible run-mass multiple")
    mass = k * RUN_MASS
    first = int(weight * mass)
    return k, mass, first, mass - first


def build_receipt(payload: dict[str, Any], payload_sha256: str) -> dict[str, Any]:
    need(
        payload.get("schema") == "oph.sim.born_context_web_payload.v1",
        "payload schema drift",
    )
    provenance = payload.get("provenance", {})
    need(
        provenance.get("run_git_commit") == EXPECTED_RUN_COMMIT,
        "pinned run commit drift",
    )

    diagonal = payload.get("context_web", {}).get("diagonal_context", {})
    need(
        diagonal.get("realized_outcome_counts") == list(RUN_COUNTS),
        "diagonal counts drift",
    )
    need(diagonal.get("block_restricted_mass") == RUN_MASS, "diagonal mass drift")
    need(decode_payload_matrix(diagonal.get("projectors")[0]) == RECORD,
         "record projector drift")

    elements = payload.get("irrep", {}).get("elements")
    need(isinstance(elements, list) and len(elements) == 6, "expected six S3 elements")
    reps: list[Mat] = []
    for index, row in enumerate(elements):
        need(row.get("index") == index, "S3 element index drift")
        rep = decode_payload_matrix(row.get("matrix"))
        need(mmul(rep, mtranspose(rep)) == IDENTITY, f"element {index} not orthogonal")
        reps.append(rep)

    conjugated = payload.get("context_web", {}).get("conjugated_contexts")
    need(isinstance(conjugated, list) and len(conjugated) == 5,
         "expected five conjugated contexts")
    projectors: list[Mat] = []
    for g in range(6):
        projector = mmul(mmul(reps[g], RECORD), mtranspose(reps[g]))
        if g >= 1:
            declared = decode_payload_matrix(conjugated[g - 1].get("projectors")[0])
            need(projector == declared, f"conjugated projector {g} drift")
        projectors.append(projector)

    commutator = msub(mmul(projectors[3], RECORD), mmul(RECORD, projectors[3]))
    operation = phase_lift(commutator)
    plus_y = mat(
        (
            (Fraction(1, 2), X3(R3(), R3.of(Fraction(-1, 2)))),
            (X3(R3(), R3.of(Fraction(1, 2))), Fraction(1, 2)),
        )
    )
    need(operation == plus_y, "declared operation is not the Pauli +Y projector")
    need(mdagger(operation) == operation, "declared operation not Hermitian")
    need(mmul(operation, operation) == operation, "declared operation not idempotent")

    state = mat(
        (
            (Fraction(RUN_COUNTS[0], RUN_MASS), 0),
            (0, Fraction(RUN_COUNTS[1], RUN_MASS)),
        )
    )
    need(mtrace(state) == X3.of(1), "committed state trace drift")

    contexts: list[dict[str, Any]] = []
    named_effects: list[tuple[str, Mat]] = [("web_diagonal", RECORD)]
    named_effects += [
        (f"web_conjugated_{g}", projectors[g]) for g in range(6)
    ]
    named_effects.append(("phase", operation))
    for name, effect in named_effects:
        weight = exact_rational_weight(state, effect)
        k, mass, first, second = deterministic_counts(weight)
        need(first > 0 and second > 0, f"context {name}: an outcome has zero mass")
        contexts.append(
            {
                "context": name,
                "effect": encode_matrix(effect),
                "born_weight": str(weight),
                "run_mass_multiple": k,
                "mass": mass,
                "counts": [first, second],
            }
        )

    diagonal_row = contexts[0]
    need(diagonal_row["counts"] == list(RUN_COUNTS), "diagonal semantics drift")

    phase_row = contexts[-1]
    a, b = phase_row["counts"]
    window_lhs = 32041 * (a - b) ** 2
    window_rhs = 30192 * (a + b) ** 2
    need(window_lhs <= window_rhs, "phase counts violate the receipt window")

    return {
        "schema": RECEIPT_SCHEMA,
        "recorded_decision": {
            "register_row": "PR-04",
            "disposition": "axiomatize",
            "date": DECISION_DATE,
            "lane_issue": 730,
            "statement": DECISION_STATEMENT,
        },
        "inputs": {
            "payload_path": "code/born_context_phase_lift/BORN_CONTEXT_WEB_PAYLOAD.v1.json",
            "payload_sha256": payload_sha256,
            "run_git_commit": EXPECTED_RUN_COMMIT,
            "run_id": provenance.get("run_id"),
            "seed": provenance.get("seed"),
            "run_mass": RUN_MASS,
            "run_counts": list(RUN_COUNTS),
        },
        "declared_operation": {
            "formula": "I/2 - (2*sqrt(3)/3)*i*(Q*P - P*Q)",
            "matrix": encode_matrix(operation),
            "equals_pauli_y_plus_projector": True,
            "scalar_encoding": "each scalar is [re_rational, re_sqrt3, im_rational, im_sqrt3] over Q",
        },
        "committed_state": {
            "matrix": encode_matrix(state),
            "reading": "record-diagonal state of the committed run literals (111, 68) at mass 179",
        },
        "semantics": SEMANTICS_DECLARATION,
        "contexts": contexts,
        "phase_receipt_window": {
            "inequality": "32041*(a - b)^2 <= 30192*(a + b)^2",
            "lhs": window_lhs,
            "rhs": window_rhs,
            "satisfied": True,
        },
        "nonclaims": NONCLAIMS,
    }


def receipt_bytes(receipt: dict[str, Any]) -> bytes:
    return (json.dumps(receipt, indent=1, sort_keys=True) + "\n").encode("utf-8")


def produce(payload_path: Path = PAYLOAD_PATH) -> bytes:
    raw = payload_path.read_bytes()
    digest = sha256_bytes(raw)
    need(digest == EXPECTED_PAYLOAD_SHA256, "pinned payload hash mismatch")
    payload = json.loads(raw.decode("utf-8"))
    return receipt_bytes(build_receipt(payload, digest))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--payload", type=Path, default=PAYLOAD_PATH)
    parser.add_argument("--output", type=Path, default=RECEIPT_PATH)
    parser.add_argument(
        "--check",
        action="store_true",
        help="fail when the committed receipt differs from the deterministic rerun",
    )
    args = parser.parse_args()
    data = produce(args.payload)
    if args.check:
        committed = args.output.read_bytes() if args.output.is_file() else b""
        need(committed == data, "committed receipt differs from deterministic rerun")
        print("phase-operation receipt: deterministic rerun matches the committed receipt")
        return
    args.output.write_bytes(data)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
