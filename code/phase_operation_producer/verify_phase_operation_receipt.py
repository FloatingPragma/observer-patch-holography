#!/usr/bin/env python3
"""Independent verifier for the PR-04 recorded-decision phase-operation receipt.

The verifier shares no arithmetic code with the producer. It reads the pinned
context-web payload and the committed receipt as adversarial inputs, rebuilds
the declared operation and the complete deterministic outcome table from the
payload alone in its own exact ``Q(sqrt(3), i)`` arithmetic (scalars as
4-tuples of ``Fraction``), and demands byte-level agreement with every
receipt field: input digests, the operation matrix, the committed state, the
per-context effects, Born weights, masses, counts, the run-mass-multiple
minimality of each mass, and the integer phase receipt window. Any mismatch
raises; there is no tolerance parameter and no floating-point arithmetic.

The verifier also demands the receipt's declaration fields: the recorded
decision naming register row PR-04 under the axiomatize disposition, the
exhaustive-deterministic-semantics declaration, and the nonclaims text. A
receipt that presents these counts as measured, sampled, or laboratory data
does not verify.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
PAYLOAD_PATH = HERE.parent / "born_context_phase_lift" / "BORN_CONTEXT_WEB_PAYLOAD.v1.json"
RECEIPT_PATH = HERE / "PHASE_OPERATION_RECEIPT.v1.json"

PINNED_PAYLOAD_SHA256 = (
    "71a06f1c15192123cd09feb2386da702b572c8ac57c9b7633f5aa60c5d404e22"
)
PINNED_RUN_COMMIT = "b39b78faf894894ebe573571e0902ccfaaeac32a"
RUN_MASS = 179
RUN_COUNTS = [111, 68]

# A scalar of Q(sqrt(3), i) is (re_rational, re_sqrt3, im_rational, im_sqrt3).
Scalar = tuple[Fraction, Fraction, Fraction, Fraction]
Matrix = list[list[Scalar]]


class VerificationError(ValueError):
    """Raised on any exact mismatch against the pinned inputs."""


def need(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def scal(re_rat: Any = 0, re_s3: Any = 0, im_rat: Any = 0, im_s3: Any = 0) -> Scalar:
    return (Fraction(re_rat), Fraction(re_s3), Fraction(im_rat), Fraction(im_s3))


def s_add(x: Scalar, y: Scalar) -> Scalar:
    return (x[0] + y[0], x[1] + y[1], x[2] + y[2], x[3] + y[3])


def s_neg(x: Scalar) -> Scalar:
    return (-x[0], -x[1], -x[2], -x[3])


def s_mul(x: Scalar, y: Scalar) -> Scalar:
    # (a + b s + i(c + d s)) * (e + f s + i(g + h s)) with s^2 = 3.
    a, b, c, d = x
    e, f, g, h = y
    re_rat = a * e + 3 * b * f - c * g - 3 * d * h
    re_s3 = a * f + b * e - c * h - d * g
    im_rat = a * g + 3 * b * h + c * e + 3 * d * f
    im_s3 = a * h + b * g + c * f + d * e
    return (re_rat, re_s3, im_rat, im_s3)


def s_conj(x: Scalar) -> Scalar:
    return (x[0], x[1], -x[2], -x[3])


def m_of(entries: list[list[Scalar]]) -> Matrix:
    need(len(entries) == 2 and all(len(r) == 2 for r in entries), "expected 2x2 matrix")
    return entries


def m_add(x: Matrix, y: Matrix) -> Matrix:
    return [[s_add(x[i][j], y[i][j]) for j in range(2)] for i in range(2)]


def m_sub(x: Matrix, y: Matrix) -> Matrix:
    return [[s_add(x[i][j], s_neg(y[i][j])) for j in range(2)] for i in range(2)]


def m_mul(x: Matrix, y: Matrix) -> Matrix:
    return [
        [
            s_add(s_mul(x[i][0], y[0][j]), s_mul(x[i][1], y[1][j]))
            for j in range(2)
        ]
        for i in range(2)
    ]


def m_scale(c: Scalar, x: Matrix) -> Matrix:
    return [[s_mul(c, x[i][j]) for j in range(2)] for i in range(2)]


def m_dagger(x: Matrix) -> Matrix:
    return [[s_conj(x[j][i]) for j in range(2)] for i in range(2)]


def m_transpose(x: Matrix) -> Matrix:
    return [[x[j][i] for j in range(2)] for i in range(2)]


def m_trace(x: Matrix) -> Scalar:
    return s_add(x[0][0], x[1][1])


M_IDENTITY = m_of([[scal(1), scal()], [scal(), scal(1)]])
M_RECORD = m_of([[scal(1), scal()], [scal(), scal()]])
M_PLUS_Y = m_of(
    [
        [scal(Fraction(1, 2)), scal(0, 0, Fraction(-1, 2))],
        [scal(0, 0, Fraction(1, 2)), scal(Fraction(1, 2))],
    ]
)


def decode_payload_matrix(value: Any) -> Matrix:
    need(isinstance(value, list) and len(value) == 2, "bad payload matrix rows")
    rows: Matrix = []
    for row in value:
        need(isinstance(row, list) and len(row) == 2, "bad payload matrix row")
        entries = []
        for cell in row:
            need(isinstance(cell, list) and len(cell) == 2, "bad payload scalar")
            entries.append(scal(Fraction(cell[0]), Fraction(cell[1])))
        rows.append(entries)
    return m_of(rows)


def decode_receipt_matrix(value: Any) -> Matrix:
    need(isinstance(value, list) and len(value) == 2, "bad receipt matrix rows")
    rows: Matrix = []
    for row in value:
        need(isinstance(row, list) and len(row) == 2, "bad receipt matrix row")
        entries = []
        for cell in row:
            need(isinstance(cell, list) and len(cell) == 4, "bad receipt scalar")
            entries.append(
                (
                    Fraction(cell[0]),
                    Fraction(cell[1]),
                    Fraction(cell[2]),
                    Fraction(cell[3]),
                )
            )
        rows.append(entries)
    return m_of(rows)


def check_window(counts: list[int]) -> tuple[int, int]:
    """The committed integer receipt window of
    ``Lean/EventAlgebra/PhaseInstrumentDetermination.lean`` on a phase count
    pair; raises when violated."""

    need(len(counts) == 2, "phase counts must be a pair")
    a, b = counts
    lhs = 32041 * (a - b) ** 2
    rhs = 30192 * (a + b) ** 2
    need(lhs <= rhs, "phase counts violate the committed receipt window")
    return lhs, rhs


def expected_table(payload: dict[str, Any]) -> tuple[Matrix, Matrix, list[tuple[str, Matrix]]]:
    """Rebuild the committed state, the declared operation, and the named
    effect table from the payload alone."""

    elements = payload.get("irrep", {}).get("elements")
    need(isinstance(elements, list) and len(elements) == 6, "expected six S3 elements")
    reps = [decode_payload_matrix(row.get("matrix")) for row in elements]
    for index, rep in enumerate(reps):
        need(m_mul(rep, m_transpose(rep)) == M_IDENTITY, f"element {index} not orthogonal")

    diagonal = payload.get("context_web", {}).get("diagonal_context", {})
    need(diagonal.get("realized_outcome_counts") == RUN_COUNTS, "diagonal counts drift")
    need(diagonal.get("block_restricted_mass") == RUN_MASS, "diagonal mass drift")
    need(
        decode_payload_matrix(diagonal.get("projectors")[0]) == M_RECORD,
        "record projector drift",
    )

    projectors = [
        m_mul(m_mul(reps[g], M_RECORD), m_transpose(reps[g])) for g in range(6)
    ]
    conjugated = payload.get("context_web", {}).get("conjugated_contexts")
    need(isinstance(conjugated, list) and len(conjugated) == 5,
         "expected five conjugated contexts")
    for g in range(1, 6):
        declared = decode_payload_matrix(conjugated[g - 1].get("projectors")[0])
        need(projectors[g] == declared, f"conjugated projector {g} drift")

    commutator = m_sub(
        m_mul(projectors[3], M_RECORD), m_mul(M_RECORD, projectors[3])
    )
    operation = m_sub(
        m_scale(scal(Fraction(1, 2)), M_IDENTITY),
        m_scale(scal(0, 0, 0, Fraction(2, 3)), commutator),
    )
    need(operation == M_PLUS_Y, "recomputed operation is not the Pauli +Y projector")
    need(m_dagger(operation) == operation, "recomputed operation not Hermitian")
    need(m_mul(operation, operation) == operation, "recomputed operation not idempotent")

    state = m_of(
        [
            [scal(Fraction(RUN_COUNTS[0], RUN_MASS)), scal()],
            [scal(), scal(Fraction(RUN_COUNTS[1], RUN_MASS))],
        ]
    )

    table: list[tuple[str, Matrix]] = [("web_diagonal", M_RECORD)]
    table += [(f"web_conjugated_{g}", projectors[g]) for g in range(6)]
    table.append(("phase", operation))
    return state, operation, table


def verify(
    receipt: dict[str, Any],
    payload: dict[str, Any],
    payload_sha256: str | None,
) -> dict[str, Any]:
    need(receipt.get("schema") == "oph.receipt.pr04_phase_operation.v1",
         "receipt schema drift")

    decision = receipt.get("recorded_decision", {})
    need(decision.get("register_row") == "PR-04", "decision register row drift")
    need(decision.get("disposition") == "axiomatize", "decision disposition drift")
    need(decision.get("date") == "2026-08-18", "decision date drift")
    need("declared architecture operation" in decision.get("statement", ""),
         "decision statement drift")

    inputs = receipt.get("inputs", {})
    need(inputs.get("run_git_commit") == PINNED_RUN_COMMIT, "run commit drift")
    need(inputs.get("payload_sha256") == PINNED_PAYLOAD_SHA256,
         "receipt payload digest drift")
    if payload_sha256 is not None:
        need(payload_sha256 == PINNED_PAYLOAD_SHA256, "payload file hash mismatch")
    need(inputs.get("run_mass") == RUN_MASS, "run mass drift")
    need(inputs.get("run_counts") == RUN_COUNTS, "run count literals drift")
    need(
        payload.get("provenance", {}).get("run_git_commit") == PINNED_RUN_COMMIT,
        "payload run commit drift",
    )

    semantics = receipt.get("semantics", "")
    need("exhaustive deterministic semantics" in semantics.lower(),
         "semantics declaration missing")
    need("no sampling" in semantics.lower(), "no-sampling declaration missing")
    nonclaims = receipt.get("nonclaims", "")
    need("never source-produced" in nonclaims, "nonclaims declaration missing")
    need("no laboratory or emergent-instrument claim" in nonclaims.lower(),
         "laboratory nonclaim missing")

    state, operation, table = expected_table(payload)

    declared_operation = decode_receipt_matrix(
        receipt.get("declared_operation", {}).get("matrix")
    )
    need(declared_operation == operation, "receipt operation matrix drift")
    need(
        receipt.get("declared_operation", {}).get("equals_pauli_y_plus_projector") is True,
        "operation projector flag drift",
    )
    declared_state = decode_receipt_matrix(receipt.get("committed_state", {}).get("matrix"))
    need(declared_state == state, "receipt committed-state drift")

    rows = receipt.get("contexts")
    need(isinstance(rows, list) and len(rows) == len(table), "context table size drift")
    for row, (name, effect) in zip(rows, table):
        where = f"context {name}"
        need(row.get("context") == name, f"{where}: name drift")
        need(decode_receipt_matrix(row.get("effect")) == effect,
             f"{where}: effect drift")
        weight = m_trace(m_mul(state, effect))
        need(weight[1] == 0 and weight[2] == 0 and weight[3] == 0,
             f"{where}: weight not a plain rational")
        value = weight[0]
        need(Fraction(row.get("born_weight")) == value, f"{where}: Born weight drift")
        k = row.get("run_mass_multiple")
        mass = row.get("mass")
        counts = row.get("counts")
        need(isinstance(k, int) and k >= 1, f"{where}: bad run-mass multiple")
        need(mass == k * RUN_MASS, f"{where}: mass is not the declared multiple")
        need(
            all((value * j * RUN_MASS).denominator != 1 for j in range(1, k)),
            f"{where}: run-mass multiple not minimal",
        )
        need(isinstance(counts, list) and len(counts) == 2, f"{where}: bad counts")
        need(all(isinstance(x, int) and x > 0 for x in counts),
             f"{where}: counts must be positive integers")
        need(counts[0] + counts[1] == mass, f"{where}: counts do not sum to the mass")
        need(Fraction(counts[0], mass) == value,
             f"{where}: counts are not the exact scaled Born weight")

    need(rows[0]["counts"] == RUN_COUNTS, "diagonal counts drift")
    phase_counts = rows[-1]["counts"]
    lhs, rhs = check_window(phase_counts)
    window = receipt.get("phase_receipt_window", {})
    need(window.get("lhs") == lhs and window.get("rhs") == rhs,
         "receipt window integers drift")
    need(window.get("satisfied") is True, "receipt window flag drift")

    return {
        "schema": "oph.verify.pr04_phase_operation_receipt.v1",
        "payload_sha256": PINNED_PAYLOAD_SHA256,
        "contexts_verified": len(rows),
        "phase_counts": phase_counts,
        "phase_window_lhs": lhs,
        "phase_window_rhs": rhs,
        "deterministic_semantics_declared": True,
        "verdict": (
            "Receipt verifies: the deterministic outcome table of the declared "
            "PR-04 operation on the committed run state is reproduced exactly "
            "from the pinned payload."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--payload", type=Path, default=PAYLOAD_PATH)
    parser.add_argument("--receipt", type=Path, default=RECEIPT_PATH)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    raw = args.payload.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    payload = json.loads(raw.decode("utf-8"))
    receipt = json.loads(args.receipt.read_text())
    report = verify(receipt, payload, digest)
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
