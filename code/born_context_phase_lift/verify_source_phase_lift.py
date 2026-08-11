#!/usr/bin/env python3
"""Independent exact verifier for the B13 source phase-lift boundary.

The verifier deliberately does not import the simulator extractor.  It reads
the vendored, hash-pinned ``BORN_CONTEXT_WEB_PAYLOAD.v1.json`` as an adversarial
input, pins its source-hash metadata, reconstructs the S3 representation and projector
orbit in exact ``Q(sqrt(3))`` arithmetic, and then works in
``Q(sqrt(3), i)`` to certify the following sharp boundary:

* the native real web has operator-span rank three and is blind to the two
  opposite Pauli-Y states;
* the phase lift built from two source-attached algebraic projector candidates
  is exactly the +Y projector and raises the operator-span rank to four;
* the payload still supplies no rotated or phase-lifted outcome receipt.

This is an operator-algebra certificate, not an instrument simulation.  No
floating-point arithmetic and no predicted outcomes are used.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable, Sequence


class VerificationError(ValueError):
    """Raised when an exact source or algebra invariant fails."""


def need(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


EXPECTED_PAYLOAD_SHA256 = (
    "71a06f1c15192123cd09feb2386da702b572c8ac57c9b7633f5aa60c5d404e22"
)
EXPECTED_RUN_COMMIT = "b39b78faf894894ebe573571e0902ccfaaeac32a"
EXPECTED_SOURCE_HASHES = {
    "receipt_sha256": "d6739274e9451295b8bb0334180231bfb1e516c03bfd6f2c80f70e4da64db749",
    "freezeout_sha256": "b962c5b80205a17d5d6bc023f5f5d487bc23c1327793978e7d1bc69494fac49e",
    "gauge_state_sha256": "5f2cd276b84b1917c4d321a532db6efd12442931fbc80fd764cbcf1918dadca2",
    "e1_payload_sha256": "005223dc4fa7442bf10e6e8b446a1e5ebba08e41796d6ac89c45251de090bf45",
}
HERE = Path(__file__).resolve().parent
VENDORED_PAYLOAD = HERE / "BORN_CONTEXT_WEB_PAYLOAD.v1.json"


@dataclass(frozen=True)
class Q3:
    """An exact element ``a + b sqrt(3)``."""

    a: Fraction = Fraction(0)
    b: Fraction = Fraction(0)

    @staticmethod
    def coerce(value: Any) -> "Q3":
        if isinstance(value, Q3):
            return value
        return Q3(Fraction(value), Fraction(0))

    def __add__(self, other: Any) -> "Q3":
        rhs = Q3.coerce(other)
        return Q3(self.a + rhs.a, self.b + rhs.b)

    __radd__ = __add__

    def __neg__(self) -> "Q3":
        return Q3(-self.a, -self.b)

    def __sub__(self, other: Any) -> "Q3":
        return self + (-Q3.coerce(other))

    def __rsub__(self, other: Any) -> "Q3":
        return Q3.coerce(other) - self

    def __mul__(self, other: Any) -> "Q3":
        rhs = Q3.coerce(other)
        return Q3(
            self.a * rhs.a + 3 * self.b * rhs.b,
            self.a * rhs.b + self.b * rhs.a,
        )

    __rmul__ = __mul__

    def inverse(self) -> "Q3":
        norm = self.a * self.a - 3 * self.b * self.b
        need(norm != 0, "division by zero in Q(sqrt(3))")
        return Q3(self.a / norm, -self.b / norm)

    def __truediv__(self, other: Any) -> "Q3":
        return self * Q3.coerce(other).inverse()

    def __pow__(self, exponent: int) -> "Q3":
        need(exponent >= 0, "negative powers are not used")
        result = Q3.coerce(1)
        base = self
        for _ in range(exponent):
            result *= base
        return result

    def text(self) -> str:
        if self.b == 0:
            return str(self.a)
        return f"{self.a}+({self.b})*sqrt(3)"


@dataclass(frozen=True)
class C3:
    """An exact element ``re + i im`` over ``Q(sqrt(3))``."""

    re: Q3 = Q3()
    im: Q3 = Q3()

    @staticmethod
    def coerce(value: Any) -> "C3":
        if isinstance(value, C3):
            return value
        return C3(Q3.coerce(value), Q3())

    def __add__(self, other: Any) -> "C3":
        rhs = C3.coerce(other)
        return C3(self.re + rhs.re, self.im + rhs.im)

    __radd__ = __add__

    def __neg__(self) -> "C3":
        return C3(-self.re, -self.im)

    def __sub__(self, other: Any) -> "C3":
        return self + (-C3.coerce(other))

    def __rsub__(self, other: Any) -> "C3":
        return C3.coerce(other) - self

    def __mul__(self, other: Any) -> "C3":
        rhs = C3.coerce(other)
        return C3(
            self.re * rhs.re - self.im * rhs.im,
            self.re * rhs.im + self.im * rhs.re,
        )

    __rmul__ = __mul__

    def conjugate(self) -> "C3":
        return C3(self.re, -self.im)

    def inverse(self) -> "C3":
        norm = self.re * self.re + self.im * self.im
        need(norm != Q3(), "division by zero in Q(sqrt(3), i)")
        return C3(self.re / norm, -self.im / norm)

    def __truediv__(self, other: Any) -> "C3":
        return self * C3.coerce(other).inverse()

    def text(self) -> str:
        if self.im == Q3():
            return self.re.text()
        return f"({self.re.text()})+i*({self.im.text()})"


ZERO = Q3()
ONE = Q3.coerce(1)
SQRT3 = Q3(Fraction(0), Fraction(1))
CZERO = C3()
CONE = C3.coerce(1)
CI = C3(Q3(), Q3.coerce(1))

QMatrix = tuple[tuple[Q3, Q3], tuple[Q3, Q3]]
CMatrix = tuple[tuple[C3, C3], tuple[C3, C3]]


def qmatrix(rows: Sequence[Sequence[Any]]) -> QMatrix:
    need(len(rows) == 2 and all(len(row) == 2 for row in rows), "expected 2x2 matrix")
    return tuple(tuple(Q3.coerce(x) for x in row) for row in rows)  # type: ignore[return-value]


def cmatrix(rows: Sequence[Sequence[Any]]) -> CMatrix:
    need(len(rows) == 2 and all(len(row) == 2 for row in rows), "expected 2x2 matrix")
    return tuple(tuple(C3.coerce(x) for x in row) for row in rows)  # type: ignore[return-value]


QZERO = qmatrix(((0, 0), (0, 0)))
QIDENTITY = qmatrix(((1, 0), (0, 1)))
QRECORD = qmatrix(((1, 0), (0, 0)))
QCOMPANION = qmatrix(((0, 0), (0, 1)))
CIDENTITY = cmatrix(((1, 0), (0, 1)))


def qadd(left: QMatrix, right: QMatrix) -> QMatrix:
    return qmatrix([[left[i][j] + right[i][j] for j in range(2)] for i in range(2)])


def qsub(left: QMatrix, right: QMatrix) -> QMatrix:
    return qmatrix([[left[i][j] - right[i][j] for j in range(2)] for i in range(2)])


def qmul(left: QMatrix, right: QMatrix) -> QMatrix:
    return qmatrix(
        [
            [sum((left[i][k] * right[k][j] for k in range(2)), ZERO) for j in range(2)]
            for i in range(2)
        ]
    )


def qtranspose(matrix: QMatrix) -> QMatrix:
    return qmatrix([[matrix[j][i] for j in range(2)] for i in range(2)])


def qtrace(matrix: QMatrix) -> Q3:
    return matrix[0][0] + matrix[1][1]


def cadd(left: CMatrix, right: CMatrix) -> CMatrix:
    return cmatrix([[left[i][j] + right[i][j] for j in range(2)] for i in range(2)])


def csub(left: CMatrix, right: CMatrix) -> CMatrix:
    return cmatrix([[left[i][j] - right[i][j] for j in range(2)] for i in range(2)])


def cmul(left: CMatrix, right: CMatrix) -> CMatrix:
    return cmatrix(
        [
            [sum((left[i][k] * right[k][j] for k in range(2)), CZERO) for j in range(2)]
            for i in range(2)
        ]
    )


def cscale(scalar: C3, matrix: CMatrix) -> CMatrix:
    return cmatrix([[scalar * matrix[i][j] for j in range(2)] for i in range(2)])


def cconj_transpose(matrix: CMatrix) -> CMatrix:
    return cmatrix([[matrix[j][i].conjugate() for j in range(2)] for i in range(2)])


def ctrace(matrix: CMatrix) -> C3:
    return matrix[0][0] + matrix[1][1]


def complexify(matrix: QMatrix) -> CMatrix:
    return cmatrix([[C3(matrix[i][j], Q3()) for j in range(2)] for i in range(2)])


def born_weight(state: CMatrix, effect: CMatrix) -> C3:
    return ctrace(cmul(state, effect))


def rank(rows: Iterable[Sequence[C3]]) -> int:
    work = [list(row) for row in rows]
    if not work:
        return 0
    width = len(work[0])
    need(all(len(row) == width for row in work), "ragged rank matrix")
    pivot_row = 0
    for column in range(width):
        pivot = next(
            (r for r in range(pivot_row, len(work)) if work[r][column] != CZERO), None
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        inverse = work[pivot_row][column].inverse()
        work[pivot_row] = [inverse * x for x in work[pivot_row]]
        for r in range(len(work)):
            if r == pivot_row or work[r][column] == CZERO:
                continue
            factor = work[r][column]
            work[r] = [work[r][c] - factor * work[pivot_row][c] for c in range(width)]
        pivot_row += 1
        if pivot_row == len(work):
            break
    return pivot_row


def flattened_pairing_row(effect: CMatrix) -> list[C3]:
    """Coefficients of ``Tr(rho E)`` in rho00,rho01,rho10,rho11."""

    return [effect[0][0], effect[1][0], effect[0][1], effect[1][1]]


def decode_q3(value: Any) -> Q3:
    need(isinstance(value, list) and len(value) == 2, "bad Q(sqrt(3)) scalar encoding")
    return Q3(Fraction(value[0]), Fraction(value[1]))


def decode_qmatrix(value: Any) -> QMatrix:
    need(isinstance(value, list) and len(value) == 2, "bad matrix row count")
    need(
        all(isinstance(row, list) and len(row) == 2 for row in value),
        "bad matrix column count",
    )
    return qmatrix([[decode_q3(value[i][j]) for j in range(2)] for i in range(2)])


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def compose_permutations(g: Sequence[int], h: Sequence[int]) -> tuple[int, int, int]:
    return tuple(g[h[i]] for i in range(3))  # type: ignore[return-value]


def phase_lift(
    commutator: QMatrix, *, include_i: bool = True, coefficient: Q3 | None = None
) -> CMatrix:
    """Construct ``I/2 - c i [Q,P]``; knobs exist for negative controls."""

    if coefficient is None:
        coefficient = Q3(Fraction(0), Fraction(2, 3))  # 2 sqrt(3) / 3
    phase = CI if include_i else CONE
    return csub(
        cscale(C3.coerce(Fraction(1, 2)), CIDENTITY),
        cscale(C3(coefficient) * phase, complexify(commutator)),
    )


def verify_payload(
    payload: dict[str, Any],
    source_root: Path | None = None,
    *,
    verify_hashes: bool = False,
    payload_path: Path | None = None,
) -> dict[str, Any]:
    need(
        payload.get("schema") == "oph.sim.born_context_web_payload.v1",
        "payload schema drift",
    )
    provenance = payload.get("provenance", {})
    need(provenance.get("run_id") == "b12_prereg_16k_20260806", "run id drift")
    need(
        provenance.get("run_git_commit") == EXPECTED_RUN_COMMIT,
        "pinned source commit drift",
    )
    need(provenance.get("seed") == 20260806, "seed drift")
    need(
        provenance.get("extraction", {}).get("script")
        == "scripts/extract_born_context_web.py",
        "extractor provenance drift",
    )
    for hash_key, expected_hash in EXPECTED_SOURCE_HASHES.items():
        need(
            provenance.get(hash_key) == expected_hash,
            f"pinned provenance hash drift: {hash_key}",
        )

    hashes_verified = False
    payload_file_sha256: str | None = None
    if payload_path is not None:
        need(payload_path.is_file(), "missing pinned context-web payload")
        payload_file_sha256 = sha256(payload_path)
        need(
            payload_file_sha256 == EXPECTED_PAYLOAD_SHA256,
            "context-web payload hash mismatch",
        )
    if verify_hashes:
        need(source_root is not None, "source root required for hash verification")
        if payload_path is None:
            authenticated_payload = source_root / "docs" / "BORN_CONTEXT_WEB_PAYLOAD.json"
            need(authenticated_payload.is_file(), "missing pinned context-web payload")
            payload_file_sha256 = sha256(authenticated_payload)
            need(
                payload_file_sha256 == EXPECTED_PAYLOAD_SHA256,
                "context-web payload hash mismatch",
            )
        for path_key, hash_key in (
            ("receipt_path", "receipt_sha256"),
            ("freezeout_path", "freezeout_sha256"),
            ("gauge_state_path", "gauge_state_sha256"),
            ("e1_payload_path", "e1_payload_sha256"),
        ):
            path = source_root / provenance[path_key]
            need(path.is_file(), f"missing pinned source file: {path}")
            need(
                sha256(path) == provenance[hash_key],
                f"source hash mismatch: {path_key}",
            )
        hashes_verified = True

    irrep = payload.get("irrep", {})
    need(irrep.get("name") == "standard_two_dimensional", "irrep name drift")
    elements = irrep.get("elements")
    need(isinstance(elements, list) and len(elements) == 6, "expected six S3 elements")
    matrices: dict[int, QMatrix] = {}
    permutations: dict[int, tuple[int, int, int]] = {}
    permutation_index: dict[tuple[int, int, int], int] = {}
    for expected_index, row in enumerate(elements):
        need(row.get("index") == expected_index, "S3 element index drift")
        permutation = tuple(row.get("permutation", ()))
        need(sorted(permutation) == [0, 1, 2], "invalid S3 permutation")
        matrix = decode_qmatrix(row.get("matrix"))
        need(
            qmul(matrix, qtranspose(matrix)) == QIDENTITY,
            f"element {expected_index} not orthogonal",
        )
        matrices[expected_index] = matrix
        permutations[expected_index] = permutation  # type: ignore[assignment]
        permutation_index[permutation] = expected_index  # type: ignore[index]
    need(len(permutation_index) == 6, "duplicate S3 permutations")
    for g in range(6):
        for h in range(6):
            product_index = permutation_index[
                compose_permutations(permutations[g], permutations[h])
            ]
            need(
                qmul(matrices[g], matrices[h]) == matrices[product_index],
                f"irrep homomorphism failure at ({g},{h})",
            )

    web = payload.get("context_web", {})
    diagonal = web.get("diagonal_context", {})
    diagonal_projectors = diagonal.get("projectors")
    need(
        isinstance(diagonal_projectors, list) and len(diagonal_projectors) == 2,
        "bad diagonal context",
    )
    need(decode_qmatrix(diagonal_projectors[0]) == QRECORD, "record projector drift")
    need(
        decode_qmatrix(diagonal_projectors[1]) == QCOMPANION,
        "companion projector drift",
    )
    need(diagonal.get("realized_outcome_counts") == [111, 68], "diagonal counts drift")
    need(diagonal.get("block_restricted_mass") == 179, "diagonal mass drift")

    conjugated = web.get("conjugated_contexts")
    need(
        isinstance(conjugated, list) and len(conjugated) == 5,
        "expected five conjugated contexts",
    )
    native_effects: list[QMatrix] = [QRECORD, QCOMPANION]
    record_orbit: dict[int, QMatrix] = {0: QRECORD}
    noncommuting_indices: list[int] = []
    for expected_g, row in enumerate(conjugated, start=1):
        need(
            row.get("gauge_element_index") == expected_g,
            "conjugated context index drift",
        )
        projectors = row.get("projectors")
        need(
            isinstance(projectors, list) and len(projectors) == 2,
            "bad conjugated binary context",
        )
        record = decode_qmatrix(projectors[0])
        companion = decode_qmatrix(projectors[1])
        expected_record = qmul(
            qmul(matrices[expected_g], QRECORD), qtranspose(matrices[expected_g])
        )
        need(
            record == expected_record, f"conjugated record projector {expected_g} drift"
        )
        need(
            record == qtranspose(record),
            f"conjugated record projector {expected_g} not symmetric",
        )
        need(
            qmul(record, record) == record,
            f"conjugated record projector {expected_g} not idempotent",
        )
        need(
            qtrace(record) == ONE,
            f"conjugated record projector {expected_g} trace drift",
        )
        need(qadd(record, companion) == QIDENTITY, f"context {expected_g} not complete")
        need(
            qmul(companion, companion) == companion,
            f"companion projector {expected_g} not idempotent",
        )
        commutator = qsub(qmul(record, QRECORD), qmul(QRECORD, record))
        need(
            decode_qmatrix(row.get("commutator_with_record_projector")) == commutator,
            f"commutator {expected_g} drift",
        )
        noncommuting = commutator != QZERO
        need(
            row.get("noncommutation_certificate") is noncommuting,
            f"noncommutation flag {expected_g} drift",
        )
        if noncommuting:
            noncommuting_indices.append(expected_g)
        need(
            row.get("realized_incident_edge_count_on_observer_support", 0) > 0,
            "unearned context on observer support",
        )
        need(
            row.get("realized_outcome_counts") is None,
            "rotated outcome counts are not source-realized",
        )
        record_orbit[expected_g] = record
        native_effects.extend((record, companion))

    need(noncommuting_indices == [2, 3, 4, 5], "noncommuting orbit drift")
    need(record_orbit[2] == record_orbit[3], "2/3 coincidence drift")
    need(record_orbit[4] == record_orbit[5], "4/5 coincidence drift")
    need(record_orbit[0] == record_orbit[1], "0/1 coincidence drift")
    need(
        web.get("context_coincidence_classes")
        == [
            ["conjugated_2", "conjugated_3"],
            ["conjugated_4", "conjugated_5"],
            ["diagonal", "conjugated_1"],
        ],
        "coincidence-class declaration drift",
    )

    boundary = payload.get("outcome_frequency_boundary", {})
    need(
        boundary.get("contexts_with_realized_frequencies") == ["diagonal"],
        "realized-frequency boundary drift",
    )
    need(
        boundary.get("contexts_without_realized_frequencies")
        == [f"conjugated_{g}" for g in range(1, 6)],
        "missing-frequency boundary drift",
    )
    need(
        "no producer" in boundary.get("named_capability_gap", ""),
        "producer capability gap missing",
    )
    need(
        "Born predictions" in boundary.get("prohibited_fill", ""),
        "prohibited-fill boundary drift",
    )

    native_complex = [complexify(effect) for effect in native_effects]
    native_rank = rank(
        [flattened_pairing_row(CIDENTITY)]
        + [flattened_pairing_row(effect) for effect in native_complex]
    )
    need(native_rank == 3, "native real web should have exact operator-span rank three")

    q = record_orbit[3]
    commutator = qsub(qmul(q, QRECORD), qmul(QRECORD, q))
    expected_commutator = qmatrix(
        ((0, Q3(Fraction(0), Fraction(1, 4))), (Q3(Fraction(0), Fraction(-1, 4)), 0))
    )
    need(commutator == expected_commutator, "exact sqrt(3)/4 source commutator drift")
    lifted = phase_lift(commutator)
    rho_plus = cmatrix(
        (
            (Fraction(1, 2), C3(Q3(), Q3.coerce(Fraction(-1, 2)))),
            (C3(Q3(), Q3.coerce(Fraction(1, 2))), Fraction(1, 2)),
        )
    )
    rho_minus = cmatrix(
        (
            (Fraction(1, 2), C3(Q3(), Q3.coerce(Fraction(1, 2)))),
            (C3(Q3(), Q3.coerce(Fraction(-1, 2))), Fraction(1, 2)),
        )
    )
    need(lifted == rho_plus, "normalized phase lift is not exactly the +Y projector")
    for name, rho in (("rho_plus", rho_plus), ("rho_minus", rho_minus)):
        need(cconj_transpose(rho) == rho, f"{name} not Hermitian")
        need(cmul(rho, rho) == rho, f"{name} not idempotent")
        need(ctrace(rho) == CONE, f"{name} trace drift")
    need(rho_plus != rho_minus, "Pauli-Y controls collapsed")
    for effect in native_complex:
        need(
            born_weight(rho_plus, effect) == born_weight(rho_minus, effect),
            "native effect unexpectedly sees Pauli-Y sign",
        )
    need(born_weight(rho_plus, lifted) == CONE, "+Y phase-lift weight should be one")
    need(born_weight(rho_minus, lifted) == CZERO, "-Y phase-lift weight should be zero")

    lifted_rank = rank(
        [
            flattened_pairing_row(CIDENTITY),
            flattened_pairing_row(complexify(QRECORD)),
            flattened_pairing_row(complexify(q)),
            flattened_pairing_row(lifted),
        ]
    )
    need(
        lifted_rank == 4, "phase-lifted frame should have exact operator-span rank four"
    )

    payload_text = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return {
        "schema": "oph.verify.b13_source_phase_lift.v1",
        "payload_file_sha256": payload_file_sha256,
        "payload_expected_sha256": EXPECTED_PAYLOAD_SHA256,
        "payload_canonical_sha256": hashlib.sha256(payload_text).hexdigest(),
        "source_run_git_commit": EXPECTED_RUN_COMMIT,
        "source_hashes_verified": hashes_verified,
        "s3_homomorphism_pairs_checked": 36,
        "native_declared_outcome_effect_count": len(native_effects),
        "native_distinct_outcome_effect_count": len(set(native_effects)),
        "native_operator_span_rank": native_rank,
        "native_complex_tomography_complete": False,
        "pauli_y_native_indistinguishability_exact": True,
        "phase_lift_formula": "I/2 - (2*sqrt(3)/3)*i*(Q*P-P*Q)",
        "phase_lift_equals_plus_y_projector": True,
        "phase_lifted_operator_span_rank": lifted_rank,
        "phase_lifted_fixed_trace_tomography_complete": True,
        "pauli_y_phase_lift_weights": ["1", "0"],
        "rotated_outcome_receipt_present": False,
        "phase_lift_instrument_receipt_present": False,
        "claim_boundary": (
            "Exact operator-algebra closure only. The source-attached algebraic web remains phase-free "
            "and has no rotated or phase-lifted outcome producer; common-preparation instrument "
            "validation and operational noncontextual additivity remain open."
        ),
    }


def default_paths() -> tuple[Path, Path]:
    workspace = Path(__file__).resolve().parents[3]
    source_root = workspace / "oph-physics-sim"
    return VENDORED_PAYLOAD, source_root


def main() -> None:
    default_payload, default_source_root = default_paths()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--payload", type=Path, default=default_payload)
    parser.add_argument("--source-root", type=Path, default=default_source_root)
    source_group = parser.add_mutually_exclusive_group()
    source_group.add_argument(
        "--verify-source-hashes",
        action="store_true",
        help="also rehash the optional sibling simulator's large raw inputs",
    )
    source_group.add_argument(
        "--skip-source-hashes",
        action="store_true",
        help="deprecated compatibility flag; hermetic verification is the default",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = json.loads(args.payload.read_text())
    report = verify_payload(
        payload,
        args.source_root,
        verify_hashes=args.verify_source_hashes,
        payload_path=args.payload,
    )
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
