#!/usr/bin/env python3
"""Exact algebra utilities for target-free source-phase enumeration.

The only accepted semantic input is the ``irrep`` table and the outcome-0
record projector in ``BORN_CONTEXT_WEB_PAYLOAD.v1.json``.  Counts, PR rows,
the historical declared phase matrix, receipt verdicts, and paper claims are
never read.  The current packet producer lives in
``source_phase_selection_family.py`` and:

1. reconstructs the complete S3 projector orbit;
2. enumerates every stable-index noncommuting pair;
3. normalizes its commutator using the absolute off-diagonal magnitude computed
   from the permitted declared irrep/projector inputs (no target coefficient);
4. obtains a rank-one phase projection for every pair; and
5. enumerates all state/event cells, treating zero-weight cells as disabled
   and creating no transition; enabled cells collapse to the event projector.

The OPH bridge is evaluated only on enabled ``SourceStep`` entries.  This
module retains a quarantined legacy totalized-stutter packet builder for
historical replay, under a distinct legacy schema.  Its command-line entry
delegates to the current enabled-domain producer and cannot write the legacy
object as v1.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable, Sequence


EXPECTED_PAYLOAD_SHA256 = (
    "71a06f1c15192123cd09feb2386da702b572c8ac57c9b7633f5aa60c5d404e22"
)
LEGACY_PACKET_SCHEMA = (
    "oph.eventalgebra.source_phase_selection.legacy_totalized_stutter.v0"
)
PACKET_SCHEMA = LEGACY_PACKET_SCHEMA
SOURCE_VIEW_SCHEMA = "oph.eventalgebra.phase_source_view.v1"


class GateError(ValueError):
    pass


def need(condition: bool, message: str) -> None:
    if not condition:
        raise GateError(message)


@dataclass(frozen=True)
class Q3:
    """Exact ``a + b*sqrt(3)``."""

    a: Fraction = Fraction(0)
    b: Fraction = Fraction(0)

    @staticmethod
    def of(value: Any) -> "Q3":
        if isinstance(value, Q3):
            return value
        return Q3(Fraction(value), Fraction(0))

    def __add__(self, other: Any) -> "Q3":
        rhs = Q3.of(other)
        return Q3(self.a + rhs.a, self.b + rhs.b)

    __radd__ = __add__

    def __neg__(self) -> "Q3":
        return Q3(-self.a, -self.b)

    def __sub__(self, other: Any) -> "Q3":
        return self + (-Q3.of(other))

    def __rsub__(self, other: Any) -> "Q3":
        return Q3.of(other) - self

    def __mul__(self, other: Any) -> "Q3":
        rhs = Q3.of(other)
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
        return self * Q3.of(other).inverse()

    def sign(self) -> int:
        """Exact sign without floating point."""

        if self.a == 0:
            return (self.b > 0) - (self.b < 0)
        if self.b == 0:
            return (self.a > 0) - (self.a < 0)
        if (self.a > 0 and self.b > 0) or (self.a < 0 and self.b < 0):
            return (self.a > 0) - (self.a < 0)
        cmp = self.a * self.a - 3 * self.b * self.b
        need(cmp != 0, "unexpected rational equality with sqrt(3)")
        if self.a > 0:  # a > 0, b < 0
            return 1 if cmp > 0 else -1
        return 1 if cmp < 0 else -1  # a < 0, b > 0

    def abs(self) -> "Q3":
        return self if self.sign() >= 0 else -self

    def encode(self) -> list[str]:
        return [str(self.a), str(self.b)]


@dataclass(frozen=True)
class C3:
    re: Q3 = Q3()
    im: Q3 = Q3()

    @staticmethod
    def of(value: Any) -> "C3":
        if isinstance(value, C3):
            return value
        return C3(Q3.of(value), Q3())

    def __add__(self, other: Any) -> "C3":
        rhs = C3.of(other)
        return C3(self.re + rhs.re, self.im + rhs.im)

    __radd__ = __add__

    def __neg__(self) -> "C3":
        return C3(-self.re, -self.im)

    def __sub__(self, other: Any) -> "C3":
        return self + (-C3.of(other))

    def __rsub__(self, other: Any) -> "C3":
        return C3.of(other) - self

    def __mul__(self, other: Any) -> "C3":
        rhs = C3.of(other)
        return C3(
            self.re * rhs.re - self.im * rhs.im,
            self.re * rhs.im + self.im * rhs.re,
        )

    __rmul__ = __mul__

    def inverse(self) -> "C3":
        norm = self.re * self.re + self.im * self.im
        need(norm != Q3(), "division by zero in Q(sqrt(3), i)")
        return C3(self.re / norm, -self.im / norm)

    def __truediv__(self, other: Any) -> "C3":
        return self * C3.of(other).inverse()

    def conjugate(self) -> "C3":
        return C3(self.re, -self.im)

    def encode(self) -> list[str]:
        return self.re.encode() + self.im.encode()


QMatrix = tuple[tuple[Q3, Q3], tuple[Q3, Q3]]
CMatrix = tuple[tuple[C3, C3], tuple[C3, C3]]
QZERO: QMatrix = ((Q3(), Q3()), (Q3(), Q3()))
QIDENTITY: QMatrix = ((Q3.of(1), Q3()), (Q3(), Q3.of(1)))
CIDENTITY: CMatrix = ((C3.of(1), C3()), (C3(), C3.of(1)))
CI = C3(Q3(), Q3.of(1))


def qmatrix(rows: Sequence[Sequence[Any]]) -> QMatrix:
    need(len(rows) == 2 and all(len(row) == 2 for row in rows), "expected 2x2")
    return tuple(tuple(Q3.of(x) for x in row) for row in rows)  # type: ignore[return-value]


def cmatrix(rows: Sequence[Sequence[Any]]) -> CMatrix:
    need(len(rows) == 2 and all(len(row) == 2 for row in rows), "expected 2x2")
    return tuple(tuple(C3.of(x) for x in row) for row in rows)  # type: ignore[return-value]


def qadd(a: QMatrix, b: QMatrix) -> QMatrix:
    return qmatrix([[a[i][j] + b[i][j] for j in range(2)] for i in range(2)])


def qsub(a: QMatrix, b: QMatrix) -> QMatrix:
    return qmatrix([[a[i][j] - b[i][j] for j in range(2)] for i in range(2)])


def qmul(a: QMatrix, b: QMatrix) -> QMatrix:
    return qmatrix([[sum((a[i][k] * b[k][j] for k in range(2)), Q3())
                     for j in range(2)] for i in range(2)])


def qtranspose(a: QMatrix) -> QMatrix:
    return qmatrix([[a[j][i] for j in range(2)] for i in range(2)])


def cadd(a: CMatrix, b: CMatrix) -> CMatrix:
    return cmatrix([[a[i][j] + b[i][j] for j in range(2)] for i in range(2)])


def csub(a: CMatrix, b: CMatrix) -> CMatrix:
    return cmatrix([[a[i][j] - b[i][j] for j in range(2)] for i in range(2)])


def cmul(a: CMatrix, b: CMatrix) -> CMatrix:
    return cmatrix([[sum((a[i][k] * b[k][j] for k in range(2)), C3())
                     for j in range(2)] for i in range(2)])


def cscale(s: C3, a: CMatrix) -> CMatrix:
    return cmatrix([[s * a[i][j] for j in range(2)] for i in range(2)])


def cdagger(a: CMatrix) -> CMatrix:
    return cmatrix([[a[j][i].conjugate() for j in range(2)] for i in range(2)])


def ctrace(a: CMatrix) -> C3:
    return a[0][0] + a[1][1]


def complexify(a: QMatrix) -> CMatrix:
    return cmatrix([[C3(a[i][j], Q3()) for j in range(2)] for i in range(2)])


def decode_q3(value: Any) -> Q3:
    need(isinstance(value, list) and len(value) == 2, "bad Q3 scalar")
    return Q3(Fraction(value[0]), Fraction(value[1]))


def decode_qmatrix(value: Any) -> QMatrix:
    need(isinstance(value, list) and len(value) == 2, "bad matrix rows")
    need(all(isinstance(row, list) and len(row) == 2 for row in value), "bad matrix cols")
    return qmatrix([[decode_q3(value[i][j]) for j in range(2)] for i in range(2)])


def encode_cmatrix(a: CMatrix) -> list[list[list[str]]]:
    return [[a[i][j].encode() for j in range(2)] for i in range(2)]


def matrix_key(a: CMatrix) -> str:
    return json.dumps(encode_cmatrix(a), sort_keys=True, separators=(",", ":"))


def matrix_sha(a: CMatrix) -> str:
    return hashlib.sha256(matrix_key(a).encode()).hexdigest()


def canonical_json_bytes(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def compose(g: Sequence[int], h: Sequence[int]) -> tuple[int, int, int]:
    return tuple(g[h[i]] for i in range(3))  # type: ignore[return-value]


def extract_source_view(payload: dict[str, Any]) -> dict[str, Any]:
    """Project the payload onto the only permitted semantic inputs."""

    irrep = payload.get("irrep", {})
    elements = irrep.get("elements")
    diagonal = payload.get("context_web", {}).get("diagonal_context", {})
    projectors = diagonal.get("projectors")
    need(isinstance(elements, list) and len(elements) == 6, "expected six irrep rows")
    need(isinstance(projectors, list) and len(projectors) == 2, "missing record projector")
    return {
        "schema": SOURCE_VIEW_SCHEMA,
        "irrep_name": irrep.get("name"),
        "elements": [
            {
                "index": row.get("index"),
                "permutation": row.get("permutation"),
                "matrix": row.get("matrix"),
            }
            for row in elements
        ],
        "record_projector": projectors[0],
    }


def validate_source_view(view: dict[str, Any]) -> tuple[list[QMatrix], QMatrix]:
    need(view.get("schema") == SOURCE_VIEW_SCHEMA, "source-view schema drift")
    need(view.get("irrep_name") == "standard_two_dimensional", "irrep name drift")
    rows = view.get("elements")
    need(isinstance(rows, list) and len(rows) == 6, "expected six elements")
    matrices: list[QMatrix] = []
    perms: list[tuple[int, int, int]] = []
    pindex: dict[tuple[int, int, int], int] = {}
    for index, row in enumerate(rows):
        need(row.get("index") == index, "unstable irrep index")
        perm = tuple(row.get("permutation", ()))
        need(sorted(perm) == [0, 1, 2], "invalid permutation")
        matrix = decode_qmatrix(row.get("matrix"))
        need(qmul(matrix, qtranspose(matrix)) == QIDENTITY, "nonorthogonal irrep row")
        matrices.append(matrix)
        perms.append(perm)  # type: ignore[arg-type]
        pindex[perm] = index  # type: ignore[index]
    need(len(pindex) == 6, "duplicate permutation")
    for g in range(6):
        for h in range(6):
            need(qmul(matrices[g], matrices[h]) == matrices[pindex[compose(perms[g], perms[h])]],
                 f"homomorphism failure {g},{h}")
    record = decode_qmatrix(view.get("record_projector"))
    need(qmul(record, record) == record, "record projector not idempotent")
    need(qtranspose(record) == record, "record projector not symmetric")
    return matrices, record


def normalized_phase_operation(commutator: QMatrix) -> tuple[CMatrix, Q3]:
    """Normalize a nonzero 2x2 skew commutator without a target constant."""

    need(qtranspose(commutator) == qmatrix([[-commutator[i][j] for j in range(2)]
                                            for i in range(2)]),
         "commutator not skew")
    c01 = commutator[0][1]
    need(c01 != Q3(), "zero commutator")
    magnitude = c01.abs()
    # 1/2 I - i [Q,P] / (2 |[Q,P]01|).
    operation = csub(
        cscale(C3.of(Fraction(1, 2)), CIDENTITY),
        cscale(CI / C3.of(2 * magnitude), complexify(commutator)),
    )
    need(cdagger(operation) == operation, "generated operation not Hermitian")
    need(cmul(operation, operation) == operation, "generated operation not idempotent")
    need(ctrace(operation) == C3.of(1), "generated operation not rank-one")
    return operation, magnitude


def born_weight(state: CMatrix, event: CMatrix) -> C3:
    return ctrace(cmul(state, event))


def lueders_enabled_transition(state: CMatrix, event: CMatrix) -> tuple[bool, CMatrix]:
    weight = born_weight(state, event)
    if weight == C3():
        return False, state
    updated = cscale(weight.inverse(), cmul(cmul(event, state), event))
    need(updated == event, "rank-one Lueders collapse failed")
    return True, updated


def build_packet(payload: dict[str, Any], payload_sha256: str | None = None) -> dict[str, Any]:
    """Build the quarantined legacy totalized-stutter packet."""
    view = extract_source_view(payload)
    matrices, record = validate_source_view(view)
    orbit = [qmul(qmul(matrix, record), qtranspose(matrix)) for matrix in matrices]
    for projector in orbit:
        need(qmul(projector, projector) == projector, "orbit projector not idempotent")
        need(qtranspose(projector) == projector, "orbit projector not symmetric")

    event_rows: list[dict[str, Any]] = []
    phase_by_key: dict[str, CMatrix] = {}
    for left in range(6):
        for right in range(left + 1, 6):
            commutator = qsub(qmul(orbit[right], orbit[left]), qmul(orbit[left], orbit[right]))
            if commutator == QZERO:
                continue
            operation, magnitude = normalized_phase_operation(commutator)
            key = matrix_key(operation)
            phase_by_key[key] = operation
            event_rows.append({
                "id": f"pair:{left:02d}-{right:02d}",
                "left_orbit_index": left,
                "right_orbit_index": right,
                "commutator_01": commutator[0][1].encode(),
                "normalization_magnitude": magnitude.encode(),
                "operation_matrix": encode_cmatrix(operation),
                "operation_matrix_sha256": matrix_sha(operation),
            })
    need(event_rows, "no noncommuting pair")

    state_by_key: dict[str, tuple[str, CMatrix, str]] = {}
    for index, projector in enumerate(orbit):
        matrix = complexify(projector)
        key = matrix_key(matrix)
        if key not in state_by_key:
            state_by_key[key] = (f"orbit:{index:02d}", matrix, "orbit_projector")
    for key in sorted(phase_by_key):
        matrix = phase_by_key[key]
        sign = "positive" if matrix[0][1].im.sign() < 0 else "negative"
        state_by_key.setdefault(key, (f"phase:{sign}", matrix, "normalized_commutator"))

    states = sorted(state_by_key.values(), key=lambda row: row[0])
    state_id = {matrix_key(matrix): sid for sid, matrix, _ in states}
    for event in event_rows:
        operation = cmatrix([[C3(Q3(Fraction(cell[0]), Fraction(cell[1])),
                                  Q3(Fraction(cell[2]), Fraction(cell[3])))
                              for cell in row] for row in event["operation_matrix"]])
        event["target_state_id"] = state_id[matrix_key(operation)]

    transitions: list[dict[str, Any]] = []
    for sid, state, _kind in states:
        for event in event_rows:
            operation = next(matrix for psid, matrix, _ in states
                             if psid == event["target_state_id"])
            enabled, nxt = lueders_enabled_transition(state, operation)
            transitions.append({
                "id": f"{sid}|{event['id']}",
                "source_state_id": sid,
                "event_id": event["id"],
                "enabled": enabled,
                "born_weight": born_weight(state, operation).encode(),
                "next_state_id": state_id[matrix_key(nxt)],
            })

    source_bytes = canonical_json_bytes(view)
    packet = {
        "schema": PACKET_SCHEMA,
        "source_binding": {
            "payload_sha256": payload_sha256,
            "source_view_sha256": hashlib.sha256(source_bytes).hexdigest(),
            "permitted_input_fields": ["irrep.elements", "context_web.diagonal_context.projectors[0]"],
            "excluded_input_classes": ["counts", "PR rows", "receipt verdicts", "declared phase target"],
        },
        "states": [
            {"id": sid, "kind": kind, "matrix": encode_cmatrix(matrix),
             "matrix_sha256": matrix_sha(matrix)} for sid, matrix, kind in states
        ],
        "events": event_rows,
        "transitions": transitions,
        "completeness": {
            "state_count": len(states),
            "event_count": len(event_rows),
            "transition_count": len(transitions),
            "expected_transition_count": len(states) * len(event_rows),
        },
    }
    return packet


def bridge_classification(packet: dict[str, Any]) -> dict[str, int]:
    """Evaluate the event-indexed canonical OPH carrier, never a row verdict."""

    targets = {event["id"]: event["target_state_id"] for event in packet["events"]}
    counts = {"accepted": 0, "equality_stutter": 0, "unmapped": 0}
    for row in packet["transitions"]:
        source = row["source_state_id"]
        target = targets[row["event_id"]]
        # Unique one-patch localRepair for carrier C_event.
        repaired = target if source != target else source
        nxt = row["next_state_id"]
        if nxt == source:
            counts["equality_stutter"] += 1
        elif nxt == repaired and repaired != source:
            counts["accepted"] += 1
        else:
            counts["unmapped"] += 1
    return counts


def verify_packet(packet: dict[str, Any]) -> dict[str, Any]:
    """Verify only the quarantined legacy totalized-stutter schema."""
    need(packet.get("schema") == PACKET_SCHEMA, "packet schema drift")
    states = {row["id"] for row in packet.get("states", [])}
    events = {row["id"]: row for row in packet.get("events", [])}
    transitions = packet.get("transitions", [])
    need(len(transitions) == len(states) * len(events), "incomplete transition cartesian product")
    need(len({row["id"] for row in transitions}) == len(transitions), "duplicate transition")
    for row in transitions:
        need(row["source_state_id"] in states, "unknown source state")
        need(row["event_id"] in events, "unknown event")
        need(row["next_state_id"] in states, "unknown next state")
        target = events[row["event_id"]]["target_state_id"]
        expected = target if row["enabled"] else row["source_state_id"]
        need(row["next_state_id"] == expected, "SOURCE_TRANSITION_MISMATCH")
    classification = bridge_classification(packet)
    need(classification["accepted"] > 0, "accepted branch vacuous")
    need(classification["equality_stutter"] > 0, "stutter branch vacuous")
    need(classification["unmapped"] == 0, "CANONICAL_OPH_UNMAPPED")
    packet_bytes = canonical_json_bytes(packet)
    return {
        "packet_sha256": hashlib.sha256(packet_bytes).hexdigest(),
        "state_count": len(states),
        "event_count": len(events),
        "transition_count": len(transitions),
        **classification,
    }


def hidden_transition_mutant(packet: dict[str, Any]) -> dict[str, Any]:
    mutant = copy.deepcopy(packet)
    row = next(row for row in mutant["transitions"]
               if row["enabled"] and row["next_state_id"] != row["source_state_id"])
    row["next_state_id"] = row["source_state_id"]
    return mutant


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    # Import lazily to avoid a cycle: the current producer imports the exact
    # scalar and matrix utilities above.  The old CLI therefore has one
    # canonical enabled-domain behavior even though the legacy replay helpers
    # remain callable under their explicitly different schema.
    import source_phase_selection_family as family

    family.main()


if __name__ == "__main__":
    main()
