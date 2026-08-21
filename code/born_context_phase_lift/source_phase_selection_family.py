#!/usr/bin/env python3
"""Enabled-domain EventAlgebra source phase-effect selection family.

This version follows the current Lean semantics exactly: zero Born weight is
outside the enabled rank-one Lueders transition relation.  It is represented
as a disabled state/event cell, not invented as an equality no-op.  The
packet contains every state/event cell and only the 48 enabled transitions.

The semantic input projection remains target-free: six exact irrep rows and
the diagonal record projector.  Counts, premise rows, receipts, manuscript
claims, and the historical declared phase target are not read.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any

import source_phase_selection_algebra as algebra


PACKET_SCHEMA = "oph.eventalgebra.source_phase_selection.v1"
EXPECTED_PAYLOAD_SHA256 = algebra.EXPECTED_PAYLOAD_SHA256


class SourceCheckError(ValueError):
    pass


def need(condition: bool, message: str) -> None:
    if not condition:
        raise SourceCheckError(message)


@dataclass(frozen=True)
class DerivedSource:
    view: dict[str, Any]
    states: tuple[tuple[str, str, algebra.CMatrix], ...]
    events: tuple[tuple[str, int, int, algebra.QMatrix, algebra.Q3, algebra.CMatrix], ...]


def decode_cmatrix(value: Any) -> algebra.CMatrix:
    need(isinstance(value, list) and len(value) == 2, "bad complex matrix rows")
    rows = []
    for row in value:
        need(isinstance(row, list) and len(row) == 2, "bad complex matrix columns")
        cells = []
        for cell in row:
            need(isinstance(cell, list) and len(cell) == 4, "bad complex scalar")
            cells.append(
                algebra.C3(
                    algebra.Q3(Fraction(cell[0]), Fraction(cell[1])),
                    algebra.Q3(Fraction(cell[2]), Fraction(cell[3])),
                )
            )
        rows.append(cells)
    return algebra.cmatrix(rows)


def derive_source(payload: dict[str, Any]) -> DerivedSource:
    view = algebra.extract_source_view(payload)
    matrices, record = algebra.validate_source_view(view)
    orbit = [algebra.qmul(algebra.qmul(matrix, record), algebra.qtranspose(matrix))
             for matrix in matrices]

    events = []
    phase_by_key: dict[str, algebra.CMatrix] = {}
    for left in range(6):
        for right in range(left + 1, 6):
            commutator = algebra.qsub(
                algebra.qmul(orbit[right], orbit[left]),
                algebra.qmul(orbit[left], orbit[right]),
            )
            if commutator == algebra.QZERO:
                continue
            operation, magnitude = algebra.normalized_phase_operation(commutator)
            event_id = f"pair:{left:02d}-{right:02d}"
            events.append((event_id, left, right, commutator, magnitude, operation))
            phase_by_key[algebra.matrix_key(operation)] = operation
    need(len(events) == 12, "EVENT_UNIVERSE_MISMATCH")

    states_by_key: dict[str, tuple[str, str, algebra.CMatrix]] = {}
    for index, projector in enumerate(orbit):
        matrix = algebra.complexify(projector)
        states_by_key.setdefault(
            algebra.matrix_key(matrix),
            (f"orbit:{index:02d}", "orbit_projector", matrix),
        )
    for key in sorted(phase_by_key):
        matrix = phase_by_key[key]
        orientation = "positive" if matrix[0][1].im.sign() < 0 else "negative"
        states_by_key.setdefault(
            key,
            (f"phase:{orientation}", "normalized_commutator", matrix),
        )
    states = tuple(sorted(states_by_key.values(), key=lambda row: row[0]))
    need(len(states) == 5, "STATE_UNIVERSE_MISMATCH")
    return DerivedSource(view, states, tuple(events))


def source_semantics(derived: DerivedSource) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    state_id = {algebra.matrix_key(matrix): sid for sid, _kind, matrix in derived.states}
    cells: list[dict[str, Any]] = []
    transitions: list[dict[str, Any]] = []
    for sid, _kind, state in derived.states:
        for event_id, _left, _right, _commutator, _magnitude, operation in derived.events:
            weight = algebra.born_weight(state, operation)
            enabled = weight != algebra.C3()
            cell = {
                "id": f"{sid}|{event_id}",
                "source_state_id": sid,
                "event_id": event_id,
                "enabled": enabled,
                "born_weight": weight.encode(),
            }
            cells.append(cell)
            if not enabled:
                continue
            updated = algebra.cscale(weight.inverse(), algebra.cmul(algebra.cmul(operation, state), operation))
            need(updated == operation, "RANK_ONE_LUEDERS_COLLAPSE_MISMATCH")
            transitions.append({
                "id": cell["id"],
                "source_state_id": sid,
                "event_id": event_id,
                "born_weight": weight.encode(),
                "next_state_id": state_id[algebra.matrix_key(updated)],
            })
    return cells, transitions


def build_packet(payload: dict[str, Any], payload_sha256: str | None = None) -> dict[str, Any]:
    derived = derive_source(payload)
    cells, transitions = source_semantics(derived)
    enabled = sum(1 for row in cells if row["enabled"])
    disabled = len(cells) - enabled
    source_view_bytes = algebra.canonical_json_bytes(derived.view)
    packet = {
        "schema": PACKET_SCHEMA,
        "source_binding": {
            "payload_sha256": payload_sha256,
            "source_view_sha256": hashlib.sha256(source_view_bytes).hexdigest(),
            "permitted_input_fields": [
                "irrep.elements",
                "context_web.diagonal_context.projectors[0]",
            ],
        },
        "states": [
            {
                "id": sid,
                "kind": kind,
                "matrix": algebra.encode_cmatrix(matrix),
                "matrix_sha256": algebra.matrix_sha(matrix),
            }
            for sid, kind, matrix in derived.states
        ],
        "events": [
            {
                "id": event_id,
                "left_orbit_index": left,
                "right_orbit_index": right,
                "commutator_01": commutator[0][1].encode(),
                "normalization_magnitude": magnitude.encode(),
                "operation_matrix": algebra.encode_cmatrix(operation),
                "operation_matrix_sha256": algebra.matrix_sha(operation),
            }
            for event_id, left, right, commutator, magnitude, operation in derived.events
        ],
        "enabledness_cells": cells,
        "transitions": transitions,
        "completeness": {
            "state_count": len(derived.states),
            "event_count": len(derived.events),
            "state_event_cell_count": len(cells),
            "enabled_cell_count": enabled,
            "disabled_cell_count": disabled,
            "transition_count": len(transitions),
        },
    }
    return packet


def verify_source_packet(packet: dict[str, Any], payload: dict[str, Any]) -> dict[str, Any]:
    """Recompute source derivation and enabled transitions, with no OPH logic."""

    need(packet.get("schema") == PACKET_SCHEMA, "PACKET_SCHEMA_MISMATCH")
    derived = derive_source(payload)
    expected_cells, expected_transitions = source_semantics(derived)

    expected_states = {
        sid: (kind, algebra.matrix_key(matrix), algebra.matrix_sha(matrix))
        for sid, kind, matrix in derived.states
    }
    actual_states = packet.get("states")
    need(isinstance(actual_states, list), "STATE_TABLE_MISSING")
    need(len(actual_states) == len(expected_states), "STATE_TABLE_SIZE_MISMATCH")
    for row in actual_states:
        sid = row.get("id")
        need(sid in expected_states, "STATE_ID_MISMATCH")
        kind, matrix_key, matrix_hash = expected_states[sid]
        need(row.get("kind") == kind, "STATE_KIND_MISMATCH")
        matrix = decode_cmatrix(row.get("matrix"))
        need(algebra.matrix_key(matrix) == matrix_key, "STATE_MATRIX_MISMATCH")
        need(row.get("matrix_sha256") == matrix_hash, "STATE_HASH_MISMATCH")

    expected_events = {
        event_id: (left, right, commutator, magnitude, operation)
        for event_id, left, right, commutator, magnitude, operation in derived.events
    }
    actual_events = packet.get("events")
    need(isinstance(actual_events, list), "EVENT_TABLE_MISSING")
    need(len(actual_events) == len(expected_events), "EVENT_TABLE_SIZE_MISMATCH")
    for row in actual_events:
        event_id = row.get("id")
        need(event_id in expected_events, "EVENT_ID_MISMATCH")
        left, right, commutator, magnitude, operation = expected_events[event_id]
        need(row.get("left_orbit_index") == left, "EVENT_LEFT_INDEX_MISMATCH")
        need(row.get("right_orbit_index") == right, "EVENT_RIGHT_INDEX_MISMATCH")
        need(row.get("commutator_01") == commutator[0][1].encode(), "EVENT_COMMUTATOR_MISMATCH")
        need(row.get("normalization_magnitude") == magnitude.encode(), "EVENT_NORMALIZATION_MISMATCH")
        need(decode_cmatrix(row.get("operation_matrix")) == operation, "EVENT_OPERATION_MISMATCH")
        need(row.get("operation_matrix_sha256") == algebra.matrix_sha(operation), "EVENT_OPERATION_HASH_MISMATCH")

    actual_cells = packet.get("enabledness_cells")
    need(isinstance(actual_cells, list), "ENABLEDNESS_TABLE_MISSING")
    need(actual_cells == expected_cells, "SOURCE_ENABLEDNESS_MISMATCH")

    actual_transitions = packet.get("transitions")
    need(isinstance(actual_transitions, list), "TRANSITION_TABLE_MISSING")
    expected_ids = {row["id"] for row in expected_transitions}
    actual_ids = {row.get("id") for row in actual_transitions}
    need(actual_ids == expected_ids, "SOURCE_TRANSITION_SET_MISMATCH")
    expected_by_id = {row["id"]: row for row in expected_transitions}
    for row in actual_transitions:
        expected = expected_by_id[row["id"]]
        need(row.get("source_state_id") == expected["source_state_id"], "SOURCE_STATE_ID_MISMATCH")
        need(row.get("event_id") == expected["event_id"], "SOURCE_EVENT_ID_MISMATCH")
        need(row.get("born_weight") == expected["born_weight"], "SOURCE_WEIGHT_MISMATCH")
        need(row.get("next_state_id") == expected["next_state_id"], "SOURCE_NEXT_STATE_MISMATCH")

    enabled_changes = sum(
        1 for row in expected_transitions
        if row["source_state_id"] != row["next_state_id"]
    )
    enabled_self = len(expected_transitions) - enabled_changes
    packet_hash = hashlib.sha256(algebra.canonical_json_bytes(packet)).hexdigest()
    return {
        "packet_sha256": packet_hash,
        "state_count": len(derived.states),
        "event_count": len(derived.events),
        "state_event_cell_count": len(expected_cells),
        "enabled_cell_count": len(expected_transitions),
        "disabled_cell_count": len(expected_cells) - len(expected_transitions),
        "transition_count": len(expected_transitions),
        "enabled_state_change_count": enabled_changes,
        "enabled_already_target_count": enabled_self,
    }


def hidden_transition_mutant(packet: dict[str, Any]) -> dict[str, Any]:
    mutant = copy.deepcopy(packet)
    row = next(
        row for row in mutant["transitions"]
        if row["source_state_id"] != row["next_state_id"]
    )
    row["next_state_id"] = row["source_state_id"]
    return mutant


def missing_transition_mutant(packet: dict[str, Any]) -> dict[str, Any]:
    mutant = copy.deepcopy(packet)
    mutant["transitions"].pop()
    return mutant


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--payload", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--verify-mutations", action="store_true")
    args = parser.parse_args()

    digest = sha256(args.payload)
    need(digest == EXPECTED_PAYLOAD_SHA256, "PAYLOAD_SHA256_MISMATCH")
    payload = json.loads(args.payload.read_text(encoding="utf-8"))
    packet = build_packet(payload, digest)
    report = verify_source_packet(packet, payload)
    if args.verify_mutations:
        for name, mutant, code in (
            ("hidden_transition_mutant", hidden_transition_mutant(packet), "SOURCE_NEXT_STATE_MISMATCH"),
            ("missing_transition_mutant", missing_transition_mutant(packet), "SOURCE_TRANSITION_SET_MISMATCH"),
        ):
            try:
                verify_source_packet(mutant, payload)
            except SourceCheckError as exc:
                need(str(exc) == code, f"{name}: unexpected failure {exc}")
                report[name] = f"REJECTED:{code}"
            else:
                raise SourceCheckError(f"{name}: mutation survived")
    if args.output:
        args.output.write_bytes(algebra.canonical_json_bytes(packet))
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
