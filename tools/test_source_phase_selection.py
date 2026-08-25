#!/usr/bin/env python3
"""Run the source-phase selector and six load-bearing mutation controls."""

from __future__ import annotations

import argparse
import copy
import json
import sys
from pathlib import Path
from typing import Any


PLUS_Y = [
    [["1/2", "0", "0", "0"], ["0", "0", "-1/2", "0"]],
    [["0", "0", "1/2", "0"], ["1/2", "0", "0", "0"]],
]


def require(condition: bool, code: str) -> None:
    if not condition:
        raise AssertionError(code)


def sign(algebra: Any, matrix: Any) -> str:
    return "positive" if matrix[0][1].im.sign() < 0 else "negative"


def empty_family_after_full_derivation(algebra: Any, family: Any,
                                       payload: dict[str, Any]) -> bool:
    view = algebra.extract_source_view(payload)
    matrices, record = algebra.validate_source_view(view)
    orbit = [
        algebra.qmul(algebra.qmul(matrix, record), algebra.qtranspose(matrix))
        for matrix in matrices
    ]
    for left in range(6):
        for right in range(left + 1, 6):
            commutator = algebra.qsub(
                algebra.qmul(orbit[right], orbit[left]),
                algebra.qmul(orbit[left], orbit[right]),
            )
            if commutator != algebra.QZERO:
                return False
    return True


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--payload", type=Path, required=True)
    parser.add_argument("--module-dir", type=Path, required=True)
    args = parser.parse_args()

    sys.path.insert(0, str(args.module_dir.resolve()))
    import source_phase_selection_algebra as algebra
    import source_phase_selection_family as family

    require(algebra.PACKET_SCHEMA != family.PACKET_SCHEMA,
            "LEGACY_AND_CURRENT_SCHEMA_COLLISION")
    require("legacy_totalized_stutter" in algebra.PACKET_SCHEMA,
            "LEGACY_SCHEMA_NOT_EXPLICIT")

    payload = json.loads(args.payload.read_text(encoding="utf-8"))
    require(algebra.sha256(args.payload) == algebra.EXPECTED_PAYLOAD_SHA256,
            "PAYLOAD_SHA256_MISMATCH")

    derived = family.derive_source(payload)
    events = sorted(derived.events, key=lambda row: (row[1], row[2]))
    packet = family.build_packet(payload, algebra.EXPECTED_PAYLOAD_SHA256)
    committed_packet_path = (
        args.module_dir / "SOURCE_PHASE_SELECTION_PACKET.v1.json"
    )
    require(committed_packet_path.is_file(), "COMMITTED_PACKET_MISSING")
    committed_packet = json.loads(
        committed_packet_path.read_text(encoding="utf-8")
    )
    require(committed_packet == packet, "COMMITTED_PACKET_MISMATCH")
    semantics_receipt_path = (
        args.module_dir / "SOURCE_PHASE_SELECTION_SEMANTICS_RECEIPT.v1.json"
    )
    require(semantics_receipt_path.is_file(), "SEMANTICS_RECEIPT_MISSING")
    semantics_receipt = json.loads(
        semantics_receipt_path.read_text(encoding="utf-8")
    )
    require(
        semantics_receipt.get("current_first_labelled_event") == "pair:00-02",
        "CURRENT_FIRST_LABELLED_EVENT_METADATA_MISMATCH",
    )
    require(
        semantics_receipt.get("orientation_status")
        == "declared_label_convention_not_source_forced",
        "ORIENTATION_STATUS_METADATA_MISMATCH",
    )
    require("selected_event" not in semantics_receipt,
            "LEGACY_SELECTED_EVENT_METADATA_PRESENT")
    report = family.verify_source_packet(packet, payload)
    selected = events[0]
    selected_bytes = algebra.encode_cmatrix(selected[5])
    require(selected[0] == "pair:00-02", "SELECTED_PAIR_MISMATCH")
    require(selected_bytes == PLUS_Y, "SELECTED_EFFECT_MISMATCH")
    reversed_first = events[-1]
    require(reversed_first[0] == "pair:03-05", "REVERSED_FIRST_PAIR_MISMATCH")
    require(algebra.encode_cmatrix(reversed_first[5]) == PLUS_Y,
            "EVENT_LIST_REVERSAL_ORIENTATION_REGRESSION")

    distinct = {algebra.matrix_key(row[5]): row[5] for row in events}
    positive = sum(sign(algebra, row[5]) == "positive" for row in events)
    negative = sum(sign(algebra, row[5]) == "negative" for row in events)
    require(len(distinct) == 2, "EFFECT_VALUE_COUNT_MISMATCH")
    require((positive, negative) == (8, 4), "EVENT_MULTIPLICITY_MISMATCH")
    require(report["state_event_cell_count"] == 60, "CELL_COUNT_MISMATCH")
    require(report["enabled_cell_count"] == 48, "ENABLED_COUNT_MISMATCH")
    require(report["disabled_cell_count"] == 12, "DISABLED_COUNT_MISMATCH")
    require(report["enabled_state_change_count"] == 36,
            "CHANGING_BRANCH_COUNT_MISMATCH")
    require(report["enabled_already_target_count"] == 12,
            "EQUALITY_BRANCH_COUNT_MISMATCH")

    runs = []

    for name, mutant, expected in (
        ("transition_target_changed", family.hidden_transition_mutant(packet),
         "SOURCE_NEXT_STATE_MISMATCH"),
        ("transition_removed", family.missing_transition_mutant(packet),
         "SOURCE_TRANSITION_SET_MISMATCH"),
    ):
        try:
            family.verify_source_packet(mutant, payload)
        except family.SourceCheckError as error:
            actual = str(error)
        else:
            actual = "SURVIVED"
        require(actual == expected, f"{name}:{actual}")
        runs.append({"case": name, "actual": actual, "expected": expected})

    zero_record = copy.deepcopy(payload)
    zero_record["context_web"]["diagonal_context"]["projectors"][0] = [
        [["0", "0"], ["0", "0"]],
        [["0", "0"], ["0", "0"]],
    ]
    actual = "SOURCE_PHASE_EFFECT_FAMILY_EMPTY" if empty_family_after_full_derivation(
        algebra, family, zero_record
    ) else "SURVIVED"
    require(actual == "SOURCE_PHASE_EFFECT_FAMILY_EMPTY", actual)
    runs.append({"case": "record_projector_zero", "actual": actual,
                 "expected": "SOURCE_PHASE_EFFECT_FAMILY_EMPTY"})

    trivial = copy.deepcopy(payload)
    identity = [[['1', '0'], ['0', '0']], [['0', '0'], ['1', '0']]]
    for row in trivial["irrep"]["elements"]:
        row["matrix"] = copy.deepcopy(identity)
    actual = "SOURCE_PHASE_EFFECT_FAMILY_EMPTY" if empty_family_after_full_derivation(
        algebra, family, trivial
    ) else "SURVIVED"
    require(actual == "SOURCE_PHASE_EFFECT_FAMILY_EMPTY", actual)
    runs.append({"case": "trivial_representation", "actual": actual,
                 "expected": "SOURCE_PHASE_EFFECT_FAMILY_EMPTY"})

    commutator = selected[3]
    reversed_commutator = algebra.qsub(algebra.QZERO, commutator)
    reversed_effect, _ = algebra.normalized_phase_operation(reversed_commutator)
    actual = "PR04_EFFECT_EXT_EQ_MISMATCH" if (
        algebra.encode_cmatrix(reversed_effect) != PLUS_Y
    ) else "SURVIVED"
    require(actual == "PR04_EFFECT_EXT_EQ_MISMATCH", actual)
    runs.append({"case": "reversed_orientation", "actual": actual,
                 "expected": "PR04_EFFECT_EXT_EQ_MISMATCH"})

    positive_only = [row for row in events if sign(algebra, row[5]) == "positive"]
    actual = "SOURCE_EFFECT_UNIVERSE_MISMATCH" if len({
        algebra.matrix_key(row[5]) for row in positive_only
    }) != 2 else "SURVIVED"
    require(actual == "SOURCE_EFFECT_UNIVERSE_MISMATCH", actual)
    runs.append({"case": "negative_effect_value_removed", "actual": actual,
                 "expected": "SOURCE_EFFECT_UNIVERSE_MISMATCH"})

    print(json.dumps({
        "verdict": "PASS",
        "selected_pair": selected[0],
        "effect_values": len(distinct),
        "event_multiplicity": {"positive": positive, "negative": negative},
        "cells": 60,
        "enabled_disabled": [48, 12],
        "accepted_equality": [36, 12],
        "mutation_count": len(runs),
        "mutations": runs,
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
