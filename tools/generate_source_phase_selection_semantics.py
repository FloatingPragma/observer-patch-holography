#!/usr/bin/env python3
"""Generate Lean phase-effect semantics from a committed source-selection packet.

The generator consumes source states, source event effect matrices, and source
event order.  It does not read a declared target effect, sourcePhaseLift, PR
rows, OPH classifications, counts as semantic premises, manuscript prose, or
an expected result class.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Any


PACKET_SCHEMA = "oph.eventalgebra.source_phase_selection.v1"
EXPECTED_PAYLOAD_SHA256 = (
    "71a06f1c15192123cd09feb2386da702b572c8ac57c9b7633f5aa60c5d404e22"
)


STATE_CTORS = {
    "orbit:00": "orbit00",
    "orbit:02": "orbit02",
    "orbit:04": "orbit04",
    "phase:negative": "phaseNegative",
    "phase:positive": "phasePositive",
}


def event_ctor(event_id: str) -> str:
    return "pair" + event_id.removeprefix("pair:").replace("-", "")


def lean_rat(value: str) -> str:
    q = Fraction(value)
    if q.denominator == 1:
        return str(q.numerator)
    return f"({q.numerator} / {q.denominator} : ℝ)"


def lean_q3(a: str, b: str) -> str:
    return f"(({lean_rat(a)} : ℝ) + ({lean_rat(b)} : ℝ) * OPH.QFT.sqrt3)"


def lean_c3(cell: list[str]) -> str:
    if len(cell) != 4:
        raise ValueError("bad complex scalar")
    re = lean_q3(cell[0], cell[1])
    im = lean_q3(cell[2], cell[3])
    return f"((({re}) : ℂ) + ((({im}) : ℂ) * Complex.I))"


def lean_matrix(matrix: list[list[list[str]]]) -> str:
    return (
        "!![" + lean_c3(matrix[0][0]) + ", " + lean_c3(matrix[0][1])
        + "; " + lean_c3(matrix[1][0]) + ", " + lean_c3(matrix[1][1]) + "]"
    )


def effect_sign(matrix: list[list[list[str]]]) -> str:
    upper = matrix[0][1]
    if Fraction(upper[3]) != 0:
        raise ValueError("unexpected sqrt3 imaginary term in effect sign")
    return "positive" if Fraction(upper[2]) < 0 else "negative"


def canonical_bytes(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--packet", type=Path, required=True)
    parser.add_argument("--lean-output", type=Path, required=True)
    parser.add_argument("--receipt-output", type=Path, required=True)
    args = parser.parse_args()

    packet_bytes = args.packet.read_bytes()
    packet_hash = hashlib.sha256(packet_bytes).hexdigest()
    packet = json.loads(packet_bytes.decode("utf-8"))
    if packet.get("schema") != PACKET_SCHEMA:
        raise ValueError("PACKET_SCHEMA_MISMATCH")
    if packet.get("source_binding", {}).get("payload_sha256") != EXPECTED_PAYLOAD_SHA256:
        raise ValueError("PAYLOAD_SHA256_MISMATCH")
    states = packet["states"]
    events = packet["events"]
    if [row["id"] for row in states] != list(STATE_CTORS):
        raise ValueError("STATE_ORDER_MISMATCH")
    if len(events) != 12:
        raise ValueError("EVENT_UNIVERSE_MISMATCH")

    effect_matrices: dict[str, list[list[list[str]]]] = {}
    event_signs: list[tuple[str, str]] = []
    for event in events:
        sign = effect_sign(event["operation_matrix"])
        event_signs.append((event["id"], sign))
        effect_matrices.setdefault(sign, event["operation_matrix"])
    if set(effect_matrices) != {"positive", "negative"}:
        raise ValueError("DISTINCT_EFFECT_VALUE_UNIVERSE_MISMATCH")

    lines = [
        "import QFT.SourcePhaseLiftBridge",
        "import EventAlgebra.SourcePhaseSelectionInput",
        "",
        "set_option autoImplicit false",
        "",
        "/-! Mechanically generated from the committed target-free source packet. -/",
        "namespace EventAlgebra.SourcePhaseSelectionSemantics",
        "",
        "open Matrix",
        "open EventAlgebra.SourcePhaseSelectionInput",
        "",
        "noncomputable section",
        "",
        "inductive EffectValue where",
        "  | positive",
        "  | negative",
        "  deriving DecidableEq, Fintype, Repr",
        "",
        "def stateMatrix : PhaseState → Matrix (Fin 2) (Fin 2) ℂ",
    ]
    for state in states:
        lines.append(
            f"  | .{STATE_CTORS[state['id']]} => {lean_matrix(state['matrix'])}"
        )
    lines.extend([
        "",
        "def effectMatrix : EffectValue → Matrix (Fin 2) (Fin 2) ℂ",
        f"  | .positive => {lean_matrix(effect_matrices['positive'])}",
        f"  | .negative => {lean_matrix(effect_matrices['negative'])}",
        "",
        "def effectValueOf : PhaseEvent → EffectValue",
    ])
    for event_id, sign in event_signs:
        lines.append(f"  | .{event_ctor(event_id)} => .{sign}")
    lines.extend([
        "",
        "def generatedEffect (event : PhaseEvent) : Matrix (Fin 2) (Fin 2) ℂ :=",
        "  effectMatrix (effectValueOf event)",
        "",
        "def sourcePairOrder : List PhaseEvent :=",
        "  [" + ", ".join(f".{event_ctor(row['id'])}" for row in events) + "]",
        "",
        "/-- The current first event under the declared stable label order.",
        "This name does not assert a source-forced orientation. -/",
        "def currentFirstLabelledEvent : PhaseEvent := sourcePairOrder[0]",
        "",
        "/-- The generated effect of the current first-labelled event. -/",
        "def currentFirstLabelledGeneratedEffect : Matrix (Fin 2) (Fin 2) ℂ :=",
        "  generatedEffect currentFirstLabelledEvent",
        "",
        "/-- Legacy compatibility alias. No source theorem forces this label",
        "or its `+Y` orientation. -/",
        "abbrev sourceSelectedEvent : PhaseEvent := currentFirstLabelledEvent",
        "",
        "/-- Legacy compatibility alias. Prefer",
        "`currentFirstLabelledGeneratedEffect`. -/",
        "abbrev sourceSelectedGeneratedEffect : Matrix (Fin 2) (Fin 2) ℂ :=",
        "  currentFirstLabelledGeneratedEffect",
        "",
        "def generatedEffectValues : List EffectValue :=",
        "  (sourcePairOrder.map effectValueOf).eraseDups",
        "",
        "def upperImaginaryNegative (value : EffectValue) : Prop :=",
        "  Complex.im (effectMatrix value 0 1) < 0",
        "",
        "def sourcePayloadSha256 : String :=",
        f"  \"{packet['source_binding']['payload_sha256']}\"",
        "def sourceSelectionPacketSha256 : String :=",
        f"  \"{packet_hash}\"",
        "",
        "#check stateMatrix",
        "#check effectMatrix",
        "#check generatedEffect",
        "#check currentFirstLabelledGeneratedEffect",
        "#check sourceSelectedGeneratedEffect",
        "",
        "end",
        "",
        "end EventAlgebra.SourcePhaseSelectionSemantics",
        "",
    ])
    lean_text = "\n".join(lines)
    args.lean_output.write_text(lean_text, encoding="utf-8", newline="\n")

    receipt = {
        "schema": "oph.eventalgebra.source_phase_selection_semantics.v1",
        "source_selection_packet_sha256": packet_hash,
        "semantic_inputs": [
            "states[*].id",
            "states[*].matrix",
            "events[*].id",
            "events[*].left_orbit_index",
            "events[*].right_orbit_index",
            "events[*].operation_matrix",
        ],
        "forbidden_inputs": {
            "declared_target_effect": 0,
            "sourcePhaseLift_lookup": 0,
            "premise_register_or_ledger": 0,
            "accepted_or_stutter_classification": 0,
            "PASS_or_result_class": 0,
            "manuscript_claim": 0,
            "completeness_count_as_semantic_input": 0,
        },
        "current_first_labelled_event": events[0]["id"],
        "orientation_status": "declared_label_convention_not_source_forced",
        "event_count": len(events),
        "distinct_effect_value_count": len(effect_matrices),
        "positive_event_witness_count": sum(sign == "positive" for _, sign in event_signs),
        "negative_event_witness_count": sum(sign == "negative" for _, sign in event_signs),
        "lean_output": {
            "bytes": args.lean_output.stat().st_size,
            "sha256": hashlib.sha256(args.lean_output.read_bytes()).hexdigest(),
        },
    }
    args.receipt_output.write_bytes(canonical_bytes(receipt))
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
