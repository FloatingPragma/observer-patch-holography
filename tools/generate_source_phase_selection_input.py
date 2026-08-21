#!/usr/bin/env python3
"""Generate the Lean source-selection input from the committed source packet.

The packet must be generated from exactly two mathematical inputs in the
committed payload: ``irrep.elements`` and
``context_web.diagonal_context.projectors[0]``.  Target effects, register
rows, manuscript prose, expected classifications, and verdicts are not
inputs to this generator.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


PACKET_SCHEMA = "oph.eventalgebra.source_phase_selection.v1"
EXPECTED_PAYLOAD_SHA256 = (
    "71a06f1c15192123cd09feb2386da702b572c8ac57c9b7633f5aa60c5d404e22"
)
EXPECTED_SOURCE_VIEW_SHA256 = (
    "06970219d5df9fc0783cf7576020f9796945a44eb433b31918587dd81a09a6ac"
)


def canonical_bytes(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def state_ctor(source_id: str) -> str:
    return {
        "orbit:00": "orbit00",
        "orbit:02": "orbit02",
        "orbit:04": "orbit04",
        "phase:negative": "phaseNegative",
        "phase:positive": "phasePositive",
    }[source_id]


def event_ctor(event_id: str) -> str:
    return "pair" + event_id.removeprefix("pair:").replace("-", "")


def nat_list(value: bytes) -> str:
    return "[" + ", ".join(str(byte) for byte in value) + "]"


def require(condition: bool, code: str) -> None:
    if not condition:
        raise ValueError(code)


def render(packet: dict[str, Any], packet_hash: str) -> str:
    states = packet["states"]
    events = packet["events"]
    cells = packet["enabledness_cells"]

    state_cases = [
        f"  | .{state_ctor(row['id'])} => {nat_list(canonical_bytes(row['matrix']))}"
        for row in states
    ]
    event_cases = [
        f"  | .{event_ctor(row['id'])} => "
        f"{nat_list(canonical_bytes(row['operation_matrix']))}"
        for row in events
    ]

    status_cases: list[str] = []
    all_cells: list[str] = []
    enabled_cells: list[str] = []
    disabled_cells: list[str] = []
    for row in cells:
        state = state_ctor(row["source_state_id"])
        event = event_ctor(row["event_id"])
        status = "enabled" if row["enabled"] else "disabled"
        status_cases.append(f"  | .{state}, .{event} => .{status}")
        pair = f"(.{state}, .{event})"
        all_cells.append(pair)
        witness = f"⟨{pair}, by rfl⟩"
        (enabled_cells if row["enabled"] else disabled_cells).append(witness)

    state_ctors = "\n".join(f"  | {state_ctor(row['id'])}" for row in states)
    event_ctors = "\n".join(f"  | {event_ctor(row['id'])}" for row in events)

    return f'''import Mathlib

set_option autoImplicit false

/-!
# Generated source phase-selection input

This module is mechanically generated from the committed target-free source
packet.  It contains source states, source operations, the complete
state/event cell table, the enabled-domain split, and custody hashes only.
-/

namespace EventAlgebra.SourcePhaseSelectionInput

inductive PhaseState where
{state_ctors}
  deriving DecidableEq, Fintype, Repr

inductive PhaseEvent where
{event_ctors}
  deriving DecidableEq, Fintype, Repr

abbrev StateBytes := List Nat
abbrev OperationBytes := List Nat

def stateBytes : PhaseState → StateBytes
{chr(10).join(state_cases)}

def operationBytes : PhaseEvent → OperationBytes
{chr(10).join(event_cases)}

inductive CellStatus where
  | enabled
  | disabled
  deriving DecidableEq, Repr

abbrev SourceCell := PhaseState × PhaseEvent

def cellStatus : PhaseState → PhaseEvent → CellStatus
{chr(10).join(status_cases)}

def allSourceCells : List SourceCell :=
  [{', '.join(all_cells)}]

abbrev EnabledCell :=
  {{cell : SourceCell // cellStatus cell.1 cell.2 = .enabled}}

abbrev DisabledCell :=
  {{cell : SourceCell // cellStatus cell.1 cell.2 = .disabled}}

def enabledCells : List EnabledCell :=
  [{', '.join(enabled_cells)}]

def disabledCells : List DisabledCell :=
  [{', '.join(disabled_cells)}]

theorem sourceCellCensus : allSourceCells.length = 60 := by decide
theorem enabledCellCensus : enabledCells.length = 48 := by decide
theorem disabledCellCensus : disabledCells.length = 12 := by decide

def sourceBoundEnabledCell : EnabledCell := enabledCells[0]
def sourceBoundDisabledCell : DisabledCell := disabledCells[0]

def sourcePayloadSha256 : String := "{EXPECTED_PAYLOAD_SHA256}"
def sourceSelectionPacketSha256 : String := "{packet_hash}"
def sourceViewSha256 : String := "{EXPECTED_SOURCE_VIEW_SHA256}"

end EventAlgebra.SourcePhaseSelectionInput
'''


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--packet", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    packet_bytes = args.packet.read_bytes()
    packet = json.loads(packet_bytes)
    require(packet.get("schema") == PACKET_SCHEMA, "PACKET_SCHEMA_MISMATCH")
    binding = packet.get("source_binding", {})
    require(binding.get("payload_sha256") == EXPECTED_PAYLOAD_SHA256,
            "PAYLOAD_SHA256_MISMATCH")
    require(binding.get("source_view_sha256") == EXPECTED_SOURCE_VIEW_SHA256,
            "SOURCE_VIEW_SHA256_MISMATCH")
    require(binding.get("permitted_input_fields") == [
        "irrep.elements", "context_web.diagonal_context.projectors[0]"
    ], "PERMITTED_FIELDS_MISMATCH")
    require(len(packet.get("states", [])) == 5, "STATE_COUNT_MISMATCH")
    require(len(packet.get("events", [])) == 12, "EVENT_COUNT_MISMATCH")
    cells = packet.get("enabledness_cells", [])
    require(len(cells) == 60, "CELL_COUNT_MISMATCH")
    require(sum(bool(row["enabled"]) for row in cells) == 48,
            "ENABLED_COUNT_MISMATCH")
    require(sum(not bool(row["enabled"]) for row in cells) == 12,
            "DISABLED_COUNT_MISMATCH")

    packet_hash = sha256_bytes(packet_bytes)
    output = render(packet, packet_hash)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(output, encoding="utf-8", newline="\n")
    print(json.dumps({
        "packet_sha256": packet_hash,
        "output_sha256": sha256_bytes(args.output.read_bytes()),
        "states": 5,
        "events": 12,
        "cells": 60,
        "enabled": 48,
        "disabled": 12,
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
