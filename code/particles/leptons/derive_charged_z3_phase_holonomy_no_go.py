#!/usr/bin/env python3
"""Certify that the current OPH data do not select the Stage-5 phase 2/9."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = (
    ROOT / "particles" / "runs" / "leptons"
    / "charged_z3_phase_holonomy_no_go.json"
)


def roots(delta: float) -> list[float]:
    return sorted(
        1.0 + math.sqrt(2.0) * math.cos(delta + 2.0 * math.pi * k / 3.0)
        for k in range(3)
    )


def trace_invariants(delta: float) -> tuple[float, float]:
    values = roots(delta)
    return sum(values), sum(value * value for value in values)


def build_artifact() -> dict[str, Any]:
    controls = {}
    for label, delta in (("identity_transport", 0.0), ("stage5", 2.0 / 9.0)):
        values = roots(delta)
        trace_one, trace_two = trace_invariants(delta)
        controls[label] = {
            "delta_radians": delta,
            "positive_roots": all(value > 0.0 for value in values),
            "ordered_roots": values,
            "trace": trace_one,
            "trace_square": trace_two,
            "Q": trace_two / (trace_one * trace_one),
        }
    return {
        "artifact": "oph_charged_z3_phase_holonomy_no_go",
        "status": "CLOSED_NO_GO_CURRENT_OPH_DATA_LEAVE_PHASE_CONTINUOUS",
        "public_phase_promotion_allowed": False,
        "countermodels": controls,
        "theorem": {
            "statement": (
                "The balanced Hermitian Z3 circulant has trace 3, squared trace 6, and Q=2/3 "
                "for every phase delta. Hence balance, Z3 symmetry, and the current invariant counts "
                "do not select delta=2/9."
            ),
            "physical_loop_invariant": "arg(C_12 C_23 C_31)=3 delta mod 2 pi",
            "family_gauge_boundary": (
                "Hypercharge acts on three generation copies as exp(i Y theta) I_3 and supplies "
                "a common phase, not the cyclic family shift."
            ),
        },
        "missing_source_objects": [
            "a physical C3 generation bundle",
            "a source-emitted oriented U(1) flavor connection and declared loop",
            "a gauge-invariant law fixing the total loop holonomy",
            "an equal-link gauge convention and orientation rule",
            "an attachment of the resulting roots to physical charged-family lines",
        ],
        "claim_boundary": (
            "The arithmetic identity (N_c+1)/(2 N_c N_g)=2/9 is not a phase-transport theorem."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(artifact["status"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
