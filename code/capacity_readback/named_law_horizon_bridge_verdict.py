#!/usr/bin/env python3
"""Issue #650 short-circuit verdict after the named-law capacity exit.

The horizon bridge can only identify an existing named-law global-capacity
object with an independently constructed horizon record.  The issue-648
consumer selected no global-capacity branch and explicitly kept the stage
gate closed.  The horizon lane is therefore not evaluable on this source
branch.  This result neither rejects an inhabited gravity carrier nor rules
out a future capacity law satisfying the recorded reopen condition.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import named_law_n_closure_verdict as capacity_verdict


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
CAPACITY_VERDICT = HERE / "runtime" / "named_law_n_closure_verdict.json"
OUTPUT_PATH = HERE / "runtime" / "named_law_horizon_bridge_verdict.json"

SCHEMA = "oph.named_law_horizon_bridge_verdict.v1"
STATUS = "NOT_EVALUABLE_NO_GLOBAL_CAPACITY_OBJECT"
ISSUE = 650
PARENT_SCHEMA = "oph.named_law_n_closure_verdict.v1"
PARENT_STATUS = "NOT_EVALUABLE_NO_GLOBAL_CAPACITY_ATTACHMENT"


def canonical(value: Any) -> bytes:
    return (
        json.dumps(value, ensure_ascii=True, separators=(",", ":"), sort_keys=True)
        + "\n"
    ).encode("ascii")


def tagged(raw: bytes) -> str:
    return "sha256:" + hashlib.sha256(raw).hexdigest()


def load(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError("parent is not an object")
    return value


def pin(path: Path) -> dict[str, Any]:
    raw = path.read_bytes()
    return {
        "path": path.relative_to(ROOT).as_posix(),
        "bytes": len(raw),
        "sha256": tagged(raw),
    }


def build() -> dict[str, Any]:
    parent = load(CAPACITY_VERDICT)
    expected_parent = capacity_verdict.build()
    if parent != expected_parent:
        raise ValueError("named-law capacity verdict failed producer replay")
    if CAPACITY_VERDICT.read_bytes() != capacity_verdict.canonical_json_bytes(
        expected_parent
    ):
        raise ValueError("named-law capacity verdict is not canonical")
    if parent.get("schema") != PARENT_SCHEMA or parent.get("status") != PARENT_STATUS:
        raise ValueError("named-law capacity verdict drift")
    boundary = parent.get("comparison_boundary", {})
    if boundary.get("issue_650_stage_gate_opened") is not False:
        raise ValueError("horizon stage gate unexpectedly opened")
    if boundary.get("cosmological_prediction") is not False:
        raise ValueError("parent promoted a cosmological prediction")
    scope = parent.get("exit_scope", {})
    if scope.get("bounded_to_declared_attachment_class") is not True:
        raise ValueError("parent lost its bounded scope")
    if scope.get("all_possible_global_capacity_laws_excluded") is not False:
        raise ValueError("parent overstates the negative result")

    value: dict[str, Any] = {
        "schema": SCHEMA,
        "issue": ISSUE,
        "status": STATUS,
        "parent_pin": pin(CAPACITY_VERDICT),
        "short_circuit": {
            "required_input": "a source-positive named-law global-capacity object",
            "input_available": False,
            "horizon_identification_attempted": False,
            "reason": (
                "issue 648 selects no global-capacity branch, so there is no "
                "typed capacity object for a horizon commuting square"
            ),
        },
        "scope": {
            "declared_named_law_branch_only": True,
            "inhabited_gravity_carrier_rejected": False,
            "all_future_capacity_laws_excluded": False,
            "finite_de_sitter_identities_affected": False,
            "cosmological_value_emitted": False,
            "comparison_permitted": False,
        },
        "unaffected_issue": {
            "issue": 503,
            "statement": (
                "the independent inhabited Einstein/de Sitter carrier campaign "
                "keeps its own theorem-or-no-go scope"
            ),
        },
        "reopen_condition": (
            "issue 648 is reopened by a target-independent source law that selects "
            "one global-capacity action; the horizon lane must then construct an "
            "independent typed horizon record and normalization-preserving commuting square"
        ),
    }
    value["verdict_sha256"] = tagged(canonical(value))
    return value


def verify() -> None:
    expected = build()
    actual = load(OUTPUT_PATH)
    if actual != expected or OUTPUT_PATH.read_bytes() != canonical(expected):
        raise SystemExit("named-law horizon bridge verdict is stale")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.write:
        OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT_PATH.write_bytes(canonical(build()))
    if args.verify:
        verify()
        print("NAMED_LAW_HORIZON_BRIDGE_VERDICT_VALID")
    if not args.write and not args.verify:
        print(canonical(build()).decode("ascii"), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
