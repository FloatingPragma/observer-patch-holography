#!/usr/bin/env python3
"""Independent boundary replay for the issue-650 short-circuit verdict."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
PARENT = HERE / "runtime" / "named_law_n_closure_verdict.json"
VERDICT = HERE / "runtime" / "named_law_horizon_bridge_verdict.json"


def canonical(value: Any) -> bytes:
    return (
        json.dumps(value, ensure_ascii=True, separators=(",", ":"), sort_keys=True)
        + "\n"
    ).encode("ascii")


def tagged(raw: bytes) -> str:
    return "sha256:" + hashlib.sha256(raw).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(message)


def verify() -> None:
    parent_raw = PARENT.read_bytes()
    parent = json.loads(parent_raw)
    value = json.loads(VERDICT.read_bytes())
    require(value.get("schema") == "oph.named_law_horizon_bridge_verdict.v1", "schema")
    require(value.get("issue") == 650, "issue")
    require(value.get("status") == "NOT_EVALUABLE_NO_GLOBAL_CAPACITY_OBJECT", "status")
    require(
        value.get("parent_pin")
        == {
            "path": PARENT.relative_to(ROOT).as_posix(),
            "bytes": len(parent_raw),
            "sha256": tagged(parent_raw),
        },
        "parent pin",
    )
    boundary = parent.get("comparison_boundary", {})
    require(boundary.get("named_law_comparison_promotable") is False, "comparison promotion")
    require(boundary.get("cosmological_prediction") is False, "parent prediction")
    short = value.get("short_circuit", {})
    require(short.get("input_available") is False, "input promotion")
    require(short.get("horizon_identification_attempted") is False, "attempt promotion")
    scope = value.get("scope", {})
    require(scope.get("declared_named_law_branch_only") is True, "scope")
    for key in (
        "inhabited_gravity_carrier_rejected",
        "all_future_capacity_laws_excluded",
        "finite_de_sitter_identities_affected",
        "cosmological_value_emitted",
        "comparison_permitted",
    ):
        require(scope.get(key) is False, f"promotion: {key}")
    expected = value.get("verdict_sha256")
    unhashed = dict(value)
    unhashed.pop("verdict_sha256", None)
    require(expected == tagged(canonical(unhashed)), "verdict hash")


if __name__ == "__main__":
    verify()
    print("NAMED_LAW_HORIZON_BRIDGE_VERDICT_INDEPENDENT_VALID")
