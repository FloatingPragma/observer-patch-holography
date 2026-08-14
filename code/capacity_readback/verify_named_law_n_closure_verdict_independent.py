#!/usr/bin/env python3
"""Independent boundary replay for the issue-648 named-law verdict."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
VERDICT = HERE / "runtime" / "named_law_n_closure_verdict.json"
BRANCH = HERE / "manifests" / "n_closure_branch_certificate.json"
ATTACHMENT = HERE / "runtime" / "global_capacity_attachment_receipt.json"


def canonical(value: Any) -> bytes:
    return (
        json.dumps(value, ensure_ascii=True, separators=(",", ":"), sort_keys=True)
        + "\n"
    ).encode("ascii")


def tagged(raw: bytes) -> str:
    return "sha256:" + hashlib.sha256(raw).hexdigest()


def pin(path: Path) -> dict[str, Any]:
    raw = path.read_bytes()
    return {
        "path": path.relative_to(ROOT).as_posix(),
        "bytes": len(raw),
        "sha256": tagged(raw),
    }


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(message)


def verify() -> None:
    value = json.loads(VERDICT.read_text(encoding="ascii"))
    branch = json.loads(BRANCH.read_text(encoding="utf-8"))
    attachment = json.loads(ATTACHMENT.read_text(encoding="ascii"))
    require(value.get("schema") == "oph.named_law_n_closure_verdict.v1", "schema")
    require(value.get("issue") == 648, "issue")
    require(
        value.get("status") == "NOT_EVALUABLE_NO_GLOBAL_CAPACITY_ATTACHMENT",
        "status",
    )
    require(value.get("parent_pins") == [pin(BRANCH), pin(ATTACHMENT)], "pins")
    require(branch.get("scope", {}).get("branch_selected") is False, "branch selected")
    require(branch.get("scope", {}).get("global_capacity_derived") is False, "capacity promoted")
    require(
        attachment.get("status")
        == "BOUNDED_NONIDENTIFIABLE_GLOBAL_COMPOSITION",
        "attachment status",
    )
    require(
        all(selected is False for selected in value.get("branch_status", {}).values()),
        "branch promotion",
    )
    boundary = value.get("comparison_boundary", {})
    require(boundary.get("retrospective_arithmetic_retained") is True, "arithmetic lost")
    for key in (
        "selective_weight",
        "cosmological_prediction",
        "prospective_forecast",
        "horizon_interpretation",
        "named_law_comparison_promotable",
    ):
        require(boundary.get(key) is False, f"promotion: {key}")
    scope = value.get("exit_scope", {})
    require(scope.get("bounded_to_declared_attachment_class") is True, "scope")
    require(scope.get("all_possible_global_capacity_laws_excluded") is False, "overclaim")
    require(scope.get("new_source_law_forbidden") is False, "reopen blocked")
    expected = value.get("verdict_sha256")
    unhashed = dict(value)
    unhashed.pop("verdict_sha256", None)
    require(expected == tagged(canonical(unhashed)), "hash")


if __name__ == "__main__":
    verify()
    print("NAMED_LAW_N_CLOSURE_VERDICT_INDEPENDENT_VALID")
