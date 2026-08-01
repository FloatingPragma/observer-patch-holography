#!/usr/bin/env python3
"""Issue #648 consumer verdict after the issue-649 attachment result.

The two reserve formulae remain exact retrospective arithmetic on their named
conditional branches.  The bounded issue-649 theorem shows that the local
finite-cut datum, even together with positive monoid composition and
refinement regrouping, does not select an action on global capacity.  The
named-law N lane therefore exits not evaluable at its first missing typed
attachment.  It emits no cosmic value and selects neither reserve branch.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import global_capacity_attachment as attachment_certificate
import n_closure_branch_certificate as branch_certificate


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
BRANCH_PACKET = HERE / "manifests" / "n_closure_branch_certificate.json"
ATTACHMENT_PACKET = HERE / "runtime" / "global_capacity_attachment_receipt.json"
OUTPUT_PATH = HERE / "runtime" / "named_law_n_closure_verdict.json"

SCHEMA = "oph.named_law_n_closure_verdict.v1"
STATUS = "NOT_EVALUABLE_NO_GLOBAL_CAPACITY_ATTACHMENT"
ISSUE = 648
BRANCH_SCHEMA = "oph.n_closure_branch_certificate.v1"
BRANCH_STATUS = "RETROSPECTIVE_UNSELECTED_CONDITIONAL_BRANCH_MENU"
ATTACHMENT_SCHEMA = "oph.global_capacity_attachment.v1"
ATTACHMENT_STATUS = "BOUNDED_NONIDENTIFIABLE_GLOBAL_COMPOSITION"


def canonical_json_bytes(value: Any) -> bytes:
    return (
        json.dumps(value, ensure_ascii=True, separators=(",", ":"), sort_keys=True)
        + "\n"
    ).encode("ascii")


def tagged(raw: bytes) -> str:
    return "sha256:" + hashlib.sha256(raw).hexdigest()


def _load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _pin(path: Path) -> dict[str, Any]:
    raw = path.read_bytes()
    return {
        "path": path.relative_to(ROOT).as_posix(),
        "bytes": len(raw),
        "sha256": tagged(raw),
    }


def build() -> dict[str, Any]:
    branch = _load(BRANCH_PACKET)
    attachment = _load(ATTACHMENT_PACKET)
    branch_certificate.validate(branch, verify_sources=True)
    attachment_certificate.validate(attachment, verify_source=True)
    if branch.get("schema") != BRANCH_SCHEMA or branch.get("status") != BRANCH_STATUS:
        raise ValueError("retrospective branch packet drift")
    if branch.get("scope", {}).get("branch_selected") is not False:
        raise ValueError("retrospective branch packet selected a branch")
    if branch.get("scope", {}).get("global_capacity_derived") is not False:
        raise ValueError("retrospective branch packet promoted global capacity")
    if (
        attachment.get("schema") != ATTACHMENT_SCHEMA
        or attachment.get("status") != ATTACHMENT_STATUS
    ):
        raise ValueError("global-capacity attachment verdict drift")
    boundary = attachment.get("comparison_boundary", {})
    if boundary.get("numeric_cosmic_capacity_emitted") is not False:
        raise ValueError("attachment packet emitted a cosmic capacity")
    if boundary.get("forecast_or_comparison_permitted") is not False:
        raise ValueError("attachment packet permitted a comparison")

    value: dict[str, Any] = {
        "schema": SCHEMA,
        "issue": ISSUE,
        "status": STATUS,
        "consumed_verdicts": {
            "reserve_branch_packet": BRANCH_STATUS,
            "global_capacity_attachment": ATTACHMENT_STATUS,
        },
        "strongest_supported_result": (
            "the finite-presence and Poisson expressions remain distinct exact "
            "retrospective conditional arithmetic; the declared source interface "
            "does not select either expression as an action on global capacity"
        ),
        "missing_antecedent": (
            "a source-derived typed action from the connected-cut survival object "
            "to global capacity, including selection of the physical blocked-event "
            "projector and its disconnected-cut/refinement composition"
        ),
        "branch_status": {
            "finite_presence_selected": False,
            "poisson_or_projective_limit_selected": False,
            "one_class_blocked_event_selected": False,
            "six_class_total_blocked_event_selected": False,
        },
        "comparison_boundary": {
            "retrospective_arithmetic_retained": True,
            "selective_weight": False,
            "cosmological_prediction": False,
            "prospective_forecast": False,
            "horizon_interpretation": False,
            "issue_650_stage_gate_opened": False,
        },
        "exit_scope": {
            "bounded_to_declared_attachment_class": True,
            "all_possible_global_capacity_laws_excluded": False,
            "new_source_law_forbidden": False,
            "reason": (
                "the countermodels establish nonselection under the declared "
                "finite-cut data plus positive monoid composition and cut-count "
                "regrouping; a stronger source-derived action can reopen the lane"
            ),
        },
        "unaffected_results": [
            "the typed self-readback equality after a same-invariant bridge",
            "the exact distinction between finite presence and Poisson factors",
            "the six-class translation-invariant uniformity theorem",
            "the bounded direct-N counterfamily and unfinished complete-source bridge",
        ],
        "reopen_condition": (
            "a target-independent source packet supplies the missing global action "
            "and blocked-event semantics, after which the E1 seam action, SR-1 seam "
            "content, inherited common-load premises, and optional Poisson carrier "
            "must be discharged on one provenance-coherent branch"
        ),
        "parent_pins": [_pin(BRANCH_PACKET), _pin(ATTACHMENT_PACKET)],
    }
    value["verdict_sha256"] = tagged(canonical_json_bytes(value))
    return value


def verify(path: Path = OUTPUT_PATH) -> None:
    expected = build()
    actual = _load(path)
    if actual != expected:
        raise SystemExit("named-law N closure verdict is stale")
    if path.read_bytes() != canonical_json_bytes(expected):
        raise SystemExit("named-law N closure verdict is not canonical")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.write:
        OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT_PATH.write_bytes(canonical_json_bytes(build()))
    if args.verify:
        verify()
        print("NAMED_LAW_N_CLOSURE_VERDICT_VALID")
    if not args.write and not args.verify:
        print(canonical_json_bytes(build()).decode("ascii"), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
