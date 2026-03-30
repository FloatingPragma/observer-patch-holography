#!/usr/bin/env python3
"""Emit the quark relative-sheet orbit scaffold.

This is a solver-extension contract, not a solved branch selector. It records
the finite orbit that must be exposed before a theorem-grade
``quark_relative_sheet_selector`` can be discussed honestly. When a caller
supplies finite candidate elements, the script packages them into the canonical
artifact and computes a compare-only debug ranking. That ranking is never
promotable to theorem grade.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_sigma_ud_orbit.json"
TARGET_THETA_12 = 0.2256
TARGET_THETA_23 = 0.0438
TARGET_THETA_13 = 0.00347


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _extract_ckm(item: dict[str, Any]) -> dict[str, Any]:
    if isinstance(item.get("ckm_invariants"), dict):
        return dict(item["ckm_invariants"])
    if isinstance(item.get("ckm"), dict):
        return dict(item["ckm"])
    raise KeyError("orbit element must contain `ckm_invariants` or `ckm`")


def _debug_log_shell_loss(ckm: dict[str, Any]) -> float:
    theta12 = float(ckm["theta_12"])
    theta23 = float(ckm["theta_23"])
    theta13 = float(ckm["theta_13"])
    return (
        math.log(theta12 / TARGET_THETA_12) ** 2
        + math.log(theta23 / TARGET_THETA_23) ** 2
        + math.log(theta13 / TARGET_THETA_13) ** 2
    )


def _rank_elements(elements: list[dict[str, Any]]) -> list[dict[str, Any]]:
    ranked: list[dict[str, Any]] = []
    for item in elements:
        ckm = _extract_ckm(item)
        theta12 = float(ckm["theta_12"])
        theta23 = float(ckm["theta_23"])
        theta13 = float(ckm["theta_13"])
        ranked.append(
            {
                "sigma_id": item.get("sigma_id"),
                "canonical_token": item.get("canonical_token"),
                "loss": _debug_log_shell_loss(ckm),
                "abs_log_error_theta13": abs(math.log(theta13 / TARGET_THETA_13)),
                "abs_log_error_theta23": abs(math.log(theta23 / TARGET_THETA_23)),
                "abs_log_error_theta12": abs(math.log(theta12 / TARGET_THETA_12)),
            }
        )
    ranked.sort(
        key=lambda item: (
            item["loss"],
            item["abs_log_error_theta13"],
            item["abs_log_error_theta23"],
            item["abs_log_error_theta12"],
            str(item["sigma_id"]),
        )
    )
    return ranked


def _load_elements(path: Path | None) -> list[dict[str, Any]]:
    if path is None:
        return []
    raw = _load_json(path)
    if isinstance(raw, dict) and isinstance(raw.get("elements"), list):
        return list(raw["elements"])
    if isinstance(raw, list):
        return list(raw)
    raise ValueError("elements-json must be a list or an object with an `elements` field")


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the quark sigma_ud orbit scaffold.")
    parser.add_argument("--elements-json", default=None)
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    elements = _load_elements(Path(args.elements_json) if args.elements_json else None)
    debug_ranking = _rank_elements(elements) if elements else []

    artifact = {
        "artifact": "oph_quark_sigma_ud_orbit",
        "generated_utc": _timestamp(),
        "status": "missing_solver_side_orbit" if not elements else "candidate_orbit_elements_supplied",
        "public_promotion_allowed": False,
        "exact_missing_object": "sigma_ud_orbit",
        "exact_missing_solver_interface": "sigma_ud_orbit_provider_interface",
        "concrete_provider_scaffold": "code/particles/flavor/sigma_ud_orbit_provider.py",
        "orbit_kind": "finite_relative_sheet_orbit",
        "branch_key": ["D12", "sigma_ud"],
        "selector_status": (
            "quark_relative_sheet_selector_not_emittable_without_orbit"
            if not elements
            else "selection_rule_still_open_target_free"
        ),
        "input_contract": {
            "must_use": [
                "forward_yukawas.json",
                "same-label transport lifts",
                "the emitted D12 reference-sheet representative",
            ],
            "must_not_use": [
                "CKM target fitting inside the theorem artifact",
                "target masses",
                "same-sheet rephasing as a repair",
            ],
        },
        "elements_schema": [
            "sigma_id",
            "canonical_token",
            "U_u_left",
            "U_d_left",
            "V_CKM",
            "ckm_invariants",
        ],
        "elements": elements,
        "theorem_grade_selection": None,
        "selection_rule_status": "open_target_free_rule_unemitted",
        "debug_compare_shell_ranking": {
            "promotable": False,
            "kind": "ckm_log_shell_loss",
            "shell": {
                "theta_12": TARGET_THETA_12,
                "theta_23": TARGET_THETA_23,
                "theta_13": TARGET_THETA_13,
            },
            "ranked": debug_ranking,
        },
        "compare_only_helper_contract": {
            "script": "score_quark_sigma_ud_orbit_against_ckm_shell.py",
            "status": "compare_only_only",
            "must_not_promote_selector": True,
        },
        "selection_gate": {
            "quark_relative_sheet_selector": None,
            "may_emit_only_if": [
                "orbit collapses to one intrinsic canonical token",
                "or an intrinsic non-target selection theorem is proved",
            ],
        },
        "notes": [
            "This scaffold exists to make the missing finite solver object explicit.",
            "The current D12 sheet is transport-closed but wrong-branch; same-sheet changes cannot move CKM invariants to the physical shell.",
            "Branch selection is discrete here. A continuous scalar cannot replace orbit exposure.",
            "On the current live corpus the more immediate implementation gap is the orbit provider itself: no same-label left-handed Sigma_ud enumerator or sigma-to-CKM evaluator is emitted yet.",
            "The provider interface has been widened to the full left-handed evaluation schema expected by this artifact; what is still missing is a real implementation that can populate those fields from emitted same-label transport data.",
            "If elements is empty, the artifact records the honest frontier rather than inventing Sigma_ud.",
            "If elements are supplied, the debug ranking remains comparison-only and cannot be promoted.",
        ],
    }

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
