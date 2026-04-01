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

from sigma_ud_orbit_provider import (
    build_emitted_reference_sheet_orbit_elements,
    build_sigma_ud_provider_frontier,
    load_already_local_diagnostic_orbit,
    load_sigma_ud_singleton_uniqueness_witness,
    load_transport_frame_diagnostic_orbit,
)


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

    elements = (
        _load_elements(Path(args.elements_json))
        if args.elements_json
        else build_emitted_reference_sheet_orbit_elements()
    )
    debug_ranking = _rank_elements(elements) if elements else []
    provider_frontier = build_sigma_ud_provider_frontier()
    already_local_diagnostic_orbit = load_already_local_diagnostic_orbit()
    diagnostic_transport_frame_orbit = load_transport_frame_diagnostic_orbit()
    singleton_uniqueness = load_sigma_ud_singleton_uniqueness_witness()
    default_singleton = bool(provider_frontier.get("emitted_reference_sheet")) and not args.elements_json
    singleton_closed = bool(singleton_uniqueness.get("theorem_grade_select")) and default_singleton

    artifact = {
        "artifact": "oph_quark_sigma_ud_orbit",
        "generated_utc": _timestamp(),
        "status": (
            "missing_solver_side_orbit"
            if not elements
            else (
                "same_label_left_handed_local_orbit_singleton_closed"
                if singleton_closed
                else (
                    "reference_sheet_singleton_emitted_orbit_incomplete"
                    if default_singleton
                    else "candidate_orbit_elements_supplied"
                )
            )
        ),
        "public_promotion_allowed": False,
        "exact_missing_object": (
            None
            if singleton_closed
            else (
                "first_distinct_same_label_left_handed_relative_sheet_evaluation"
                if default_singleton
                else "sigma_ud_orbit"
            )
        ),
        "orbit_completion_target": "sigma_ud_orbit",
        "exact_missing_solver_interface": "sigma_ud_orbit_provider_interface",
        "exact_missing_provider_methods": [
            "enumerate_relative_sheets_d12()",
            "evaluate_relative_sheet(token)",
        ],
        "concrete_provider_scaffold": "code/particles/flavor/sigma_ud_orbit_provider.py",
        "provider_frontier": provider_frontier,
        "orbit_kind": "finite_relative_sheet_orbit",
        "branch_key": ["D12", "sigma_ud"],
        "elements_origin": (
            "reference_sheet_singleton_provider"
            if default_singleton
            else ("external_elements_json" if args.elements_json else "no_provider_output")
        ),
        "selector_status": (
            "quark_relative_sheet_selector_not_emittable_without_orbit"
            if not elements
            else (
                "quark_relative_sheet_selector_closed_to_reference_singleton"
                if singleton_closed
                else (
                    "quark_relative_sheet_selector_not_emittable_without_distinct_relative_sheet"
                    if default_singleton
                    else "selection_rule_still_open_target_free"
                )
            )
        ),
        "next_exact_object_after_orbit_closure": (
            "quark_d12_t1_value_law" if singleton_closed else None
        ),
        "selected_sigma": (
            singleton_uniqueness.get("selected_sigma") if singleton_closed else None
        ),
        "singleton_uniqueness_theorem": singleton_uniqueness,
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
        "already_local_diagnostic_orbit": already_local_diagnostic_orbit,
        "diagnostic_transport_frame_orbit": diagnostic_transport_frame_orbit,
        "theorem_grade_selection": (
            singleton_uniqueness.get("selected_sigma") if singleton_closed else None
        ),
        "selection_rule_status": (
            "closed_by_singleton_uniqueness_theorem"
            if singleton_closed
            else "open_target_free_rule_unemitted"
        ),
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
            "quark_relative_sheet_selector": (
                singleton_uniqueness.get("selected_sigma") if singleton_closed else None
            ),
            "may_emit_only_if": [
                "orbit collapses to one intrinsic canonical token",
                "or an intrinsic non-target selection theorem is proved",
            ],
        },
        "notes": [
            "This scaffold exists to make the missing finite solver object explicit.",
            "The current D12 sheet is transport-closed but wrong-branch; same-sheet changes cannot move CKM invariants to the physical shell.",
            "Branch selection is discrete here. A continuous scalar cannot replace orbit exposure.",
            "The common-refinement line-lift now feeds a derived transport-frame diagnostic orbit artifact rather than a hardcoded note. Its self-overlap F0^dagger F1 is a real already-local compare-only witness, but it is still not a sector-attached Sigma_ud element.",
            (
                "The emitted local same-label left-handed orbit now closes to the singleton sigma_ref because only L/L survives the chirality admissibility check and the published five-anchor standard CKM gauge fixes the remaining diagonal rephasing to the trivial global phase."
                if singleton_closed
                else (
                    "The current live corpus now emits one real same-label left-handed D12 reference-sheet evaluation in the full orbit-element schema."
                    if default_singleton
                    else "On the current live corpus the more immediate implementation gap is the first non-empty provider output itself: no same-label left-handed Sigma_ud enumerator or sigma-to-CKM evaluator is emitted yet."
                )
            ),
            (
                "That closes the solver-side sigma_ud orbit and shifts the exact next object to the selected-branch intrinsic scale law on D12_ud_mass_ray."
                if singleton_closed
                else (
                    "That honest singleton does not close Sigma_ud: the smaller exact blocker is now one additional non-reference same-label left-handed sheet evaluation, or an intrinsic uniqueness theorem proving the singleton is the full orbit."
                    if default_singleton
                    else "The provider interface has been widened to the full left-handed evaluation schema expected by this artifact; what is still missing is a real implementation that can populate those fields from emitted same-label transport data."
                )
            ),
            "The already-local chirality-basis orbit is now threaded in explicitly as a diagnostic exclusion surface, so future solver work cannot mistake it for Sigma_ud.",
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
