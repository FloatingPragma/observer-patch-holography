#!/usr/bin/env python3
"""Emit the theorem-side quark relative-sheet selector object.

This script packages the exact next discrete theorem object even when the
current corpus cannot yet emit a selector value. A value is emitted only if one
orbit element carries an explicit theorem-grade selection witness.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
REPAIR_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_physical_branch_repair_theorem.json"
ORBIT_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_sigma_ud_orbit.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_relative_sheet_selector.json"
CANONICAL_SINGLETON_BRANCH_KEY = ["D12", "sigma_ref"]


def _timestamp() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _theorem_selected(elements: list[dict[str, Any]]) -> list[dict[str, Any]]:
    selected: list[dict[str, Any]] = []
    for item in elements:
        proof = dict(item.get("selection_proof") or {})
        if bool(proof.get("theorem_grade_select", False)):
            selected.append(item)
    return selected


def _closed_singleton_selector_surface(repair: dict[str, Any], orbit: dict[str, Any]) -> bool:
    del repair
    return orbit.get("selector_status") == "quark_relative_sheet_selector_closed_to_reference_singleton"


def _assert_closed_singleton_branch_key_alignment(
    repair: dict[str, Any],
    orbit: dict[str, Any],
    sigma_value: dict[str, Any] | None,
) -> None:
    if not _closed_singleton_selector_surface(repair, orbit):
        return

    expected = list(CANONICAL_SINGLETON_BRANCH_KEY)
    checks = {
        "orbit.selected_sigma.branch_key": ((orbit.get("selected_sigma") or {}).get("branch_key")),
        "orbit.provider_frontier.emitted_reference_sheet.branch_key": (
            (((orbit.get("provider_frontier") or {}).get("emitted_reference_sheet") or {}).get("branch_key"))
        ),
        "repair.current_d12_sheet.branch_key": ((repair.get("current_d12_sheet") or {}).get("branch_key")),
        "repair.minimal_branch_shift_repair_theorem.branch_key_after_repair": (
            ((repair.get("minimal_branch_shift_repair_theorem") or {}).get("branch_key_after_repair"))
        ),
        "repair.minimal_branch_shift_repair_theorem.selected_value.branch_key": (
            (((repair.get("minimal_branch_shift_repair_theorem") or {}).get("selected_value") or {}).get("branch_key"))
        ),
        "repair.relative_sheet_scan.selector_value.branch_key": (
            (((repair.get("relative_sheet_scan") or {}).get("selector_value") or {}).get("branch_key"))
        ),
        "selector.value.branch_key": ((sigma_value or {}).get("branch_key")),
    }
    mismatched = {label: value for label, value in checks.items() if value != expected}
    if mismatched:
        details = ", ".join(f"{label}={value!r}" for label, value in mismatched.items())
        raise ValueError(
            "closed singleton quark selector surface has inconsistent branch keys; "
            f"expected {expected!r} across orbit/repair/selector artifacts, got {details}"
        )


def build_artifact(repair: dict[str, Any], orbit: dict[str, Any]) -> dict[str, Any]:
    elements = list(orbit.get("elements") or [])
    selected = _theorem_selected(elements)
    debug_ranked = ((orbit.get("debug_compare_shell_ranking") or {}).get("ranked") or [])
    theorem = dict(repair["minimal_branch_shift_repair_theorem"])

    if len(selected) == 1:
        winner = selected[0]
        sigma_value = {
            "sigma_id": winner["sigma_id"],
            "canonical_token": winner["canonical_token"],
            "branch_key": winner.get("branch_key") or ["D12", winner["sigma_id"]],
        }
        selection_rule_status = "theorem_grade_rule_emitted_on_current_solver_surface"
        selection_status = "theorem_grade_value_emitted"
        reason = "Exactly one orbit element carries a theorem-grade selection witness."
    elif len(selected) > 1:
        sigma_value = None
        selection_rule_status = "ambiguous_theorem_rule_surface"
        selection_status = "ambiguous_theorem_witnesses"
        reason = "More than one orbit element claims theorem-grade selection; no value emitted."
    elif len(elements) == 1:
        sigma_value = None
        selection_rule_status = "open_target_free_rule_unemitted"
        selection_status = "singleton_reference_sheet_not_promoted"
        reason = (
            "The current solver emits one D12 reference-sheet orbit element, but no theorem-grade uniqueness or "
            "selection witness identifies that singleton as sigma_ud."
        )
    else:
        sigma_value = None
        selection_rule_status = "open_target_free_rule_unemitted"
        selection_status = "not_emitted_from_current_corpus"
        reason = (
            "No orbit element carries a theorem-grade selection witness. The current corpus therefore emits only the selector object, "
            "not a selector value."
        )

    _assert_closed_singleton_branch_key_alignment(repair, orbit, sigma_value)

    return {
        "artifact": "oph_quark_relative_sheet_selector",
        "generated_utc": _timestamp(),
        "scope": "D12_relative_sheet_selection",
        "id": theorem["id"],
        "definition": theorem["definition"],
        "selector_domain": theorem.get("selector_domain", "left_handed_same_label_relative_sheet_orbit_only"),
        "branch_key_after_repair": theorem["branch_key_after_repair"],
        "must_not_use_compare_fit_masses": theorem["must_not_use_compare_fit_masses"],
        "must_not_use_same_sheet_rephasing": theorem["must_not_use_same_sheet_rephasing"],
        "selection_rule_status": selection_rule_status,
        "selection_status": selection_status,
        "reason": reason,
        "quark_relative_sheet_selector": {
            "symbol": "sigma_ud",
            "carrier": "Sigma_ud",
            "value": sigma_value,
        },
        "debug_best_candidate": debug_ranked[0] if debug_ranked else None,
        "debug_best_candidate_promotable": False,
        "next_after_selection": {
            "id": "quark_d12_t1_value_law",
            "must_emit": "quark_d12_t1_value_law",
            "scalar_name": "t1",
            "identifies": "ray_modulus = t1",
            "unique_intersection_with": "D12_ud_mass_ray",
            "then_emits": ["t1", "ray_modulus", "Delta_ud_overlap", "eta_Q_centered"],
            "derived_wrapper": "intrinsic_scale_law_D12",
            "must_not_use_target_masses": True,
            "must_not_use_ckm_cp": True,
        },
        "dependencies": {
            "repair_theorem_artifact": "quark_physical_branch_repair_theorem.json",
            "solver_artifact": "quark_sigma_ud_orbit.json",
        },
        "notes": [
            (
                "The local same-label left-handed selector is now theorem-emitted on the current solver surface and its unique value is sigma_ref."
                if sigma_value is not None
                else "This is the exact next theorem-side object even when the selector value remains open."
            ),
            "This script refuses to convert compare-only CKM-shell ranking into a theorem-grade selection.",
            (
                "That emitted value is a negative closure: sigma_ref is the current D12 reference sheet, so the physical CKM-shell no-go survives unchanged on the selected branch."
                if sigma_value is not None
                else "A singleton reference-sheet orbit element is still not enough unless the solver also emits an intrinsic uniqueness witness."
            ),
            (
                "The next exact theorem-side quark object is therefore quark_d12_t1_value_law on the emitted D12_ud_mass_ray, with intrinsic_scale_law_D12 retained only as the derived wrapper."
                if sigma_value is not None
                else "Once a theorem-grade orbit witness exists, rerunning this script will emit sigma_ud directly."
            ),
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the quark relative-sheet selector artifact.")
    parser.add_argument("--repair-theorem", default=str(REPAIR_JSON))
    parser.add_argument("--orbit", default=str(ORBIT_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    repair = _load_json(Path(args.repair_theorem))
    orbit = _load_json(Path(args.orbit))
    payload = build_artifact(repair, orbit)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
