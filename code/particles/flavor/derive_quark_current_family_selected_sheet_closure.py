#!/usr/bin/env python3
"""Emit the exact same-family selected-sheet quark closure artifact."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SELECTOR_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_relative_sheet_selector.json"
EXACT_READOUT_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_exact_readout.json"
QUADRATIC_THEOREM_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_quadratic_readout_theorem.json"
PHYSICAL_BRANCH_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_physical_branch_repair_theorem.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_selected_sheet_closure.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(selector: dict, exact_readout: dict, quadratic_theorem: dict, physical_branch: dict) -> dict:
    selected_sheet_payload = dict(selector["quark_relative_sheet_selector"])
    selected_value = dict(selected_sheet_payload["value"])
    current_sheet = dict(physical_branch["current_d12_sheet"])
    return {
        "artifact": "oph_quark_current_family_selected_sheet_exact_closure",
        "generated_utc": _timestamp(),
        "proof_status": "closed_current_family_selected_sheet_exact_readout",
        "theorem_scope": "current_family_only",
        "public_promotion_allowed": False,
        "supporting_selector_artifact": selector["artifact"],
        "supporting_exact_readout_artifact": exact_readout["artifact"],
        "supporting_quadratic_readout_theorem": quadratic_theorem["artifact"],
        "physical_branch_boundary_artifact": physical_branch["artifact"],
        "theorem_statement": (
            "On the emitted same-label left-handed selected sheet sigma_ref and the fixed ordered three-point family, "
            "the target-anchored quadratic readout together with the exact sector geometric means g_u and g_d "
            "emit one exact same-family quark sextet. This closes the selected-sheet exact readout chain on "
            "current_family_only."
        ),
        "selected_sheet": selected_value,
        "selected_sheet_payload": selected_sheet_payload,
        "current_d12_sheet": current_sheet,
        "exact_sector_geometric_means": {
            "g_u": float(exact_readout["g_u"]),
            "g_d": float(exact_readout["g_d"]),
        },
        "exact_outputs_gev": {
            "u": float(exact_readout["predicted_singular_values_u"][0]),
            "c": float(exact_readout["predicted_singular_values_u"][1]),
            "t": float(exact_readout["predicted_singular_values_u"][2]),
            "d": float(exact_readout["predicted_singular_values_d"][0]),
            "s": float(exact_readout["predicted_singular_values_d"][1]),
            "b": float(exact_readout["predicted_singular_values_d"][2]),
        },
        "exact_fit_residuals_u": exact_readout["exact_fit_residuals_u"],
        "exact_fit_residuals_d": exact_readout["exact_fit_residuals_d"],
        "do_not_claim_now": [
            "a theorem-grade quark_d12_t1_value_law",
            "a branch change away from the selected wrong-branch sigma_ref sheet",
            "a theorem-grade physical CKM closure on the selected D12 sheet",
        ],
        "notes": [
            "This closes the exact same-family quark readout on the selected local sheet.",
            "It does not remove the separate physical-branch no-go on sigma_ref and does not emit quark_d12_t1_value_law.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the exact same-family selected-sheet quark closure artifact.")
    parser.add_argument("--selector", default=str(SELECTOR_JSON))
    parser.add_argument("--exact-readout", default=str(EXACT_READOUT_JSON))
    parser.add_argument("--quadratic-theorem", default=str(QUADRATIC_THEOREM_JSON))
    parser.add_argument("--physical-branch", default=str(PHYSICAL_BRANCH_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    selector = json.loads(Path(args.selector).read_text(encoding="utf-8"))
    exact_readout = json.loads(Path(args.exact_readout).read_text(encoding="utf-8"))
    quadratic_theorem = json.loads(Path(args.quadratic_theorem).read_text(encoding="utf-8"))
    physical_branch = json.loads(Path(args.physical_branch).read_text(encoding="utf-8"))
    artifact = build_artifact(selector, exact_readout, quadratic_theorem, physical_branch)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
