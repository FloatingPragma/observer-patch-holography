#!/usr/bin/env python3
"""Emit the charged determinant-line section extension scaffold.

This does not close the charged lane. It sharpens the absolute-anchor problem:
the smallest constructive breaker beyond the current common-shift quotient is a
section of the charged determinant line. Once such a section exists, the
charged absolute anchor is read out as ``A_ch = (1/3) log det(Y_e)`` on the
promoted charged operator surface.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
ANCHOR_SECTION_JSON = ROOT / "particles" / "runs" / "leptons" / "charged_absolute_anchor_section.json"
CENTRAL_SPLIT_EXTENSION_JSON = (
    ROOT / "particles" / "runs" / "flavor" / "charged_central_split_transfer_extension.json"
)
DEFAULT_OUT = ROOT / "particles" / "runs" / "leptons" / "charged_determinant_line_section_extension.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def build_artifact(anchor_section: dict[str, Any], transfer_extension: dict[str, Any]) -> dict[str, Any]:
    return {
        "artifact": "oph_charged_determinant_line_section_extension",
        "generated_utc": _timestamp(),
        "status": "minimal_constructive_extension",
        "public_promotion_allowed": False,
        "extension_kind": "determinant_line_section",
        "exact_missing_object": "charged_determinant_line_section",
        "possible_smaller_missing_object": "charged_determinant_line_trivialization_or_metric_compatibility",
        "upstream_prerequisites": {
            "promotion_theorem": "oph_generation_bundle_branch_generator_splitting",
            "promotion_effect": "theorem_grade_C_hat_e",
            "required_operator_surface": "theorem_grade_physical_Y_e_or_equivalent_determinant_line",
            "transfer_extension_artifact": transfer_extension.get("artifact"),
        },
        "section_contract": {
            "kind": "determinant_torsor_section",
            "required_covariance": "s_det(exp(3c) * detY) = s_det(detY) + 3c",
            "charged_anchor_readout": "A_ch = (1/3) * s_det(det Y_e)",
            "derived_quantities": {
                "g_e": "exp(A_ch)",
                "Delta_e_abs": "log(g_ch_shared) - A_ch",
                "masses": "exp(A_ch + centered_log_i)",
            },
        },
        "canonical_formula_on_fill": {
            "A_ch": "(1/3) * log(det(Y_e))",
            "g_e": "exp((1/3) * log(det(Y_e)))",
            "determinant_rule": "det(Y_e) = exp(3 * A_ch)",
        },
        "why_this_is_sharper_than_bare_A_ch": [
            "The common-shift quotient acts on the determinant line by det(Y_e) -> exp(3c) det(Y_e).",
            "A section of that line is exactly the affine object needed to break the quotient symmetry.",
            "This packages the charged absolute anchor as a determinant-torsor coordinate rather than an abstract free scalar.",
        ],
        "input_contract": anchor_section.get("input_contract"),
        "hard_rejections": anchor_section.get("hard_rejections"),
        "notes": [
            "This is a constructive extension route, not a theorem hidden in the current corpus.",
            "Promoting C_hat_e^cand is still upstream and necessary, but the determinant-line section is the first object that can realize the affine common-shift breaker geometrically.",
            "If a determinant-line section still fails to become canonical on the live refinement family, the smaller exact missing object is a refinement-stable determinant-line trivialization or equivalent metric-compatibility theorem for that line bundle.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the charged determinant-line section extension scaffold.")
    parser.add_argument("--anchor-section", default=str(ANCHOR_SECTION_JSON))
    parser.add_argument("--transfer-extension", default=str(CENTRAL_SPLIT_EXTENSION_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    anchor_section = _load_json(Path(args.anchor_section))
    transfer_extension = _load_json(Path(args.transfer_extension))
    payload = build_artifact(anchor_section, transfer_extension)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
