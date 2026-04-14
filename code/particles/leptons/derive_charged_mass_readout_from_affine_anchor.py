#!/usr/bin/env python3
"""Emit the algebraic charged-mass readout theorem from theorem-grade A_ch(P)."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

from charged_absolute_route_common import (
    ANCHOR_SECTION_JSON,
    MASS_READOUT_FROM_AFFINE_JSON,
    PHYSICAL_CANONICAL_LIFT_JSON,
    PHYSICAL_CLASS_AFFINE_SCALAR_REDUCTION_JSON,
    artifact_ref,
    load_json,
)


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(anchor: dict, canonical_lift: dict, physical_reduction: dict) -> dict:
    return {
        "artifact": "oph_charged_mass_readout_from_affine_anchor",
        "generated_utc": _timestamp(),
        "proof_status": "closed_algebraic_readout_theorem",
        "public_promotion_allowed": False,
        "theorem_scope": "conditional_on_theorem_grade_A_ch_on_physical_charged_data",
        "theorem_id": "charged_mass_readout_from_affine_anchor",
        "theorem_statement": (
            "Assume theorem-grade physical charged data Y_e(P) on which the affine charged anchor "
            "A_ch(P) is theorem-grade. Then the charged absolute scale and charged mass triple are "
            "algebraic readouts: g_e(P) = exp(A_ch(P)), Delta_e_abs(P) = log(g_ch_shared) - A_ch(P), "
            "and m_i(P) = exp(A_ch(P) + ell_i_centered(P)). Equivalently, once the determinant-line "
            "datum has fixed A_ch(P) = (1/3) * log(det(Y_e(P))), no extra theorem remains between "
            "the affine anchor and the emitted charged masses."
        ),
        "dependencies": {
            "absolute_anchor_section": {
                "artifact": anchor.get("artifact"),
                "artifact_ref": artifact_ref(ANCHOR_SECTION_JSON),
                "exact_missing_object": anchor.get("exact_missing_object"),
                "induced_formula_on_fill": anchor.get("induced_formula_on_fill"),
            },
            "canonical_lift": {
                "artifact": canonical_lift.get("artifact"),
                "artifact_ref": artifact_ref(PHYSICAL_CANONICAL_LIFT_JSON),
                "theorem_id": canonical_lift.get("theorem_id"),
            },
            "physical_class_reduction": {
                "artifact": physical_reduction.get("artifact"),
                "artifact_ref": artifact_ref(PHYSICAL_CLASS_AFFINE_SCALAR_REDUCTION_JSON),
                "theorem_id": physical_reduction.get("theorem_id"),
            },
        },
        "exact_readout": {
            "affine_anchor": "A_ch(P)",
            "absolute_scale": "g_e(P) = exp(A_ch(P))",
            "absolute_gap": "Delta_e_abs(P) = log(g_ch_shared) - A_ch(P)",
            "charged_masses": "m_i(P) = exp(A_ch(P) + ell_i_centered(P))",
            "determinant_line_form": "A_ch(P) = (1/3) * log(det(Y_e(P)))",
        },
        "why_this_closes_the_downstream_lane": [
            "The common-shift breaking datum is exactly the affine anchor A_ch(P).",
            "Once A_ch(P) is theorem-grade, the absolute scale g_e(P) is a single exponential readout.",
            "The charged masses are the same centered charged logs shifted by the theorem-grade affine anchor.",
            "No further bridge-sized theorem sits between A_ch(P) and the charged mass triple.",
        ],
        "scope_exclusions": [
            "No theorem-grade source landing B_D10(P) -> physical_Y_e(P) is proved by this artifact.",
            "No theorem-grade source landing B_D10(P) -> charged_determinant_line_section(P) is proved by this artifact.",
            "No promotion C_hat_e^{cand} -> C_hat_e is proved by this artifact.",
        ],
        "issue_repair_interpretation": {
            "closeable_issue_form": (
                "Once theorem-grade A_ch(P) exists, emit g_e(P) and the charged mass triple algebraically."
            ),
            "separate_open_upstream_math": "d10_to_charged_determinant_line_bridge",
        },
        "notes": [
            "This is the exact downstream theorem beneath the charged source-landing frontier.",
            "It is algebraic and does not use target anchoring, compare-only normalization, or same-family substitution.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the charged mass readout theorem from theorem-grade A_ch(P).")
    parser.add_argument("--anchor", default=str(ANCHOR_SECTION_JSON))
    parser.add_argument("--canonical-lift", default=str(PHYSICAL_CANONICAL_LIFT_JSON))
    parser.add_argument("--physical-reduction", default=str(PHYSICAL_CLASS_AFFINE_SCALAR_REDUCTION_JSON))
    parser.add_argument("--output", default=str(MASS_READOUT_FROM_AFFINE_JSON))
    args = parser.parse_args()

    artifact = build_artifact(
        load_json(Path(args.anchor)),
        load_json(Path(args.canonical_lift)),
        load_json(Path(args.physical_reduction)),
    )

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
