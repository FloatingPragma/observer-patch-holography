#!/usr/bin/env python3
"""Emit the strongest honest charged selected-surface affine-scalar reduction theorem."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

from charged_absolute_route_common import (
    P_TO_AFFINE_REDUCTION_JSON,
    PHYSICAL_CANONICAL_LIFT_JSON,
    PHYSICAL_CLASS_AFFINE_SCALAR_REDUCTION_JSON,
    artifact_ref,
    load_json,
)


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(bridge_reduction: dict, canonical_lift: dict) -> dict:
    return {
        "artifact": "oph_charged_physical_class_affine_scalar_reduction",
        "generated_utc": _timestamp(),
        "proof_status": "closed_selected_surface_reduction_theorem",
        "public_promotion_allowed": False,
        "theorem_scope": "theorem_grade_physical_charged_class_only",
        "theorem_id": "charged_physical_class_affine_scalar_reduction",
        "theorem_statement": (
            "Let B_D10(P) denote the theorem-grade D10 calibration descendant package emitted from P. "
            "Assume there is a theorem-grade landing of B_D10(P) on theorem-grade physical charged data Y_e(P), "
            "or equivalently on a charged determinant-line section s_det(P). Assume further that on that same "
            "theorem-grade physical charged surface the promoted charged operator C_hat_e exists. Then the "
            "canonical determinant-line lift theorem defines mu_phys(Y_e(P)) = (1/3) * log(det(Y_e(P))) and "
            "C_tilde_e(Y_e(P)) = C_hat_e(Y_e(P)) + mu_phys(Y_e(P)) I, so the physical identity-mode equalizer, "
            "determinant-line datum, and affine anchor are fixed canonically by s_det(Y_e(P)) = 3 * mu_phys(Y_e(P)) "
            "and A_ch(P) = mu_phys(Y_e(P)). Therefore g_e(P), Delta_e_abs(P), and the charged mass triple are "
            "algebraic descendants of that same physical-class datum."
        ),
        "scope_exclusions": [
            "No public selected-class charged source descent analogous to the quark theorem is claimed.",
            "No theorem-grade landing B_D10(P) -> Y_e(P) is proved by this artifact.",
            "No promotion C_hat_e^{cand} -> C_hat_e is proved by this artifact.",
            "No global classification of charged physical classes is claimed.",
        ],
        "dependencies": {
            "bridge_reduction": {
                "artifact": bridge_reduction.get("artifact"),
                "artifact_ref": artifact_ref(P_TO_AFFINE_REDUCTION_JSON),
                "theorem_id": bridge_reduction.get("reduction_theorem", {}).get("id"),
                "exact_smallest_bridge_target": bridge_reduction.get("exact_smallest_bridge_target", {}).get("id"),
                "P_threaded_source_exactness_reformulation": bridge_reduction.get(
                    "P_threaded_source_exactness_reformulation", {}
                ).get("id"),
            },
            "canonical_lift": {
                "artifact": canonical_lift.get("artifact"),
                "artifact_ref": artifact_ref(PHYSICAL_CANONICAL_LIFT_JSON),
                "theorem_id": canonical_lift.get("theorem_id"),
            },
        },
        "canonical_consequences_on_fill": {
            "descended_scalar": "mu_phys(Y_e(P))",
            "determinant_line": "s_det(Y_e(P)) = 3 * mu_phys(Y_e(P))",
            "affine_anchor": "A_ch(P) = mu_phys(Y_e(P)) = (1/3) * log(det(Y_e(P)))",
            "absolute_scale": "g_e(P) = exp(A_ch(P))",
            "absolute_gap": "Delta_e_abs(P) = log(g_ch_shared) - A_ch(P)",
            "charged_masses": "m_i(P) = exp(A_ch(P) + ell_i_centered(P))",
        },
        "why_this_is_the_strongest_honest_selected_surface_theorem": [
            "The charged lane does not contain a public selected-class source descent theorem analogous to the quark theorem.",
            "The common-shift quotient no-go excludes any centered-only route to determinant-line or affine data.",
            "Once theorem-grade physical charged data and theorem-grade C_hat_e are granted, the determinant line fixes one canonical physical affine scalar and therefore the uncentered lift.",
            "The D10 bridge reduction proves that any P-driven charged absolute law reduces to a theorem-grade landing on physical Y_e(P) or on the charged determinant line.",
        ],
        "notes": [
            "This is a reduction theorem on theorem-grade physical charged classes, not a closure theorem for the present charged issue set.",
            "It packages the strongest honest selected-surface charged statement supported by the local theorem surface.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build the charged physical-class affine-scalar reduction theorem artifact."
    )
    parser.add_argument("--bridge-reduction", default=str(P_TO_AFFINE_REDUCTION_JSON))
    parser.add_argument("--canonical-lift", default=str(PHYSICAL_CANONICAL_LIFT_JSON))
    parser.add_argument("--output", default=str(PHYSICAL_CLASS_AFFINE_SCALAR_REDUCTION_JSON))
    args = parser.parse_args()

    artifact = build_artifact(
        load_json(Path(args.bridge_reduction)),
        load_json(Path(args.canonical_lift)),
    )

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
