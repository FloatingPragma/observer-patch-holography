#!/usr/bin/env python3
"""Emit the current-corpus no-go for the D10-to-charged affine bridge.

This artifact isolates the third charged blocker explicitly: the current corpus
does not yet contain a theorem-grade bridge from the D10 calibration family of
``P`` to the charged determinant line or affine anchor ``A_ch``.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from charged_absolute_route_common import (
    ANCHOR_SECTION_JSON,
    P_TO_AFFINE_REDUCTION_JSON,
    POST_PROMOTION_ROUTE_JSON,
    artifact_ref,
    load_json,
)


ROOT = Path(__file__).resolve().parents[2]
D10_FORWARD_JSON = ROOT / "particles" / "runs" / "calibration" / "d10_ew_forward_transmutation_certificate.json"
D10_REPAIR_JSON = ROOT / "particles" / "runs" / "calibration" / "d10_ew_target_free_repair_value_law.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "leptons" / "charged_p_to_affine_anchor_bridge_no_go.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _as_dict(payload: dict[str, Any], key: str) -> dict[str, Any]:
    value = payload.get(key)
    return dict(value) if isinstance(value, dict) else {}


def build_artifact(
    d10_forward: dict[str, Any],
    d10_repair: dict[str, Any],
    anchor: dict[str, Any],
    post_promotion_route: dict[str, Any],
    reduction: dict[str, Any],
) -> dict[str, Any]:
    forward_theorem = _as_dict(d10_forward, "theorem")
    repair_theorem = _as_dict(d10_repair, "theorem")
    post_promotion_slot = _as_dict(post_promotion_route, "post_promotion_single_slot")
    exact_descended_scalar = _as_dict(post_promotion_slot, "exact_descended_scalar")

    return {
        "artifact": "oph_charged_p_to_affine_anchor_bridge_no_go",
        "generated_utc": _timestamp(),
        "verdict": "no_current_corpus_P_to_A_ch_bridge",
        "status": "exact_bridge_gap_certificate",
        "public_promotion_allowed": False,
        "theorem_statement": (
            "On the current corpus, the D10 calibration family of P closes only to the electroweak "
            "source-only and repaired electroweak observable surfaces. No theorem-grade map from those "
            "D10 descendants to theorem-grade physical Y_e, the charged determinant line, or the affine "
            "charged anchor A_ch is presently derived. Therefore the current corpus does not emit A_ch(P), "
            "g_e(P), or charged masses from P."
        ),
        "closed_d10_surface_from_P": {
            "input": "P",
            "forward_transmutation_certificate": {
                "artifact": d10_forward.get("artifact"),
                "artifact_ref": artifact_ref(D10_FORWARD_JSON),
                "object_id": d10_forward.get("object_id"),
                "theorem_name": forward_theorem.get("name"),
                "emitted_from_P": [
                    "alpha_U(P)",
                    "t_unified(P)",
                    "t_transmutation(P)",
                    "v_report(P)",
                ],
            },
            "target_free_repair_theorem": {
                "artifact": d10_repair.get("artifact"),
                "artifact_ref": artifact_ref(D10_REPAIR_JSON),
                "object_id": d10_repair.get("object_id"),
                "theorem_name": repair_theorem.get("name"),
                "emitted_from_P_branch": [
                    "alpha2_prime(P)",
                    "alphaY_prime(P)",
                    "MW_pole(P)",
                    "MZ_pole(P)",
                    "alpha_em_eff_inv(P)",
                    "sin2w_eff(P)",
                ],
            },
        },
        "charged_target_surface": {
            "required_affine_object": anchor.get("exact_missing_object"),
            "post_promotion_slot": post_promotion_slot.get("id"),
            "equivalent_descended_scalar": exact_descended_scalar.get("id"),
            "induced_formula_on_fill": anchor.get("induced_formula_on_fill"),
            "upstream_operator_gate": _as_dict(post_promotion_route, "operator_promotion_gate").get(
                "missing_theorem"
            ),
        },
        "bridge_absence_certificate": {
            "no_hidden_landing_on_physical_Y_e_from_P": True,
            "no_hidden_landing_on_charged_determinant_line_from_P": True,
            "no_hidden_landing_on_A_ch_from_P": True,
            "why_more_digits_of_P_do_not_help": (
                "The missing object is theorem-sized bridge content, not numerical precision inside the "
                "already-closed D10 electroweak lane."
            ),
        },
        "exact_missing_object": {
            "id": "d10_to_charged_affine_anchor_bridge",
            "smallest_honest_contract": (
                "one theorem-grade map from the D10 calibration descendants of P to the charged determinant "
                "line, or to a smaller intermediate object that canonically forces A_ch(P)"
            ),
            "exact_reduction_artifact": reduction.get("artifact"),
            "exact_reduction_artifact_ref": artifact_ref(P_TO_AFFINE_REDUCTION_JSON),
            "exact_smallest_bridge_target": reduction.get("exact_smallest_bridge_target", {}).get("id"),
            "admissible_landings": [
                "theorem_grade_physical_Y_e(P)",
                "charged_determinant_line_section(P)",
                "charged_absolute_anchor_A_ch(P)",
            ],
        },
        "forbidden_current_outputs": [
            "A_ch(P)",
            "g_e(P)",
            "m_e(P)",
            "m_mu(P)",
            "m_tau(P)",
        ],
        "future_symbolic_surface_if_bridge_closes": {
            "A_ch(P)": "mu_phys(Y_e(P)) = (1/3) * log(det(Y_e(P)))",
            "g_e(P)": "exp(A_ch(P))",
        },
        "hard_rejections": [
            "measured charged masses",
            "same-family affine anchor reuse",
            "common-scale tuning against centered charged logs",
            "assuming D10 electroweak closure already implies charged determinant data",
        ],
        "notes": [
            "This artifact does not deny that a future universal P-driven charged-mass law may exist; it records that no such bridge is presently derived.",
            "The charged affine anchor remains downstream of the post-promotion uncentered trace lift rather than a visible D10 readout today.",
            "This is the exact local result needed before spending Pro cycles on a Step-3 bridge prompt.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the current-corpus D10-to-charged affine bridge no-go.")
    parser.add_argument("--d10-forward", default=str(D10_FORWARD_JSON))
    parser.add_argument("--d10-repair", default=str(D10_REPAIR_JSON))
    parser.add_argument("--anchor-section", default=str(ANCHOR_SECTION_JSON))
    parser.add_argument("--post-promotion-route", default=str(POST_PROMOTION_ROUTE_JSON))
    parser.add_argument("--reduction", default=str(P_TO_AFFINE_REDUCTION_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    artifact = build_artifact(
        load_json(Path(args.d10_forward)),
        load_json(Path(args.d10_repair)),
        load_json(Path(args.anchor_section)),
        load_json(Path(args.post_promotion_route)),
        load_json(Path(args.reduction)),
    )

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
