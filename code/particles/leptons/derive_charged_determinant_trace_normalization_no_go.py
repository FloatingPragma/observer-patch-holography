#!/usr/bin/env python3
"""Prove the current-corpus additive normalization no-go for charged determinant descent."""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path

from charged_absolute_route_common import (
    DETERMINANT_CHARACTER_FRONTIER_JSON,
    DETERMINANT_TRACE_NORMALIZATION_NO_GO_JSON,
    PHYSICAL_EQUALIZER_JSON,
    TRACE_LIFT_PHYSICAL_DESCENT_JSON,
    artifact_ref,
    load_json,
)


DEFAULT_OUT = DETERMINANT_TRACE_NORMALIZATION_NO_GO_JSON


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(frontier: dict, physical_equalizer: dict, physical_descent: dict) -> dict:
    q_values = frontier.get("live_same_label_q_e", {})
    log_q = {key: math.log(float(value)) for key, value in q_values.items()}
    return {
        "artifact": "oph_charged_determinant_trace_normalization_no_go",
        "generated_utc": _timestamp(),
        "proof_status": "source_character_defined_additive_normalization_open",
        "public_promotion_allowed": False,
        "scope": "current_corpus_relative_to_populated_q_e_and_post_promotion_trace_lift_scaffold",
        "source_side_character": {
            "formula": "S_M(P) = sum_e M_e^ch log q_e(P)",
            "live_same_label_q_e": q_values,
            "live_same_label_log_q_e": log_q,
            "readback_artifact": frontier.get("artifact"),
            "readback_artifact_ref": artifact_ref(DETERMINANT_CHARACTER_FRONTIER_JSON),
        },
        "post_promotion_trace_lift_scaffold": {
            "equalizer_artifact": physical_equalizer.get("artifact"),
            "equalizer_artifact_ref": artifact_ref(PHYSICAL_EQUALIZER_JSON),
            "equalizer_contract": "mu(r') - mu(r) = 0 on a fixed physical charged fiber",
            "physical_descent_artifact": physical_descent.get("artifact"),
            "physical_descent_artifact_ref": artifact_ref(TRACE_LIFT_PHYSICAL_DESCENT_JSON),
            "trace_lift_formula_on_fill": "C_tilde_e(r) = C_hat_e(r) + mu(r) I",
            "determinant_formula_on_fill": "s_det(r) = 3 * mu(r)",
        },
        "unresolved_scalar_defect": {
            "id": "charged_determinant_normalization_defect",
            "physical_formula": "N_det(P) = s_det(P) - sum_e M_e^ch log q_e(P)",
            "trace_lift_formula": "Theta(r) = 3 * mu(r) - sum_e M_e^ch log q_e(r)",
            "meaning": (
                "The current corpus fixes the source-side determinant character and the "
                "post-promotion identity-mode scalar only up to their relative additive origin."
            ),
        },
        "affine_shift_witness": {
            "parameter": "kappa in R",
            "shift_family": [
                "mu -> mu + kappa",
                "C_tilde_e -> C_tilde_e + kappa I",
                "s_det -> s_det + 3 kappa",
                "A_ch -> A_ch + kappa",
            ],
            "preserved_current_contracts": [
                "physical identity-mode equalizer on a fixed physical charged fiber",
                "common-shift covariance of the charged affine scalar",
                "populated same-label q_e readback",
                "source-side determinant character S_M",
            ],
            "effect_on_defect": "Theta -> Theta + 3 kappa",
        },
        "no_go_theorem": {
            "id": "charged_determinant_trace_normalization_not_forced_by_current_corpus",
            "statement": (
                "The populated same-label readback fixes the source-side determinant character "
                "S_M = sum_e M_e^ch log q_e, and the post-promotion identity-mode equalizer fixes "
                "the fiberwise differences of mu. But the charged corpus does not force the "
                "cross-normalization 3 mu(r) = S_M(r), because the additive affine shift mu -> mu + "
                "kappa preserves the existing cocycle/equalizer contracts while leaving S_M fixed. "
                "Therefore one additive determinant normalization constant stays unresolved."
            ),
            "proof_skeleton": [
                "The live same-label q_e readback fixes S_M = sum_e M_e^ch log q_e on the realized support.",
                "The post-promotion equalizer constrains only pairwise differences mu(r') - mu(r) on one physical fiber.",
                "For any constant kappa, the shifted scalar mu^kappa = mu + kappa preserves those pairwise differences.",
                "The induced shifts s_det -> s_det + 3 kappa and A_ch -> A_ch + kappa preserve the required affine covariance.",
                "The source-side q_e readback and therefore S_M are unchanged under this shift.",
                "Hence the defect Theta = 3 mu - S_M is not fixed by current theorem-grade data.",
            ],
        },
        "exact_remaining_theorem": {
            "id": "charged_determinant_trace_lift_attachment",
            "statement": (
                "Attach the fixed determinant exponent vector M^ch sector-isolatedly to the charged "
                "L-H-E determinant channel and prove 3 mu(r) = sum_e M_e^ch log q_e(r) on that "
                "channel. Equivalently, force the determinant normalization defect N_det to vanish."
            ),
            "equivalent_fiber_form": (
                "On one physical charged fiber, this is equivalent to determinant-cocycle vanishing "
                "Delta_M(r,r') = 0 together with one basepoint normalization S_M(r_0) = 3 mu(r_0)."
            ),
        },
        "notes": [
            "This theorem states a relative underdetermination result for the charged determinant normalization.",
            "It isolates the missing scalar normalization theorem.",
            "The live same-label q_e readback is fully populated and theorem-relevant on the source side.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build the charged determinant trace-normalization no-go artifact."
    )
    parser.add_argument("--frontier", default=str(DETERMINANT_CHARACTER_FRONTIER_JSON))
    parser.add_argument("--physical-equalizer", default=str(PHYSICAL_EQUALIZER_JSON))
    parser.add_argument("--physical-descent", default=str(TRACE_LIFT_PHYSICAL_DESCENT_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    artifact = build_artifact(
        load_json(Path(args.frontier)),
        load_json(Path(args.physical_equalizer)),
        load_json(Path(args.physical_descent)),
    )

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
