#!/usr/bin/env python3
"""Shrink the charged post-promotion scalar slot to one physical affine scalar.

This emitter also records the exact smaller forcing object beneath the
descended scalar: the physical identity-mode equalizer on common physical
fibers.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

from charged_absolute_route_common import (
    ANCHOR_SECTION_JSON,
    DETERMINANT_LINE_JSON,
    PHYSICAL_EQUALIZER_JSON,
    TRACE_LIFT_COCYCLE_JSON,
    TRACE_LIFT_JSON,
    TRACE_LIFT_PHYSICAL_DESCENT_JSON,
    artifact_ref,
    load_json,
    trace_lift_physical_descent_contract,
)


DEFAULT_OUT = TRACE_LIFT_PHYSICAL_DESCENT_JSON


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(
    trace_lift: dict,
    cocycle: dict,
    determinant_line: dict,
    anchor: dict,
    physical_equalizer: dict,
) -> dict:
    descent_contract = trace_lift_physical_descent_contract()
    descent_contract["smaller_forcing_object"] = physical_equalizer.get(
        "exact_smaller_forcing_object",
        "charged_physical_identity_mode_equalizer",
    )
    descent_contract["smaller_forcing_rule"] = physical_equalizer.get(
        "equalizer_contract",
        {},
    ).get("fiber_rule")

    return {
        "artifact": "oph_charged_mu_physical_descent_reduction",
        "generated_utc": _timestamp(),
        "status": "exact_descended_scalar_reduction",
        "public_promotion_allowed": False,
        "larger_missing_object": trace_lift.get("exact_missing_object"),
        "larger_missing_object_artifact_ref": artifact_ref(TRACE_LIFT_JSON),
        "exact_smaller_missing_object": descent_contract["exact_smaller_missing_object"],
        "exact_smaller_missing_object_kind": "single_affine_scalar_on_theorem_grade_physical_Y_e",
        "exact_smaller_forcing_object": {
            "id": physical_equalizer.get(
                "exact_smaller_forcing_object",
                "charged_physical_identity_mode_equalizer",
            ),
            "artifact": physical_equalizer.get("artifact"),
            "artifact_ref": artifact_ref(PHYSICAL_EQUALIZER_JSON),
            "kind": physical_equalizer.get(
                "exact_smaller_forcing_object_kind",
                "fiberwise_zero_cocycle_certificate_on_theorem_grade_physical_Y_e",
            ),
        },
        "descent_contract": descent_contract,
        "input_scalarization_artifact_ref": artifact_ref(TRACE_LIFT_COCYCLE_JSON),
        "reduction_theorem": {
            "id": "charged_refinement_stable_mu_descends_to_physical_scalar",
            "statement": (
                "Assume theorem-grade C_hat_e is promoted on theorem-grade physical Y_e and that an "
                "uncentered lift C_tilde_e = C_hat_e + mu I is refinement-stable on that physical surface. "
                "Then the refinement identity-mode cocycle vanishes on every pair of refinement representatives "
                "of the same physical Y_e, so mu descends to a unique physical scalar mu_phys(Y_e). Conversely, "
                "any such mu_phys defines the uncentered lift, determinant-line section, and affine anchor by "
                "C_tilde_e = C_hat_e + mu_phys I, s_det = 3 mu_phys, and A_ch = mu_phys. This fiberwise "
                "zero-cocycle certificate is the exact smaller forcing object beneath mu_phys(Y_e)."
            ),
        },
        "P_threaded_source_exactness_reformulation": {
            "id": "d10_selected_charged_source_exactness_forces_physical_affine_scalar_mu",
            "statement": (
                "Assume OPH axioms, the realized D9 branch, the D10 calibration sector with external input P, "
                "and a theorem-grade source landing P -> Y_e(P) on theorem-grade physical charged data. "
                "Assume theorem-grade C_hat_e is promoted on that physical Y_e(P) surface and that an admissible "
                "refinement-stable uncentered lift C_tilde_e(r;P) = C_hat_e(Y_e(P)) + mu(r;P) I exists on "
                "refinement representatives of the same physical Y_e(P). Then the physical identity-mode equalizer "
                "forces delta(r,r';P) = 0 on each same-fiber pair, so mu(r;P) descends uniquely to "
                "mu_phys(Y_e(P)). Equivalently, once the D10-selected charged source is exact on theorem-grade "
                "physical Y_e(P), the affine scalar is forced by the existing descent package rather than remaining "
                "an extra independent bridge-sized object."
            ),
            "bridge_role": (
                "This reformulates the charged P-bridge burden as a theorem-grade landing on physical Y_e(P), "
                "after which the post-promotion affine scalar follows from the descent reduction."
            ),
        },
        "forced_vanishing": {
            "refinement_pair_rule_before_descent": cocycle.get("scalar_cocycle_contract", {}).get(
                "pairwise_difference_rule"
            ),
            "on_same_physical_Y_e": descent_contract["forced_refinement_identity_mode"],
            "reason": "refinement stability identifies all refinement representatives of the same physical operator",
        },
        "equivalent_presentations_on_fill": {
            "descended_scalar": "mu_phys(Y_e)",
            "uncentered_trace_lift": {
                "artifact": trace_lift.get("artifact"),
                "artifact_ref": artifact_ref(TRACE_LIFT_JSON),
                "formula": "C_tilde_e(Y_e) = C_hat_e(Y_e) + mu_phys(Y_e) I",
            },
            "determinant_line_section": {
                "artifact": determinant_line.get("artifact"),
                "artifact_ref": artifact_ref(DETERMINANT_LINE_JSON),
                "formula": "s_det(Y_e) = 3 * mu_phys(Y_e)",
            },
            "affine_anchor": {
                "artifact": anchor.get("artifact"),
                "artifact_ref": artifact_ref(ANCHOR_SECTION_JSON),
                "formula": "A_ch(Y_e) = mu_phys(Y_e)",
            },
        },
        "why_this_is_smaller": [
            "The family-wise scalar primitive mu(r) is only needed before refinement stability is imposed.",
            "The exact forcing object beneath mu_phys(Y_e) is the physical identity-mode equalizer delta(r,r') = 0 on each physical fiber.",
            "Once that equalizer holds, mu_phys(Y_e) is only the common value of mu(r) on the fiber, not an extra independent theorem slot.",
        ],
        "do_not_claim_now": [
            "theorem-grade mu_phys(Y_e) on the live corpus",
            "theorem-grade determinant-line section on the live corpus",
            "theorem-grade A_ch on the live corpus",
        ],
        "notes": [
            "This sharpens the post-promotion frontier only under the same refinement-stability contract already required by the lift scaffold.",
            "It does not promote current-corpus closure or bypass the upstream C_hat_e promotion theorem.",
            "The exact smaller forcing object beneath the descended scalar is emitted separately as the physical identity-mode equalizer.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build the charged physical-scalar descent reduction beneath the uncentered trace lift."
    )
    parser.add_argument("--trace-lift", default=str(TRACE_LIFT_JSON))
    parser.add_argument("--cocycle-reduction", default=str(TRACE_LIFT_COCYCLE_JSON))
    parser.add_argument("--determinant-line", default=str(DETERMINANT_LINE_JSON))
    parser.add_argument("--anchor-section", default=str(ANCHOR_SECTION_JSON))
    parser.add_argument("--physical-equalizer", default=str(PHYSICAL_EQUALIZER_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    artifact = build_artifact(
        load_json(Path(args.trace_lift)),
        load_json(Path(args.cocycle_reduction)),
        load_json(Path(args.determinant_line)),
        load_json(Path(args.anchor_section)),
        load_json(Path(args.physical_equalizer)),
    )

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
