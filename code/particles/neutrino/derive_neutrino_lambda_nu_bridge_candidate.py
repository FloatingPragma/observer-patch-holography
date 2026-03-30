#!/usr/bin/env python3
"""Record the sharpest current local bridge candidate for lambda_nu."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
SCALE_ANCHOR = ROOT / "particles" / "runs" / "neutrino" / "neutrino_scale_anchor.json"
COMPARE_FIT = ROOT / "particles" / "runs" / "neutrino" / "neutrino_compare_only_scale_fit.json"
SCALAR_EVALUATOR = ROOT / "particles" / "runs" / "neutrino" / "majorana_overlap_defect_scalar_evaluator.json"
THEOREM_OBJECT = ROOT / "particles" / "runs" / "neutrino" / "neutrino_weighted_cycle_theorem_object.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "neutrino" / "neutrino_lambda_nu_bridge_candidate.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def build_payload(
    *,
    scale_anchor: dict[str, Any],
    compare_fit: dict[str, Any],
    scalar_evaluator: dict[str, Any],
    theorem_object: dict[str, Any],
) -> dict[str, Any]:
    m_star_eV = float(scale_anchor["anchors"]["m_star_gev"]) * 1.0e9
    lambda_cmp = float(compare_fit["fits"]["weighted_least_squares"]["lambda_nu"])
    bridge_factor = lambda_cmp / m_star_eV
    gamma = float(theorem_object["live_inputs"]["gamma"])
    ratio_hat = float(theorem_object["live_outputs"]["dimensionless_ratio_dm21_over_dm32"])
    lambda_closed_form = gamma / (ratio_hat ** 0.5)
    delta_hat = compare_fit["scale_free_branch"]["delta_hat_m_sq_eV2"]
    masses_hat = compare_fit["scale_free_branch"]["m_hat_eV"]
    residuals = {
        "21": (
            lambda_closed_form**2 * float(delta_hat["21"])
            - float(compare_fit["reference_central_values"]["delta_m21_sq_eV2"])
        )
        / float(compare_fit["reference_central_values"]["delta_m21_sq_sigma_eV2"]),
        "32": (
            lambda_closed_form**2 * float(delta_hat["32"])
            - float(compare_fit["reference_central_values"]["delta_m32_sq_eV2"])
        )
        / float(compare_fit["reference_central_values"]["delta_m32_sq_sigma_eV2"]),
    }
    phase_clause = scalar_evaluator["strictly_smaller_missing_clause_if_not_closed"]
    normalizer_id = scalar_evaluator["attachment_normalizer_candidate_id"]
    return {
        "artifact": "oph_neutrino_lambda_nu_bridge_candidate",
        "generated_utc": _timestamp(),
        "status": "candidate_interface_not_promoted",
        "public_promotion_allowed": False,
        "current_dimensionless_theorem_object": theorem_object["theorem_object"]["name"],
        "dimensionful_anchor": {
            "name": "m_star",
            "formula": scale_anchor["anchors"]["m_star_formula"],
            "m_star_gev": scale_anchor["anchors"]["m_star_gev"],
            "m_star_eV": m_star_eV,
        },
        "current_candidate_interface_artifact": "oph_majorana_overlap_defect_scalar_evaluator",
        "exact_next_theorem_object": normalizer_id,
        "strictly_smaller_missing_clause": phase_clause,
        "bridge_ansatz": "lambda_nu = m_star_eV * F_nu",
        "where_F_nu_should_come_from": (
            "A normalized same-label overlap-defect weight section built on top of the finite-angle "
            "Majorana overlap-defect scalar evaluator on the selected weighted-cycle / midpoint branch."
        ),
        "bridge_interface_theorem_stack": [
            {
                "id": phase_clause,
                "status": "open_clause",
                "role": "selector-centered phase-cocycle triviality gate beneath the attachment route",
            },
            {
                "id": normalizer_id,
                "status": "candidate_only",
                "role": "normalized attachment section that converts the centered edge-norm scalar data into the positive bridge factor F_nu",
            },
            {
                "id": scalar_evaluator["remaining_theorem_object"],
                "status": scalar_evaluator["proof_status"],
                "role": "parent finite-angle centered edge-norm theorem furnishing the scalar side of the attachment route",
            },
            {
                "id": "neutrino_weighted_cycle_absolute_attachment",
                "status": "open",
                "role": "absolute-attachment theorem upgrading the weighted-cycle normal form to absolute masses and splittings",
            },
        ],
        "why_this_is_the_sharpest_local_candidate": [
            "The weighted-cycle theorem object closes the dimensionless branch only.",
            "The scale anchor closes only the isotropic dimensionful amplitude m_star.",
            "The scalar evaluator is the only current local object family already sitting at the interface between those two emitted sectors.",
        ],
        "compare_only_bridge_factor": {
            "lambda_nu_weighted_fit": lambda_cmp,
            "F_nu_star": bridge_factor,
            "interpretation": "If the bridge ansatz is correct, the theorem-grade evaluator would have to emit this positive dimensionless factor on the current branch.",
        },
        "target_free_closed_form_candidates": [
            {
                "name": "gamma_over_sqrt_ratio_hat",
                "status": "exactly_refuted_as_theorem_grade_absolute_scale_law",
                "formula": "lambda_nu = gamma / sqrt(Delta_hat_21 / Delta_hat_32)",
                "equivalent_statement": "Delta m21^2 = gamma^2 * Delta_hat_32",
                "lambda_nu": lambda_closed_form,
                "masses_eV": [lambda_closed_form * float(val) for val in masses_hat],
                "delta_m_sq_eV2": {
                    key: lambda_closed_form**2 * float(val) for key, val in delta_hat.items()
                },
                "residual_sigma": residuals,
                "proof_obstruction": "positive_rescaling_nonidentifiability",
                "why_refuted": (
                    "gamma and Delta_hat_21 / Delta_hat_32 are orbit invariants on the positive-rescaling family, so no theorem-grade absolute-scale law can be a function of those invariants alone."
                ),
            }
        ],
        "ruled_out_shortcuts": {
            "lambda_nu_equals_m_star_eV": {
                "status": "ruled_out",
                "source_artifact": "neutrino_weighted_cycle_absolute_amplitude_bridge.json",
            },
            "trace_determinant_minor_repackages_of_isotropic_majorana_matrix": {
                "status": "not_weighted_cycle_sensitive",
                "reason": "They reconstruct the isotropic amplitude sector again rather than the weighted-cycle branch attachment.",
            },
            "edge_amplitude_only_formulas": {
                "status": "ruled_out_local_shape_only",
                "reason": "Current edge-amplitude invariants are far too small and do not carry the required bridge factor.",
            },
        },
        "next_theorem_if_this_route_is_right": {
            "id": normalizer_id,
            "current_status": "candidate_only",
            "promotion_gate": scalar_evaluator["phase_cocycle_triviality_candidate_id"],
            "smallest_missing_clause": phase_clause,
        },
        "notes": [
            "This bridge candidate does not claim lambda_nu is already emitted.",
            "It packages the strongest current local interface between the emitted D10 amplitude scale and the emitted weighted-cycle theorem object.",
            "The closed-form gamma-over-sqrt-ratio numerology is retained only as a refuted compare-only audit target; it is incompatible with the exact positive-rescaling no-go.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the neutrino lambda_nu bridge candidate artifact.")
    parser.add_argument("--scale-anchor", default=str(SCALE_ANCHOR))
    parser.add_argument("--compare-fit", default=str(COMPARE_FIT))
    parser.add_argument("--scalar-evaluator", default=str(SCALAR_EVALUATOR))
    parser.add_argument("--theorem-object", default=str(THEOREM_OBJECT))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    payload = build_payload(
        scale_anchor=_load_json(Path(args.scale_anchor)),
        compare_fit=_load_json(Path(args.compare_fit)),
        scalar_evaluator=_load_json(Path(args.scalar_evaluator)),
        theorem_object=_load_json(Path(args.theorem_object)),
    )
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
