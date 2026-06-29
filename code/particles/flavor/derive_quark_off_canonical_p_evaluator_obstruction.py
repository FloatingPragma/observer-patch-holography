#!/usr/bin/env python3
"""Emit the closure certificate for theorem-grade off-canonical quark P transport.

Chain role: make the current failure to close issue #212 executable. The
selected public quark theorem is closed, but the arbitrary-P moving evaluator
still needs a theorem-grade sigma lift and odd scale law.

Mathematics: the ordered even family reduces exactly to a free positive pair
``(sigma_u, sigma_d)``. Current edge statistics provide candidate coefficients
for that pair, but the artifacts do not derive those coefficients. Likewise,
the odd response closes a direction ``Xi_q`` but leaves its scalar
normalization ``kappa`` as a ray coordinate. The final obstruction is stronger
than "not found": the present corpus admits canonical-preserving off-canonical
counterfamilies with different masses and odd payloads.

Output: a machine-readable non-closure witness with explicit numerical
nonuniqueness data.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_EDGE = ROOT / "particles" / "runs" / "flavor" / "quark_edge_statistics_spread_candidate.json"
DEFAULT_SPREAD = ROOT / "particles" / "runs" / "flavor" / "quark_spread_map.json"
DEFAULT_ODD = ROOT / "particles" / "runs" / "flavor" / "quark_odd_response_law.json"
DEFAULT_PUBLIC_YUKAWA = ROOT / "particles" / "runs" / "flavor" / "quark_public_exact_yukawa_end_to_end_theorem.json"
DEFAULT_PUBLIC_SOURCE = ROOT / "particles" / "runs" / "flavor" / "quark_d12_public_source_payload.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_off_canonical_p_evaluator_obstruction.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _scale_vector(values: Any, scale: float) -> Any:
    if isinstance(values, list) and values and isinstance(values[0], list):
        return [[float(item) * scale for item in row] for row in values]
    return [float(item) * scale for item in values]


def build_artifact(
    edge_candidate: dict[str, Any],
    spread_map: dict[str, Any],
    odd_response: dict[str, Any],
    public_yukawa: dict[str, Any],
    public_source_payload: dict[str, Any] | None = None,
) -> dict[str, Any]:
    edge_inputs = dict(edge_candidate["ordered_family_inputs"])
    edge_stats = dict(edge_candidate["edge_statistics_inputs"]["shared_charged_suppression_seed"])
    rho = float(edge_inputs["rho_ord"])
    x2 = float(edge_inputs["x2"])
    delta21 = float(edge_inputs["delta21"])
    s13 = float(edge_stats["S_13"])
    s23 = float(edge_stats["S_23"])
    candidate_sigmas = dict(edge_candidate["candidate_sigmas"])
    active_sigma_u = float(spread_map["sigma_u_total_log_per_side"])
    active_sigma_d = float(spread_map["sigma_d_total_log_per_side"])
    used_coeff_u = rho / (1.0 + rho)
    used_coeff_d = 1.0 / (2.0 * (1.0 + rho - x2 * x2))
    fitted_coeff_u = (active_sigma_u - s13) / delta21
    fitted_coeff_d = (active_sigma_d - s23) / delta21

    kappa_unit = float(odd_response.get("kappa", 1.0))
    kappa_alt = 0.5 * min(float(odd_response["kappa_max_admissible"]), max(kappa_unit, 1.0))
    kappa_witnesses = []
    for kappa in (0.0, kappa_alt, kappa_unit):
        kappa_witnesses.append(
            {
                "kappa": kappa,
                "Delta_logD_left_q": _scale_vector(odd_response["Delta_logD_left_q"], kappa / kappa_unit),
                "Delta_logD_right_q": _scale_vector(odd_response["Delta_logD_right_q"], kappa / kappa_unit),
                "Delta_S_q": _scale_vector(odd_response["Delta_S_q"], kappa / kappa_unit),
                "Delta_Phi_q": _scale_vector(odd_response["Delta_Phi_q"], kappa / kappa_unit),
            }
        )

    selected_sigma_datum = public_yukawa["descended_physical_sigma_datum"]
    sigma_u_selected = float(selected_sigma_datum["sigma_u"])
    sigma_d_selected = float(selected_sigma_datum["sigma_d"])
    canonical_even_countermodels = [
        {
            "name": "constant_selected_class_extension",
            "sigma_u(P)": "sigma_u_star",
            "sigma_d(P)": "sigma_d_star",
            "canonical_preservation_condition": "trivial",
            "interpretation": (
                "Extends the selected-class exact sigma pair off-canonical without using edge statistics. "
                "It preserves the canonical point but carries no source-emitted off-canonical response."
            ),
        },
        {
            "name": "edge_candidate_response_extension",
            "sigma_u(P)": "S_13(P) + (rho(P)/(1 + rho(P))) * delta21(P)",
            "sigma_d(P)": "S_23(P) + (1/(2*(1 + rho(P) - x2(P)^2))) * delta21(P)",
            "canonical_preservation_condition": (
                "requires an extra normalization/correction because the emitted candidate does not equal "
                "the active selected-class sigma pair at the current canonical point"
            ),
            "canonical_residuals": {
                "sigma_u_candidate_minus_selected": float(candidate_sigmas["sigma_u_total_log_per_side"])
                - sigma_u_selected,
                "sigma_d_candidate_minus_selected": float(candidate_sigmas["sigma_d_total_log_per_side"])
                - sigma_d_selected,
            },
            "interpretation": (
                "Uses the live edge-shaped candidate coefficients. It is source-shaped but not theorem-selected, "
                "and it fails selected-class exactness without an additional correction law."
            ),
        },
        {
            "name": "canonical_vanishing_free_perturbation_family",
            "sigma_u(P)": "sigma_u_star + epsilon_u * h(P)",
            "sigma_d(P)": "sigma_d_star + epsilon_d * h(P)",
            "canonical_preservation_condition": "h(P_star)=0 and sigmas remain positive",
            "interpretation": (
                "For any nonzero source-unfixed h and small epsilons, this preserves every selected-class "
                "artifact while producing different off-canonical masses. Current artifacts provide no rule "
                "that excludes it."
            ),
        },
    ]
    odd_countermodels = [
        {
            "name": "zero_odd_extension",
            "kappa(P)": "0",
            "canonical_preservation_condition": "admissible because 0 <= kappa <= kappa_max",
            "interpretation": "Keeps the odd ray identically silent off-canonical.",
        },
        {
            "name": "canonical_vanishing_odd_response_family",
            "kappa(P)": "epsilon * h(P)",
            "canonical_preservation_condition": "h(P_star)=0 and 0 <= epsilon*h(P) <= kappa_max(P)",
            "interpretation": (
                "Turns on a nonzero odd payload away from the selected point while preserving the canonical "
                "selected-class artifacts. Current artifacts do not emit epsilon or h."
            ),
        },
        {
            "name": "unit_representative_extension",
            "kappa(P)": "1",
            "canonical_preservation_condition": (
                "only a coordinate representative unless a physical normalization theorem identifies it"
            ),
            "interpretation": "The current odd artifact records this as a unit representative, not a value law.",
        },
    ]

    return {
        "artifact": "oph_quark_off_canonical_p_evaluator_obstruction",
        "generated_utc": _timestamp(),
        "scope": "arbitrary_P_off_canonical_quark_evaluator",
        "proof_status": "off_canonical_p_evaluator_underdetermined_current_corpus",
        "theorem_grade_closure": False,
        "public_promotion_allowed": False,
        "input_artifacts": {
            "edge_statistics_candidate": edge_candidate.get("artifact"),
            "spread_map": spread_map.get("artifact"),
            "odd_response": odd_response.get("artifact"),
            "public_exact_yukawa_theorem": public_yukawa.get("artifact"),
            "public_source_payload": None if public_source_payload is None else public_source_payload.get("artifact"),
        },
        "selected_public_exact_surface": {
            "closed": public_yukawa.get("proof_status") == "closed_source_only_public_exact_yukawa_end_to_end_theorem",
            "conditional_downstream_closed": bool(public_yukawa.get("downstream_algebra_closed")),
            "source_only_sigma_emitted": bool(public_yukawa.get("source_only_sigma_emitted")),
            "proof_status": public_yukawa.get("proof_status"),
            "theorem_scope": public_yukawa.get("theorem_scope"),
            "minimal_exact_blocker_set": public_yukawa.get("minimal_exact_blocker_set"),
            "selected_by": public_yukawa["selected_public_physical_frame_class"].get("selected_by"),
        },
        "selected_public_pure_B_payload": {
            "closed": public_source_payload is not None
            and public_source_payload.get("proof_status") == "closed_public_selected_class_pure_B_source_payload",
            "proof_status": None if public_source_payload is None else public_source_payload.get("proof_status"),
            "off_canonical_promotion_allowed": None
            if public_source_payload is None
            else public_source_payload.get("off_canonical_promotion_allowed"),
        },
        "lane_closure_verdict": {
            "issue_212_acceptance_met": False,
            "closure_kind": "hard_no_go_current_corpus",
            "closed_theorem_grade_surface": "selected_public_physical_quark_frame_class_only",
            "blocked_surface": "arbitrary_P_off_canonical_quark_evaluator",
            "repo_runtime_browser_contract_status": (
                "shared candidate/evidence surface only; theorem-grade off-canonical promotion disallowed"
            ),
            "why_this_is_final_for_current_corpus": (
                "The attached source artifacts admit distinct canonical-preserving off-canonical completions "
                "for both the even sigma pair and the odd kappa payload. A direct arbitrary-P evaluator would "
                "therefore choose extra data not emitted by the present theorem corpus."
            ),
        },
        "edge_statistics_sigma_lift_obstruction": {
            "missing_object": "oph_edge_statistics_sigma_lift",
            "bridge_status": edge_candidate.get("bridge_status"),
            "free_sigma_branch_formulas": {
                "E_u_log": "sigma_u * v_u",
                "E_d_log": "sigma_d * v_d",
                "gamma21_u": "rho * sigma_u / (1 + rho)",
                "gamma32_u": "sigma_u / (1 + rho)",
                "gamma21_d": "sigma_d / (1 + rho)",
                "gamma32_d": "rho * sigma_d / (1 + rho)",
            },
            "candidate_formula_coefficients": {
                "sigma_u": "S_13 + c_u * delta21",
                "sigma_d": "S_23 + c_d * delta21",
                "c_u_used": used_coeff_u,
                "c_d_used": used_coeff_d,
            },
            "coefficients_that_match_active_closed_pair": {
                "c_u_fit": fitted_coeff_u,
                "c_d_fit": fitted_coeff_d,
            },
            "candidate_sigmas": {
                "sigma_u_total_log_per_side": float(candidate_sigmas["sigma_u_total_log_per_side"]),
                "sigma_d_total_log_per_side": float(candidate_sigmas["sigma_d_total_log_per_side"]),
            },
            "active_closed_sigmas": {
                "sigma_u_total_log_per_side": active_sigma_u,
                "sigma_d_total_log_per_side": active_sigma_d,
            },
            "residuals": {
                "sigma_u_candidate_minus_active": float(candidate_sigmas["sigma_u_total_log_per_side"]) - active_sigma_u,
                "sigma_d_candidate_minus_active": float(candidate_sigmas["sigma_d_total_log_per_side"]) - active_sigma_d,
                "c_u_used_minus_fit": used_coeff_u - fitted_coeff_u,
                "c_d_used_minus_fit": used_coeff_d - fitted_coeff_d,
            },
            "underdetermination_statement": (
                "The ordered even family is exact once sigma_u and sigma_d are fixed, but current edge statistics "
                "do not derive the coefficients selecting sigma_u and sigma_d from S_13, S_23, delta21, rho, and x2."
            ),
        },
        "odd_response_scale_obstruction": {
            "missing_object": "oph_off_canonical_odd_response_kappa_value_law",
            "proof_status": odd_response.get("proof_status"),
            "lift_parameterization_kind": odd_response.get("lift_parameterization_kind"),
            "kappa_unit_representative": kappa_unit,
            "kappa_max_admissible": odd_response.get("kappa_max_admissible"),
            "admissible_kappa_witnesses": kappa_witnesses,
            "underdetermination_statement": (
                "The current odd response closes a signed projector direction but leaves a one-scalar kappa ray. "
                "Different admissible kappa values change the odd payload while satisfying the present structural guards."
            ),
        },
        "formal_countermodel_witness": {
            "canonical_point": {
                "name": "P_star",
                "sigma_u_star": sigma_u_selected,
                "sigma_d_star": sigma_d_selected,
            },
            "even_sigma_countermodels": canonical_even_countermodels,
            "odd_kappa_countermodels": odd_countermodels,
            "arbitrary_frame_transport_countermodel": {
                "name": "selected_class_identity_vs_free_off_canonical_frame_lift",
                "frame_A": (
                    "Keep the selected public physical frame class fixed at P_star and extend no additional "
                    "off-canonical labels."
                ),
                "frame_B": (
                    "Attach any nonzero off-canonical frame label f(P)=h(P)L with h(P_star)=0 and with L "
                    "in a source-unclassified public quark-frame direction."
                ),
                "why_both_fit_current_corpus": (
                    "All selected-class artifacts inspect only P_star and the selected public class. The corpus "
                    "does not emit a classifier that transports or rejects L for arbitrary P."
                ),
            },
            "nonidentifiability_result": (
                "Because these countermodels agree on every selected-class exact artifact but disagree away from "
                "P_star, sigma_u(P), sigma_d(P), kappa(P), and public frame transport are not identifiable from "
                "the current source basis."
            ),
        },
        "closure_blockers": [
            "edge_statistics_sigma_lift_not_derived",
            "off_canonical_odd_response_kappa_value_law_not_derived",
            "arbitrary_P_public_quark_frame_classification_not_derived",
        ],
        "required_constructive_objects": [
            "oph_edge_statistics_sigma_lift",
            "oph_off_canonical_odd_response_kappa_value_law",
            "oph_arbitrary_P_public_quark_frame_transport_classification",
        ],
        "notes": [
            "The selected public exact quark theorem is closed and can be used as a theorem-grade static evaluator.",
            "The selected-public-class pure-B source payload is now data-bearing when present.",
            "Those selected-class closures do not supply an arbitrary-P off-canonical sigma or odd-scale family.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the quark off-canonical P evaluator obstruction artifact.")
    parser.add_argument("--edge", default=str(DEFAULT_EDGE))
    parser.add_argument("--spread", default=str(DEFAULT_SPREAD))
    parser.add_argument("--odd", default=str(DEFAULT_ODD))
    parser.add_argument("--public-yukawa", default=str(DEFAULT_PUBLIC_YUKAWA))
    parser.add_argument("--public-source-payload", default=str(DEFAULT_PUBLIC_SOURCE))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    public_source_path = Path(args.public_source_payload)
    payload = build_artifact(
        _load_json(Path(args.edge)),
        _load_json(Path(args.spread)),
        _load_json(Path(args.odd)),
        _load_json(Path(args.public_yukawa)),
        _load_json(public_source_path) if public_source_path.exists() else None,
    )
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
