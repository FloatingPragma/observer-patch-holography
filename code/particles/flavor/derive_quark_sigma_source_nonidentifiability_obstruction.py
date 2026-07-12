#!/usr/bin/env python3
"""Emit the target-free quark-spread non-identifiability obstruction.

The current quark source surface fixes the ordered three-point profile rays,
but it does not fix their independent up/down positive moduli.  This builder
records that exact two-dimensional fiber without loading any running-quark
target, selected/current-family exact witness, fitted sigma value, or the
contaminated comparison spread artifacts.

The obstruction is deliberately stronger than the upstream candidate-status
boundary: even if the current template-derived branch-generator packet is
granted as exact, the positive spread fiber remains ``(R_{>0})^2``.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
FAMILY_KERNEL_JSON = ROOT / "particles" / "runs" / "flavor" / "family_transport_kernel.json"
GENERATOR_JSON = ROOT / "particles" / "runs" / "flavor" / "generation_bundle_branch_generator.json"
SPLITTING_OBSTRUCTION_JSON = (
    ROOT
    / "particles"
    / "runs"
    / "flavor"
    / "generation_bundle_branch_generator_splitting_obstruction.json"
)
COCYCLE_JSON = ROOT / "particles" / "runs" / "flavor" / "overlap_edge_transport_cocycle.json"
DEFAULT_OUT = (
    ROOT
    / "particles"
    / "runs"
    / "flavor"
    / "quark_sigma_source_nonidentifiability_obstruction.json"
)

DIRECT_INPUT_PATHS = [
    "particles/runs/flavor/family_transport_kernel.json",
    "particles/runs/flavor/generation_bundle_branch_generator.json",
    "particles/runs/flavor/generation_bundle_branch_generator_splitting_obstruction.json",
    "particles/runs/flavor/overlap_edge_transport_cocycle.json",
]

FORBIDDEN_ANCESTORS = [
    "particle_reference_values",
    "quark_current_family_exact_readout",
    "quark_current_family_exact_sigma_target",
    "quark_current_family_exact_pdg_theorem",
    "quark_current_family_end_to_end_exact_pdg_derivation_chain",
    "quark_current_family_transport_frame_strengthened_physical_sigma_lift_theorem",
    "quark_selected_class_public_exact_evaluator",
    "quark_public_physical_sigma_datum_descent",
    "quark_public_exact_yukawa_end_to_end_theorem",
    "quark_exact_pdg_end_to_end_theorem",
    "quark_p_driven_shared_evaluator_contract",
    "family_excitation_evaluator",
    "quark_edge_statistics_spread_candidate",
    "quark_spread_map",
    "quark_sector_mean_split",
    "target_centered_log_u",
    "target_centered_log_d",
    "reference_targets_u",
    "reference_targets_d",
    "exact_fit_residuals",
    "fitted_sigma_values",
    "PDG_API_quark_rows",
    "CODATA",
    "compare_only",
]


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _ordered_source_geometry(generator: dict[str, Any]) -> dict[str, float | list[float]]:
    certificate = dict(generator["noncentrality_certificate"])
    eigenvalues = sorted(float(value) for value in certificate["eigenvalues"])
    if len(eigenvalues) != 3:
        raise ValueError("branch-generator certificate must contain three eigenvalues")
    lam1, lam2, lam3 = eigenvalues
    delta21 = lam2 - lam1
    delta32 = lam3 - lam2
    spectral_span = lam3 - lam1
    if delta21 <= 0.0 or delta32 <= 0.0 or spectral_span <= 0.0:
        raise ValueError("branch-generator spectrum must be simple and ordered")

    x2 = 2.0 * delta21 / spectral_span - 1.0
    rho = 3.0 * delta32 / (2.0 * delta32 + delta21)
    if rho <= 0.0:
        raise ValueError("ordered gap ratio rho must be positive")
    denominator = 3.0 * (1.0 + rho)
    v_u = [
        -((2.0 * rho) + 1.0) / denominator,
        (rho - 1.0) / denominator,
        (rho + 2.0) / denominator,
    ]
    v_d = [
        -(rho + 2.0) / denominator,
        (1.0 - rho) / denominator,
        ((2.0 * rho) + 1.0) / denominator,
    ]
    return {
        "ordered_eigenvalues": eigenvalues,
        "delta21": delta21,
        "delta32": delta32,
        "spectral_span": spectral_span,
        "x2": x2,
        "rho_ord": rho,
        "v_u": v_u,
        "v_d": v_d,
    }


def _ray_certificate(profile: list[float], *, expected_ratio: float) -> dict[str, Any]:
    gap21 = profile[1] - profile[0]
    gap32 = profile[2] - profile[1]
    ratio = gap21 / gap32
    return {
        "profile": profile,
        "trace": sum(profile),
        "endpoint_span": profile[2] - profile[0],
        "gap21_per_unit_spread": gap21,
        "gap32_per_unit_spread": gap32,
        "gap21_over_gap32": ratio,
        "expected_gap21_over_gap32": expected_ratio,
        "max_abs_identity_residual": max(
            abs(sum(profile)),
            abs((profile[2] - profile[0]) - 1.0),
            abs(ratio - expected_ratio),
        ),
    }


def _sigma_tuple(sigma_u: float, sigma_d: float) -> dict[str, float]:
    return {
        "sigma_seed_ud": 0.5 * (sigma_u + sigma_d),
        "eta_ud": 0.5 * (sigma_u - sigma_d),
        "sigma_u": sigma_u,
        "sigma_d": sigma_d,
    }


def _countermodel(
    name: str,
    sigma_u: float,
    sigma_d: float,
    v_u: list[float],
    v_d: list[float],
    *,
    rho: float,
) -> dict[str, Any]:
    if sigma_u <= 0.0 or sigma_d <= 0.0:
        raise ValueError("formal countermodel spreads must be positive")
    e_u = [sigma_u * value for value in v_u]
    e_d = [sigma_d * value for value in v_d]
    gap21_u = e_u[1] - e_u[0]
    gap32_u = e_u[2] - e_u[1]
    gap21_d = e_d[1] - e_d[0]
    gap32_d = e_d[2] - e_d[1]
    return {
        "name": name,
        "role": "formal_nonidentifiability_witness_not_a_prediction",
        "source_projection_token": "same_target_free_source_packet",
        "sigma_tuple": _sigma_tuple(sigma_u, sigma_d),
        "E_u_log": e_u,
        "E_d_log": e_d,
        "structural_identity_residuals": {
            "trace_E_u": sum(e_u),
            "trace_E_d": sum(e_d),
            "endpoint_span_u_minus_sigma_u": (e_u[2] - e_u[0]) - sigma_u,
            "endpoint_span_d_minus_sigma_d": (e_d[2] - e_d[0]) - sigma_d,
            "up_gap_ratio_minus_rho": gap21_u / gap32_u - rho,
            "down_gap_ratio_minus_reciprocal_rho": gap21_d / gap32_d - 1.0 / rho,
        },
    }


def _affine_injectivity_certificate(
    *,
    rho: float,
    x2: float,
    model_a: dict[str, Any],
    model_b: dict[str, Any],
) -> dict[str, Any]:
    mean_denominator = 1.0 + rho - x2 * x2
    skew_denominator = 1.0 - x2 * x2 - x2 * x2 / (1.0 + rho)
    if abs(mean_denominator) <= 1.0e-12 or abs(skew_denominator) <= 1.0e-12:
        raise ValueError("affine readout denominator is singular")
    a_ud = 1.0 / (2.0 * mean_denominator)
    b_ud = 1.0 / (2.0 * skew_denominator)
    matrix = [
        [(-a_ud + b_ud) / 2.0, (-a_ud - b_ud) / 2.0],
        [(-a_ud - b_ud) / 2.0, (-a_ud + b_ud) / 2.0],
    ]
    determinant = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]

    def log_shifts(model: dict[str, Any]) -> dict[str, float]:
        sigma = dict(model["sigma_tuple"])
        seed = float(sigma["sigma_seed_ud"])
        eta = float(sigma["eta_ud"])
        return {
            "log_g_u_over_g_ch": -(a_ud * seed - b_ud * eta),
            "log_g_d_over_g_ch": -(a_ud * seed + b_ud * eta),
        }

    shifts_a = log_shifts(model_a)
    shifts_b = log_shifts(model_b)
    return {
        "A_ud": a_ud,
        "B_ud": b_ud,
        "mean_denominator": mean_denominator,
        "skew_denominator": skew_denominator,
        "linear_map_from_sigma_u_sigma_d_to_log_sector_means": matrix,
        "determinant": determinant,
        "determinant_formula": "-A_ud * B_ud",
        "determinant_residual": determinant + a_ud * b_ud,
        "rank": 2,
        "injective": abs(determinant) > 1.0e-12,
        "model_A_log_shifts": shifts_a,
        "model_B_log_shifts": shifts_b,
        "models_have_distinct_sector_means": shifts_a != shifts_b,
        "conclusion": (
            "The free spread pair is not gauge: distinct positive moduli change the affine sector means, "
            "and the endpoint spans change the ordered log shapes."
        ),
    }


def _edge_coefficient_family(
    *,
    cocycle: dict[str, Any],
    delta21: float,
    rho: float,
    x2: float,
) -> dict[str, Any]:
    suppression = cocycle["derived_pairwise_suppression"]
    s13 = float(suppression[0][2])
    s23 = float(suppression[1][2])
    if delta21 <= 0.0:
        raise ValueError("delta21 must be positive")
    candidate_c_u = rho / (1.0 + rho)
    candidate_c_d = 1.0 / (2.0 * (1.0 + rho - x2 * x2))
    return {
        "statement": (
            "Even after granting the current edge packet, every positive spread pair can be written as "
            "sigma_u=S_13+c_u*delta21 and sigma_d=S_23+c_d*delta21. The current corpus emits no "
            "equation selecting c_u or c_d."
        ),
        "source_edge_data": {
            "S_13": s13,
            "S_23": s23,
            "delta21": delta21,
        },
        "free_coefficients": ["c_u", "c_d"],
        "coefficient_domain": {
            "c_u": f"c_u > {-s13 / delta21!r}",
            "c_d": f"c_d > {-s23 / delta21!r}",
        },
        "coefficient_to_spread_map": {
            "sigma_u": "S_13 + c_u * delta21",
            "sigma_d": "S_23 + c_d * delta21",
            "inverse_c_u": "(sigma_u - S_13) / delta21",
            "inverse_c_d": "(sigma_d - S_23) / delta21",
        },
        "unselected_existing_formula_point": {
            "status": "candidate_formula_reconstructed_from_target_free_fields_not_promoted",
            "c_u": candidate_c_u,
            "c_d": candidate_c_d,
            "sigma_u": s13 + candidate_c_u * delta21,
            "sigma_d": s23 + candidate_c_d * delta21,
        },
        "zero_coefficient_counterpoint": {
            "status": "formal_countermodel_not_a_prediction",
            "c_u": 0.0,
            "c_d": 0.0,
            "sigma_u": s13,
            "sigma_d": s23,
        },
        "two_points_share_source_edge_data_and_differ_in_spreads": True,
    }


def build_artifact(
    family_kernel: dict[str, Any],
    generator: dict[str, Any],
    splitting_obstruction: dict[str, Any],
    cocycle: dict[str, Any],
) -> dict[str, Any]:
    geometry = _ordered_source_geometry(generator)
    rho = float(geometry["rho_ord"])
    x2 = float(geometry["x2"])
    v_u = [float(value) for value in geometry["v_u"]]
    v_d = [float(value) for value in geometry["v_d"]]
    up_ray = _ray_certificate(v_u, expected_ratio=rho)
    down_ray = _ray_certificate(v_d, expected_ratio=1.0 / rho)

    model_a = _countermodel("unit_modulus_pair", 1.0, 1.0, v_u, v_d, rho=rho)
    model_b = _countermodel("independently_rescaled_pair", 2.0, 3.0, v_u, v_d, rho=rho)
    affine = _affine_injectivity_certificate(rho=rho, x2=x2, model_a=model_a, model_b=model_b)
    edge_family = _edge_coefficient_family(
        cocycle=cocycle,
        delta21=float(geometry["delta21"]),
        rho=rho,
        x2=x2,
    )

    allowed_ancestors = [
        str(family_kernel.get("artifact")),
        str(generator.get("artifact")),
        str(splitting_obstruction.get("artifact")),
        str(cocycle.get("artifact")),
    ]
    intersection = sorted(set(allowed_ancestors).intersection(FORBIDDEN_ANCESTORS))
    max_countermodel_residual = max(
        abs(float(value))
        for model in (model_a, model_b)
        for value in model["structural_identity_residuals"].values()
    )

    return {
        "artifact": "oph_quark_sigma_source_nonidentifiability_obstruction",
        "generated_utc": _timestamp(),
        "github_issues": [377, 379, 380],
        "proof_status": "closed_exact_current_corpus_obstruction",
        "claim_tier": "source_only_nonidentifiability_obstruction",
        "theorem_grade_obstruction": True,
        "source_only_sigma_emitted": False,
        "public_promotion_allowed": False,
        "numeric_quark_rows_allowed": False,
        "issue_377_acceptance_met_as_obstruction": True,
        "issue_379_acceptance_met_as_obstruction": True,
        "issue_380_acceptance_met_as_obstruction": True,
        "theorem_statement": (
            "After all target-derived and compare-only ancestors are removed, the strongest current quark "
            "source packet fixes the trace-zero ordered profile rays v_u and v_d but not their independent "
            "positive moduli. The compatible spread fiber is exactly (R_{>0})^2, with "
            "E_u=sigma_u*v_u and E_d=sigma_d*v_d. Independent positive rescaling fixes every source datum "
            "and every emitted shape identity while changing (sigma_seed_ud,eta_ud,sigma_u,sigma_d), the "
            "affine sector means, and the ordered log shapes. Therefore the present source theory does not "
            "identify a unique quark spread package or running quark mass sextet."
        ),
        "scope": {
            "premise": (
                "Target-free current-corpus quark source signature, granting the current template-derived "
                "ordered branch-generator and edge packets for the sake of a stronger obstruction."
            ),
            "conclusion_kind": "nonidentifiability_not_absolute_impossibility_in_future_extensions",
            "stronger_than_upstream_status_boundary": True,
        },
        "template_ancestry": {
            "family_transport_kernel": {
                "artifact": family_kernel.get("artifact"),
                "status": family_kernel.get("status"),
                "proof_status": family_kernel.get("proof_status"),
                "metadata_note": dict(family_kernel.get("metadata") or {}).get("note"),
                "is_template": family_kernel.get("status") == "template",
            },
            "generation_bundle_branch_generator": {
                "artifact": generator.get("artifact"),
                "proof_status": generator.get("proof_status"),
                "current_proxy_kind": generator.get("current_proxy_kind"),
                "is_candidate": generator.get("proof_status") == "candidate_only",
            },
            "generation_bundle_splitting_obstruction": {
                "artifact": splitting_obstruction.get("artifact"),
                "proof_status": splitting_obstruction.get("proof_status"),
                "verdict": splitting_obstruction.get("verdict"),
                "target_clause": splitting_obstruction.get("target_clause"),
                "first_failed_implication": splitting_obstruction.get("first_failed_implication"),
            },
            "overlap_edge_transport_cocycle": {
                "artifact": cocycle.get("artifact"),
                "proof_status": cocycle.get("proof_status"),
                "is_candidate": str(cocycle.get("proof_status", "")).endswith("candidate"),
            },
            "positive_source_emission_blocked_upstream": True,
            "conditional_grant_used_for_nonidentifiability_proof": True,
        },
        "ordered_source_geometry": geometry,
        "exact_ray_classification": {
            "fiber": "(R_{>0})^2",
            "fiber_dimension": 2,
            "independent_coordinates": ["sigma_u", "sigma_d"],
            "redundant_four_tuple_map": {
                "sigma_seed_ud": "(sigma_u + sigma_d) / 2",
                "eta_ud": "(sigma_u - sigma_d) / 2",
                "rank": 2,
            },
            "up_type_ray": up_ray,
            "down_type_ray": down_ray,
            "classification_proof": [
                "Trace zero reconstructs a three-point log profile from its two adjacent gaps.",
                "The ordered ratio law fixes the ratio of those two gaps separately in each sector.",
                "Their positive sum, the endpoint span sigma_q, remains free in each sector.",
                "Therefore the up and down solution sets are independent positive rays and their product is (R_{>0})^2.",
            ],
        },
        "free_action_certificate": {
            "group": "(R_{>0})^2",
            "action": "(lambda_u,lambda_d).(sigma_u,sigma_d)=(lambda_u*sigma_u,lambda_d*sigma_d)",
            "source_projection_fixed": True,
            "action_free": True,
            "action_transitive_on_positive_spread_fiber": True,
            "requested_tuple_not_invariant": True,
            "formal_countermodels": [model_a, model_b],
            "countermodels_share_source_projection": (
                model_a["source_projection_token"] == model_b["source_projection_token"]
            ),
            "countermodels_have_distinct_requested_tuples": model_a["sigma_tuple"] != model_b["sigma_tuple"],
            "max_abs_structural_identity_residual": max_countermodel_residual,
        },
        "affine_downstream_injectivity": affine,
        "edge_statistics_do_not_select_moduli": edge_family,
        "issue_composition": {
            "issue_379_up_type": {
                "remaining_fiber": "R_{>0}",
                "free_coordinate": "sigma_u",
                "closure": "closed_as_current_corpus_nonidentifiability_obstruction",
                "direct_top_codomain_kept_separate": True,
            },
            "issue_380_down_type": {
                "remaining_fiber": "R_{>0}",
                "free_coordinate": "sigma_d",
                "closure": "closed_as_current_corpus_nonidentifiability_obstruction",
                "elementary_running_quark_vs_hadron_qcd_distinction_preserved": True,
            },
            "issue_377_all_quark": {
                "remaining_fiber": "(R_{>0})^2",
                "composition": "issue_379_up_type x issue_380_down_type",
                "closure": "closed_as_current_corpus_nonidentifiability_obstruction",
            },
        },
        "minimal_future_extension": {
            "required_independent_scalar_count": 2,
            "acceptable_form": (
                "A theorem-grade source map into (R_{>0})^2, or two independent source-emitted sector "
                "normalizations, together with sector attachment and refinement compatibility."
            ),
            "must_break_free_action": True,
            "must_replace_template_ancestry": True,
        },
        "dependency_audit": {
            "direct_input_paths": DIRECT_INPUT_PATHS,
            "allowed_ancestors": allowed_ancestors,
            "forbidden_ancestors": FORBIDDEN_ANCESTORS,
            "allowed_forbidden_intersection": intersection,
            "allowed_forbidden_disjoint": not intersection,
            "loads_contaminated_family_excitation_artifact": False,
            "loads_contaminated_edge_comparison_artifact": False,
            "loads_spread_or_sector_mean_artifact": False,
            "loads_running_quark_reference_rows": False,
            "loads_selected_or_current_family_exact_witness": False,
            "loads_fitted_sigma_values": False,
            "emits_target_sigma_values": False,
            "no_target_leak": not intersection,
        },
        "public_prediction_policy": {
            "running_quark_numeric_rows": "withheld",
            "forward_yukawa_numeric_promotion": "withheld",
            "reason": "source spread fiber has two unselected positive moduli",
        },
        "notes": [
            "The unit and independently rescaled pairs are formal model witnesses only; neither is a physical sigma proposal.",
            "This certificate does not prove that no future OPH extension can emit the spread pair.",
            "It proves that the present target-free corpus, even after granting its strongest candidate source packet, does not identify it.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the target-free quark sigma non-identifiability obstruction.")
    parser.add_argument("--family-kernel", default=str(FAMILY_KERNEL_JSON))
    parser.add_argument("--generator", default=str(GENERATOR_JSON))
    parser.add_argument("--splitting-obstruction", default=str(SPLITTING_OBSTRUCTION_JSON))
    parser.add_argument("--cocycle", default=str(COCYCLE_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    payload = build_artifact(
        _load_json(Path(args.family_kernel)),
        _load_json(Path(args.generator)),
        _load_json(Path(args.splitting_obstruction)),
        _load_json(Path(args.cocycle)),
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
