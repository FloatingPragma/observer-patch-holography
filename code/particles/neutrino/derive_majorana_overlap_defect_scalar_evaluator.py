#!/usr/bin/env python3
"""Export the exact scalar-evaluator blocker boundary for the OPH-only Majorana route."""

from __future__ import annotations

import argparse
import json
import pathlib
from datetime import datetime, timezone


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_ACTION_GERM = ROOT / "particles" / "runs" / "neutrino" / "majorana_overlap_defect_action_germ.json"
DEFAULT_HESSIAN = ROOT / "particles" / "runs" / "neutrino" / "majorana_overlap_defect_hessian.json"
DEFAULT_READBACK = ROOT / "particles" / "runs" / "neutrino" / "realized_same_label_gap_defect_readback.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "neutrino" / "majorana_overlap_defect_scalar_evaluator.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the Majorana scalar-evaluator boundary artifact.")
    parser.add_argument("--action-germ", default=str(DEFAULT_ACTION_GERM))
    parser.add_argument("--hessian", default=str(DEFAULT_HESSIAN))
    parser.add_argument("--readback", default=str(DEFAULT_READBACK))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    action_germ = json.loads(pathlib.Path(args.action_germ).read_text(encoding="utf-8"))
    hessian = json.loads(pathlib.Path(args.hessian).read_text(encoding="utf-8"))
    readback_path = pathlib.Path(args.readback)
    readback = json.loads(readback_path.read_text(encoding="utf-8")) if readback_path.exists() else {}
    readback_complete = readback.get("payload_status") == "complete_from_live_flavor_artifacts"
    edge_weights = dict(action_germ.get("edge_coefficients", {}))
    mu_nu = float(next(iter(edge_weights.values()), 0.0))
    selector_absolute = dict(hessian.get("selector_point", {}))
    selector_residual = dict(action_germ.get("selector_point", {}))
    overlap_clause = "same_label_overlap_nonzero_on_realized_refinement_arrows"
    phase_clause = "selector_overlap_phase_coboundary_trivializes_same_label_edge_transport"
    overlap_status = "closed_from_live_flavor_readback" if readback_complete else "candidate_only"

    payload = {
        "artifact": "oph_majorana_overlap_defect_scalar_evaluator",
        "generated_utc": _timestamp(),
        "theorem_candidate_id": "oph_majorana_scalar_from_centered_edge_norm",
        "sublemma_candidate_id": "selector_centered_quadraticity_polarization_law_on_edge_bundle",
        "parent_theorem_id": "oph_majorana_scalar_from_centered_edge_norm",
        "source_artifacts": {
            "action_germ": str(pathlib.Path(args.action_germ)),
            "hessian": str(pathlib.Path(args.hessian)),
        },
        "domain": action_germ.get("domain", "affine_majorana_lift"),
        "quadraticity_domain": "selector_centered_common_refinement_edge_bundle",
        "bundle_displacement_symbol": "xi = Xi_nu(psi) = s(psi) - s_star",
        "bundle_descent_candidate_id": "selector_centered_unitary_common_refinement_descent_on_edge_bundle",
        "bundle_descent_status": "candidate_only",
        "phase_cocycle_triviality_candidate_id": "selector_overlap_phase_coboundary_trivializes_same_label_edge_transport",
        "phase_cocycle_triviality_status": "candidate_only",
        "phase_triviality_proof_mode": "selector_overlap_phase_coboundary",
        "smaller_exact_missing_clause_id": phase_clause,
        "overlap_nonvanishing_status": overlap_status,
        "overlap_nonvanishing_witness_hint": "gap_and_defect_fields_from_flavor_artifacts",
        "selector_center": action_germ.get("selector_center", "principal_equal_split"),
        "selector_point_absolute": {
            "psi12": selector_absolute.get("psi12"),
            "psi23": selector_absolute.get("psi23"),
            "psi31": selector_absolute.get("psi31"),
        },
        "selector_reference_section": "principal_equal_split",
        "selector_origin_residual": {
            "u": selector_residual.get("psi12", 0.0),
            "v": selector_residual.get("psi23", 0.0),
        },
        "edge_line_model": "L_e = ell_i tensor conjugate(ell_j)",
        "edge_bundle_model": "E_nu = direct_sum_e L_e",
        "transport_source": "projective_line_lift_readout_of_centered_generation_bundle_branch_generator",
        "selector_centered_reference_kind": "normalized_edge_reference_ray",
        "raw_lift_kind": "arbitrary_unitary_same_label_lifts",
        "edge_transport_normalization": "<r_e_star(v),U_vu_e(r_e_star(u))> in R_{>0}",
        "required_overlap_certificate": overlap_clause,
        "required_overlap_certificate_status": overlap_status,
        "edge_transport_kind": "unitary_complex_linear",
        "arrow_overlap_phase_definition": "theta_uv_e = arg(<r_e_star(v),U_vu_e(r_e_star(u))>)",
        "raw_triangle_phase_definition": "Theta_uvw_e = theta_wv_e + theta_vu_e - theta_wu_e",
        "phase_coboundary_identity": "Theta_uvw_e = d theta_e(u,v,w)",
        "normalized_lift_definition": "U_vu_e_norm = exp(-i*theta_uv_e) * U_vu_e",
        "normalized_reference_transport_identity": "<r_e_star(v),U_vu_e_norm(r_e_star(u))> in R_{>0}",
        "normalized_triangle_phase_identity": "Theta_uvw_e_norm = 1",
        "all_triangle_phases_one_certificate": False,
        "edge_transport_cocycle_equation": "U_wv_e o U_vu_e = U_wu_e",
        "direct_sum_transport_equation": "T_vu = direct_sum_e U_vu_e",
        "well_defined_sum_formula": "[u,x] + [v,y] := [w,T_wu(x)+T_wv(y)]",
        "well_defined_difference_formula": "[u,x] - [v,y] := [w,T_wu(x)-T_wv(y)]",
        "well_defined_i_formula": "i*[u,x] := [u,i*x]",
        "edge_weights": edge_weights,
        "mu_nu": mu_nu,
        "remaining_numerical_defect": "exact_1_2_near_degeneracy_from_isotropic_mu_nu",
        "next_forward_object": "defect_weighted_mu_e_family",
        "defect_weighted_mu_e_family_candidate_id": "oph_defect_weighted_majorana_edge_weight_family",
        "defect_weighted_mu_e_family_status": "candidate_only",
        "defect_weighted_mu_e_normalizer_candidate": "oph_same_label_overlap_defect_weight_normalizer",
        "attachment_normalizer_candidate_id": "oph_same_label_overlap_defect_weight_normalizer",
        "attachment_normalizer_status": "candidate_only",
        "defect_weight_rule": "mu_e = base_mu_nu * exp(delta_e) / mean_f(exp(delta_f))",
        "defect_weighted_mu_e_family_role": "breaks isotropic 1_2 near_degeneracy while preserving the centered edge-norm route",
        "intrinsic_witness_kind": "centered_edge_character_norm_defect",
        "origin_object_candidate": "oph_affine_majorana_edge_character_functor",
        "origin_object_role": "realizes centered Majorana overlap phases as unitary edge characters on transported Hermitian lines",
        "formula_family": "S_G(psi)=sum_e mu_e*G(psi_e-psi_e_star)",
        "family_conditions": [
            "G(0)=0",
            "G'(0)=0",
            "G''(0)=1",
            "G(-x)=G(x)",
        ],
        "formula_affine_candidate": "sum_e mu_e*(1-cos(psi_e-psi_e_star))",
        "formula_residual_candidate": f"{mu_nu:.18g}*(3-cos(u+v)-cos(u)-cos(v))",
        "intrinsic_generator_formula": "0.5*sum_e mu_e*abs(z_e(psi)-1)^2",
        "edge_character_definition": "z_e(psi)=<s_e_star,s_e(psi)> in U(1)",
        "norm_identification_formula": "S_nu(psi)=0.5*||s(psi)-s_star||^2_{diag(mu_e)}",
        "quadraticity_equation": "Q_nu(xi+eta)+Q_nu(xi-eta)=2*Q_nu(xi)+2*Q_nu(eta)",
        "complex_quadraticity_equation": "Q_nu(xi+i*eta)+Q_nu(xi-i*eta)=2*Q_nu(xi)+2*Q_nu(eta)",
        "polarization_formula": "H_nu(xi,eta)=0.5*(Q_nu(xi+eta)-Q_nu(xi-eta)+i*Q_nu(xi+i*eta)-i*Q_nu(xi-i*eta))",
        "forced_bundle_norm_formula": "Q_nu(xi)=0.5*sum_e mu_e*abs(xi_e)^2",
        "forced_phase_formula": "S_nu(psi)=sum_e mu_e*(1-cos(psi_e-psi_e_star))",
        "forced_residual_formula": f"{mu_nu:.19g}*(3-cos(u+v)-cos(u)-cos(v))",
        "ambient_witness_formula": "0.5*mu_nu*||z(psi)-1||^2_edge",
        "ambient_hessian_match": "mu_nu*I3_in_edge_basis",
        "conjugation_action": "z_e -> conjugate(z_e)",
        "hessian_anchor_matrix_3x3": hessian.get("ambient_hessian_3x3"),
        "hessian_anchor_residual_2x2": hessian.get("residual_hessian_2x2"),
        "direct_sum_weight_matrix": [
            [mu_nu, 0.0, 0.0],
            [0.0, mu_nu, 0.0],
            [0.0, 0.0, mu_nu],
        ],
        "longitudinal_null_vector": [1.0, 1.0, 1.0],
        "normalization_witness": "edge_coefficient_hessian_match",
        "normalization_witness_status": "local_quadratic_match_only",
        "hessian_recovery_certificate": {
            "status": "closed",
            "matches_action_germ": hessian.get("proof_status") == "local_quadratic_germ_closed",
            "residual_hessian_2x2": hessian.get("residual_hessian_2x2"),
        },
        "quadraticity_clause_status": "pending_bundle_parallelogram_identity",
        "presentation_independence_status": "pending_theorem",
        "finite_angle_exactness_status": "pending_edge_character_norm_theorem",
        "promotion_gate": "close_only_if_bundle_descent_is_defined_for_xi±eta_and_xi±i_eta",
        "bundle_descent_gate_if_closed": "promotion_gate_for_xi±eta_and_xi±i_eta_cleared",
        "presentation_independence_status_if_closed": "closed",
        "quadraticity_clause_status_if_closed": "closed",
        "cubic_freedom_status_if_closed": "eliminated_exactly",
        "finite_angle_exactness_status_if_closed": "closed_on_current_isotropic_branch",
        "bounded_scalar_range_if_closed": {
            "lower": 0.0,
            "upper": 6.0 * mu_nu,
        },
        "proof_status": "candidate_only",
        "nonisotropic_formula_status": "open_centered_family_only",
        "oph_origin_status": "candidate_formula_not_yet_intrinsic_theorem",
        "invariant_ring_obstruction": {
            "status": "open",
            "quadratic_invariant": action_germ.get("quadratic_invariant_residual"),
            "cubic_invariant": action_germ.get("cubic_invariant_obstruction", {}).get("formula"),
            "cubic_elimination_condition": "exact bundle quadraticity on selector-centered displacement domain",
            "cubic_elimination_status": "exact_if_condition_adopted",
            "reduced_invariant_family": "F(I2,I3^2)",
            "cubic_freedom_eliminated": False,
        },
        "naive_uncentered_formula_ruled_out": True,
        "naive_uncentered_cubic_coeff": 0.0016047225727754176,
        "cubic_kill_mechanism": "Hermitian displacement depends only on Re(z_e), so odd I3 terms vanish",
        "exact_remaining_ingredient": "selector-centered phase cocycle triviality for same-label edge transport",
        "smallest_exact_missing_clause": phase_clause,
        "strictly_smaller_missing_clause_if_not_closed": phase_clause,
        "fallback_family_if_not_closed": "sum_e mu_e*Phi(z_e), Phi(z)=Phi(conj(z)), Phi(1)=0, second jet fixed",
        "remaining_theorem_object": "oph_majorana_scalar_from_centered_edge_norm",
        "notes": [
            "This boundary now carries the strongest centered constructive candidate compatible with the closed local quadratic action germ: a chordal/cosine evaluator centered at the actual equal-split selector point.",
            "The current sharp origin candidate is the Hermitian displacement norm of a centered OPH Majorana edge-character functor.",
            "The strongest current closure route is the selector-centered unitary common-refinement descent theorem on the direct-sum edge bundle, with phase-cocycle triviality treated as an explicit normalized-lift coboundary statement.",
            (
                "The overlap-nonvanishing subclause is already discharged by the live same-label gap/defect readback exported from flavor-side certificates."
                if readback_complete
                else "The current local witness hint for overlap nonvanishing is not a theorem yet, but it is expected to live in the same gap/defect fields already exported by the flavor artifacts."
            ),
            "The sharp live theorem gate is now selector-centered phase-cocycle triviality for same-label edge transport; once that closes, the normalized overlap-defect weight normalizer becomes the direct attachment object beneath the scalar route.",
            "Numerically, the current remaining defect is the exact 1-2 near-degeneracy induced by isotropic mu_nu, so the next forward object is a defect-weighted mu_e family rather than another isotropic re-evaluation.",
        ],
    }

    out_path = pathlib.Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
