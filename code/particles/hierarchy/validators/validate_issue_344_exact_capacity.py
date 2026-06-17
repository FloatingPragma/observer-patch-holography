#!/usr/bin/env python3
"""Validate the OPH issue #344 exact-capacity certificate."""

from __future__ import annotations

import json
import pathlib
import sys
from decimal import Decimal


def D(value: str) -> Decimal:
    return Decimal(value)


def main(path: str = "certificates/R_EW_global_capacity_certificate.json") -> int:
    cert_path = pathlib.Path(path)
    cert = json.loads(cert_path.read_text(encoding="utf-8"))
    exact = cert.get("exact_capacity_fixed_point", {})
    contraction = cert.get("contraction_certificate", {})
    rounded = cert.get("rounded_capacity_diagnostic", {})
    forbidden = cert.get("forbidden_calibrations", [])
    boundary = cert.get("claim_boundary", {})

    derivation = cert.get("derivation_chain", [])
    factors = cert.get("factor_origins", {})
    branch_scope = cert.get("branch_scope", {})
    acceptance = cert.get("acceptance_criteria_status", {})
    deps = cert.get("dependency_artifacts", {})
    consumers = cert.get("consumer_artifacts", {})
    acyclic = cert.get("dependency_acyclicity_note", {})

    derivation_steps = {item.get("step"): item for item in derivation if isinstance(item, dict)}

    def _step_conclusion(step: int) -> str:
        node = derivation_steps.get(step, {})
        return node.get("conclusion", "") if isinstance(node, dict) else ""

    checks = {
        "issue_is_344": cert.get("issue") == 344,
        "certificate_id_v2": cert.get("certificate_id") == "issue-344-exact-ew-refined-global-capacity-v2",
        "accepted": cert.get("accepted") is True,
        "status_closed_exact_capacity": (
            cert.get("status") == "closed_bridge_refined_global_capacity_fixed_point_certificate"
        ),
        "target_relation_states_zero_bridge_residual": (
            "B_EW(P_star, N_CRC^EW)" in cert.get("target_relation", "")
            and "= 0" in cert.get("target_relation", "")
            and "Pi_EW(P_star, N_CRC^EW) = 4*P_star" in cert.get("target_relation", "")
        ),
        "bridge_residual_zero": D(exact.get("bridge_residual", "1")) == 0,
        "fixed_point_residual_zero": D(exact.get("fixed_point_residual_x", "1")) == 0,
        "projection_exponent_matches_4P": (
            abs(D(exact.get("projection_exponent", "0")) - D(exact.get("target_exponent_4P", "1")))
            <= Decimal("1e-40")
        ),
        "v_identity_matches": abs(D(exact.get("v_identity_error", "1"))) <= Decimal("1e-40"),
        "contraction_factor_half": D(contraction.get("lipschitz_constant", "0")) == Decimal("0.5"),
        "residual_contracts": (
            abs(
                D(contraction.get("sample_residual_ratio", "0"))
                - D(contraction.get("residual_contracts_by", "1"))
            )
            <= Decimal("1e-40")
        ),
        "rounded_capacity_is_diagnostic": (
            rounded.get("status") == "diagnostic_only_not_exact_bridge_certificate"
        ),
        "rounded_capacity_fails_bridge": abs(D(rounded.get("bridge_residual", "0"))) > Decimal("1e-6"),
        "weak_scale_forbidden": any("measured weak scale" in item for item in forbidden),
        "lambda_forbidden": any("measured Lambda" in item for item in forbidden),
        "rounded_display_forbidden": any("3.31e122" in item for item in forbidden),
        "finite_readback_recorded_elsewhere": any(
            "R_readback_resolution_certificate" in item
            for item in boundary.get("closed_elsewhere", [])
        ),
        "round_count_recorded_elsewhere": any(
            "R_m_rep_24_certificate" in item
            for item in boundary.get("closed_elsewhere", [])
        ),
        "no_remaining_boundary": boundary.get("not_closed_here", []) == [],
        "boundary_scope_present": "scope" in boundary and "N_CRC^EW" in boundary.get("scope", ""),
        "derivation_chain_has_eight_steps": (
            len(derivation) == 8 and set(derivation_steps.keys()) == set(range(1, 9))
        ),
        "step_1_imports_pixel_fixed_point": (
            "P_star" in _step_conclusion(1) and "alpha_U" in _step_conclusion(1)
        ),
        "step_2_imports_d10_beta_ew": "beta_EW" in _step_conclusion(2) and "= 4" in _step_conclusion(2),
        "step_3_imports_m_rep_24": (
            "m_rep = 24" in _step_conclusion(3)
            and ("6 = m_rep / beta_EW" in _step_conclusion(3) or "6 = m_rep/beta_EW" in _step_conclusion(3))
        ),
        "step_4_imports_pi_ew_resonance_target": (
            "Pi_EW(P_star, N_CRC^EW) = 4*P_star" in _step_conclusion(4)
        ),
        "step_5_equates_resonance_to_bridge_residual": (
            "B_EW(P_star,N)" in _step_conclusion(5) and "= 0" in _step_conclusion(5)
        ),
        "step_6_solves_closed_form_capacity": (
            "N_CRC^EW(P_star) = pi*exp[6*pi/(P_star*alpha_U(P_star))]" in _step_conclusion(6)
        ),
        "step_7_certifies_banach_contraction": (
            "Banach" in _step_conclusion(7) and "1 - lambda = 1/2" in _step_conclusion(7)
        ),
        "step_8_records_numerical_witness_and_rejects_rounded": (
            "B_EW(P_star, N_CRC^EW) = 0" in _step_conclusion(8) and "3.31e122" in _step_conclusion(8)
        ),
        "factor_origin_p_star_recorded": "P_star_pixel_fixed_point" in factors,
        "factor_origin_alpha_u_recorded": "alpha_U_unification_width" in factors,
        "factor_origin_beta_ew_value_4": (
            factors.get("beta_EW_transmutation_multiplicity", {}).get("value") == "4"
        ),
        "factor_origin_m_rep_value_24": (
            factors.get("m_rep_doubled_sm_adjoint_round_count", {}).get("value") == "24"
        ),
        "factor_origin_six_recorded_as_m_rep_over_beta_ew": (
            factors.get("factor_six_in_bridge_residual", {}).get("value") == "6"
            and "m_rep / beta_EW" in factors.get("factor_six_in_bridge_residual", {}).get("expression", "")
        ),
        "factor_origin_twenty_four_recorded": (
            factors.get("factor_twenty_four_in_projection_map", {}).get("value") == "24"
        ),
        "factor_origin_lambda_one_half_recorded": (
            factors.get("banach_contraction_lambda_one_half", {}).get("value") == "1/2"
        ),
        "factor_origin_pi_recorded": "factor_pi_in_capacity_normalisation" in factors,
        "factor_origin_resonance_target_recorded": (
            "factor_four_P_star_resonance_target" in factors
            and factors["factor_four_P_star_resonance_target"].get("value") == "4 * P_star"
        ),
        "branch_scope_pixel_branch_present": "oph_local_pixel_branch" in branch_scope,
        "branch_scope_d10_branch_present": "d10_transmutation_branch" in branch_scope,
        "branch_scope_rep_branch_present": "representation_to_spectrum_branch" in branch_scope,
        "branch_scope_pi_ew_branch_present": "ew_tick_projection_branch" in branch_scope,
        "branch_scope_banach_branch_present": "banach_contraction_branch" in branch_scope,
        "branch_scope_note_present": (
            "scope_note" in branch_scope
            and "N_CRC^EW(P_star) = pi*exp[6*pi/(P_star*alpha_U(P_star))]" in branch_scope.get("scope_note", "")
        ),
        "acceptance_exact_global_capacity_fixed_point_defined": (
            acceptance.get("exact_global_capacity_fixed_point_defined") is True
        ),
        "acceptance_bridge_residual_zero": (
            acceptance.get("bridge_residual_zero_on_source_side") is True
        ),
        "acceptance_banach_contraction_certified": (
            acceptance.get("banach_contraction_certified_with_explicit_lipschitz_constant") is True
        ),
        "acceptance_rounded_display_rejected": (
            acceptance.get("rounded_capacity_display_rejected_as_exact_witness") is True
        ),
        "acceptance_forbidden_calibrations_listed": (
            acceptance.get("forbidden_calibrations_listed_and_unused") is True
        ),
        "acceptance_machine_readable_published": (
            acceptance.get("machine_readable_verifier_and_certificate_published") is True
        ),
        "acceptance_downstream_unchanged": (
            acceptance.get("downstream_surfaces_unchanged_because_status_unchanged") is True
        ),
        "dependency_artifacts_pixel_present": "oph_local_pixel_fixed_point" in deps,
        "dependency_artifacts_alpha_u_present": "alpha_u_unification_width" in deps,
        "dependency_artifacts_d10_present": "d10_transmutation_theorem_beta_ew" in deps,
        "dependency_artifacts_m_rep_present": "representation_to_spectrum_m_rep_24" in deps,
        "dependency_artifacts_pi_ew_present": "ew_tick_projection_pi_ew_definition" in deps,
        "dependency_artifacts_banach_present": "banach_fixed_point_theorem" in deps,
        "consumer_artifacts_pi_ew_specialisation_present": (
            "ew_tick_projection_specialisation" in consumers
        ),
        "consumer_artifacts_readback_dependency_present": (
            "finite_readback_resolution_dependency" in consumers
        ),
        "consumer_artifacts_umbrella_present": (
            "local_global_hierarchy_resonance_umbrella" in consumers
        ),
        "acyclicity_summary_states_peer_cross_reference": (
            "peer cross-reference" in acyclic.get("summary", "")
            and "acyclic" in acyclic.get("summary", "")
        ),
        "acyclicity_primary_theorems_independent": (
            "primary_theorems_are_independent" in acyclic
            and "ew_tick_projection_primary" in acyclic.get("primary_theorems_are_independent", {})
            and "exact_capacity_primary" in acyclic.get("primary_theorems_are_independent", {})
        ),
        "acyclicity_specialised_corollary_explained": (
            "specialised_corollary_is_a_composition_not_a_circle" in acyclic
        ),
        "acyclicity_umbrella_resolves_composition": (
            "umbrella_certificate_resolves_the_composition" in acyclic
        ),
        "acyclicity_other_branches_upstream_only": (
            "other_remaining_branches_are_upstream_only" in acyclic
        ),
    }
    payload = {"checks": checks, "pass": all(checks.values())}
    print(json.dumps(payload, indent=2))
    return 0 if payload["pass"] else 1


if __name__ == "__main__":
    sys.exit(
        main(sys.argv[1] if len(sys.argv) > 1 else "certificates/R_EW_global_capacity_certificate.json")
    )
