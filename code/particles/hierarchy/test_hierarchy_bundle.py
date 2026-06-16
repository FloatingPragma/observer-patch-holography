#!/usr/bin/env python3
"""Guard the frozen OPH hierarchy proof bundle."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parent


def _run(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, *args],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=True,
    )


def test_hierarchy_bundle_validators_pass() -> None:
    result = _run("validators/validate_bundle.py")
    payload = json.loads(result.stdout)

    assert len(payload) == 11
    assert all(entry["returncode"] == 0 for entry in payload)
    validator_outputs = [json.loads(entry["stdout"]) for entry in payload]
    assert all(output["pass"] is True for output in validator_outputs)


def test_ru_krawczyk_certificate_is_unique_root_witness() -> None:
    result = _run(
        "validators/validate_ru_interval_certificate.py",
        "certificates/R_U_krawczyk_certificate.json",
    )
    payload = json.loads(result.stdout)

    assert payload["K_subset_interior_I"] is True
    assert payload["derivative_excludes_zero"] is True
    assert payload["pass"] is True

    cert = json.loads((ROOT / "certificates/R_U_krawczyk_certificate.json").read_text())
    assert cert["status"] == "krawczyk_inclusion_witness_supplied"
    assert cert["center_c"] == "0.041124336195630495"
    assert cert["inclusion"]["K_I_subset_interior_I_U"] is True


def test_hierarchy_numeric_witness_keeps_public_and_source_audit_branches_separate() -> None:
    witness = json.loads((ROOT / "computations/hierarchy_numeric_witness.json").read_text())

    public = witness["public_endpoint_branch"]
    source_audit = witness["source_audit_branch"]
    assert public["P_C"] == "1.630968209403959324879279847782648941"
    assert public["alpha_U_display"] == "0.041124336195630495"
    assert public["v_over_E_star"] == "2.0199803239725553e-17"
    assert source_audit["P_cand"] == "1.63097209569432901817967892561191884270169"
    assert source_audit["alpha_U"] == "0.04112424744557487"
    assert source_audit["v_over_E_star"] == "2.0198114150099223e-17"


def test_global_repair_tick_lemma_is_closed_with_derived_round_count() -> None:
    result = _run(
        "validators/validate_global_repair_tick_certificate.py",
        "certificates/R_N_global_repair_tick_certificate.json",
    )
    payload = json.loads(result.stdout)

    assert payload["pass"] is True
    assert payload["status_is_closed_theorem_with_derived_round_count"] is True
    assert payload["theorem_kind_is_theorem_with_derived_round_count"] is True
    assert payload["tick_exponent_is_minus_one_over_48"] is True
    assert payload["tick_times_rounds_is_full_cycle"] is True
    assert payload["full_cycle_map_recorded"] is True
    assert payload["one_tick_map_recorded"] is True
    assert payload["numeric_display_matches_formula"] is True
    assert payload["full_cycle_multiplier_is_derived_from_closure"] is True
    assert payload["f_interface_equivalence_derived"] is True
    assert payload["round_count_is_derived"] is True
    assert payload["round_count_source_recorded"] is True
    assert payload["no_open_round_count_boundary"] is True
    assert payload["ew_inputs_excluded_from_derived_uses"] is True

    cert = json.loads((ROOT / "certificates/R_N_global_repair_tick_certificate.json").read_text())
    assert cert["status"] == "closed_global_repair_tick_theorem_with_derived_round_count"
    assert cert["theorem_kind"] == "theorem_with_derived_round_count"
    assert cert["normalization"]["abs_g_star_prime"] == "(N_CRC/pi)^(-1/48)"
    assert cert["exponent_law"]["per_tick_exponent_for_m_ticks"] == "-1/(2m)"
    assert cert["normalization"]["round_count_source"] == "R_m_rep_24_certificate.json"
    acceptance = cert["acceptance_criteria_status"]
    assert acceptance["proves_declared_screen_capacity_fixed_point_emits_tick_contraction"] is True
    assert acceptance["round_count_derived_from_first_principles"] is True
    assert acceptance["round_count_certificate_recorded"] is True
    assert acceptance["closure_transport_derived_from_F_interface"] is True
    assert acceptance["readback_counting_model_is_modeling_identification"] is True
    assert acceptance["concrete_finite_machinery_verification_open"] is False
    assert acceptance["finite_readback_resolution_certificate_recorded"] is True
    declared = cert["claim_boundary"]["declared_not_derived"]
    assert any("modeling identification" in item for item in declared)
    boundary = cert["claim_boundary"]["not_closed_by_certificate"]
    assert boundary == []
    closed_elsewhere = cert["claim_boundary"]["closed_elsewhere"]
    assert any("R_m_rep_24_certificate" in item for item in closed_elsewhere)
    assert any("R_EW_tick_projection_certificate" in item for item in closed_elsewhere)
    assert any("R_EW_global_capacity_certificate" in item for item in closed_elsewhere)
    assert any("R_readback_resolution_certificate" in item for item in closed_elsewhere)


def test_joint_pn_fixed_point_certificate_records_product_closure_and_coupled_boundary() -> None:
    result = _run(
        "validators/validate_joint_fixed_point_certificate.py",
        "certificates/R_PN_joint_fixed_point_certificate_report.json",
    )
    payload = json.loads(result.stdout)

    assert payload["pass"] is True
    checks = payload["checks"]
    assert checks["product_theorem_status"] is True
    assert checks["backsolve_is_diagnostic_only"] is True
    assert checks["coupled_residual_boundary_recorded"] is True

    cert = json.loads((ROOT / "certificates/R_PN_joint_fixed_point_certificate_report.json").read_text())
    assert cert["status"] == "closed_product_branch_theorem_with_explicit_coupled_branch_boundary"
    assert cert["product_contraction_certificate"]["status"] == "conditional_on_component_contractions"
    assert cert["coupled_contraction_certificate"]["status"] == "residual_coupled_branch_boundary"
    assert "CIRCULAR_DIAGNOSTIC_ONLY" in cert["N_backsolved_warning"]


def test_issue_337_electroweak_projection_certificate_records_exact_bridge_condition() -> None:
    result = _run(
        "validators/validate_issue_337_electroweak_projection.py",
        "certificates/R_EW_tick_projection_certificate.json",
    )
    payload = json.loads(result.stdout)

    assert payload["pass"] is True
    checks = payload["checks"]
    assert checks["projection_exponent_matches_4P"] is True
    assert checks["rounded_N_is_diagnostic"] is True
    assert checks["rounded_N_fails_exact_bridge"] is True
    assert checks["derivation_chain_has_seven_steps"] is True
    assert checks["step_3_derives_projection_formula"] is True
    assert checks["step_4_records_resonance_target_scope"] is True
    assert checks["step_6_cites_capacity_certificate"] is True
    assert checks["step_7_derives_P_over_12_form"] is True
    assert checks["factor_origin_beta_EW_recorded"] is True
    assert checks["factor_origin_m_rep_recorded"] is True
    assert checks["factor_origin_48_recorded"] is True
    assert checks["factor_origin_12_recorded"] is True
    assert checks["acceptance_projection_map_defined"] is True
    assert checks["acceptance_4P_proved_under_resonance_target"] is True
    assert checks["acceptance_factor_4_origin_documented"] is True
    assert checks["acceptance_factor_12_origin_documented"] is True
    assert checks["acceptance_compatible_with_local_D10"] is True
    assert checks["acceptance_no_measured_weak_inputs"] is True
    assert checks["acceptance_resonance_target_scoped_as_oph_condition"] is True
    assert checks["boundary_records_closed_elsewhere"] is True
    assert checks["boundary_includes_scope_note"] is True

    cert = json.loads((ROOT / "certificates/R_EW_tick_projection_certificate.json").read_text())
    assert cert["accepted"] is True
    assert cert["status"] == "closed_projection_map_with_exact_bridge_condition"
    assert cert["exact_bridge"]["bridge_residual"] == "0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"
    assert cert["rounded_capacity_diagnostic"]["status"] == "diagnostic_only_not_exact_bridge_certificate"

    chain = cert["derivation_chain"]
    assert len(chain) == 7
    assert chain[0]["premise"] == "D10 forward transmutation theorem"
    assert "beta_EW" in str(chain[0]["uses"]) or "beta_EW" in str(chain[0].get("conclusion", ""))
    assert "m_rep" in str(chain[1]["uses"]) or "m_rep" in str(chain[1].get("conclusion", ""))
    assert "Pi_EW" in chain[2]["conclusion"] and "24*pi" in chain[2]["conclusion"]
    assert "scope_note" in chain[3]
    assert "B_EW" in chain[4]["conclusion"]
    assert "R_EW_global_capacity" in chain[5]["source"]
    assert "(P_star/12)" in chain[6]["conclusion"] and "48/4" in chain[6]["conclusion"]

    factors = cert["factor_origins"]
    assert factors["beta_EW"]["value"] == "4"
    assert "D10" in factors["beta_EW"]["source_theorem"]
    assert factors["m_rep"]["value"] == "24"
    assert "representation-to-spectrum" in factors["m_rep"]["source_theorem"]
    assert factors["tick_exponent_denominator_48"]["value"] == "48"
    assert "global repair-tick" in factors["tick_exponent_denominator_48"]["source_theorem"]
    assert factors["projection_target_factor_4_in_4P"]["identification"] == "beta_EW"
    assert factors["projection_target_denominator_12_in_P_over_12"]["value"] == "12"
    assert "48 / 4" in factors["projection_target_denominator_12_in_P_over_12"]["definition"]

    acceptance = cert["acceptance_criteria_status"]
    assert acceptance["projection_map_defined"] is True
    assert acceptance["sampling_exponent_4P_proved_under_resonance_target"] is True
    assert acceptance["factor_4_origin_documented"] is True
    assert acceptance["factor_12_origin_documented"] is True
    assert acceptance["compatibility_with_local_D10_transmutation_certificate"] is True
    assert acceptance["no_measured_weak_scale_inputs"] is True
    assert acceptance["no_measured_higgs_top_W_Z_inputs"] is True
    assert acceptance["no_measured_gravity_inputs"] is True
    assert acceptance["rounded_N_display_rejected_as_high_precision_bridge"] is True
    assert acceptance["resonance_target_scoped_as_oph_condition"] is True

    boundary = cert["claim_boundary"]
    assert "derivation chain" in boundary["closed_here"]
    closed_elsewhere = boundary["closed_elsewhere"]
    assert any("D10" in item for item in closed_elsewhere)
    assert any("representation-to-spectrum" in item or "R_m_rep_24" in item for item in closed_elsewhere)
    assert any("global repair-tick" in item or "R_N_global_repair_tick" in item for item in closed_elsewhere)
    assert any("R_EW_global_capacity" in item or "EW-refined" in item for item in closed_elsewhere)
    scope = boundary["scope"]
    assert "resonance target" in scope
    assert "N_CRC^EW" in scope
    assert "EW-refined" in scope


def test_issue_344_exact_capacity_certificate_is_fixed_point_source_record() -> None:
    result = _run(
        "validators/validate_issue_344_exact_capacity.py",
        "certificates/R_EW_global_capacity_certificate.json",
    )
    payload = json.loads(result.stdout)

    assert payload["pass"] is True
    checks = payload["checks"]
    assert checks["bridge_residual_zero"] is True
    assert checks["fixed_point_residual_zero"] is True
    assert checks["residual_contracts"] is True
    assert checks["rounded_capacity_fails_bridge"] is True

    cert = json.loads((ROOT / "certificates/R_EW_global_capacity_certificate.json").read_text())
    assert cert["accepted"] is True
    assert cert["status"] == "closed_bridge_refined_global_capacity_fixed_point_certificate"
    assert cert["contraction_certificate"]["lipschitz_constant"] == "0.5"
    assert cert["exact_capacity_fixed_point"]["bridge_residual"] == "0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"


def test_issue_342_readback_resolution_certificate_is_singleton_resolution() -> None:
    result = _run(
        "validators/validate_issue_342_readback_resolution.py",
        "certificates/R_readback_resolution_certificate.json",
    )
    payload = json.loads(result.stdout)

    assert payload["pass"] is True
    checks = payload["checks"]
    assert checks["one_selected_atom"] is True
    assert checks["positive_root_extractor"] is True
    assert checks["refinement_bound_recorded"] is True
    assert checks["no_remaining_boundary"] is True
    assert checks["round_count_recorded_elsewhere"] is True
    assert checks["exact_capacity_recorded_elsewhere"] is True

    cert = json.loads((ROOT / "certificates/R_readback_resolution_certificate.json").read_text())
    assert cert["accepted"] is True
    assert cert["status"] == "closed_finite_readback_resolution_certificate"
    assert cert["source_status"]["remaining_for_full_hierarchy_resonance"] == []
    assert cert["claim_boundary"]["not_closed_here"] == []
    assert cert["capacity_register"]["selected_variance"] == "0"


def test_issue_343_m_rep_certificate_derives_twenty_four_rounds() -> None:
    result = _run(
        "validators/validate_issue_343_m_rep_24.py",
        "certificates/R_m_rep_24_certificate.json",
    )
    payload = json.loads(result.stdout)

    assert payload["pass"] is True
    checks = payload["checks"]
    assert checks["component_dimensions_are_8_3_1"] is True
    assert checks["oriented_support_is_24"] is True
    assert checks["spectral_period_is_24"] is True
    assert checks["su5_negative_control_recorded"] is True
    assert checks["no_forbidden_inputs_used"] is True

    cert = json.loads((ROOT / "certificates/R_m_rep_24_certificate.json").read_text())
    assert cert["accepted"] is True
    assert cert["status"] == "closed_representation_to_spectrum_round_count"
    assert cert["representation_sector"]["unoriented_adjoint_dimension"] == 12
    assert cert["representation_sector"]["oriented_support_dimension"] == 24
    assert cert["result"]["specialized_exponent"] == "-1/48"
    assert cert["claim_boundary"]["not_closed_here"] == []


def test_issue_332_rg_higgs_naturality_certificate_is_zero_defect() -> None:
    result = _run(
        "validators/validate_issue_332_rg_naturality.py",
        "issue_332_rg_naturality_certificate.json",
    )
    payload = json.loads(result.stdout)

    assert payload["pass"] is True
    checks = payload["checks"]
    assert checks["epsilon_H_zero"] is True
    assert checks["weak_scale_forbidden"] is True
    assert checks["higgs_mass_forbidden"] is True

    cert = json.loads((ROOT / "issue_332_rg_naturality_certificate.json").read_text())
    assert cert["accepted"] is True
    assert cert["epsilon_H_interval"] == ["0", "0"]
    assert cert["optional_upstream_resonance_check"]["strict_resonance"] is False
    forbidden = cert["forbidden_calibrations"]
    assert any("measured weak scale" in item for item in forbidden)
    assert any("measured Higgs" in item for item in forbidden)


def test_issue_335_local_global_resonance_closes_as_full_selected_branch_statement() -> None:
    result = _run(
        "validators/validate_issue_335_local_global_resonance.py",
        "certificates/R_local_global_hierarchy_resonance_closeout_335.json",
    )
    payload = json.loads(result.stdout)

    assert payload["pass"] is True
    checks = payload["checks"]
    assert checks["full_theorem_promoted"] is True
    assert checks["exact_capacity_supplied"] is True
    assert checks["finite_readback_supplied"] is True
    assert checks["round_count_supplied"] is True
    assert checks["screen_sieve_supplied"] is True
    assert checks["screen_sieve_dependency_present"] is True
    assert checks["rounded_capacity_rejected"] is True
    assert checks["no_promotion_gates_remain"] is True
    assert checks["derivation_chain_has_nine_steps"] is True
    assert checks["step_5_is_screen_sieve_geometric_strengthening"] is True
    assert checks["step_8_composes_target_relation"] is True
    assert checks["factor_origin_icosahedral_orbit_recorded"] is True
    assert checks["factor_origin_cell_entropy_disclosed"] is True
    assert checks["branch_scope_records_screen_branch"] is True
    assert checks["honest_residual_disclosed_in_acceptance"] is True

    cert = json.loads((ROOT / "certificates/R_local_global_hierarchy_resonance_closeout_335.json").read_text())
    assert cert["accepted"] is True
    assert cert["status"] == "closed_full_local_global_hierarchy_resonance"
    assert cert["full_theorem_grade_resonance_promoted"] is True
    acceptance = cert["acceptance_criteria_status"]
    assert acceptance["prerequisite_steps_accounted_for"] is True
    assert acceptance["full_theorem_grade_resonance_proved"] is True
    assert acceptance["exact_capacity_source_certificate_supplied"] is True
    assert acceptance["finite_readback_resolution_supplied"] is True
    assert acceptance["round_count_derivation_supplied"] is True
    assert acceptance["screen_sieve_geometric_strengthening_supplied"] is True
    assert acceptance["honest_residual_definitional_residue_disclosed"] is True
    assert "P/beta_EW" in acceptance["honest_residual_disclosure"]
    assert cert["remaining_promotion_gates"] == []

    chain = cert["derivation_chain"]
    assert len(chain) == 9
    assert chain[0]["premise"] == "D10 forward transmutation theorem"
    assert "icosahedral screen-sieve" in chain[4]["premise"]
    assert "issue-#337" in chain[4]["geometric_strengthening_note"]
    assert "(P/12)" in chain[4]["conclusion"]
    assert "electroweak tick-projection bridge" in chain[5]["premise"]
    assert "EW-refined exact-capacity" in chain[6]["premise"]
    assert "(P_*/12)*log(N_CRC^EW/pi)" in chain[7]["conclusion"]
    assert "RG/Higgs naturality" in chain[8]["premise"]

    factors = cert["factor_origins"]
    assert factors["beta_EW"]["value"] == "4"
    assert factors["m_rep"]["value"] == "24"
    assert factors["icosahedral_orbit_size_12"]["value"] == "12"
    assert factors["icosahedral_orbit_size_12"]["definition"] == "|A5| / |C5| = 60 / 5"
    assert factors["total_curvature_charge_12"]["value"] == "12"
    assert factors["cell_entropy_factor_one_over_four"]["value"] == "1/4"
    assert "honest_disclosure" in factors["cell_entropy_factor_one_over_four"]
    assert (
        "icosahedral"
        in factors["projection_target_denominator_12_in_P_over_12"]["source_theorem"]
    )

    branch_scope = cert["branch_scope"]
    assert "triangulated S^2" in branch_scope["screen_branch"]
    assert "product adjoint" in branch_scope["oph_product_gauge_branch"]
    assert "cell-entropy" in branch_scope["honest_disclosure"]

    deps = cert["dependencies"]
    assert deps["screen_sieve_icosahedral_geometric_strengthening"] is True

    screen_sieve_summary = cert["screen_sieve_certificate"]
    assert screen_sieve_summary["orbit_size"] == 12
    assert screen_sieve_summary["gamma_EW"] == "(P/12)*log(N/pi)"
    assert (
        screen_sieve_summary["status"] == "closed_on_declared_triangulated_screen_branch"
    )
