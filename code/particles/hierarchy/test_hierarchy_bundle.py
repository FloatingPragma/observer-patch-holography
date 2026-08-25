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

    assert len(payload) == 14
    assert all(entry["returncode"] == 0 for entry in payload)
    validator_outputs = [json.loads(entry["stdout"]) for entry in payload]
    assert all(output["pass"] is True for output in validator_outputs)


def test_screen_sieve_derives_one_normalized_seidel_graph_class() -> None:
    cert = json.loads(
        (ROOT / "certificates/R_screen_sieve_icosahedral_certificate.json").read_text()
    )
    seidel = cert["d_optimal_tomography"]["seidel_uniqueness"]
    graph_class = seidel["normalized_negative_graph_class"]

    assert graph_class["distinct_degree_profiles"] == [[2, 2, 2, 2, 2]]
    assert graph_class["every_solution_connected"] is True
    assert graph_class["derived_graph_type"] == "C5"
    assert graph_class["solutions_equal_all_labelled_c5_graphs"] is True
    assert graph_class["derived_labelled_c5_graph_count"] == 12
    assert len(graph_class["canonical_edge_codes"]) == 1
    assert graph_class["isomorphism_class_count"] == 1
    assert seidel["switching_isomorphism_classes"] == graph_class[
        "isomorphism_class_count"
    ]

    arithmetic = cert["screen_load_arithmetic"]
    assert arithmetic["local_port_read"] == "X/12"
    assert arithmetic["gamma_screen_simplified"] == "(P/12)*log(N/pi)"
    assert "hierarchy_readout" not in arithmetic
    gate = cert["hierarchy_screen_readout_gate"]
    assert gate["premise_id"] == "HIERARCHY-SCREEN-READOUT"
    assert gate["supplied_by_screen_sieve"] is False
    assert gate["required_identification"] == "log(E_cell/v)=Gamma_screen"


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


def test_ru_outward_rounded_log_reproduces_witness_from_directed_rounding() -> None:
    result = _run(
        "validators/verify_ru_outward_rounded_log.py",
        "certificates/R_U_outward_rounded_interval_log.json",
    )
    payload = json.loads(result.stdout)

    assert payload["pass"] is True
    checks = payload["checks"]
    assert checks["backend_declaration_matches_module"] is True
    assert checks["declared_input_keys_match_allowlist"] is True
    assert checks["computed_block_bit_exact"] is True
    assert checks["endpoint_signs_strict"] is True
    assert checks["derivative_enclosure_strictly_negative"] is True
    assert checks["krawczyk_image_strictly_inside_interior"] is True
    assert checks["high_precision_root_in_outward_witness"] is True
    assert checks["no_measured_endpoint_constants"] is True

    log = json.loads(
        (ROOT / "certificates/R_U_outward_rounded_interval_log.json").read_text()
    )
    assert log["status"] == "outward_rounded_interval_log_supplied"
    assert log["backend"]["number_format"] == "IEEE-754 binary64"
    declared = log["declared_inputs"]["declared"]
    assert declared["P_fwd_decimal"] == "1.630968209403959324879279847782648941"
    assert declared["N2"] == 128 and declared["N3"] == 64
    assert log["checks"]["pass"] is True


def test_ht_interval_box_log_reproduces_declared_surface_enclosure() -> None:
    result = _run(
        "validators/verify_ht_interval_box.py",
        "certificates/R_HT_interval_input_box_log.json",
    )
    payload = json.loads(result.stdout)

    assert payload["pass"] is True
    checks = payload["checks"]
    assert checks["backend_declaration_matches_module"] is True
    assert checks["declared_input_keys_match_allowlist"] is True
    assert checks["resolved_input_blocks_bit_exact"] is True
    assert checks["interval_extension_matches_module"] is True
    assert checks["computed_block_bit_exact"] is True
    assert checks["declared_outputs_inside_box_enclosures"] is True
    assert checks["center_enclosures_inside_box_enclosures"] is True
    assert checks["top_diagonal_strictly_positive"] is True
    assert checks["higgs_diagonal_strictly_negative"] is True
    assert checks["determinant_recomputed_bit_exact"] is True
    assert checks["determinant_excludes_zero"] is True
    assert checks["chart_diagonals_reproduce_from_input_boxes"] is True
    assert checks["uniqueness_statement_scoped"] is True
    assert checks["provenance_classes_allowed"] is True
    assert checks["declared_branch_inputs_match_whitelist"] is True
    assert checks["no_measured_endpoint_constants"] is True

    log = json.loads(
        (ROOT / "certificates/R_HT_interval_input_box_log.json").read_text()
    )
    assert log["status"] == "raw_interval_input_box_and_jacobian_enclosure_supplied"
    assert log["backend"]["number_format"] == "IEEE-754 binary64"
    declared = log["declared_inputs"]["declared"]
    assert declared["box_relative_half_width_decimal"] == "1e-9"
    resolved = log["declared_inputs"]["resolved"]
    assert len(resolved) == 11
    assert all("provenance_class" in entry for entry in resolved.values())
    assert log["checks"]["pass"] is True

    surface = json.loads(
        (ROOT / "certificates/R_HT_declared_surface_certificate.json").read_text()
    )
    assert surface["status"] == "declared_surface_output_recorded_with_raw_interval_input_box"
    cert_block = surface["interval_certification"]
    assert cert_block["raw_interval_input_box_log"] == (
        "certificates/R_HT_interval_input_box_log.json"
    )
    assert cert_block["verifier"] == "validators/verify_ht_interval_box.py"


def test_hierarchy_numeric_witness_keeps_public_and_source_audit_branches_separate() -> None:
    witness = json.loads((ROOT / "computations/hierarchy_numeric_witness.json").read_text())

    public = witness["public_endpoint_branch"]
    source_audit = witness["source_audit_branch"]
    assert public["P_C"] == "1.630968209403959324879279847782648941"
    assert public["alpha_U_display"] == "0.041124336195630495"
    assert public["v_over_E_star"] == "2.0199803239725553e-17"
    assert source_audit["P_cand"] == "1.63097209585889737696451390350695562847912625483895268486516"
    assert source_audit["alpha_U"].startswith("0.041124247441816685")
    assert source_audit["v_over_E_star"].startswith("2.019811407857633059")


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
    assert checks["step_4_scopes_resonance_target"] is True
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
    assert checks["acceptance_resonance_scoped_as_oph_condition"] is True
    assert checks["boundary_records_closed_elsewhere"] is True
    assert checks["boundary_includes_scope_note"] is True
    assert checks["boundary_records_open_strengthening"] is True

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
    assert checks["target_relation_states_zero_bridge_residual"] is True
    assert checks["boundary_avoids_unmatured_physical_gates"] is True
    assert checks["boundary_scope_present"] is True
    assert checks["derivation_chain_has_eight_steps"] is True
    assert checks["step_1_selects_public_endpoint_branch"] is True
    assert checks["step_1_records_branch_locator_a_t_public"] is True
    assert checks["step_1_cites_public_endpoint_pixel_artifact"] is True
    assert checks["step_1_cites_full_precision_joint_artifact"] is True
    assert checks["step_1_records_parallel_source_audit_branch_witness"] is True
    assert checks["step_1_distinguishes_p_public_from_p_source_audit"] is True
    assert checks["step_3_imports_m_rep_24"] is True
    assert checks["step_5_equates_resonance_to_bridge_residual"] is True
    assert checks["step_6_solves_closed_form_capacity"] is True
    assert checks["step_7_certifies_banach_contraction"] is True
    assert checks["step_8_records_numerical_witness_and_rejects_rounded"] is True
    assert checks["factor_origin_p_star_branch_is_public_endpoint"] is True
    assert checks["factor_origin_p_star_source_is_public_endpoint_pixel_cert"] is True
    assert checks["factor_origin_p_star_full_precision_source_is_joint_cert"] is True
    assert checks["factor_origin_p_star_records_parallel_source_audit_witness"] is True
    assert checks["factor_origin_p_star_value_matches_p_public"] is True
    assert checks["factor_origin_alpha_u_source_is_krawczyk_cert"] is True
    assert checks["factor_origin_six_recorded_as_m_rep_over_beta_ew"] is True
    assert checks["factor_origin_lambda_one_half_recorded"] is True
    assert checks["branch_scope_public_endpoint_pixel_branch_present"] is True
    assert checks["branch_scope_krawczyk_unification_width_branch_present"] is True
    assert checks["branch_scope_records_parallel_source_audit_branch"] is True
    assert checks["branch_scope_note_present"] is True
    assert checks["source_values_record_p_star_branch"] is True
    assert checks["branch_selection_block_present_and_public_endpoint"] is True
    assert checks["allowed_inputs_record_a_t_public_as_branch_locator"] is True
    assert checks["forbidden_calibrations_block_a_t_public_as_upstream_input"] is True
    assert checks["dependency_artifacts_public_endpoint_pixel_closure_present"] is True
    assert checks["dependency_artifacts_full_precision_joint_present"] is True
    assert checks["dependency_artifacts_parallel_source_audit_present"] is True
    assert checks["dependency_artifacts_m_rep_present"] is True
    assert checks["dependency_artifacts_pi_ew_present"] is True
    assert checks["consumer_artifacts_umbrella_present"] is True
    assert checks["acyclicity_summary_states_peer_cross_reference"] is True
    assert checks["acyclicity_primary_theorems_independent"] is True

    cert = json.loads((ROOT / "certificates/R_EW_global_capacity_certificate.json").read_text())
    assert cert["accepted"] is True
    assert cert["status"] == "closed_bridge_refined_global_capacity_fixed_point_certificate"
    assert cert["certificate_id"] == "issue-344-exact-ew-refined-global-capacity-v2"
    assert cert["contraction_certificate"]["lipschitz_constant"] == "0.5"
    assert cert["exact_capacity_fixed_point"]["bridge_residual"] == "0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"

    chain = cert["derivation_chain"]
    assert len(chain) == 8
    assert {step["step"] for step in chain} == set(range(1, 9))
    assert chain[0]["branch_selection"] == "public_endpoint_branch"
    assert "A_T_public = 137.035999177" in chain[0]["branch_locator"]
    assert chain[0]["source_artifact"] == "certificates/R_P_public_pixel_certificate.json"
    assert chain[0]["full_precision_source_artifact"] == "certificates/R_PN_joint_fixed_point_certificate_report.json"
    assert chain[0]["parallel_source_audit_branch_witness"] == "certificates/R_P_source_audit_pixel_certificate.json"
    assert "P_public" in chain[0]["conclusion"]
    assert "P_cand = 1.63097209585889737696451390350695562847912625483895268486516" in chain[0]["conclusion"]
    assert "beta_EW" in chain[1]["conclusion"]
    assert "m_rep = 24" in chain[2]["conclusion"]
    assert "OPH local/global EW resonance condition" in chain[3]["conclusion"]
    assert "source-side criterion" in chain[3]["conclusion"]
    assert "B_EW(P_star,N) := alpha_U*log(N/pi) - 6*pi/P_star" in chain[4]["conclusion"]
    assert "N_CRC^EW(P_star) = pi*exp[6*pi/(P_star*alpha_U(P_star))]" in chain[5]["conclusion"]
    assert "Banach" in chain[6]["conclusion"] and "1 - lambda = 1/2" in chain[6]["conclusion"]
    assert "3.31e122" in chain[7]["conclusion"]

    factors = cert["factor_origins"]
    p_star_factor = factors["P_star_pixel_fixed_point"]
    assert p_star_factor["branch"] == "public_endpoint_branch"
    assert p_star_factor["source_artifact"] == "certificates/R_P_public_pixel_certificate.json"
    assert p_star_factor["full_precision_source_artifact"] == "certificates/R_PN_joint_fixed_point_certificate_report.json"
    assert p_star_factor["parallel_source_audit_witness"] == "certificates/R_P_source_audit_pixel_certificate.json"
    assert p_star_factor["parallel_source_audit_value"] == "1.63097209585889737696451390350695562847912625483895268486516"
    assert p_star_factor["value"] == "1.6309682094039593248792798477826489413359828516279250606661507533907793398933432"
    assert factors["alpha_U_unification_width"]["source_artifact"] == "certificates/R_U_krawczyk_certificate.json"
    assert factors["beta_EW_transmutation_multiplicity"]["value"] == "4"
    assert factors["m_rep_doubled_sm_adjoint_round_count"]["value"] == "24"
    assert factors["factor_six_in_bridge_residual"]["value"] == "6"
    assert "m_rep / beta_EW" in factors["factor_six_in_bridge_residual"]["expression"]
    assert factors["banach_contraction_lambda_one_half"]["value"] == "1/2"

    branch_scope = cert["branch_scope"]
    assert "scope_note" in branch_scope
    assert "N_CRC^EW(P_star) = pi*exp[6*pi/(P_star*alpha_U(P_star))]" in branch_scope["scope_note"]
    assert "public-endpoint pixel" in branch_scope["scope_note"]
    assert "A_T_public = 137.035999177" in branch_scope["scope_note"]
    assert "P_public" in branch_scope["public_endpoint_pixel_branch"]
    assert "R_P_public_pixel_certificate.json" in branch_scope["public_endpoint_pixel_branch"]
    assert "R_U_krawczyk_certificate.json" in branch_scope["krawczyk_unification_width_branch"]
    assert "R_P_source_audit_pixel_certificate.json" in branch_scope["parallel_source_audit_branch_note"]
    assert "lambda = 1/2" in branch_scope["banach_contraction_branch"]

    branch_selection = cert["branch_selection"]
    assert branch_selection["selected_branch"] == "public_endpoint_pixel_branch"
    assert "A_T_public = 137.035999177" in branch_selection["branch_locator"]
    assert "branch locator" in branch_selection["branch_locator"]
    assert "R_P_source_audit_pixel_certificate.json" in branch_selection["parallel_source_audit_branch"]

    source_values = cert["source_values"]
    assert source_values["P_star_branch"] == "public_endpoint_branch"
    assert source_values["P_star_source_artifact"] == "certificates/R_P_public_pixel_certificate.json"
    assert source_values["P_star_full_precision_source_artifact"] == "certificates/R_PN_joint_fixed_point_certificate_report.json"
    assert source_values["alpha_U_source_artifact"] == "certificates/R_U_krawczyk_certificate.json"

    allowed = cert["allowed_inputs"]
    assert any("A_T_public = 137.035999177" in item and "branch locator" in item for item in allowed)
    forbidden = cert["forbidden_calibrations"]
    assert any("A_T_public" in item and "upstream source-map input" in item for item in forbidden)
    boundary = cert["claim_boundary"]
    assert boundary["not_closed_here"] == []
    assert "diagnostic-only" in boundary["scope"]

    acceptance = cert["acceptance_criteria_status"]
    assert all(acceptance.values())
    assert acceptance["bridge_residual_zero_on_source_side"] is True
    assert acceptance["banach_contraction_certified_with_explicit_lipschitz_constant"] is True
    assert acceptance["rounded_capacity_display_rejected_as_exact_witness"] is True

    deps = cert["dependency_artifacts"]
    assert deps["public_endpoint_local_pixel_closure"] == "certificates/R_P_public_pixel_certificate.json"
    assert deps["public_endpoint_local_pixel_full_precision_record"] == "certificates/R_PN_joint_fixed_point_certificate_report.json"
    assert deps["alpha_u_unification_width"] == "certificates/R_U_krawczyk_certificate.json"
    assert deps["parallel_source_audit_pixel_branch_witness"] == "certificates/R_P_source_audit_pixel_certificate.json"
    assert "R_m_rep_24_certificate" in deps["representation_to_spectrum_m_rep_24"]
    assert "R_EW_tick_projection_certificate" in deps["ew_tick_projection_pi_ew_definition"]

    consumers = cert["consumer_artifacts"]
    assert "R_EW_tick_projection_certificate" in consumers["ew_tick_projection_specialisation"]
    assert "R_readback_resolution_certificate" in consumers["finite_readback_resolution_dependency"]
    assert "R_local_global_hierarchy_resonance_closeout_335" in consumers["local_global_hierarchy_resonance_umbrella"]

    acyclic = cert["dependency_acyclicity_note"]
    assert "peer cross-reference" in acyclic["summary"]
    assert "ew_tick_projection_primary" in acyclic["primary_theorems_are_independent"]
    assert "exact_capacity_primary" in acyclic["primary_theorems_are_independent"]
    assert "specialised_corollary_is_a_composition_not_a_circle" in acyclic
    assert "umbrella_certificate_resolves_the_composition" in acyclic


def test_issue_342_readback_resolution_certificate_is_scope_limited() -> None:
    result = _run(
        "validators/validate_issue_342_readback_resolution.py",
        "certificates/R_readback_resolution_certificate.json",
    )
    payload = json.loads(result.stdout)

    assert payload["pass"] is True
    checks = payload["checks"]
    assert checks["status_is_closed_finite_readback"] is True
    assert checks["normal_form_is_unique_at_fixed_cutoff"] is True
    assert checks["observer_sector_is_stable_and_central"] is True
    assert checks["one_positive_sharp_capacity_atom"] is True
    assert checks["positive_root_extractor"] is True
    assert checks["exact_capacity_readback_matches"] is True
    assert checks["refinement_statement_is_explicit"] is True
    assert checks["rounded_display_is_diagnostic_only"] is True
    assert checks["scope_is_conjunction_limited"] is True
    assert checks["does_not_promote_full_hierarchy"] is True

    cert = json.loads((ROOT / "certificates/R_readback_resolution_certificate.json").read_text())
    assert cert["accepted"] is True
    assert cert["status"] == "closed_finite_readback_resolution_certificate"
    assert cert["claim_boundary"]["not_closed_here"] == []
    assert cert["capacity_register"]["selected_variance"] == "0"

    assert cert["dependencies"]["ew_refined_exact_capacity"] is True
    assert (
        "R_EW_global_capacity_certificate.json"
        in cert["dependency_artifacts"]["ew_refined_exact_capacity"]
    )


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
    assert checks["doubled_su5_rejected"] is True
    assert checks["graviton_excluded"] is True
    assert checks["no_forbidden_inputs_used"] is True
    assert checks["target_relation_recorded"] is True
    assert checks["full_cycle_decomposition_recorded"] is True
    assert checks["claim_boundary_records_scope"] is True
    assert checks["derivation_chain_has_eight_steps"] is True
    assert checks["derivation_step1_product_adjoint_branch"] is True
    assert checks["derivation_step4_orientation_doubling"] is True
    assert checks["derivation_step7_specializes_tick_law"] is True
    assert checks["every_derivation_step_has_source_artifact"] is True
    assert checks["tex_sources_use_stable_sections_or_labels_not_line_numbers"] is True
    assert checks["factor_origins_keys_complete"] is True
    assert checks["factor_origin_unoriented_cites_corpus"] is True
    assert checks["factor_origin_orientation_cites_corpus"] is True
    assert checks["factor_origin_exponent_cites_global_repair_tick_cert"] is True
    assert checks["branch_scope_keys_complete"] is True
    assert checks["branch_scope_records_patch_carrier_pipeline"] is True
    assert checks["branch_scope_includes_scope_note"] is True
    assert checks["dependency_artifacts_keys_complete"] is True
    assert checks["consumer_artifacts_keys_complete"] is True
    assert checks["acceptance_criteria_keys_complete"] is True
    assert checks["acceptance_criteria_all_satisfied"] is True
    assert checks["used_inputs_cite_compact_proof_corpus"] is True
    assert checks["certificate_id_v2"] is True

    cert = json.loads((ROOT / "certificates/R_m_rep_24_certificate.json").read_text())
    assert cert["accepted"] is True
    assert cert["status"] == "closed_representation_to_spectrum_round_count"
    assert cert["certificate_id"] == "issue-343-m-rep-24-doubled-sm-adjoint-v2"
    assert cert["representation_sector"]["unoriented_adjoint_dimension"] == 12
    assert cert["representation_sector"]["oriented_support_dimension"] == 24
    assert cert["result"]["specialized_exponent"] == "-1/48"
    assert cert["result"]["full_cycle_decomposition"] == "G_N = g_N^m_rep = g_N^24"
    assert cert["claim_boundary"]["not_closed_here"] == []
    assert "(N_CRC/pi)^(-1/48)" in cert["claim_boundary"]["scope"]

    chain = cert["derivation_chain"]
    assert len(chain) == 8
    assert chain[0]["premise"].startswith(
        "Observer-visible product-adjoint branch"
    )
    assert "compact_proof_of_oph.tex" in chain[0]["source_artifact"]
    assert "section 'The compact-gauge branch'" in chain[0]["source_artifact"]
    assert "line 482" not in chain[0]["source_artifact"]
    assert "orientation-doubling grammar" in chain[3]["premise"]
    assert "compact_proof_of_oph.tex" in chain[3]["source_artifact"]
    assert "section 'The QCD-free hierarchy witness'" in chain[3]["source_artifact"]
    assert (
        chain[6]["premise"]
        == "Specialization of the parametric global repair-tick law"
    )
    assert (
        "R_N_global_repair_tick_certificate.json" in chain[6]["source_artifact"]
    )

    factor_origins = cert["factor_origins"]
    assert factor_origins["m_rep"]["value"] == 24
    assert factor_origins["exponent_denominator"]["value"] == 48
    assert (
        factor_origins["orientation_multiplier"]["value"] == 2
    )
    assert (
        "compact_proof_of_oph.tex"
        in factor_origins["orientation_multiplier"]["source_artifact"]
    )

    branch_scope = cert["branch_scope"]
    assert "su(3)+su(2)+u(1)" in branch_scope["product_adjoint_branch"]
    assert "global quotient" in branch_scope["product_adjoint_branch"]
    assert "patch-carrier" in branch_scope["reversible_repair_orientation_branch"]
    assert "m_rep=24" in branch_scope["scope_note"]

    consumers = cert["consumer_artifacts"]
    assert (
        "R_local_global_hierarchy_resonance_closeout_335.json"
        in consumers["local_global_resonance_closeout"]
    )

    acyclicity = cert["dependency_acyclicity_note"]
    assert "peer cross-reference" in acyclicity["summary"]
    assert "not a circular dependency" in acyclicity["summary"]
    primary = acyclicity["primary_theorems_are_independent"]
    assert "m is a free parameter" in primary["global_repair_tick_lemma_primary"]
    assert "do not use the tick law" in primary["m_rep_24_primary"]
    assert (
        "composition of the two independent primary theorems"
        in acyclicity["specialized_corollary_is_a_composition_not_a_circle"]
    )
    assert (
        "R_local_global_hierarchy_resonance_closeout_335.json"
        in acyclicity["umbrella_certificate_resolves_the_composition"]
    )
    assert (
        "strictly upstream"
        in acyclicity["other_branches_are_upstream_only"]
    )

    acceptance = cert["acceptance_criteria_status"]
    assert acceptance["all_acceptance_criteria_satisfied"] is True
    assert (
        acceptance[
            "parametric_tick_law_specializes_to_minus_one_over_48_at_m_rep_24"
        ]
        is True
    )
    assert acceptance["m_rep_24_proved_on_source_side_oph_data"] is True
    assert acceptance["negative_controls_for_nearby_round_counts_supplied"] is True
    assert (
        acceptance[
            "no_measured_weak_higgs_g_planck_area_lambda_or_hierarchy_ratio_inputs_used"
        ]
        is True
    )
    assert (
        acceptance["factor_origins_documented_for_every_numerical_factor"] is True
    )


def test_issue_332_declared_map_identity_is_zero_defect_and_not_promoted() -> None:
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
    assert cert["mode"] == "exact_declared_comparison_map_identity_not_source_naturality"
    assert cert["claim_boundary"]["receipt_class"] == "conditional_identity"
    assert "Declared-map" in cert["theorem"]
    assert cert["source_comparison_maps_derived"] is False
    assert cert["physical_naturality_result"] is False
    assert cert["predictive_promotion_allowed"] is False
    assert cert["claim_boundary"]["antecedent_only_countermodel_receipt"] == (
        "certificates/antecedent_only_nonidentifiability_receipt.json"
    )
    forbidden = cert["forbidden_calibrations"]
    assert any("measured weak scale" in item for item in forbidden)
    assert any("measured Higgs" in item for item in forbidden)


def test_naturality_boundary_is_consistent_across_paper_registry_and_status_sources() -> None:
    repo = ROOT.parents[2]
    paper = (repo / "paper" / "deriving_the_particle_zoo_from_observer_consistency.tex").read_text(
        encoding="utf-8"
    )
    registry = (repo / "claims" / "claim_registry.yaml").read_text(encoding="utf-8")
    ledger = (ROOT.parent / "ledger.yaml").read_text(encoding="utf-8")

    assert "declared-map packaging identity" in paper
    assert "source-derived comparison maps" in paper.lower()
    assert "The current packet supplies no OPH result about Higgs naturalness" in paper
    assert "zero Higgs naturality defect in its declared chart" not in paper
    assert "naturality bridge is closed on its selected branch" not in paper
    assert "declared_naturality_identity_quarantined" in registry
    assert "same-source-antecedent" in ledger
    assert "not a source-derived Higgs-naturalness result" in ledger


def test_issue_335_local_global_resonance_is_exact_conditional_statement() -> None:
    result = _run(
        "validators/validate_issue_335_local_global_resonance.py",
        "certificates/R_local_global_hierarchy_resonance_closeout_335.json",
    )
    payload = json.loads(result.stdout)

    assert payload["pass"] is True
    checks = payload["checks"]
    assert checks["full_theorem_not_promoted"] is True
    assert checks["exact_capacity_supplied"] is True
    assert checks["finite_readback_supplied"] is True
    assert checks["round_count_supplied"] is True
    assert checks["screen_sieve_supplied"] is True
    assert checks["screen_sieve_dependency_present"] is True
    assert checks["rounded_capacity_rejected"] is True
    assert checks["screen_source_receipt_is_work_in_progress"] is True
    assert checks["derivation_chain_has_nine_steps"] is True
    assert checks["step_5_is_screen_sieve_arithmetic_only"] is True
    assert checks["step_8_composes_conditional_target_relation"] is True
    assert checks["factor_origin_icosahedral_orbit_recorded"] is True
    assert checks["factor_origin_cell_entropy_scoped"] is True
    assert checks["branch_scope_records_screen_branch"] is True
    assert checks["branch_scope_records_hierarchy_readout_premise"] is True
    assert checks["residual_residue_scoped_in_acceptance"] is True

    cert = json.loads((ROOT / "certificates/R_local_global_hierarchy_resonance_closeout_335.json").read_text())
    assert cert["accepted"] is True
    assert cert["status"] == "exact_conditional_local_global_hierarchy_resonance"
    assert cert["full_theorem_grade_resonance_promoted"] is False
    acceptance = cert["acceptance_criteria_status"]
    assert acceptance["prerequisite_steps_accounted_for"] is True
    assert acceptance["exact_conditional_resonance_proved"] is True
    assert acceptance["full_theorem_grade_resonance_proved"] is False
    assert acceptance["exact_capacity_source_certificate_supplied"] is True
    assert acceptance["finite_readback_resolution_supplied"] is True
    assert acceptance["round_count_derivation_supplied"] is True
    assert acceptance["screen_sieve_geometric_strengthening_supplied"] is True
    assert acceptance["hierarchy_screen_readout_premise_declared"] is True
    assert acceptance["hierarchy_screen_readout_closed"] is False
    assert acceptance["residual_definitional_residue_scoped_as_oph_identification"] is True
    assert "P/beta_EW" in acceptance["residual_definitional_residue_scope_note"]
    assert cert["work_in_progress_receipts"] == [
        "source production of the unit cost, inverse pairing, and D-optimal selector",
        "HIERARCHY-SCREEN-READOUT: identify log(E_cell/v)=Gamma_screen and attach it to the alpha_U/B_EW branch",
    ]

    chain = cert["derivation_chain"]
    assert len(chain) == 9
    assert chain[0]["premise"] == "D10 forward transmutation theorem"
    assert "icosahedral screen-sieve" in chain[4]["premise"]
    assert "does not identify log(E_cell/v)" in chain[4]["scope_note"]
    assert "(P/12)" in chain[4]["conclusion"]
    assert "HIERARCHY-SCREEN-READOUT" in chain[5]["premise"]
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
    assert "scope_note" in factors["cell_entropy_factor_one_over_four"]
    assert (
        "icosahedral"
        in factors["projection_target_denominator_12_in_P_over_12"]["source_theorem"]
    )

    branch_scope = cert["branch_scope"]
    assert "triangulated S^2" in branch_scope["screen_branch"]
    assert "HIERARCHY-SCREEN-READOUT" in branch_scope["hierarchy_screen_readout_branch"]
    assert "product adjoint" in branch_scope["oph_product_gauge_branch"]
    assert "cell-entropy" in branch_scope["scope_note"]

    deps = cert["dependencies"]
    assert deps["screen_sieve_icosahedral_geometric_strengthening"] is True

    screen_sieve_summary = cert["screen_sieve_certificate"]
    assert screen_sieve_summary["orbit_size"] == 12
    assert screen_sieve_summary["local_port_read"] == "X/12"
    assert screen_sieve_summary["gamma_screen_algebra"] == "(P/12)*log(N/pi)"
    assert screen_sieve_summary["hierarchy_readout_premise"]["premise_id"] == (
        "HIERARCHY-SCREEN-READOUT"
    )
    assert (
        screen_sieve_summary["status"] == "conditional_finite_selector_theorem"
    )


def test_pixel_screen_resonance_summary_records_tile_and_cosmology_identities() -> None:
    result = _run("verify_pixel_screen_resonance_summary.py", "--check")
    payload = json.loads(result.stdout)

    assert payload["accepted"] is True
    assert payload["status"] == "closed_receipt_summary_from_existing_hierarchy_certificates"
    assert payload["source_pair"]["P_star_branch"] == "public_endpoint_branch"
    assert payload["source_pair"]["N_CRC_EW_source_artifact"] == (
        "certificates/R_EW_global_capacity_certificate.json"
    )

    tile = payload["pixel_screen_tile_identity"]
    assert tile["cell_entropy_formula"] == "ell_cell = P_star/4"
    assert tile["cell_count_formula"] == "K_cell = N_CRC^EW / (P_star/4) = 4*N_CRC^EW/P_star"
    assert tile["reconstruction_formula"] == "K_cell*(P_star/4) = N_CRC^EW"
    assert tile["relative_reconstruction_error"] == "0"

    alignment = payload["shared_12_24_port_alignment"]
    assert alignment["screen_ports"] == 12
    assert alignment["screen_oriented_slots"] == 24
    assert alignment["product_adjoint_dimension"] == 12
    assert alignment["product_adjoint_rounds_m_rep"] == 24
    assert alignment["equal_cardinality_status"] == (
        "arithmetic_alignment_without_physical_identification"
    )
    assert "screen_microphysics" in alignment["screen_register_source"]
    assert "R_m_rep_24_certificate" in alignment["product_adjoint_source"]
    assert alignment["screen_load_split"] == "X/12"
    assert alignment["identification_premise"] == "HIERARCHY-SCREEN-READOUT"
    assert alignment["identification_status"] == "work_in_progress"
    assert alignment["conditional_hierarchy_readout"] == (
        "v/E_cell = (N_CRC^EW/pi)^(-P_star/12)"
    )

    ds = payload["dimensionless_de_sitter_coordinate"]
    assert ds["Lambda_lstar2_formula"] == "Lambda_CRC*l_star^2 = 3*pi/N_CRC^EW"
    assert ds["Lambda_a_cell_formula"] == "Lambda_CRC*a_cell = 3*pi*P_star/N_CRC^EW"
    assert ds["Lambda_a_cell_count_formula"] == "Lambda_CRC*a_cell = 12*pi/K_cell"
    assert ds["relative_cell_coordinate_error"] == "0"

    boundary = payload["claim_boundary"]
    assert "primitive one-cell" in str(boundary["not_closed_here"])
    assert "SI Lambda" in str(boundary["not_closed_here"])
