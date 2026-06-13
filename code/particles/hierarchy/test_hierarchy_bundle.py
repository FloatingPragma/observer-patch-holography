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

    assert len(payload) == 6
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


def test_global_repair_tick_lemma_is_closed_on_declared_rounds_with_claim_boundary() -> None:
    result = _run(
        "validators/validate_global_repair_tick_certificate.py",
        "certificates/R_N_global_repair_tick_certificate.json",
    )
    payload = json.loads(result.stdout)

    assert payload["pass"] is True
    assert payload["status_is_closed_lemma_on_declared_rounds"] is True
    assert payload["tick_exponent_is_minus_one_over_48"] is True
    assert payload["tick_times_rounds_is_full_cycle"] is True
    assert payload["full_cycle_map_recorded"] is True
    assert payload["one_tick_map_recorded"] is True
    assert payload["numeric_display_matches_formula"] is True
    assert payload["full_cycle_multiplier_is_derived_from_closure"] is True
    assert payload["f_interface_equivalence_derived"] is True
    assert payload["round_count_is_declared_not_derived"] is True
    assert payload["ew_inputs_excluded_from_derived_uses"] is True

    cert = json.loads((ROOT / "certificates/R_N_global_repair_tick_certificate.json").read_text())
    assert cert["status"] == "closed_global_repair_tick_lemma_on_declared_round_structure"
    assert cert["normalization"]["abs_g_star_prime"] == "(N_CRC/pi)^(-1/48)"
    assert cert["exponent_law"]["per_tick_exponent_for_m_ticks"] == "-1/(2m)"
    acceptance = cert["acceptance_criteria_status"]
    assert acceptance["proves_declared_screen_capacity_fixed_point_emits_tick_contraction"] is True
    assert acceptance["round_count_derived_from_first_principles"] is False
    assert acceptance["closure_transport_derived_from_F_interface"] is True
    assert acceptance["readback_counting_model_is_modeling_identification"] is True
    assert acceptance["concrete_finite_machinery_verification_open"] is True
    declared = cert["claim_boundary"]["declared_not_derived"]
    assert any("modeling identification" in item for item in declared)
    boundary = cert["claim_boundary"]["not_closed_by_certificate"]
    assert any("round count" in item for item in boundary)
    assert any("electroweak tick-projection lemma" in item for item in boundary)
    assert any("finite repair machinery" in item for item in boundary)


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
