"""Tests for the fail-closed clock-to-source-energy closure audit."""

from __future__ import annotations

from decimal import Decimal, localcontext

import derive_clock_source_energy_closure_audit as lane


def _synthetic_certificate(epsilon: str) -> dict:
    return {
        "status": "synthetic_in_memory_certificate",
        "candidate_values_from_sources": {"epsilon_Cs": epsilon},
    }


def test_conversion_arithmetic_on_synthetic_epsilon():
    artifact = lane.build_artifact(certificate=_synthetic_certificate("1e-33"))
    receipt = artifact["clock_certificate_receipt"]
    with localcontext() as ctx:
        ctx.prec = lane.WORKING_PREC
        expected_j = lane.H_SI * lane.NU_CS_HZ / Decimal("1e-33")
        expected_gev = lane.H_SI * lane.NU_CS_HZ / (Decimal("1e-33") * lane.GEV_J)
    assert Decimal(receipt["candidate_E_star_J"]) == expected_j
    assert Decimal(receipt["candidate_E_star_GeV"]) == expected_gev
    assert receipt["path"] == "in_memory_certificate"
    assert artifact["checks"]["conversion_round_trip_residual_below_1e_60"] is True


def test_candidate_scale_helper_matches_exact_formulas():
    with localcontext() as ctx:
        ctx.prec = lane.WORKING_PREC
        e_j, e_gev = lane.candidate_scale_from_gap(Decimal("1e-33"))
        assert e_j == lane.H_SI * lane.NU_CS_HZ / Decimal("1e-33")
        assert e_gev == lane.H_SI * lane.NU_CS_HZ / (Decimal("1e-33") * lane.GEV_J)


def test_component_ledger_against_live_repository():
    artifact = lane.build_artifact()
    ledger = artifact["component_ledger"]
    assert ledger["R_U"] == "supplied"
    for name in ("R_alpha", "R_e_abs", "R_QCD_nuc_133Cs", "R_atom_133Cs"):
        assert ledger[name] == "absent"
    assert artifact["component_producer_files"]["R_U"] == [
        "R_U_interval_certificate.json",
        "R_U_krawczyk_certificate.json",
    ]
    assert artifact["components_absent"] == [
        "R_QCD_nuc_133Cs",
        "R_alpha",
        "R_atom_133Cs",
        "R_e_abs",
    ]


def test_fail_closed_verdict_and_promotion_gate():
    artifact = lane.build_artifact()
    assert artifact["verdict"] == "NOT_FINAL_SOURCE_ONLY_PACKET"
    assert artifact["final_packet_ready"] is False
    assert artifact["promotion_allowed"] is False
    assert artifact["application_status"] == "OPEN_BECAUSE_EPSILON_CS_IS_NOT_SOURCE_EMITTED"
    receipt = artifact["clock_certificate_receipt"]
    assert receipt["epsilon_role"] == "calibration_checksum_reproducing_displayed_G"
    assert receipt["displayed_G_comparison_role"] == "taint_diagnosis_never_a_prediction"
    assert artifact["checks_pass"] is True


def test_theorem_block_names_and_required_next_source_objects():
    artifact = lane.build_artifact()
    theorem = artifact["exact_conversion_theorem"]
    assert theorem["theorem"] == "CLOCK_TO_SOURCE_ENERGY_THEOREM"
    assert theorem["alias"] == "CLOCK_YO_SOURCE_ENERGY_THEOREM"
    assert theorem["conclusion"] == "E_star = h*nu_clk/epsilon_clk"
    assert "h*nu_clk/eps_hi" in theorem["interval_corollary"]
    assert "r_q*h*nu_clk/(epsilon_clk*J_GeV)" in theorem["quark_corollary"]
    assert theorem["application_status"] == "OPEN_BECAUSE_EPSILON_CS_IS_NOT_SOURCE_EMITTED"
    assert artifact["required_next_source_objects"] == list(lane.REQUIRED_NEXT_SOURCE_OBJECTS)
    assert "R_ALPHA_ATOMIC_SCHEME_PACKET" in artifact["required_next_source_objects"]


def test_candidate_e_star_and_g_identity_against_live_certificate():
    artifact = lane.build_artifact()
    receipt = artifact["clock_certificate_receipt"]
    assert receipt["candidate_epsilon_Cs"].lower().startswith("3.113930513416012823901050434132426029496608693041876")
    assert receipt["candidate_E_star_GeV"].startswith("12208901289579269")
    g_value = Decimal(receipt["G_algebraically_inferred_from_candidate_epsilon"])
    g_display = Decimal(lane.G_DISPLAY_SI)
    assert abs(g_value - g_display) / g_display < Decimal("1e-6")
    assert receipt["reproduces_displayed_G_within_1e_6"] is True
