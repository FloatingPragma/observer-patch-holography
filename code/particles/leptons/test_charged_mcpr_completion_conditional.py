"""Tests for the conditional MCPR charged-response completion lane."""

from __future__ import annotations

import json

import mpmath as mp

import derive_charged_mcpr_completion_conditional as lane


import pytest


@pytest.fixture(autouse=True)
def _scoped_precision():
    """Set a high mpmath precision for each test and restore it on teardown.

    The tight-tolerance comparisons in this module require extended precision.
    Scoping keeps the global mpmath precision unchanged for other test modules.
    """

    with mp.workdps(120):
        yield


def test_icosahedral_incidence_solve_is_unique():
    geom = lane.icosahedral_incidence_solve()
    assert geom["unique"] is True
    assert (geom["V"], geom["E"], geom["F"]) == (12, 30, 20)


def test_register_dimensions_match_declared_architecture():
    geom = lane.icosahedral_incidence_solve()
    dims = lane.register_dimensions(geom)
    assert list(dims.values()) == lane.EXPECTED_REGISTER_DIMENSIONS


def test_fixed_point_matches_declared_coefficients():
    geom = lane.icosahedral_incidence_solve()
    dims = lane.register_dimensions(geom)
    source = json.loads(lane.SOURCE_CERTIFICATE.read_text(encoding="utf-8"))
    p_value = mp.mpf(str(source["P_cand"]))
    alpha_u = mp.mpf(str(source["alpha_U_P_cand"]))
    fp = lane.response_fixed_point(p_value, alpha_u, dims)
    assert abs(fp["kappa"] - mp.mpf("0.0008213908200217955948997382298787794")) < mp.mpf("1e-37")
    assert fp["contraction"] < mp.mpf("0.0016")
    assert max(abs(x) for x in fp["residual"]) < mp.mpf("1e-100")


def test_artifact_reproduces_conditional_triple_and_stays_fail_closed():
    artifact = lane.build()
    assert artifact["checks_pass"] is True
    assert artifact["source_only"] is False
    assert artifact["public_charged_mass_promotion_allowed"] is False
    assert artifact["charged_reference_data_consumed"] is False
    assert artifact["charged_reference_data_consumed_scope"] == "runtime_execution_only"
    assert artifact["row_class"] == "conditional_on_declared_mcpr_response_architecture"
    claim_class = artifact["claim_class"]
    assert claim_class["runtime_continuous_fit_parameter_count"] == 0
    assert claim_class["declared_architecture_contains_stipulated_continuous_choices"] is True
    assert claim_class["historical_target_informed_architecture"] is True
    assert claim_class["historically_blind_prospective_prediction"] is False

    provenance = artifact["provenance_boundary"]
    assert provenance["runtime_reference_separation"]["charged_reference_packet_consumed"] is False
    assert provenance["historical_selection"]["charged_target_exposure"] is True
    assert provenance["historical_selection"]["architecture_target_informed"] is True
    assert (
        provenance["historical_selection"]["prospective_freeze_before_charged_evaluation"]
        is False
    )
    assert len(provenance["historical_selection"]["evidence"]) == 2
    assert len(provenance["stipulated_architecture_inputs"]) >= 8
    assert "at runtime" in artifact["runtime_no_target_dependency_statement"]
    assert "historical" in artifact["runtime_no_target_dependency_statement"]

    ratios = artifact["conditional_prediction"]["ratios"]
    assert abs(mp.mpf(ratios["m_mu_over_m_e"]) - mp.mpf("206.76830059518097174")) < mp.mpf("1e-12")
    assert abs(mp.mpf(ratios["m_tau_over_m_e"]) - mp.mpf("3477.3655365960225913")) < mp.mpf("1e-9")
    assert abs(mp.mpf(ratios["m_tau_over_m_mu"]) - mp.mpf("16.81769171863604088")) < mp.mpf("1e-12")

    masses = [mp.mpf(x) for x in artifact["conditional_prediction"]["masses_over_E_star"]]
    expected = [
        mp.mpf("4.185110648223775e-23"),
        mp.mpf("8.653482165360261e-21"),
        mp.mpf("1.455315953497439e-19"),
    ]
    for got, want in zip(masses, expected, strict=True):
        assert abs(got / want - 1) < mp.mpf("1e-14")


def test_audit_boundary_gates_are_open():
    artifact = lane.build()
    route = artifact["audit_boundary"]["completion_route"]
    assert set(route) == {"W5_ORB", "A5_FAM", "DET_CAN", "QFT_POLE", "SCALE", "PROV"}
    assert all(entry["status"] == "open" for entry in route.values())
    assert artifact["audit_boundary"]["verdict"] == "SOURCE_ONLY_PHYSICAL_PREDICTION_OPEN"


def test_display_column_is_annotated_as_checksum_scale():
    artifact = lane.build()
    display = artifact["optional_scale_display"]
    assert "calibration checksum" in display["display_scale_status"]
    assert len(display["masses_MeV"]) == 3
    assert abs(mp.mpf(display["masses_MeV"][0]) - mp.mpf("0.510955973930992")) < mp.mpf("1e-12")
