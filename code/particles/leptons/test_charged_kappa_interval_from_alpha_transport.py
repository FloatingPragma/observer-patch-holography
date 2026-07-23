"""Focused tests for the charged kappa-interval transport inversion."""

from __future__ import annotations

import json
import math

import pytest

import derive_charged_kappa_interval_from_alpha_transport as lane


@pytest.fixture(scope="module")
def result(tmp_path_factory):
    out = tmp_path_factory.mktemp("kappa") / "charged_kappa_interval.json"
    return lane.build(out)


def test_packet_inversion_round_trip():
    ratios = (206.76828270846718, 3477.3652619063414)
    m_e = 0.00051099895069
    packet = lane.lepton_packet_asymptotic(m_e, ratios)
    assert lane.invert_packet_for_m_e(packet, ratios) == pytest.approx(m_e, rel=1e-12)


def test_packet_derivative_is_minus_two_over_pi():
    ratios = (206.76828270846718, 3477.3652619063414)
    m_e = 0.00051099895069
    kappa = 1e-6
    forward = lane.lepton_packet_asymptotic(m_e * math.exp(kappa), ratios)
    base = lane.lepton_packet_asymptotic(m_e, ratios)
    assert (forward - base) / kappa == pytest.approx(-2.0 / math.pi, rel=1e-6)


def test_kappa_interval_is_finite_and_ordered(result):
    lo, hi = result["kappa_interval"]["interval"]
    assert math.isfinite(lo) and math.isfinite(hi)
    assert lo < result["kappa_interval"]["central_gap_midpoint"] < hi
    assert hi - lo < 1.0, "interval no longer excludes the continuum meaningfully"


def test_witness_masses_inside_certified_intervals(result):
    assert result["compare_only"]["witness_inside_certified_intervals"] is True
    for row, ref in zip(
        result["conditional_mass_rows"],
        result["compare_only"]["witness_masses_gev"],
        strict=True,
    ):
        lo, hi = row["mass_interval"]
        assert lo < ref < hi


def test_reference_deficit_point_inside_interval(result):
    lo, hi = result["kappa_interval"]["interval"]
    assert lo < result["kappa_interval"]["reference_deficit_point"]["kappa"] < hi


def test_guards_fail_closed(result):
    guards = result["guards"]
    assert guards["source_only"] is False
    assert guards["new_axiom_introduced"] is False
    assert guards["promotable_as_oph_source_theorem"] is False
    assert guards["measured_alpha_in_solve_path"] is True
    assert guards["measured_lepton_masses_directly_supplied_to_inversion"] is False
    assert guards["target_anchored_lepton_ratios_in_solve_path"] is True
    assert guards["measured_lepton_triple_used_to_calibrate_higher_order_remainder"] is True
    assert guards["charged_mass_information_in_solve_path"] is True
    assert guards["blind_normalization_prediction"] is False
    assert guards["usable_for_public_final_values"] is False
    assert guards["usable_as_diagnostic_route_finder"] is True
    assert guards["satisfies_production_constructive_next_artifact"] is False


def test_ratios_come_from_centered_family_functional(result):
    r1, r2 = result["inputs"]["ratio_provenance"]["ratios_mu_over_e_tau_over_e"]
    assert r1 == pytest.approx(206.76828270846718, rel=1e-12)
    assert r2 == pytest.approx(3477.3652619063414, rel=1e-9)


def test_stage5_scale_consistency_is_small(result):
    assert abs(result["stage5_scale_consistency"]["implied_kappa_offset"]) < 0.05


def test_source_only_no_go_unchanged(result):
    lemma = result["kappa_symmetry_breaking_lemma"]
    assert lemma["source_only_no_go_status"] == "unchanged"
    assert lemma["packet_derivative_in_kappa"] == pytest.approx(-2.0 / math.pi)


def test_blind_route_starts_with_source_only_ratios(result):
    route = result["blind_prediction_route"]
    assert route["stage_1_shape"]["current_status"] == "open_target_anchored_shape_only"
    assert route["stage_2_normalization"]["current_status"] == "open_empirical_transport_only"
    assert route["stage_3_freeze_then_compare"]["current_status"] == "not_ready"


def test_hadronic_contract_parent_cited(result):
    parent = result["inputs"]["hadronic_contract_parent"]
    assert parent["issue"] == 317
    assert parent["artifact"] == "oph_ward_projected_spectral_measure_proof_packet"
    assert parent["contract_certified"] is True
    assert parent["physical_source_payload_available"] is False
    assert parent["verifier_command"]
    assert "hadronic_contract_parent" in (
        result["interval_width_attribution_kappa_units"]["reduction"]
    )


def _doctored_packet(tmp_path, mutate):
    packet = json.loads(
        lane.HADRONIC_PROOF_PACKET_JSON.read_text(encoding="utf-8")
    )
    mutate(packet)
    doctored = tmp_path / "doctored_packet.json"
    doctored.write_text(json.dumps(packet), encoding="utf-8")
    return doctored


def test_fail_closed_on_missing_packet(tmp_path):
    with pytest.raises(SystemExit, match="proof packet missing"):
        lane.build(
            tmp_path / "out.json",
            hadronic_packet_path=tmp_path / "absent_packet.json",
        )


def test_fail_closed_on_unaccepted_packet(tmp_path):
    doctored = _doctored_packet(
        tmp_path, lambda p: p.update({"accepted": False})
    )
    with pytest.raises(SystemExit, match="not accepted"):
        lane.build(tmp_path / "out.json", hadronic_packet_path=doctored)


def test_fail_closed_on_uncertified_contract(tmp_path):
    def mutate(p):
        p["verdicts"]["contract_certified"]["value"] = False

    doctored = _doctored_packet(tmp_path, mutate)
    with pytest.raises(SystemExit, match="contract_certified"):
        lane.build(tmp_path / "out.json", hadronic_packet_path=doctored)


def test_fail_closed_once_source_payload_reported_available(tmp_path):
    def mutate(p):
        p["verdicts"]["physical_source_payload_available"]["value"] = True

    doctored = _doctored_packet(tmp_path, mutate)
    with pytest.raises(SystemExit, match="re-derive this lane"):
        lane.build(tmp_path / "out.json", hadronic_packet_path=doctored)
