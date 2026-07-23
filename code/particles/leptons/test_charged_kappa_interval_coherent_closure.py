"""Tests for the payload-coherent charged kappa-interval closure."""

from __future__ import annotations

import json
import math

import pytest

import derive_charged_kappa_interval_coherent_closure as lane
import derive_charged_kappa_interval_from_alpha_transport as rectangle_lane


@pytest.fixture(scope="module")
def result(tmp_path_factory):
    out = tmp_path_factory.mktemp("coherent") / "coherent_closure.json"
    return lane.build(out)


def test_identity_gate_passes_on_recorded_endpoint(result):
    gate = result["coherence_premise"]["identity_gate"]
    assert gate["passed"] is True
    assert gate["breaking_term"] is None


def test_coherent_interval_ordered_and_contains_witness(result):
    lo, hi = result["kappa_interval"]["interval"]
    assert math.isfinite(lo) and math.isfinite(hi)
    assert lo < result["kappa_interval"]["central"] < hi
    assert result["compare_only"]["witness_inside_certified_intervals"] is True
    for row, ref in zip(
        result["conditional_mass_rows"],
        result["compare_only"]["witness_masses_gev"],
        strict=True,
    ):
        row_lo, row_hi = row["mass_interval"]
        assert row_lo < ref < row_hi


def test_width_reduction_is_material(result):
    # the factor is relative to the rectangle lane at the current payload
    # uncertainty; with the published-compilation payload the rectangle is
    # itself tight, so the coherent gain is bounded but stays material
    assert result["kappa_interval"]["width_reduction_factor"] > 2.0
    rect_lo, rect_hi = result["kappa_interval"]["rectangle_interval"]
    lo, hi = result["kappa_interval"]["interval"]
    assert rect_lo < lo < hi < rect_hi


def test_coherent_central_equals_rectangle_central(result):
    match = result["kappa_interval"]["central_match"]
    assert match["defect"] < 1.0e-9
    assert result["kappa_interval"]["central"] == pytest.approx(
        match["rectangle_central"], abs=1.0e-9
    )


def test_payload_terms_cancel_in_attribution(result):
    attribution = result["interval_width_attribution_kappa_units"]
    assert attribution["hadronic_payload_uncertainty"] == 0.0
    assert attribution["anchor_gap_half_width"] == 0.0
    assert attribution["higher_order_lepton_budget"] > 0.0


def test_guards_fail_closed(result):
    guards = result["guards"]
    assert guards["source_only"] is False
    assert guards["new_axiom_introduced"] is False
    assert guards["promotable_as_oph_source_theorem"] is False
    assert guards["payload_double_count_removed"] is True
    assert guards["conditional_on_payload_coherent_anchor_gap"] is True
    assert guards["usable_for_public_final_values"] is False
    assert guards["satisfies_production_constructive_next_artifact"] is False


def test_parents_cited(result):
    inputs = result["inputs"]
    assert inputs["hadronic_contract_parent"]["issue"] == 317
    assert inputs["hadronic_contract_parent"]["contract_certified"] is True
    assert inputs["rectangle_lane_ref"].endswith(
        "charged_kappa_interval_from_alpha_transport.json"
    )


def test_row_class_is_coherent_closure(result):
    assert (
        result["row_class"]
        == "target_shape_plus_empirical_transport_coherent_closure"
    )
    for row in result["conditional_mass_rows"]:
        assert row["status"] == "certified_empirical_closure_interval_coherent"


def test_identity_gate_trips_on_doctored_endpoint(tmp_path):
    endpoint = json.loads(lane.ENDPOINT_JSON.read_text(encoding="utf-8"))
    recorded = endpoint["compare_only"]["same_scheme_anchor_gap_interval_inv_alpha"]
    doctored_lo = str(float(recorded[0]) + 1.0e-4)
    endpoint["compare_only"]["same_scheme_anchor_gap_interval_inv_alpha"] = [
        doctored_lo,
        recorded[1],
    ]
    doctored = tmp_path / "doctored_endpoint.json"
    doctored.write_text(json.dumps(endpoint), encoding="utf-8")
    out = tmp_path / "out.json"
    with pytest.raises(SystemExit, match="coherence identity gate failed"):
        lane.build(out, endpoint_path=doctored)
    verdict = json.loads(out.read_text(encoding="utf-8"))
    assert verdict["proof_status"] == "no_improvement_identity_gate_failed"
    assert (
        verdict["no_improvement_verdict"]["breaking_term"]
        == "anchor_gap_lower_endpoint"
    )


def test_fail_closed_on_missing_rectangle_artifact(tmp_path):
    with pytest.raises(SystemExit, match="rectangle-lane artifact missing"):
        lane.build(
            tmp_path / "out.json",
            rectangle_path=tmp_path / "absent_rectangle.json",
        )


def test_fail_closed_on_doctored_hadronic_packet(tmp_path):
    packet = json.loads(
        rectangle_lane.HADRONIC_PROOF_PACKET_JSON.read_text(encoding="utf-8")
    )
    packet["accepted"] = False
    doctored = tmp_path / "doctored_packet.json"
    doctored.write_text(json.dumps(packet), encoding="utf-8")
    with pytest.raises(SystemExit, match="not accepted"):
        lane.build(tmp_path / "out.json", hadronic_packet_path=doctored)
