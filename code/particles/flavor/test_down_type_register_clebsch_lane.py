"""Tests for the down-type register Clebsch lane."""

from __future__ import annotations

import functools
import copy
import hashlib

import pytest

import derive_down_type_register_clebsch_lane as lane


@functools.lru_cache(maxsize=1)
def _artifact() -> dict:
    return lane.build()


def test_clebsch_boundary_is_the_declared_pattern():
    artifact = _artifact()
    factors = artifact["clebsch_boundary"]["factors"]
    assert factors == {"b_over_tau": 1.0, "s_over_mu": 1.0 / 3.0, "d_over_e": 3.0}
    assert artifact["promotion_allowed"] is False


def test_retrospective_comparison_rejects_all_six_orders():
    artifact = _artifact()
    scan = artifact["permutation_scan"]
    assert scan["permutation_count"] == 6
    assert scan["exhaustive"] is True
    assert scan["all_permutations_rejected_by_conservative_flag_gate"] is True
    assert scan["retrospective_metric"]["target_informed"] is True
    assert scan["retrospective_metric"]["preregistered"] is False
    assert scan["retrospective_metric"]["physical_order_selected"] is False
    assert (
        scan["retrospective_metric"][
            "current_assumed_order_uniquely_least_discrepant"
        ]
        is True
    )
    assert (
        artifact["status"]
        == "CONDITIONAL_DECLARED_ROUTE_RETROSPECTIVELY_REJECTED"
    )


def test_normalization_tension_is_carried_openly():
    artifact = _artifact()
    compare = artifact["compare_only"]
    assert 0.25 < compare["mb_relative"] < 0.60
    tension = artifact["normalization_tension"]
    assert "THIRD_GENERATION_REGISTER_FACTOR" in tension["open_objects"]
    assert "proximity criterion" in tension["statement"]


def test_predictions_are_positive_and_ordered():
    artifact = _artifact()
    predictions = artifact["predictions"]
    assert (
        predictions["mb_mb_gev"]
        > predictions["ms_2gev_gev"]
        > predictions["md_2gev_gev"]
        > 0
    )
    assert artifact["arithmetic_checks_pass"] is True
    assert artifact["physical_premises_discharged"] is False
    assert not any(artifact["physical_premises"].values())


def test_dependency_artifacts_are_hash_bound():
    artifact = _artifact()
    hashes = artifact["dependency_audit"]["input_sha256"]
    assert hashes == {
        "mcpr": hashlib.sha256(lane.MCPR.read_bytes()).hexdigest(),
        "clebsch_selection": hashlib.sha256(
            lane.SELECTION.read_bytes()
        ).hexdigest(),
        "flag_2024_fixture": hashlib.sha256(
            lane.FLAG_FIXTURE.read_bytes()
        ).hexdigest(),
        "declared_boundary_calibration_producer": hashlib.sha256(
            lane.CALIBRATION_PRODUCER.read_bytes()
        ).hexdigest(),
    }


def test_ratio_only_core_does_not_read_optional_scale_display():
    mcpr = lane._load_json(lane.MCPR)
    selection = lane._load_json(lane.SELECTION)
    flag = lane._load_json(lane.FLAG_FIXTURE)
    baseline = lane.build_artifact(mcpr, selection, flag)
    mutant = copy.deepcopy(mcpr)
    mutant["optional_scale_display"]["masses_MeV"] = [
        str(float(value) * 2.0)
        for value in mutant["optional_scale_display"]["masses_MeV"]
    ]
    changed = lane.build_artifact(mutant, selection, flag)
    assert baseline["ratio_only_core"] == changed["ratio_only_core"]
    assert (
        baseline["predictions"]["mb_mb_gev"]
        != changed["predictions"]["mb_mb_gev"]
    )


def test_selection_and_flag_inputs_fail_closed_when_doctored():
    mcpr = lane._load_json(lane.MCPR)
    selection = lane._load_json(lane.SELECTION)
    flag = lane._load_json(lane.FLAG_FIXTURE)

    doctored_selection = copy.deepcopy(selection)
    doctored_selection["pairing_theorem"][
        "independent_yukawa_coefficients_equated"
    ] = True
    with pytest.raises(ValueError, match="conditional channel compatibility"):
        lane.build_artifact(mcpr, doctored_selection, flag)

    doctored_flag = copy.deepcopy(flag)
    doctored_flag["claim_boundary"]["significance_gate_preregistered"] = True
    with pytest.raises(ValueError, match="compare-only claim boundary"):
        lane.build_artifact(mcpr, selection, doctored_flag)


def test_common_transport_identity_is_scale_qualified():
    transport = _artifact()["common_transport"]
    assert transport["generation_dependent_light_threshold_present"] is False
    assert transport["relative_identity_gap"] < 1e-13
    assert "five-Yukawa" in transport["approximation"]


def test_quantization_audit_excludes_impossible_ratios_and_ckm_relabel():
    audit = _artifact()["quantization_audit"]
    assert audit["permitted_distinct_assignment_factor_ratios"] == [
        "1/9",
        "1/3",
        "3",
        "9",
    ]
    assert "27" in audit["excluded_values_from_prior_informal_note"]
    assert audit["nearest_permitted_factor_ratio"] == "1/9"
    assert 1.15 < audit["predicted_over_flag_nf_2+1+1"] < 1.16
    assert audit["diagonal_ansatz_ckm_consequence"]["ckm_matrix"] == "identity"
    assert (
        audit["diagonal_ansatz_ckm_consequence"][
            "nonzero_cabibbo_angle_derived"
        ]
        is False
    )
