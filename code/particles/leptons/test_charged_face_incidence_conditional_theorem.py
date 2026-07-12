"""Tests for the declared-kernel conditional face-incidence theorem."""

from __future__ import annotations

import json

import mpmath as mp

import derive_charged_face_incidence_conditional_theorem as lane


def artifact() -> dict:
    raw = lane.FACE_RECEIPT.read_bytes()
    return lane.build_artifact(json.loads(raw), lane.sha256(raw))


def test_declared_map_is_strictly_contractive_with_one_fixed_point():
    result = artifact()
    contraction = mp.mpf(result["repair_map"]["contraction_constant"])
    residual = [
        abs(mp.mpf(value))
        for value in result["repair_map"]["fixed_point_residual"]
    ]
    assert contraction < 1
    assert max(residual) < mp.mpf("1e-90")
    assert result["checks"]["strict_contraction"] is True


def test_conditional_mass_readout_reproduces_frozen_candidate_without_promotion():
    result = artifact()
    masses = [mp.mpf(value) for value in result["conditional_spectrum"]["masses_mev"]]
    expected = [
        mp.mpf("0.5109989508433937143"),
        mp.mpf("105.6583755014812777"),
        mp.mpf("1776.930000014351498"),
    ]
    assert max(abs(left - right) for left, right in zip(masses, expected, strict=True)) < mp.mpf("1e-15")
    assert result["runtime_charged_reference_consumed"] is False
    assert result["historical_charged_target_informed"] is True
    assert result["branch_tuple_coherent"] is False
    assert result["public_prediction_promotion_allowed"] is False


def test_explicit_multiplier_witnesses_prove_scoped_non_rigidity():
    result = artifact()
    witnesses = result["scoped_source_normalization_non_rigidity"]["witnesses"]
    assert [row["source_multipliers"] for row in witnesses] == [
        [2, 1, 1],
        [1, 2, 1],
        [1, 1, 2],
    ]
    assert result["checks"]["explicit_source_multiplier_witnesses_change_masses"] is True
    assert "not a non-entailment theorem for every OPH axiom" in (
        result["scoped_source_normalization_non_rigidity"]["statement"]
    )


def test_bare_to_endpoint_stage_gap_is_explicit():
    result = artifact()
    boundary = result["bare_to_endpoint_stage_boundary"]
    assert mp.mpf(boundary["endpoint_over_bare"]) < 1
    assert result["checks"]["bare_and_endpoint_amplitude_ratios_are_distinct"] is True
    assert result["checks_pass"] is True
