"""Tests for the model-level D11 carrier census."""

from __future__ import annotations

import functools

import derive_d11_carrier_census as lane


@functools.lru_cache(maxsize=1)
def _artifact() -> dict:
    return lane.build()


def test_cf1_census_finds_exactly_two_anchor_ports():
    cf1 = _artifact()["cf1_dependency_census"]
    assert cf1["verified"] is True
    assert cf1["port_count_is_two"] is True
    assert sorted(cf1["anchor_rooted_ports"]) == [
        "transmutation_register",
        "unification_register",
    ]
    assert cf1["non_anchor_ports"] == []


def test_cf2_census_certifies_class_depth_and_parents():
    cf2 = _artifact()["cf2_anchor_class_census"]
    assert cf2["verified"] is True
    assert cf2["equal_depth"] is True
    assert cf2["same_parent_set"] is True
    assert cf2["midpoint_exponents"][
        "matches_E_star_exp_minus_pi_P_minus_sixth"
    ] is True


def test_faithfulness_gate_stays_open():
    artifact = _artifact()
    assert artifact["carrier_model_faithfulness"]["status"] == "open"
    assert artifact["promotion_allowed"] is False
    assert artifact["status"] == "CENSUS_VERIFIED_AT_MODEL_LEVEL_FAITHFULNESS_OPEN"


def test_higgs_derivation_trace_emits_the_conditional_value():
    trace = _artifact()["higgs_derivation"]
    assert abs(trace["m_H_GeV"] - 125.77) < 0.05
    assert abs(trace["m_t_GeV"] - 172.63) < 0.05
    assert trace["bands_GeV"]["loop_truncation_half_shift"] > 1.0
    assert "CARRIER_MODEL_FAITHFULNESS" in trace["conditions"]


def test_anchor_depth_is_uniform():
    depths = _artifact()["cf2_anchor_class_census"]["anchor_depths"]
    assert len(set(depths.values())) == 1
