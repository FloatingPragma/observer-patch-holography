"""Tests for the anchor-reconciliation reduction theorems."""

from __future__ import annotations

import functools

import derive_d11_anchor_reconciliation_reduction_theorems as lane


@functools.lru_cache(maxsize=1)
def _artifact() -> dict:
    return lane.build()


def test_chart_affinity_is_exact():
    cert = _artifact()["theorem_B_quadratic_cost"]["chart_affinity"]
    assert cert["exact"] is True
    assert cert["worst_residual"] < 1.0e-9


def test_gaussian_quadratic_and_barycenter_are_exact():
    artifact = _artifact()
    assert artifact["theorem_B_quadratic_cost"]["gaussian_kl_quadratic"]["exact"] is True
    assert artifact["theorem_A_barycenter_placement"]["exact"] is True


def test_premise_ledger_reflects_the_reduction():
    ledger = _artifact()["premise_ledger"]
    assert ledger["AR2_quadratic_cost"]["status"] == "discharged_under_RM"
    assert ledger["AR1_reconciliation"]["status"] == "mechanism_discharged_parentage_open"
    assert ledger["AR3_equal_capacity"]["status"] == "reduced_to_CF2"


def test_carrier_facts_stay_open_and_fail_closed():
    artifact = _artifact()
    facts = artifact["carrier_facts"]
    assert facts["CF1_two_parent_ports"]["status"] == "open"
    assert facts["CF2_same_register_class_equal_depth"]["status"] == "open"
    assert artifact["promotion_allowed"] is False
    assert artifact["status"] == "PREMISES_REDUCED_TO_CARRIER_FACTS"


def test_consequence_matches_the_scan():
    consequence = _artifact()["two_loop_consequence"]
    assert abs(consequence["mt_pole_gev"] - 172.63) < 0.05
    assert abs(consequence["mh_tree_gev"] - 125.77) < 0.05
