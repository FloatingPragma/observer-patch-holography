"""Tests for the W/Z/H residual-elimination boundary lane."""

from __future__ import annotations

import mpmath as mp

from fractions import Fraction

import derive_wzh_residual_elimination_boundary as lane


import pytest


@pytest.fixture(autouse=True)
def _scoped_precision():
    """Set a high mpmath precision for each test and restore it on teardown.

    The tight-tolerance comparisons in this module require extended precision.
    Scoping keeps the global mpmath precision unchanged for other test modules.
    """

    with mp.workdps(120):
        yield


def test_j10_factorization_is_an_exact_polynomial_identity():
    checks = lane.j10_rigidity_checks(sample_count=25)
    assert checks["identity_holds"] is True
    assert checks["factorization_residual_exactly_zero"] == 25


def test_j10_vanishes_exactly_at_the_origin_only():
    eta = Fraction(3, 137)
    kappa = Fraction(13, 3)
    assert lane.j10_value(Fraction(0), Fraction(0), eta, kappa) == 0
    assert lane.j10_value(Fraction(1, 1000), Fraction(0), eta, kappa) > 0
    assert lane.j10_value(Fraction(0), Fraction(1, 1000), eta, kappa) > 0


def test_two_law_boundary_reproduces_research_coordinates():
    artifact = lane.build()
    assert artifact["checks_pass"] is True
    two_law = artifact["strict_branch_two_law_boundary"]

    zero = two_law["zero_selector_law"]
    nonzero = two_law["nonzero_carrier_law"]
    assert abs(zero["MW_over_E_star"] / 6.5796308665770565e-18 - 1) < 1e-12
    assert abs(zero["MZ_over_E_star"] / 7.46333486303666e-18 - 1) < 1e-12
    assert abs(nonzero["MW_over_E_star"] / 6.578870367942109e-18 - 1) < 1e-12
    assert abs(nonzero["MZ_over_E_star"] / 7.463750125274474e-18 - 1) < 1e-12

    assert zero["J10"] == 0.0
    assert nonzero["J10"] > 1e-7
    assert abs(nonzero["J10"] / 3.1206348013466805e-07 - 1) < 1e-3


def test_carrier_reproductions_and_exact_readout_pass():
    artifact = lane.build()
    assert artifact["prospective_c10_carrier"]["checks_pass"] is True
    assert artifact["prospective_c11_carrier"]["checks_pass"] is True
    assert artifact["d11_exact_readout"]["checks_pass"] is True
    assert (
        artifact["prospective_c10_carrier"]["carrier_status"]
        == "prospective_post_exposure_candidate_not_current_source_evidence"
    )


def test_literal_ivp_endpoints_are_algebraically_consistent():
    artifact = lane.build()
    ivp = artifact["literal_source_ivp"]
    assert ivp["checks_pass"] is True
    assert ivp["provenance"]["ivp_reproduced_in_repo"] is False
    display = ivp["display_GeV_using_unclosed_clock"]
    assert abs(display["mH_tree"] - 115.101) < 0.01
    assert abs(display["mt_MSbar"] - 154.779) < 0.01
    assert abs(display["mt_QCD_converted"] - 164.131) < 0.01


def test_lane_is_fail_closed():
    artifact = lane.build()
    assert artifact["promotion_allowed"] is False
    assert (
        artifact["remaining_discrete_selection"]["selection_principle_status"]
        == "SOURCE_LAW_SELECTION_PRINCIPLE_ABSENT"
    )
    width = artifact["remaining_discrete_selection"]["certified_ambiguity_width_GeV"]
    assert 0.001 < width["MW"] < 0.02
    assert 0.001 < width["MZ"] < 0.02
    assert artifact["d11_synchronization_provenance"][
        "declared_surface_target_ancestral"
    ] is True
