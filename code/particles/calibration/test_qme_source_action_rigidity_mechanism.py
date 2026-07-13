"""Tests for the quotient-MaxEnt source-action rigidity mechanism receipt."""

from __future__ import annotations

import functools

import mpmath as mp
import pytest

import derive_qme_source_action_rigidity_mechanism as lane


@pytest.fixture(autouse=True)
def _scoped_precision():
    """Set a high mpmath precision for each test and restore it on teardown.

    The tight-tolerance comparisons in this module require extended precision.
    Scoping keeps the global mpmath precision unchanged for other test modules.
    """

    with mp.workdps(120):
        yield


@functools.lru_cache(maxsize=1)
def _artifact() -> dict:
    return lane.build_artifact(lane.load_certificate())


def test_bernoulli_moment_matches_the_expected_value():
    certificate = lane.load_certificate()
    p = lane.bernoulli_moment(certificate["P_cand"])
    assert abs(p - mp.mpf("0.007299545168876648")) < mp.mpf("1e-15")
    branch = _artifact()["branch"]
    assert branch["p_agrees_to_1e-15"] is True
    assert mp.mpf(branch["p_agreement_residual"]) < mp.mpf("1e-15")


def test_maxent_laws_sum_to_one_and_meet_the_moments():
    checks = _artifact()["mechanism_checks"]
    assert checks["laws_normalized"] is True
    assert mp.mpf(checks["mu0_sum_residual"]) < mp.mpf("1e-40")
    assert mp.mpf(checks["mu1_sum_residual"]) < mp.mpf("1e-40")
    assert checks["moment_constraints_hold"] is True
    assert checks["gibbs_reconstruction_matches"] is True


def test_projective_compatibility_residual_is_zero():
    checks = _artifact()["mechanism_checks"]
    assert checks["projective_compatibility_exact"] is True
    assert mp.mpf(checks["projective_residual"]) < mp.mpf("1e-30")


def test_pythagorean_identity_residual_is_small():
    checks = _artifact()["mechanism_checks"]
    assert checks["pythagorean_identity_holds"] is True
    assert mp.mpf(checks["pythagorean_residual"]) < mp.mpf("1e-14")


def test_selector_gap_meets_the_pinsker_lower_bound():
    checks = _artifact()["mechanism_checks"]
    assert checks["selector_gap_exceeds_pinsker_bound"] is True
    gap = mp.mpf(checks["selector_gap"])
    bound = mp.mpf(checks["pinsker_lower_bound"])
    assert bound > 0
    assert gap >= bound


def test_legendre_hessian_matches_the_declared_value():
    checks = _artifact()["mechanism_checks"]
    hessian = mp.mpf(checks["Gamma_hessian"])
    assert abs(hessian - mp.mpf("138.00219014027388")) < mp.mpf("1e-9")


def test_moment_gate_is_absent_and_promotion_is_disallowed():
    artifact = _artifact()
    gate = artifact["physical_standard_model_moment_packet"]
    assert gate["status"] == "ABSENT"
    assert gate["gate"] == "fail_closed"
    assert gate["open_proofs"] == [
        "operator completeness modulo BRST",
        "c_r = c_r(P_star, N_star) source emission",
        "reflection positivity and continuum nontriviality",
        "source-complete matter spectral cuts",
        "physical 1PI two-point blocks with analytic continuation",
        "source clock closure",
    ]
    assert artifact["promotion_allowed"] is False
    assert artifact["status"] == "SELECTION_MECHANISM_CLOSED_PHYSICAL_MOMENT_VECTOR_OPEN"
    assert artifact["checks_pass"] is True
    registry = artifact["theorem_registry"]
    assert "SOURCE_ACTION_RIGIDITY_THEOREM" in registry
    assert "PINSKER_SELECTOR_GAP" in registry
    assert "REFINEMENT_PUSHFORWARD_THEOREM" in registry
    assert "BRST_CLASS_UNIQUENESS" in registry
