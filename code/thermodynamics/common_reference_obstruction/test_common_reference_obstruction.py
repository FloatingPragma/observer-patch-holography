"""Mutation tests for the exact B12 common-reference obstruction."""

from __future__ import annotations

from fractions import Fraction

import verify_common_reference_obstruction as verifier


def test_canonical_packet_passes_every_check():
    result = verifier.verify_packet()
    assert result["ok"]
    assert all(result["checks"].values())
    assert result["empirical_bracket"] == ["1905/16384", "953/8192"]


def test_count_mutation_breaks_earned_stationarity_and_eigenpair():
    result = verifier.verify_packet(counts=((1325, 106), (5, 503)))
    assert not result["ok"]
    assert not result["checks"]["stationary_exact"]
    assert not result["checks"]["eigenpair_exact"]


def test_mode_mutation_breaks_exact_eigenpair():
    result = verifier.verify_packet(mode=(Fraction(54357), Fraction(-7155)))
    assert not result["ok"]
    assert not result["checks"]["eigenpair_exact"]


def test_zero_eigenvalue_removes_intertwiner_obstruction():
    result = verifier.verify_packet(eigenvalue=Fraction(0))
    assert not result["ok"]
    assert not result["checks"]["genuinely_mixing_mode"]
    assert not result["checks"]["intertwiner_obstruction"]


def test_unit_eigenvalue_removes_intertwiner_obstruction():
    result = verifier.verify_packet(eigenvalue=Fraction(1))
    assert not result["ok"]
    assert not result["checks"]["genuinely_mixing_mode"]
    assert not result["checks"]["intertwiner_obstruction"]


def test_matching_sample_denominator_removes_empirical_no_go():
    result = verifier.verify_packet(sample_total=61511)
    assert not result["ok"]
    assert not result["checks"]["adjacent_grid_bracket"]
    assert not result["checks"]["coprime_denominators"]
    assert not result["checks"]["deterministic_pushforward_obstruction"]


def test_stationary_mutation_is_detected_independently():
    result = verifier.verify_packet(
        stationary=(Fraction(7156, 61511), Fraction(54355, 61511))
    )
    assert not result["ok"]
    assert result["checks"]["stationary_normalized"]
    assert not result["checks"]["stationary_exact"]
