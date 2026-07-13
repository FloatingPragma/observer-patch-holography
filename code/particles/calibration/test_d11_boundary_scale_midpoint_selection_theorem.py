"""Tests for the conditional midpoint selection theorem."""

from __future__ import annotations

import functools
from fractions import Fraction

import derive_d11_boundary_scale_midpoint_selection_theorem as lane


@functools.lru_cache(maxsize=1)
def _artifact() -> dict:
    return lane.build()


def test_weighted_minimizer_is_exact():
    t = lane.weighted_minimizer_exact(
        Fraction(1), Fraction(9), Fraction(1), Fraction(3)
    )
    assert t == Fraction(28, 4)
    equal = lane.weighted_minimizer_exact(
        Fraction(2), Fraction(10), Fraction(7), Fraction(7)
    )
    assert equal == Fraction(6)


def test_implication_certificates_all_pass():
    certs = _artifact()["implication_certificates"]
    assert certs["implication_holds"] is True
    assert certs["gradient_zero_at_minimizer"] == certs["sampled_rational_points"]


def test_closed_form_identity_is_exact():
    identity = _artifact()["closed_form_identity"]
    assert identity["exact"] is True
    assert "exp(-pi)" in identity["identity"]


def test_premises_stay_open_and_promotion_blocked():
    artifact = _artifact()
    assert artifact["status"] == "CONDITIONAL_SELECTION_THEOREM_PREMISES_OPEN"
    assert artifact["promotion_allowed"] is False
    assert all(
        gate["status"] == "open" for gate in artifact["premises"].values()
    )


def test_two_loop_consequence_matches_the_scan():
    consequence = _artifact()["two_loop_consequence"]
    assert abs(consequence["mt_pole_gev"] - 172.63) < 0.05
    assert abs(consequence["mh_tree_gev"] - 125.77) < 0.05
