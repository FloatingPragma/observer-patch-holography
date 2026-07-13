"""Tests for the conditional nucleon-mass lane."""

from __future__ import annotations

import functools

import derive_nucleon_mass_external_qcd_ratio as lane


@functools.lru_cache(maxsize=1)
def _artifact() -> dict:
    return lane.build()


def test_row_class_declares_the_external_qcd_condition():
    artifact = _artifact()
    assert artifact["row_class"] == "oph_plus_external_qcd_theory"
    assert artifact["promotion_allowed"] is False
    assert "lattice" in artifact["external_theory_factor"]["source"]


def test_nucleon_central_lands_in_the_expected_window():
    prediction = _artifact()["prediction"]
    assert 0.85 < prediction["m_nucleon_gev_display"] < 1.00
    lo, hi = prediction["m_nucleon_interval_gev_display"]
    assert lo < prediction["m_nucleon_gev_display"] < hi


def test_compare_block_is_labeled_and_interval_contains_measurement():
    compare = _artifact()["compare_only"]
    assert compare["role"] == "comparison only, outside the prediction ancestry"
    assert compare["interval_contains_measured"] is True
    assert abs(compare["central_relative_difference"]) < 0.05


def test_checks_pass():
    assert _artifact()["checks_pass"] is True
