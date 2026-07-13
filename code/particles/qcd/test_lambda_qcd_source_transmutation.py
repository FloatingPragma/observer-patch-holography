"""Tests for the Lambda_QCD source transmutation lane."""

from __future__ import annotations

import functools

import derive_lambda_qcd_source_transmutation as lane


@functools.lru_cache(maxsize=1)
def _artifact() -> dict:
    return lane.build()


def test_machinery_is_certified_against_published_world_values():
    validation = _artifact()["machinery_validation"]
    assert all(validation["checks"].values())
    assert validation["role"] == "code_certification_only_never_a_source_row"


def test_source_lambda3_lands_in_the_expected_window():
    artifact = _artifact()
    lam3 = artifact["central"]["lambda3_gev"]
    assert 0.30 < lam3 < 0.37
    lo, hi = artifact["lambda3_interval_gev"]
    assert lo < lam3 < hi


def test_threshold_brackets_are_swept_and_fail_closed():
    artifact = _artifact()
    assert artifact["checks"]["threshold_brackets_swept"] is True
    assert artifact["promotion_allowed"] is False
    assert artifact["row_class"] == "conditional_on_P_and_declared_threshold_inputs"


def test_over_e_star_ratio_is_consistent_with_display():
    artifact = _artifact()
    ratio = artifact["over_E_star"]["lambda3"]
    display = artifact["central"]["lambda3_gev"]
    assert abs(ratio * lane.E_STAR_DISPLAY_GEV / display - 1.0) < 1.0e-12
