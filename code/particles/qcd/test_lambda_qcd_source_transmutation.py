"""Tests for the Lambda_QCD source transmutation lane."""

from __future__ import annotations

import functools
import math

import pytest

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
    assert artifact['range_scope']['full_threshold_rectangle_enclosed'] is False
    assert artifact['range_scope']['legacy_interval_fields']=='sampled_threshold_corner_ranges'


def test_over_e_star_ratio_is_consistent_with_display():
    artifact = _artifact()
    ratio = artifact["over_E_star"]["lambda3"]
    display = artifact["central"]["lambda3_gev"]
    assert abs(ratio * lane.E_STAR_DISPLAY_GEV / display - 1.0) < 1.0e-12


@pytest.mark.parametrize('n_light', [3, 4, 5])
def test_down_matching_has_the_inverse_direction_of_pdg_equation(n_light):
    # PDG 9.4 gives alpha_high = alpha_low*(1-11*a_low^2/72-...).
    # Composition with the down map cancels through absolute alpha^4.
    def residual(alpha):
        low = lane.decouple_down(alpha, n_light)
        a = low/math.pi
        c2 = -11/72
        c3 = -564731/124416+82043*lane.ZETA3/27648+2633*n_light/31104
        return abs(low*(1+c2*a*a+c3*a**3)-alpha)
    assert lane.decouple_down(0.2, n_light) > 0.2
    assert residual(0.1)/residual(0.05) > 28


@pytest.mark.parametrize('n_flavors', [3, 4, 5, 6])
def test_lambda_is_invariant_under_the_same_four_loop_flow(n_flavors):
    alpha = 0.15
    evolved = lane.run_alpha_s(alpha, 20, 30, n_flavors)
    first = lane.lambda_msbar(alpha, 20, n_flavors)
    second = lane.lambda_msbar(evolved, 30, n_flavors)
    assert abs(second/first-1) < 2e-12
