"""Tests for the D11 double-criticality boundary-scale scan."""

from __future__ import annotations

import functools

import derive_d11_criticality_boundary_scan as lane


@functools.lru_cache(maxsize=1)
def _artifact() -> dict:
    return lane.build()


def test_archived_endpoints_are_reproduced_to_1e_6():
    validation = _artifact()["archived_endpoint_validation"]
    assert validation["passes_1e_6"] is True
    assert validation["max_residual"] < 1.0e-6


def test_named_source_scales_are_ordered_and_present():
    scales = _artifact()["source_scales_gev"]
    assert (
        scales["mu_U_gauge_unification"]
        < scales["E_cell_pixel_energy"]
        < scales["E_star"]
    )
    assert abs(scales["v_transmutation_gev"] - 246.6) < 0.5


def test_boundary_scale_moves_the_pair_toward_the_upper_window():
    named = _artifact()["one_loop_named_boundaries"]
    mu_u = named["mu_U_gauge_unification"]
    e_cell = named["E_cell_pixel_energy"]
    assert abs(mu_u["mt_pole_gev"] - 164.13) < 0.02
    assert abs(mu_u["mh_tree_gev"] - 115.10) < 0.02
    assert 170.0 < e_cell["mt_pole_gev"] < 171.0
    assert 127.0 < e_cell["mh_tree_gev"] < 128.5


def test_curve_is_monotone_and_fail_closed():
    artifact = _artifact()
    assert artifact["checks"]["curve_mh_monotone_in_boundary_scale"] is True
    assert artifact["promotion_allowed"] is False
    assert artifact["boundary_law"]["free_continuous_parameters"] == 0
    assert (
        artifact["single_anchor_coherence"]["single_anchor_branch"]
        == "E_cell_pixel_energy"
    )


def test_double_criticality_yukawa_matches_beta_lambda_zero():
    g_y, g2 = 0.36, 0.52
    y = lane.critical_surface_yukawa(g_y, g2)
    residual = lane.beta_lambda_1loop(0.0, y, g_y, g2, 1.0)
    assert abs(residual) < 1.0e-18
