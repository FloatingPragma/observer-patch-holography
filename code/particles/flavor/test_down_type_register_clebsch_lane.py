"""Tests for the down-type register Clebsch lane."""

from __future__ import annotations

import functools

import derive_down_type_register_clebsch_lane as lane


@functools.lru_cache(maxsize=1)
def _artifact() -> dict:
    return lane.build()


def test_clebsch_boundary_is_the_declared_pattern():
    artifact = _artifact()
    factors = artifact["clebsch_boundary"]["factors"]
    assert factors == {"b_over_tau": 1.0, "s_over_mu": 1.0 / 3.0, "d_over_e": 3.0}
    assert artifact["promotion_allowed"] is False


def test_ratio_outputs_land_at_the_ten_percent_scale():
    compare = _artifact()["compare_only"]
    assert abs(compare["cabibbo_relative"]) < 0.12
    assert abs(compare["ms_over_md_relative"]) < 0.25


def test_normalization_tension_is_carried_openly():
    artifact = _artifact()
    compare = artifact["compare_only"]
    assert 0.25 < compare["mb_relative"] < 0.60
    tension = artifact["normalization_tension"]
    assert "THIRD_GENERATION_REGISTER_FACTOR" in tension["open_objects"]


def test_predictions_are_positive_and_ordered():
    predictions = _artifact()["predictions"]
    assert (
        predictions["mb_mb_gev"]
        > predictions["ms_2gev_gev"]
        > predictions["md_2gev_gev"]
        > 0
    )
    assert 0.18 < predictions["cabibbo_gst_sqrt_md_over_ms"] < 0.24
