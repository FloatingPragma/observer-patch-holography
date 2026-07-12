#!/usr/bin/env python3
"""Tests for the color amplitude/loop split derivation (#521 A2/A3)."""

from __future__ import annotations

import math
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from derive_color_amplitude_loop_split import build  # noqa: E402


def test_guards_no_source_theorem_claimed() -> None:
    report = build()
    assert report["guards"]["measured_values_in_any_oph_solve_path"] is False
    assert report["guards"]["source_only_theorem_emitted"] is False
    assert report["guards"]["public_promotion_allowed"] is False
    assert report["guards"]["complete_archived_value_law_tested_as_constant_c_competitor"] is False


def test_mW_selects_sqrt_Nc_over_2() -> None:
    report = build()
    ds = report["data_selection"]
    # the clean W-mass optimum sits on sqrt(3)/2 to better than 0.5 percent
    assert abs(ds["MW_best_c"] - math.sqrt(3) / 2) < 0.005


def test_mZ_selects_Nc_over_2() -> None:
    report = build()
    ds = report["data_selection"]
    assert abs(ds["MZ_best_d"] - 1.5) < 0.06


def test_competitors_excluded_prediction_inside() -> None:
    report = build()
    comps = report["data_selection"]["competitors"]
    assert comps["color_amplitude_loop  c=sqrt(Nc)/2"]["inside_MWMZ_1sigma"] is True
    assert comps["loop_symmetric        c=Nc/2"]["inside_MWMZ_1sigma"] is False
    assert comps["loop_symmetric        c=Nc"]["inside_MWMZ_1sigma"] is False
    assert comps["leading_constant_proxy c=1/(4 beta)"]["inside_MWMZ_1sigma"] is False
    assert "not excluded" in report["verdict"]["complete_archived_value_law"]


def test_channel_identification_remains_open() -> None:
    report = build()
    assert "channel_identification_open" in report["verdict"]["status"]
