#!/usr/bin/env python3
"""Tests for the repair-tuple selection theorem and conditional predictions."""

from __future__ import annotations

import json
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
RUNS = HERE.parent / "runs" / "calibration"
sys.path.insert(0, str(HERE))

from derive_d10_repair_tuple_selection_theorem import (  # noqa: E402
    build,
    candidate_a_tuple,
    candidate_b_tuple,
    forward_ht,
    forward_wz,
    load_basis,
)
from emit_conditional_ew_predictions import build as build_predictions  # noqa: E402


def test_candidate_a_reproduces_value_law_artifact() -> None:
    basis = load_basis()
    with open(RUNS / "d10_ew_target_free_repair_value_law.json") as f:
        law = json.load(f)
    cand = candidate_a_tuple(basis)
    assert abs(cand["tau2_exact"] - law["repair_chart"]["tau2_tree_exact"]) < 1e-12
    assert abs(cand["delta_n_exact"] - law["repair_chart"]["delta_n_tree_exact"]) < 1e-12
    wz = forward_wz(basis, cand["tau2_exact"], cand["delta_n_exact"])
    quintet = law["coherent_emitted_quintet"]
    assert abs(wz["MW_pole_gev"] - quintet["MW_pole"]) < 1e-6
    assert abs(wz["MZ_pole_gev"] - quintet["MZ_pole"]) < 1e-6


def test_candidate_a_reproduces_d11_split_artifact() -> None:
    basis = load_basis()
    with open(RUNS / "d11_live_exact_split_pair_theorem.json") as f:
        split = json.load(f)
    with open(RUNS / "d11_declared_calibration_surface.json") as f:
        surface = json.load(f)
    cand = candidate_a_tuple(basis)
    ht = forward_ht(basis, cand["tau2_exact"], cand["delta_n_exact"], surface)
    pair = split["exact_split_pair"]
    assert abs(ht["mH_gev"] - pair["mH_gev"]) < 1e-6
    assert abs(ht["mt_pole_gev"] - pair["mt_pole_gev"]) < 1e-6


def test_candidate_b_reproduces_minimal_conditional_artifact() -> None:
    basis = load_basis()
    with open(RUNS / "d10_ew_minimal_conditional_theorem.json") as f:
        minimal = json.load(f)
    cand = candidate_b_tuple(basis)
    wz = forward_wz(basis, cand["tau2_exact"], cand["delta_n_exact"])
    quintet = minimal["n_c_3_specialization"]["coherent_quintet"]
    assert abs(wz["MW_pole_gev"] - quintet["MW_pole_gev"]) < 1e-6
    assert abs(wz["MZ_pole_gev"] - quintet["MZ_pole_gev"]) < 1e-6


def test_selection_theorem_report_structure() -> None:
    report = build()
    assert report["promotion_allowed"] is False
    assert "A2_charged_contraction" in report["selection_axioms"]
    spread = report["selection_spread"]
    # the candidate spread sits below current experimental resolution
    assert spread["MW_pole_gev"] < 0.0133
    assert spread["mH_gev"] < 0.11
    assert spread["mt_pole_gev"] < 0.6


def test_conditional_predictions_guards_and_envelope() -> None:
    report = build_predictions()
    assert report["row_class"] == "conditional_on_P_and_repair_selection"
    assert report["guards"]["public_promotion_allowed"] is False
    cmp_rows = report["comparison_compare_only"]
    for obs in ("mH_gev", "mt_pole_gev", "MW_pole_gev", "MZ_pole_gev"):
        lo, hi = cmp_rows[obs]["conditional_envelope"]
        assert lo < hi
        assert cmp_rows[obs]["envelope_overlaps_one_sigma_band"]
