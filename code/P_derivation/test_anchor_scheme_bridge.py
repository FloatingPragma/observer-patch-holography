#!/usr/bin/env python3
"""Tests for the anchor scheme-bridge analysis (#545)."""

from __future__ import annotations

import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from anchor_scheme_bridge import build  # noqa: E402


def test_no_measured_value_in_solve_path_and_no_source_bridge_claimed() -> None:
    report = build()
    assert report["guards"]["measured_values_in_any_oph_solve_path"] is False
    assert report["guards"]["source_only_scheme_bridge_emitted"] is False
    assert report["guards"]["public_promotion_allowed"] is False


def test_reference_deficit_inside_certified_gap() -> None:
    report = build()
    ref = report["reference_decomposition_compare_only"]
    gap = report["certified_gap_from_endpoint"]["same_scheme_anchor_gap_interval"]
    # with the published-compilation payload the standard hadronic+HO deficit
    # sits strictly inside the certified interval
    assert gap[0] <= ref["gap_phys_minus_oph"] <= gap[1]
    assert report["verdict"]["reference_deficit_inside_certified_gap"] is True


def test_verdict_is_route_B_false_and_reduces_to_425() -> None:
    report = build()
    v = report["verdict"]
    assert v["route_B_anchor_exact_no_go"] == "false"
    assert v["route_A_scheme_bridge"] == "available_empirical_class_only"
    assert "425" in v["source_only_reduction"]
    assert v["issue_545_status"].startswith("structure_resolved")
