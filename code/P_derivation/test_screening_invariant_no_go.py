#!/usr/bin/env python3
"""Tests for the screening-invariant non-identifiability certificate."""

from __future__ import annotations

from copy import deepcopy
from decimal import Decimal
import json
from pathlib import Path

from screening_invariant_no_go import build_no_go


PACKAGE = Path(__file__).resolve().parent / "runtime" / "thomson_endpoint_package_current.json"


def _package() -> dict:
    return json.loads(PACKAGE.read_text(encoding="utf-8"))


def test_target_equivalences_hold() -> None:
    payload = build_no_go(_package())
    assert "generated_utc" not in payload
    assert "github_issue" not in payload
    assert "successor_github_issue" not in payload["required_source_object"]
    eq = payload["target_equivalences"]

    assert eq["all_equal_to_precision"] is True
    tolerance = Decimal(eq["tolerance"])
    assert abs(Decimal(eq["missing_delta"]) - Decimal(eq["quark_naive_times_screen_gap"])) <= tolerance
    assert abs(Decimal(eq["missing_delta"]) - Decimal(eq["quark_naive_times_x2_cq"])) <= tolerance


def test_source_hash_excludes_compare_targets() -> None:
    package = _package()
    payload = build_no_go(package)

    changed = deepcopy(package)
    packet = changed["codata_mapped_endpoint_packet"]["exact_one_loop_package"]
    packet["missing_source_transport_delta_alpha_inv"] = "99"
    packet["required_transport_delta_alpha_inv"] = "100"
    packet["screening_scalar"]["required_screening_factor"] = "101"
    packet["screening_scalar"]["residual_second_order_coefficient"] = "102"
    changed["codata_mapped_endpoint_packet"]["compare_alpha_inv"] = "103"
    changed_payload = build_no_go(changed)

    assert payload["source_packet_hash"] == changed_payload["source_packet_hash"]


def test_beta_candidate_and_alpha_u_candidate_fail() -> None:
    payload = build_no_go(_package())
    failures = payload["candidate_failures"]

    assert failures["c_Q_beta_EW_over_2Nc"]["passes"] is False
    assert Decimal(failures["c_Q_beta_EW_over_2Nc"]["inverse_alpha_overshoot"]) > 0
    assert failures["missing_delta_equals_alpha_U"]["passes"] is False
    assert Decimal(failures["missing_delta_equals_alpha_U"]["missing_minus_candidate"]) != 0


def test_lambda_nonidentifiability_changes_endpoint() -> None:
    payload = build_no_go(_package())
    witness = payload["nonidentifiability_witness"]

    assert witness["source_hash_independent_of_lambda"] is True
    assert Decimal(witness["endpoint_delta_difference"]) != 0
    assert payload["promotion_allowed"] is False
