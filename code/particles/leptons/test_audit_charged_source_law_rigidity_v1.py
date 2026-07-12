"""Tests for the compare-only source-law rigidity archive audit."""

from __future__ import annotations

import json

import mpmath as mp

import audit_charged_source_law_rigidity_v1 as lane


def review() -> dict:
    return lane.build_review(lane.DEFAULT_ARCHIVE, lane.CANONICAL_CFQ, lane.DECLARED_MAP)


def test_archive_is_safe_hash_bound_and_crc_clean():
    result = review()
    assert result["archive"]["checks_pass"] is True
    assert all(result["archive"]["checks"].values())
    assert result["archive"]["sha256"] == lane.EXPECTED_ARCHIVE_SHA256
    assert result["archive"]["submitted_receipt_sha256"] == lane.EXPECTED_RECEIPT_SHA256


def test_submitted_packet_matches_canonical_conditional_theorem():
    result = review()
    agreement = result["independent_agreement"]
    assert agreement["checks_pass"] is True
    assert all(agreement["registers"].values())
    assert all(agreement["paths"].values())
    assert all(agreement["graph_counts"].values())
    assert max(
        abs(mp.mpf(value)) for value in agreement["coefficient_differences"].values()
    ) < mp.mpf("1e-75")
    assert mp.mpf(agreement["fixed_point_max_absolute_difference"]) == 0
    assert mp.mpf(agreement["conditional_mass_max_absolute_difference_mev"]) == 0


def test_archive_remains_compare_only_nonpromoting_and_historically_exposed():
    result = review()
    assert result["compare_only"] is True
    assert result["forbidden_as_candidate_ancestor"] is True
    assert result["runtime_archive_contains_downstream_charged_mass_values"] is True
    assert result["historical_charged_target_informed"] is True
    assert result["global_source_only"] is False
    assert result["branch_tuple_coherent"] is False
    assert result["mass_scheme_certified"] is False
    assert result["public_prediction_promotion_allowed"] is False
    assert "current OPH axioms imply CFQ1-CFQ8" in result["not_established"]
    assert result["checks_pass"] is True


def test_committed_review_is_byte_reproducible():
    result = review()
    expected = (json.dumps(result, indent=2, sort_keys=True) + "\n").encode()
    assert lane.DEFAULT_OUT.read_bytes() == expected
