"""Tests for the submitted conditional-theorem bundle review."""

from __future__ import annotations

import mpmath as mp

import audit_charged_face_incidence_theorem_bundle as lane


def test_bundle_is_safe_hash_bound_and_byte_reproduced():
    review = lane.build_review(lane.DEFAULT_ARCHIVE, lane.CANONICAL)
    archive = review["archive"]
    assert archive["archive_hash_matches"] is True
    assert archive["safe_non_symlink_members"] is True
    assert archive["no_duplicate_member_names"] is True
    assert archive["exact_member_set"] is True
    assert archive["crc_check_pass"] is True
    assert archive["declared_hash_checks_pass"] is True
    assert archive["submitted_receipt_hash_matches"] is True


def test_submitted_and_canonical_conditional_numbers_agree():
    review = lane.build_review(lane.DEFAULT_ARCHIVE, lane.CANONICAL)
    agreement = review["independent_numerical_agreement"]
    assert mp.mpf(agreement["contraction_absolute_difference"]) == 0
    assert mp.mpf(agreement["fixed_point_max_absolute_difference"]) == 0
    assert mp.mpf(agreement["mass_max_absolute_difference_mev"]) == 0
    assert all(
        mp.mpf(value) == 0
        for value in agreement["countermodel_mass_max_absolute_differences_mev"]
    )


def test_bundle_remains_compare_only_and_nonpromoting():
    review = lane.build_review(lane.DEFAULT_ARCHIVE, lane.CANONICAL)
    assert review["compare_only"] is True
    assert review["runtime_archive_contains_charged_references"] is True
    assert review["historical_charged_target_informed"] is True
    assert review["global_source_only"] is False
    assert review["public_prediction_promotion_allowed"] is False
    assert review["status"] == (
        "FROZEN_CONDITIONAL_SUPPLEMENT_VERIFIED_SOURCE_ENTAILMENT_OPEN"
    )
