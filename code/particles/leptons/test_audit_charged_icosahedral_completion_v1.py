"""Tests for the compare-only Pro face-incidence candidate review."""

from __future__ import annotations

import math

import audit_charged_icosahedral_completion_v1 as lane


def test_face_candidate_rebuilds_submitted_masses():
    row = lane.face_candidate(1.6309681897, 246.6174823334856, 0.041124336195630495)
    expected = [
        0.0005109989508433937,
        0.10565837550148128,
        1.7769300000143515,
    ]
    assert all(math.isclose(left, right, rel_tol=0.0, abs_tol=2.0e-15) for left, right in zip(row["masses_gev"], expected, strict=True))


def test_inverse_map_recovers_all_three_spectral_coordinates():
    masses = [0.00051099895069, 0.1056583755, 1.77693]
    coordinates = lane.inverse_coordinates(masses, 246.6174823334856)
    assert math.isclose(coordinates["kappa"], 0.000821392482942, rel_tol=0.0, abs_tol=5.0e-15)
    assert math.isclose(coordinates["chi"], 3.30490696744e-6, rel_tol=0.0, abs_tol=5.0e-15)
    assert math.isclose(coordinates["zeta"], 2.5396723704e-6, rel_tol=0.0, abs_tol=5.0e-15)


def test_review_is_compare_only_and_nonpromoting():
    review = lane.build_review(lane.DEFAULT_ARCHIVE)
    assert review["archive"]["safe_member_paths"] is True
    assert review["archive"]["crc_check_pass"] is True
    assert review["archive"]["no_duplicate_member_names"] is True
    assert review["archive"]["exact_manifest_member_set"] is True
    assert review["archive"]["path"] == "correspondence/oph_charged_icosahedral_completion_v1.zip"
    assert review["archive"]["manifest_checks_pass"] is True
    assert len(review["upstream_receipts"]) == 5
    assert review["archive"]["independent_rebuild_max_abs_mass_difference_gev"] < 2.0e-15
    assert review["compare_only"] is True
    assert review["historical_charged_target_informed"] is True
    assert review["branch_tuple_coherent"] is False
    assert review["public_prediction_promotion_allowed"] is False
    assert review["status"] == "FROZEN_RETROSPECTIVE_TARGET_INFORMED_FACE_CANDIDATE_NOT_COMPLETION"


def test_submitted_branch_mixes_p_and_alpha_u():
    review = lane.build_review(lane.DEFAULT_ARCHIVE)
    provenance = review["branch_provenance"]
    assert provenance["submitted_alpha_U"] != provenance["alpha_U_at_submitted_P_from_D10_recompute"]
    assert provenance["submitted_P_matches_compare_only_D10_probe"] is True
    assert provenance["submitted_v_matches_D10_recompute"] is True
    assert provenance["public_alpha_U_matches_Krawczyk_center"] is True
    assert provenance["source_alpha_U_matches_source_certificate"] is True
    assert provenance["d10_probe_artifact_status"] == "compare_only_diagnostic"


def test_coherent_branch_and_tau_precision_audits_are_frozen():
    review = lane.build_review(lane.DEFAULT_ARCHIVE)
    branches = review["branch_provenance"]["branch_rows"]
    assert math.isclose(
        branches["coherent_canonical_public"]["max_absolute_residual_ppm"],
        0.4239531735672486,
        rel_tol=0.0,
        abs_tol=1.0e-6,
    )
    assert math.isclose(
        branches["coherent_source_audit"]["max_absolute_residual_ppm"],
        84.09987608726243,
        rel_tol=0.0,
        abs_tol=1.0e-6,
    )
    assert math.isclose(
        review["tau_precision_audit"]["residual_ppm_against_packet_api_value"],
        -1.3872894929489732,
        rel_tol=0.0,
        abs_tol=1.0e-6,
    )
