#!/usr/bin/env python3
"""Smoke tests for the particle derivation gap ledger."""

from __future__ import annotations

import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "scripts"))

from build_derivation_gap_ledger import build_ledger  # noqa: E402


def test_gap_ledger_keeps_compressed_trunk_claim_safe() -> None:
    ledger = build_ledger()

    assert ledger["promotion_policy"]["compressed_p_trunk_is_certified_prediction_root"] is False
    assert ledger["promotion_policy"]["torus_mode_language_allowed_in_pipeline"] is False
    assert ledger["promotion_policy"]["address_remaining_blockers_one_by_one"] is False
    assert ledger["promotion_policy"]["obstruction_only_worker_result_allowed"] is True
    assert ledger["promotion_policy"]["hadron_backend_in_current_local_scope"] is False
    gap_ids = {row["id"] for row in ledger["rows"]}
    assert "d10.ward-projected-thomson-endpoint" in gap_ids
    assert "d10.source-residual-map-and-interval-certificate" in gap_ids
    assert "pclosure.certified-codepath-adoption" in gap_ids
    assert "calibration.direct-top-bridge" in gap_ids
    row_statuses = {row["id"]: row["status"] for row in ledger["rows"]}
    assert row_statuses["hadron.production-backend-systematics"] == "closed_out_of_scope_computationally_blocked"
    assert row_statuses["charged.determinant-normalization-transport"] == (
        "closed_current_corpus_charged_end_to_end_no_go"
    )
    assert row_statuses["calibration.direct-top-bridge"] == "closed_current_corpus_codomain_no_go"
    assert row_statuses["d10.ward-projected-thomson-endpoint"] == "closed_blocker_isolated_endpoint_package"
    assert row_statuses["d10.source-residual-map-and-interval-certificate"] == (
        "open_source_residual_map_and_interval_certificate"
    )
    assert row_statuses["quark.selected-class-vs-global-classification"] == (
        "selected_class_closed_global_classification_no_go"
    )
    bundle_ids = {bundle["id"] for bundle in ledger["bundles"]}
    assert "electroweak-root-closure-bundle" in bundle_ids
    assert "spectrum-source-bundle" in bundle_ids
    assert "qcd-thomson-backend-bundle" in bundle_ids
    assert "top-codomain-bridge-bundle" in bundle_ids
    assert "particle-root-integration-gate" in bundle_ids
    bundle_statuses = {bundle["id"]: bundle["status"] for bundle in ledger["bundles"]}
    assert bundle_statuses["qcd-thomson-backend-bundle"] == "closed_out_of_scope_scope_lock_emitted"
    assert bundle_statuses["top-codomain-bridge-bundle"] == "closed_current_corpus_codomain_no_go"
    assert bundle_statuses["spectrum-source-bundle"] == "closed_current_corpus_source_boundaries_emitted"
    assert bundle_statuses["particle-root-integration-gate"] == "keep_candidate_with_constructive_next_artifacts"
