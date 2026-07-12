#!/usr/bin/env python3
"""Smoke tests for the particle derivation gap ledger."""

from __future__ import annotations

import json
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "scripts"))

import build_derivation_gap_ledger as gap_ledger  # noqa: E402


build_ledger = gap_ledger.build_ledger


def test_gap_ledger_keeps_compressed_trunk_claim_safe() -> None:
    ledger = build_ledger()

    assert ledger["promotion_policy"]["compressed_p_trunk_is_certified_prediction_root"] is False
    assert ledger["promotion_policy"]["p_trunk_issue_closed_as_guarded_candidate_adoption"] is True
    assert ledger["promotion_policy"]["torus_mode_language_allowed_in_pipeline"] is False
    assert ledger["promotion_policy"]["address_remaining_blockers_one_by_one"] is False
    assert ledger["promotion_policy"]["obstruction_only_worker_result_allowed"] is True
    assert ledger["promotion_policy"]["hadron_backend_in_current_local_scope"] is False
    assert ledger["electroweak_hierarchy"]["claim_status"] == (
        "closed_local_global_hierarchy_and_naturality_certificate"
    )
    assert ledger["electroweak_hierarchy"]["full_theorem_grade_resonance_promoted"] is True
    assert ledger["electroweak_hierarchy"]["remaining_promotion_gates"] == []
    assert ledger["electroweak_hierarchy"]["bridge_residual"].strip("0.") == ""
    assert ledger["electroweak_hierarchy"]["epsilon_H"] == "0"
    gap_ids = {row["id"] for row in ledger["rows"]}
    assert "d10.ward-projected-thomson-endpoint" in gap_ids
    assert "d10.source-residual-map-and-interval-certificate" in gap_ids
    assert "pclosure.certified-codepath-adoption" in gap_ids
    assert "qcd.strong-cp-angle" in gap_ids
    assert "calibration.direct-top-bridge" in gap_ids
    assert "hadron.empirical-ee-spectral-closure" in gap_ids
    row_statuses = {row["id"]: row["status"] for row in ledger["rows"]}
    assert row_statuses["hadron.production-backend-systematics"] == (
        "source_backend_absent_empirical_policy_emitted"
    )
    assert row_statuses["hadron.empirical-ee-spectral-closure"] == (
        "payload_populated_endpoint_evaluated_gap_anchor_localized"
    )
    assert row_statuses["charged.determinant-normalization-transport"] == (
        "closed_current_corpus_charged_end_to_end_no_go"
    )
    charged_row = next(
        row for row in ledger["rows"] if row["id"] == "charged.determinant-normalization-transport"
    )
    assert charged_row["github_issue"] == 546
    assert charged_row["closed_issue_refs"] == [201]
    assert charged_row["trace_lift_claim_label"] == "no_go_confirmed_new_source_needed"
    assert row_statuses["qcd.strong-cp-angle"] == "open_theta_qcd_bar_theta_vanishing_gap"
    assert row_statuses["calibration.direct-top-bridge"] == "closed_current_corpus_codomain_no_go"
    assert row_statuses["d10.ward-projected-thomson-endpoint"] == "closed_blocker_isolated_endpoint_package"
    assert row_statuses["d10.source-residual-map-and-interval-certificate"] == (
        "closed_blocker_isolated_source_residual_no_go"
    )
    assert row_statuses["d10.rg-matching-threshold-scheme"] == "closed_declared_convention_contract"
    assert row_statuses["pclosure.compressed-trunk-artifact"] == "closed_canonical_guarded_candidate_artifact"
    assert row_statuses["pclosure.certified-codepath-adoption"] == "closed_guarded_codepath_adoption"
    assert row_statuses["quark.selected-class-vs-global-classification"] == (
        "selected_class_descent_closed_global_classification_no_go"
    )
    assert row_statuses["quark.source-spread-identifiability"] == (
        "closed_current_corpus_two_modulus_nonidentifiability_obstruction"
    )
    quark_spread_row = next(
        row for row in ledger["rows"] if row["id"] == "quark.source-spread-identifiability"
    )
    assert quark_spread_row["axiom_level_obstruction_certified"] is True
    assert quark_spread_row["axiom_level_obstruction_artifact"].endswith(
        "quark_axiom_level_yukawa_moduli_nonidentifiability.json"
    )
    assert row_statuses["quark.running-scheme-and-yukawa-normalization"] == (
        "closed_structural_scheme_nonidentifiability_obstruction"
    )
    bundle_ids = {bundle["id"] for bundle in ledger["bundles"]}
    assert "electroweak-root-closure-bundle" in bundle_ids
    assert "spectrum-source-bundle" in bundle_ids
    assert "strong-cp-closure-bundle" in bundle_ids
    assert "qcd-thomson-backend-bundle" in bundle_ids
    assert "top-codomain-bridge-bundle" in bundle_ids
    assert "particle-root-integration-gate" in bundle_ids
    bundle_statuses = {bundle["id"]: bundle["status"] for bundle in ledger["bundles"]}
    assert bundle_statuses["qcd-thomson-backend-bundle"] == (
        "source_backend_boundary_empirical_policy_emitted"
    )
    assert bundle_statuses["electroweak-root-closure-bundle"] == "endpoint_package_closed_source_measure_payload_absent"
    assert bundle_statuses["top-codomain-bridge-bundle"] == "closed_current_corpus_codomain_no_go"
    assert bundle_statuses["spectrum-source-bundle"] == "closed_current_corpus_source_boundaries_emitted"
    assert bundle_statuses["strong-cp-closure-bundle"] == "open_physical_invariant_gap"
    assert bundle_statuses["particle-root-integration-gate"] == "keep_candidate_with_constructive_next_artifacts"
    assert ledger["promotion_policy"]["empirical_hadron_closure_class_declared"] is True
    assert ledger["promotion_policy"]["empirical_hadron_closure_source_only"] is False


def _certified_trace_lift_payload() -> dict:
    checks = {
        "certificate_present": True,
        "artifact_type": True,
        "theorem_grade": True,
        "source_only": True,
        "no_target_leak_flag": True,
        "no_target_leak_ancestry": True,
        "numeric_sector_isolated_M_ch": True,
        "source_closed_stage_indexed_q": True,
        "physical_label_map": True,
        "charged_central_projector": True,
        "factorization_certified": True,
        "zero_leakage": True,
        "reference_stage_named": True,
        "lift_constant_source_named": True,
        "lift_constant_numeric": True,
        "lift_constant_interval": True,
        "determinant_scalar_interval": True,
        "P_interval": True,
        "mass_space_affine_anchor_interval": True,
        "mass_space_anchor_source_named": True,
        "residual_certified": True,
        "residual_singleton_zero": True,
    }
    return {
        "artifact": "oph_charged_trace_lift_theorem",
        "claim_label": "trace_lift_certified",
        "source_only": True,
        "factorization_lemma": {
            "status": "certified",
            "leakage_bound": {"interval": ["0", "0"], "certified_zero": True},
        },
        "uncentered_lift_constant": {
            "source_object_name": "oph_charged_determinant_reference_receipt",
            "value": "-3",
        },
        "attachment_identity_residual": {
            "computable": True,
            "certified_zero": True,
            "interval": ["0", "0"],
        },
        "source_certificate": {"checks": checks},
    }


def test_charged_gate_flips_only_for_a_certified_singleton_zero(
    tmp_path: pathlib.Path, monkeypatch
) -> None:
    artifact = tmp_path / "charged_trace_lift_theorem.json"
    monkeypatch.setattr(gap_ledger, "CHARGED_TRACE_LIFT_THEOREM", artifact)

    certified = _certified_trace_lift_payload()
    artifact.write_text(json.dumps(certified), encoding="utf-8")
    assert gap_ledger._charged_trace_lift_gate()["status"] == (
        "trace_lift_certified_conditional_on_P"
    )

    certified["attachment_identity_residual"]["interval"] = ["-1e-12", "1e-12"]
    artifact.write_text(json.dumps(certified), encoding="utf-8")
    assert gap_ledger._charged_trace_lift_gate()["status"] == (
        "closed_current_corpus_charged_end_to_end_no_go"
    )

    certified["attachment_identity_residual"]["interval"] = ["0", "1e-99999"]
    artifact.write_text(json.dumps(certified), encoding="utf-8")
    assert gap_ledger._charged_trace_lift_gate()["status"] == (
        "closed_current_corpus_charged_end_to_end_no_go"
    )

    certified["attachment_identity_residual"]["interval"] = ["0", "0"]
    certified["uncentered_lift_constant"].pop("source_object_name")
    artifact.write_text(json.dumps(certified), encoding="utf-8")
    assert gap_ledger._charged_trace_lift_gate()["status"] == (
        "closed_current_corpus_charged_end_to_end_no_go"
    )
