#!/usr/bin/env python3
"""Tests for the Ward-projected source-spectral theorem artifact."""

from __future__ import annotations

from source_spectral_theorem import build_source_spectral_theorem


def _budget() -> dict:
    return {"bound_interval": {"lo": "0", "hi": "1e-30"}}


def _valid_source_payload() -> dict:
    return {
        "artifact": "oph_source_ward_projected_thomson_transport",
        "source_only": True,
        "source_family_id": "d10_running_tree",
        "current": "U1_Q",
        "scheme": {
            "same_subtraction_as_a0": True,
            "scheme_id": "d10_u1q_thomson_v1",
            "normalization_convention": "Q=T3+Y",
        },
        "interval_backend": {
            "theorem_grade": True,
            "directed_outward_rounding": True,
            "library": "arb",
            "certificate": "arb-directed-outward-proof",
        },
        "p_interval": {"lo": "1.6308", "hi": "1.6311"},
        "endpoint_map": {
            "alpha_inv_image": {"lo": "136", "hi": "138"},
            "derivative_abs_bound": "5",
            "transport_error_bound": "1e-30",
            "components": {
                "a0_image": {"lo": "128.3", "hi": "128.4"},
                "Delta_lep_image": {"lo": "4.30", "hi": "4.31"},
                "Delta_had_image": {"lo": "4.30", "hi": "4.40"},
                "Delta_EW_image": {"lo": "0", "hi": "0.01"},
                "R_Q_image": {"lo": "0.03", "hi": "0.05"},
            },
        },
        "source_measure": {
            "artifact": "oph_qcd_ward_projected_hadronic_spectral_measure",
            "finite_volume_levels": [
                {
                    "ensemble_id": "ens0",
                    "channel": "U1_Q",
                    "levels": [{"level_id": "rho0", "s": "1.0", "energy": "1.0", "weight": "1.0"}],
                }
            ],
            "ward_projected_residues": [
                {"level_id": "rho0", "residue": "0.25", "current_normalization": "Q=T3+Y"}
            ],
            "current_normalization": "Q=T3+Y",
            "rho_had_or_measure": {
                "representation": "primitive_spectral_measure",
                "support_variable": "s",
                "pushforward_rule": "certified finite-volume continuum pushforward",
                "positivity_status": "certified_positive",
            },
            "transport_moment_certificate": {
                "kernel": "mZ(P)^2/(3*pi*s*(s+mZ(P)^2))",
                "Delta_had_image": {"lo": "4.30", "hi": "4.40"},
                "quadrature_error_bound": "1e-40",
                "tail_bound": "1e-40",
            },
            "systematics": {
                "statistical_budget": _budget(),
                "continuum_budget": _budget(),
                "finite_volume_budget": _budget(),
                "chiral_budget": _budget(),
                "current_matching_budget": _budget(),
                "quadrature_budget": _budget(),
                "endpoint_remainder_budget": _budget(),
            },
            "guards": {
                "stable_channel_only": False,
                "surrogate_hadron_artifact": False,
                "compare_only_external_endpoint": False,
            },
        },
        "delta_EW": {"source_bound": {"lo": "0", "hi": "0"}},
        "fixed_point_certificate": {
            "self_map_pass": True,
            "uniqueness_pass": True,
            "G_image": {"lo": "1.6308", "hi": "1.6311"},
            "contraction_kappa": "0.001",
        },
    }


def test_source_spectral_theorem_blocks_without_measure_payload() -> None:
    payload = build_source_spectral_theorem()

    assert payload["artifact"] == "oph_ward_projected_source_spectral_theorem"
    assert "generated_utc" not in payload
    assert "github_issue" not in payload
    assert payload["theorem_id"] == "WardProjectedHadronicSpectralEmission_Q"
    assert payload["status"] == "source_spectral_reduction_theorem_emitted_measure_payload_absent"
    assert payload["promotion_allowed"] is False
    assert payload["external_inputs_used"] is False
    assert payload["source_payload_validation"]["status"] == "blocked_source_spectral_measure_missing"
    assert payload["current_corpus_obstruction"]["production_dump_present"] is False
    assert payload["current_corpus_obstruction"]["hadron_contract_classification"] == (
        "scientific_contract_without_production_data"
    )
    assert payload["current_corpus_obstruction"]["finite_volume_levels_populated"] is False
    assert payload["conclusion"]["source_spectral_reduction_closed"] is True
    assert payload["conclusion"]["exact_alpha_emitted_from_current_corpus"] is False


def test_source_spectral_theorem_rejects_fitted_scalars() -> None:
    payload = build_source_spectral_theorem()

    assert payload["source_only_guard"]["residual_fit_allowed"] is False
    assert payload["nonidentifiability_corollary"]["status"] == (
        "constructive_no_external_input_no_go_closed_for_current_source_packet"
    )
    rejected = set(payload["nonidentifiability_corollary"]["fitted_scalars_rejected"])
    assert "c_Q" in rejected
    assert "S_required" in rejected
    assert "missing_source_transport_delta_alpha_inv" in rejected


def test_source_spectral_theorem_constructs_unequal_thomson_moments() -> None:
    payload = build_source_spectral_theorem()
    witness = payload["nonidentifiability_corollary"]["constructive_witness"]

    assert witness["status"] == "constructive_counterexample_emitted"
    assert witness["external_inputs_used"] is False
    assert witness["common_current_corpus_projection"]["source_measure_payload_populated"] is False
    assert witness["measure_a"]["positive"] is True
    assert witness["measure_b"]["positive"] is True
    assert witness["measure_a"]["atoms"][0]["weight"] == witness["measure_b"]["atoms"][0]["weight"]
    assert witness["dimensionless_thomson_moment_a"] == "1/6"
    assert witness["dimensionless_thomson_moment_b"] == "1/12"
    assert witness["dimensionless_thomson_moment_difference"] == "1/12"
    assert witness["moments_equal"] is False


def test_source_spectral_theorem_accepts_complete_source_payload_contract() -> None:
    payload = build_source_spectral_theorem(source_transport_payload=_valid_source_payload())

    assert payload["status"] == "source_spectral_payload_contract_satisfied"
    assert payload["promotion_allowed"] is True
    assert payload["source_payload_validation"]["status"] == "source_transport_interval_certificate_satisfied"
    assert payload["conclusion"]["source_transport_payload_accepted"] is True
    assert payload["conclusion"]["exact_alpha_promotion_allowed_for_supplied_payload"] is True
