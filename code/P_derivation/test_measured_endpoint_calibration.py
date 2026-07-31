#!/usr/bin/env python3
"""Tests for the measured endpoint calibration artifact."""

from __future__ import annotations

from decimal import Decimal

from measured_endpoint_calibration import build_measured_endpoint_calibration


def _synthetic_report() -> dict[str, object]:
    return {
        "mode": "synthetic",
        "precision": 10,
        "alpha_inv": "10",
        "alpha": "0.1",
        "source_anchor_alpha_inv": "8",
        "p": "1.2",
        "phi": "1",
        "sqrt_pi": "2",
        "structured_running": {"total_delta_alpha_inv": "2"},
    }


def test_measured_endpoint_calibration_is_external_and_nonpromoting() -> None:
    payload = build_measured_endpoint_calibration(
        _synthetic_report(),
        source_spectral_theorem={
            "current_corpus_obstruction": {
                "production_dump_present": False,
                "finite_volume_levels_populated": False,
            }
        },
        compare_alpha_inv=Decimal("11"),
        compare_alpha_inv_uncertainty=Decimal("0.1"),
        precision=20,
    )

    assert payload["artifact"] == "oph_measured_fine_structure_endpoint_calibration"
    assert payload["status"] == "oph_plus_empirical_hadron_closure_endpoint"
    assert payload["row_class"] == "oph_plus_empirical_hadron_closure"
    assert payload["promotion_allowed"] is False
    assert payload["exact_alpha_promoted"] is False
    assert payload["external_input_used"] is True
    assert payload["empirical_hadron_closure"]["measured_thomson_endpoint_used"] is True
    assert (
        payload["empirical_hadron_closure"][
            "raw_cross_section_data_used_in_this_calibration"
        ]
        is False
    )
    assert (
        payload["empirical_hadron_closure"]["external_input_scope"]
        == "measured_thomson_endpoint_only"
    )
    assert payload["empirical_hadron_closure"]["source_only_theorem_status"] == "not_promoted"
    assert payload["source_only_guard"]["codata_enters_solver"] is False
    assert payload["source_only_guard"]["may_satisfy_source_spectral_payload_gate"] is False
    assert payload["consumer_policy"]["may_feed_numeric_prediction_tables"] is True
    assert payload["consumer_policy"]["may_feed_source_theorem_claim"] is False
    assert payload["current_source_candidate"]["missing_inverse_alpha_units"] == "1"
    assert payload["calibrated_values"]["alpha_inv_0"] == "11"
