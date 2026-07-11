from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import numpy as np
import pytest


MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "programs"
    / "cayley_blind_likelihood_analysis.py"
)
SPEC = importlib.util.spec_from_file_location("cayley_blind_likelihood_analysis", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
analysis = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = analysis
SPEC.loader.exec_module(analysis)


def _probability_table(valid_codes: list[int], weights: list[float]) -> dict[str, float]:
    triples = (
        (0, 1, 0),
        (1, 1, 0),
        (1, 0, 1),
        (2, 0, 2),
        (3, 1, 3),
        (4, 0, 4),
    )
    assert math.isclose(sum(weights), 1.0, abs_tol=1e-12)
    return {
        analysis.outcome_key(heated, decision, final, valid_codes): probability
        for (heated, decision, final), probability in zip(triples, weights)
    }


def _candidates(valid_codes: list[int]) -> dict[str, dict[str, float]]:
    return {
        "record_gated_repair": _probability_table(
            valid_codes,
            [0.58, 0.22, 0.10, 0.05, 0.03, 0.02],
        ),
        "lazy_heat": _probability_table(
            valid_codes,
            [0.12, 0.16, 0.18, 0.20, 0.18, 0.16],
        ),
        "delayed_record": _probability_table(
            valid_codes,
            [0.22, 0.20, 0.18, 0.16, 0.13, 0.11],
        ),
        "shuffled_record": _probability_table(
            valid_codes,
            [0.26, 0.18, 0.16, 0.14, 0.14, 0.12],
        ),
        "inverted_record": _probability_table(
            valid_codes,
            [0.02, 0.03, 0.05, 0.10, 0.30, 0.50],
        ),
        "label_only": _probability_table(
            valid_codes,
            [0.10, 0.10, 0.20, 0.20, 0.20, 0.20],
        ),
        "calibrated_noise": _probability_table(
            valid_codes,
            [1.0 / 6.0] * 6,
        ),
    }


def _build_lock() -> dict:
    s3_codes = [0, 1, 2, 3, 4, 5]
    z5_codes = [0, 1, 2, 3, 4]
    calibrations = {
        "cal_backend_a": analysis.factorized_calibration_packet(
            state_error=0.015,
            decision_error=0.01,
            receipt_sha256="1" * 64,
        ),
        "cal_backend_b": analysis.factorized_calibration_packet(
            state_error=0.018,
            decision_error=0.012,
            receipt_sha256="2" * 64,
        ),
    }
    expected_rows = [
        {
            "row_id": "opaque_primary_a",
            "endpoint": "primary_s3",
            "backend_id": "heldout_backend_a",
            "layout_id": "opaque_layout_1",
            "calibration_id": "cal_backend_a",
            "shots": 5000,
            "valid_codes": s3_codes,
            "candidate_probabilities": _candidates(s3_codes),
        },
        {
            "row_id": "opaque_primary_b",
            "endpoint": "primary_s3",
            "backend_id": "heldout_backend_b",
            "layout_id": "opaque_layout_2",
            "calibration_id": "cal_backend_b",
            "shots": 5000,
            "valid_codes": s3_codes,
            "candidate_probabilities": _candidates(s3_codes),
        },
        {
            "row_id": "opaque_secondary_z5",
            "endpoint": "secondary_z5",
            "backend_id": "heldout_backend_a",
            "layout_id": "opaque_layout_1",
            "calibration_id": "cal_backend_a",
            "shots": 2500,
            "valid_codes": z5_codes,
            "candidate_probabilities": _candidates(z5_codes),
        },
    ]
    return analysis.build_analysis_lock(
        expected_rows=expected_rows,
        calibrations=calibrations,
        secondary_tests=[
            {
                "test_id": "secondary_z5_goodness_of_fit",
                "row_ids": ["opaque_secondary_z5"],
                "model": "record_gated_repair",
            }
        ],
        blind_manifest_commitment="a" * 64,
        created_utc="2026-07-11T00:00:00+00:00",
    )


def _reseal(lock: dict, rows: list[dict]) -> dict:
    return analysis.seal_data_packet(
        analysis_lock_sha256=lock["analysis_lock_sha256"],
        rows=rows,
        created_utc="synthetic-resealed",
    )


def test_factorized_calibration_convolution_is_normalized_and_exposes_leakage() -> None:
    valid_codes = [0, 1, 2, 3, 4, 5]
    latent = analysis.probability_mapping_to_vector(
        _candidates(valid_codes)["record_gated_repair"],
        valid_codes,
    )
    calibration = analysis.factorized_calibration_packet(
        state_error=0.02,
        decision_error=0.01,
        receipt_sha256="3" * 64,
    )
    observed = analysis.convolve_calibration(latent, calibration)
    assert observed.shape == (128,)
    assert np.all(observed > 0.0)
    assert math.isclose(float(observed.sum()), 1.0, abs_tol=1e-12)
    leakage_mass = sum(
        observed[index]
        for index in range(analysis.JOINT_CARDINALITY)
        if analysis.leakage_bit(
            analysis.joint_tuple(index)[0],
            analysis.joint_tuple(index)[2],
            valid_codes,
        )
    )
    assert leakage_mass > 0.0


def test_synthetic_repair_process_clears_all_frozen_primary_gates() -> None:
    lock = _build_lock()
    data = analysis.simulate_data_packet(
        lock,
        generating_model="record_gated_repair",
        seed=509,
    )
    report = analysis.run_blind_analysis(lock, data)

    assert report["decision"]["verdict"] == "passes_frozen_reduced_repair_kernel_gate"
    assert report["decision"]["primary_pass"] is True
    primary = report["primary_endpoint"]
    assert primary["backend_count"] == 2
    assert primary["global_99_percent_simultaneous_envelope"]["pass"] is True
    assert set(primary["pooled_bayes_factors"]) == set(analysis.REQUIRED_NULL_MODELS)
    assert primary["pooled_bayes_factors"]["delayed_record"][
        "passes_pooled_threshold"
    ] is True
    assert all(
        item["passes_pooled_threshold"]
        for item in primary["pooled_bayes_factors"].values()
    )
    assert all(
        item["passes_per_backend_threshold"]
        for backend in primary["per_backend_bayes_factors"].values()
        for item in backend.values()
    )
    assert report["secondary_family"]["correction"] == "Holm"
    assert len(report["secondary_family"]["tests"]) == 1
    assert any(row["leakage_shots"] > 0 for row in report["shot_audit"])
    assert all(
        row["declared_submitted_retrieved_counted_shots"] > 0
        for row in report["shot_audit"]
    )

    report_body = dict(report)
    claimed_report_hash = report_body.pop("blind_report_sha256")
    assert claimed_report_hash == analysis.sha256_json(report_body)
    assert report["analysis_lock_sha256"] == lock["analysis_lock_sha256"]
    assert report["heldout_data_packet_sha256"] == data["data_packet_sha256"]


def test_synthetic_lazy_heat_process_favors_null_and_rejects_repair() -> None:
    lock = _build_lock()
    data = analysis.simulate_data_packet(
        lock,
        generating_model="lazy_heat",
        seed=510,
    )
    report = analysis.run_blind_analysis(lock, data)

    assert report["decision"]["primary_pass"] is False
    assert report["decision"]["kernel_failure"] is True
    assert report["decision"]["verdict"] == "fails_frozen_reduced_repair_kernel"
    assert (
        report["primary_endpoint"]["pooled_bayes_factors"]["lazy_heat"][
            "log_bayes_factor_repair_over_null"
        ]
        < -math.log(100.0)
    )
    assert all(
        backend["lazy_heat"]["log_bayes_factor_repair_over_null"]
        < -math.log(100.0)
        for backend in report["primary_endpoint"]["per_backend_bayes_factors"].values()
    )


def test_missing_expected_row_fails_closed_even_when_packet_is_resealed() -> None:
    lock = _build_lock()
    data = analysis.simulate_data_packet(
        lock,
        generating_model="record_gated_repair",
        seed=511,
    )
    missing = _reseal(lock, list(data["rows"][:-1]))
    with pytest.raises(analysis.AnalysisValidationError, match="held-out row mismatch"):
        analysis.run_blind_analysis(lock, missing)


def test_postselected_flag_fails_closed_even_when_packet_is_resealed() -> None:
    lock = _build_lock()
    data = analysis.simulate_data_packet(
        lock,
        generating_model="record_gated_repair",
        seed=512,
    )
    rows = analysis._json_copy(data["rows"])
    rows[0]["postselected"] = True
    postselected = _reseal(lock, rows)
    with pytest.raises(analysis.AnalysisValidationError, match="postselected"):
        analysis.run_blind_analysis(lock, postselected)


def test_dropped_leakage_or_other_shots_fails_shot_conservation() -> None:
    lock = _build_lock()
    data = analysis.simulate_data_packet(
        lock,
        generating_model="record_gated_repair",
        seed=513,
    )
    rows = analysis._json_copy(data["rows"])
    first_key = next(iter(rows[0]["counts"]))
    rows[0]["counts"][first_key] -= 1
    if rows[0]["counts"][first_key] == 0:
        del rows[0]["counts"][first_key]
    incomplete = _reseal(lock, rows)
    with pytest.raises(analysis.AnalysisValidationError, match="counted .* shots"):
        analysis.run_blind_analysis(lock, incomplete)


def test_lock_and_data_hash_mutations_fail_closed() -> None:
    lock = _build_lock()
    changed_lock = analysis._json_copy(lock)
    changed_lock["thresholds"]["pooled_bayes_factor"] = 99.0
    with pytest.raises(analysis.AnalysisValidationError, match="lock hash mismatch"):
        analysis.validate_analysis_lock(changed_lock)

    data = analysis.simulate_data_packet(
        lock,
        generating_model="record_gated_repair",
        seed=514,
    )
    changed_data = analysis._json_copy(data)
    changed_data["rows"][0]["all_outcomes_included"] = False
    with pytest.raises(analysis.AnalysisValidationError, match="data packet hash mismatch"):
        analysis.run_blind_analysis(lock, changed_data)


def test_lock_without_delayed_record_table_is_rejected() -> None:
    valid_codes = [0, 1, 2, 3, 4, 5]
    candidates = _candidates(valid_codes)
    candidates.pop("delayed_record")
    calibration = analysis.factorized_calibration_packet(
        state_error=0.01,
        decision_error=0.01,
        receipt_sha256="4" * 64,
    )
    rows = [
        {
            "row_id": f"primary_{backend}",
            "endpoint": "primary_s3",
            "backend_id": backend,
            "layout_id": f"layout_{backend}",
            "calibration_id": "cal",
            "shots": 100,
            "valid_codes": valid_codes,
            "candidate_probabilities": candidates,
        }
        for backend in ("a", "b")
    ]
    rows.append(
        {
            "row_id": "secondary",
            "endpoint": "secondary_z5",
            "backend_id": "a",
            "layout_id": "layout_a",
            "calibration_id": "cal",
            "shots": 100,
            "valid_codes": valid_codes,
            "candidate_probabilities": candidates,
        }
    )
    with pytest.raises(analysis.AnalysisValidationError, match="candidate set"):
        analysis.build_analysis_lock(
            expected_rows=rows,
            calibrations={"cal": calibration},
            secondary_tests=[
                {"test_id": "secondary", "row_ids": ["secondary"], "model": "lazy_heat"}
            ],
            blind_manifest_commitment="b" * 64,
            created_utc="2026-07-11T00:00:00+00:00",
        )


def test_holm_adjustment_is_monotone_and_familywise() -> None:
    results = analysis.holm_adjust(
        [("a", 0.01), ("b", 0.04), ("c", 0.03)],
        alpha=0.05,
    )
    by_id = {result["test_id"]: result for result in results}
    assert math.isclose(by_id["a"]["holm_adjusted_p_value"], 0.03)
    assert math.isclose(by_id["b"]["holm_adjusted_p_value"], 0.06)
    assert math.isclose(by_id["c"]["holm_adjusted_p_value"], 0.06)
    assert by_id["a"]["holm_reject_at_family_alpha"] is True
    assert by_id["b"]["holm_reject_at_family_alpha"] is False
    assert by_id["c"]["holm_reject_at_family_alpha"] is False
