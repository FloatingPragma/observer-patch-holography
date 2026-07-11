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
    repair = _probability_table(
        valid_codes,
        [0.58, 0.22, 0.10, 0.05, 0.03, 0.02],
    )
    candidates = {
        "record_gated_repair": repair,
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
        "calibrated_noise": _probability_table(
            valid_codes,
            [1.0 / 6.0] * 6,
        ),
    }
    candidates["state_preparation_only"] = analysis.state_preparation_only_joint_null(
        repair, valid_codes
    )
    return candidates


def _row(
    *,
    row_id: str,
    endpoint: str,
    backend_role: str,
    backend_name: str,
    layout_id: str,
    physical_layout: list[int],
    calibration_id: str,
    shots: int,
    valid_codes: list[int],
    family: str = "cayley",
) -> dict:
    candidates = _candidates(valid_codes)
    logical_circuit_sha256 = analysis.sha256_json({"opaque_id": row_id})
    return {
        "row_id": row_id,
        "opaque_id": row_id,
        "logical_circuit_sha256": logical_circuit_sha256,
        "family": family,
        "state_preparation_stratum_id": f"stratum_{row_id}",
        "endpoint": endpoint,
        "backend_role": backend_role,
        "backend_name": backend_name,
        "layout_id": layout_id,
        "physical_layout": physical_layout,
        "calibration_id": calibration_id,
        "shots": shots,
        "max_leakage_fraction": 0.25,
        "valid_codes": valid_codes,
        "candidate_probabilities": candidates,
        "candidate_provenance": analysis.build_candidate_provenance(
            candidates,
            {model: analysis.sha256_json({"model": model}) for model in candidates},
            logical_circuit_sha256,
        ),
    }


def _build_lock(*, factorized: bool = False, shots: int = 5000) -> dict:
    s3_codes = [0, 1, 2, 3, 4, 5]
    z5_codes = [0, 1, 2, 3, 4]
    if factorized:
        calibrations = {
            "cal_backend_a": analysis.factorized_calibration_packet(
                state_error=0.015,
                decision_error=0.01,
                receipt_sha256="1" * 64,
                diagnostic_opaque_ids=("synthetic-calibration-a",),
            ),
            "cal_backend_b": analysis.factorized_calibration_packet(
                state_error=0.018,
                decision_error=0.012,
                receipt_sha256="2" * 64,
                diagnostic_opaque_ids=("synthetic-calibration-b",),
            ),
        }
    else:
        calibrations = {
            "cal_backend_a": analysis.contamination_calibration_packet(
                contamination_probability=0.03,
                sensitivity_probabilities=[0.02, 0.05],
                receipt_sha256="1" * 64,
                derivation_sha256="6" * 64,
                diagnostic_opaque_ids=("synthetic-calibration-a",),
            ),
            "cal_backend_b": analysis.contamination_calibration_packet(
                contamination_probability=0.035,
                sensitivity_probabilities=[0.02, 0.055],
                receipt_sha256="2" * 64,
                derivation_sha256="7" * 64,
                diagnostic_opaque_ids=("synthetic-calibration-b",),
            ),
        }
    expected_rows = [
        _row(
            row_id="opaque_primary_a",
            endpoint="primary_s3",
            backend_role="backend_role_a",
            backend_name="heldout_backend_a",
            layout_id="opaque_layout_1",
            physical_layout=[1, 2, 3, 4],
            calibration_id="cal_backend_a",
            shots=shots,
            valid_codes=s3_codes,
        ),
        _row(
            row_id="opaque_primary_b",
            endpoint="primary_s3",
            backend_role="backend_role_b",
            backend_name="heldout_backend_b",
            layout_id="opaque_layout_2",
            physical_layout=[5, 6, 7, 8],
            calibration_id="cal_backend_b",
            shots=shots,
            valid_codes=s3_codes,
        ),
        _row(
            row_id="opaque_secondary_z5",
            endpoint="secondary_z5",
            backend_role="backend_role_a",
            backend_name="heldout_backend_a",
            layout_id="opaque_layout_1",
            physical_layout=[1, 2, 3, 4],
            calibration_id="cal_backend_a",
            shots=shots,
            valid_codes=z5_codes,
        ),
    ]
    label_components = []
    label_weights = (
        [0.10, 0.10, 0.20, 0.20, 0.20, 0.20],
        [0.20, 0.20, 0.10, 0.10, 0.20, 0.20],
    )
    for index, weights in enumerate(label_weights):
        row_probabilities = {
            row["row_id"]: _probability_table(row["valid_codes"], list(weights))
            for row in expected_rows
            if row["endpoint"] == "primary_s3"
        }
        label_components.append(
            analysis.build_label_layout_component(
                component_id=f"opaque_label_{index}",
                prior_weight=0.5,
                row_probabilities=row_probabilities,
                component_derivation_sha256=str(8 + index) * 64,
                row_derivation_sha256={
                    row_id: str(8 + index) * 64 for row_id in row_probabilities
                },
            )
        )
    label_layout_model = {
        "mapping_scope": "global_shared_across_primary_rows",
        "component_set_derivation_sha256": "a" * 64,
        "components": label_components,
        "reference_component_id": label_components[1]["component_id"],
    }
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
        catalog_precommitment_sha256="b" * 64,
        label_layout_model=label_layout_model,
        created_utc="2026-07-11T00:00:00+00:00",
    )


def _reseal(lock: dict, rows: list[dict]) -> dict:
    calibration_results = {
        calibration_id: analysis.basis_calibration_control_result(
            diagnostic_counts_by_opaque_id={
                opaque_id: {f"{basis_code:04b}": 512}
                for opaque_id, basis_code in calibration["control_rule"][
                    "expected_basis_code_by_opaque_id"
                ].items()
            },
            control_rule=calibration["control_rule"],
            provider_job_ids=[f"synthetic-calibration-job-{calibration_id}"],
            calibration_receipt_sha256=analysis.sha256_json(
                {"synthetic_calibration": calibration_id}
            ),
        )
        for calibration_id, calibration in lock["calibrations"].items()
    }
    return analysis.seal_data_packet(
        analysis_lock_sha256=lock["analysis_lock_sha256"],
        rows=rows,
        manifest_sha256="e" * 64,
        submission_journal_sha256="c" * 64,
        harvest_journal_sha256="d" * 64,
        source_kind="synthetic_preflight",
        calibration_results=calibration_results,
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
    assert calibration["primary_eligible"] is False
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
    assert set(primary["pooled_conditional_likelihood_ratios"]) == set(
        analysis.REQUIRED_NULL_MODELS
    )
    assert primary["pooled_conditional_likelihood_ratios"]["delayed_record"][
        "passes_pooled_threshold"
    ] is True
    assert all(
        item["passes_pooled_threshold"]
        for item in primary["pooled_conditional_likelihood_ratios"].values()
    )
    assert all(
        item["passes_per_backend_threshold"]
        for backend in primary["per_backend_conditional_likelihood_ratios"].values()
        for item in backend.values()
    )
    multiplicity = primary["label_layout_multiplicity"]
    assert multiplicity["component_count"] == 2
    assert multiplicity["reference_component_is_unique_top"] is True
    assert math.isclose(
        sum(
            item["posterior_within_label_model"]
            for item in multiplicity["ranked_components"]
        ),
        1.0,
        abs_tol=1e-12,
    )
    assert report["secondary_family"]["correction"] == "Holm"
    assert len(report["secondary_family"]["tests"]) == 1
    assert any(row["leakage_shots"] > 0 for row in report["shot_audit"])
    assert all(
        row["declared_submitted_retrieved_counted_shots"] > 0
        for row in report["shot_audit"]
    )
    assert all(
        "logical_circuit_sha256" in row and "logical_qpy_sha256" not in row
        for row in report["shot_audit"]
    )
    assert all(
        candidate["logical_circuit_sha256"] == lock_row["logical_circuit_sha256"]
        for lock_row in lock["expected_rows"]
        for candidate in report["candidate_probability_hashes"][lock_row["row_id"]].values()
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
        report["primary_endpoint"]["pooled_conditional_likelihood_ratios"][
            "lazy_heat"
        ][
            "conditional_log_likelihood_ratio_repair_over_null"
        ]
        < -math.log(100.0)
    )
    assert all(
        backend["lazy_heat"]["sensitivity_log_likelihood_ratio_bounds"][1]
        < -math.log(100.0)
        for backend in report["primary_endpoint"][
            "per_backend_conditional_likelihood_ratios"
        ].values()
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
    changed_lock["thresholds"]["pooled_likelihood_ratio"] = 99.0
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
    valid = _build_lock()
    rows = analysis._json_copy(valid["expected_rows"])
    rows[0]["candidate_probabilities"].pop("delayed_record")
    rows[0]["candidate_provenance"].pop("delayed_record")
    with pytest.raises(analysis.AnalysisValidationError, match="candidate set"):
        analysis.build_analysis_lock(
            expected_rows=rows,
            calibrations=valid["calibrations"],
            secondary_tests=valid["secondary_tests"],
            catalog_precommitment_sha256="b" * 64,
            label_layout_model=valid["label_layout_model"],
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


def test_qiskit_joined_bit_converter_is_exhaustive_and_strict() -> None:
    valid_codes = [0, 1, 2, 3, 4, 5]
    joined_counts = {}
    for heated in range(8):
        for decision in range(2):
            for final in range(8):
                joined = f"{final:03b}{decision}{heated:03b}"
                expected = analysis.outcome_key(heated, decision, final, valid_codes)
                assert (
                    analysis.qiskit_joined_key_to_outcome_key(joined, valid_codes)
                    == expected
                )
                joined_counts[joined] = 1
    converted = analysis.qiskit_joined_counts_to_analysis_counts(
        joined_counts,
        valid_codes,
        expected_shots=128,
    )
    assert set(converted) == set(analysis.all_outcome_keys(valid_codes))
    assert sum(converted.values()) == 128
    with pytest.raises(analysis.AnalysisValidationError, match="exactly seven"):
        analysis.qiskit_joined_key_to_outcome_key("011 1 100", valid_codes)
    with pytest.raises(analysis.AnalysisValidationError, match="contain 1 of 2"):
        analysis.qiskit_joined_counts_to_analysis_counts(
            {"0111100": 1}, valid_codes, expected_shots=2
        )


def test_basis_calibration_control_uses_exact_binomial_bonferroni_rule() -> None:
    calibration = analysis.contamination_calibration_packet(
        contamination_probability=0.08,
        sensitivity_probabilities=[0.04, 0.15],
        receipt_sha256="a" * 64,
        derivation_sha256="b" * 64,
        diagnostic_opaque_ids=("cal_0", "cal_1"),
        diagnostic_expected_codes={"cal_0": 0, "cal_1": 15},
    )
    passing = analysis.basis_calibration_control_result(
        diagnostic_counts_by_opaque_id={
            "cal_0": {"0000": 500, "0001": 12},
            "cal_1": {"1111": 496, "1110": 16},
        },
        control_rule=calibration["control_rule"],
        provider_job_ids=["job-cal"],
        calibration_receipt_sha256="c" * 64,
    )
    assert passing["gof_p_value"] >= analysis.CALIBRATION_CONTROL_MIN_P_VALUE
    assert passing["minimum_count_per_prepared_state"] == 512
    failing = analysis.basis_calibration_control_result(
        diagnostic_counts_by_opaque_id={
            "cal_0": {"0000": 300, "0001": 212},
            "cal_1": {"1111": 496, "1110": 16},
        },
        control_rule=calibration["control_rule"],
        provider_job_ids=["job-cal"],
        calibration_receipt_sha256="d" * 64,
    )
    assert failing["gof_p_value"] < analysis.CALIBRATION_CONTROL_MIN_P_VALUE


def test_factorized_dynamic_calibration_can_only_produce_invalid_verdict() -> None:
    lock = _build_lock(factorized=True)
    data = analysis.simulate_data_packet(
        lock,
        generating_model="record_gated_repair",
        seed=515,
    )
    report = analysis.run_blind_analysis(lock, data)
    assert report["validity"]["calibration_gate_pass"] is False
    assert report["decision"]["primary_pass"] is False
    assert report["decision"]["kernel_failure"] is False
    assert report["decision"]["verdict"] == "invalid_calibration_or_leakage_gate"


def test_excess_leakage_is_retained_and_invalidates_without_kernel_failure() -> None:
    lock = _build_lock()
    data = analysis.simulate_data_packet(
        lock,
        generating_model="record_gated_repair",
        seed=516,
    )
    rows = analysis._json_copy(data["rows"])
    shots = rows[0]["declared_shots"]
    rows[0]["counts"] = {"h7|d0|f7|l1": shots}
    rows[0]["raw_joined_counts_sha256"] = analysis.sha256_json(rows[0]["counts"])
    leaked = _reseal(lock, rows)
    report = analysis.run_blind_analysis(lock, leaked)
    assert report["validity"]["leakage_gate_pass"] is False
    assert report["decision"]["kernel_failure"] is False
    assert report["decision"]["verdict"] == "invalid_calibration_or_leakage_gate"
    audited = {row["row_id"]: row for row in report["shot_audit"]}
    assert audited["opaque_primary_a"]["leakage_shots"] == shots


def test_data_identity_must_match_locked_opaque_and_logical_circuit() -> None:
    lock = _build_lock(shots=192)
    assert all(row["shots"] == 192 for row in lock["expected_rows"])
    data = analysis.simulate_data_packet(
        lock,
        generating_model="record_gated_repair",
        seed=517,
    )
    rows = analysis._json_copy(data["rows"])
    rows[0]["logical_circuit_sha256"] = "f" * 64
    changed = _reseal(lock, rows)
    with pytest.raises(
        analysis.AnalysisValidationError,
        match="locked logical_circuit_sha256",
    ):
        analysis.run_blind_analysis(lock, changed)


def test_removed_logical_qpy_identity_is_rejected_in_lock_and_data() -> None:
    valid = _build_lock(shots=192)
    lock_rows = analysis._json_copy(valid["expected_rows"])
    lock_rows[0]["logical_qpy_sha256"] = lock_rows[0]["logical_circuit_sha256"]
    with pytest.raises(analysis.AnalysisValidationError, match="removed logical_qpy_sha256"):
        analysis.build_analysis_lock(
            expected_rows=lock_rows,
            calibrations=valid["calibrations"],
            secondary_tests=valid["secondary_tests"],
            catalog_precommitment_sha256=valid["catalog_precommitment_sha256"],
            label_layout_model=valid["label_layout_model"],
            created_utc="synthetic-legacy-lock",
        )

    data = analysis.simulate_data_packet(
        valid,
        generating_model="record_gated_repair",
        seed=518,
    )
    data_rows = analysis._json_copy(data["rows"])
    data_rows[0]["logical_qpy_sha256"] = data_rows[0]["logical_circuit_sha256"]
    changed = _reseal(valid, data_rows)
    with pytest.raises(analysis.AnalysisValidationError, match="removed logical_qpy_sha256"):
        analysis.run_blind_analysis(valid, changed)


def test_candidate_provenance_binds_stable_logical_circuit_hash() -> None:
    valid = _build_lock()
    rows = analysis._json_copy(valid["expected_rows"])
    first_model = next(iter(rows[0]["candidate_provenance"]))
    rows[0]["candidate_provenance"][first_model]["logical_circuit_sha256"] = "f" * 64
    with pytest.raises(analysis.AnalysisValidationError, match="different logical circuit"):
        analysis.build_analysis_lock(
            expected_rows=rows,
            calibrations=valid["calibrations"],
            secondary_tests=valid["secondary_tests"],
            catalog_precommitment_sha256=valid["catalog_precommitment_sha256"],
            label_layout_model=valid["label_layout_model"],
            created_utc="synthetic-provenance-tamper",
        )


def test_candidate_probability_table_hash_is_enforced_at_lock_build() -> None:
    valid = _build_lock()
    rows = analysis._json_copy(valid["expected_rows"])
    table = rows[0]["candidate_probabilities"]["lazy_heat"]
    first, second = list(table)[:2]
    table[first], table[second] = table[second], table[first]
    with pytest.raises(analysis.AnalysisValidationError, match="probability-table hash mismatch"):
        analysis.build_analysis_lock(
            expected_rows=rows,
            calibrations=valid["calibrations"],
            secondary_tests=valid["secondary_tests"],
            catalog_precommitment_sha256="b" * 64,
            label_layout_model=valid["label_layout_model"],
            created_utc="2026-07-11T00:00:00+00:00",
        )


def test_label_marginal_is_one_global_shared_mapping_not_rowwise_mixture() -> None:
    lock = _build_lock()
    data = analysis.simulate_data_packet(
        lock,
        generating_model="record_gated_repair",
        seed=518,
    )
    report = analysis.run_blind_analysis(lock, data)
    label = report["primary_endpoint"]["label_layout_multiplicity"]
    manual = analysis._logsumexp(
        [item["log_weighted_likelihood"] for item in label["ranked_components"]]
    )
    assert math.isclose(label["log_marginal_likelihood"], manual, abs_tol=1e-12)
    assert label["mapping_scope"] == "global_shared_across_primary_rows"
