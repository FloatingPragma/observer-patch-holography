from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path


PROGRAMS = Path(__file__).resolve().parents[1] / "programs"
sys.path.insert(0, str(PROGRAMS))

import blind_preregister  # noqa: E402
import cayley_blind_likelihood_analysis as analysis  # noqa: E402
from build_blinded_analysis_lock import build_locked_analysis  # noqa: E402


def test_full_catalog_lock_coverage_components_and_ideal_power() -> None:
    public, reveal = blind_preregister.build_blinded_preregistration(
        test_seed="full-analysis-lock-integration-test"
    )
    lock = build_locked_analysis(public, reveal)

    assert len(public["circuits"]) == 3104
    assert len(lock["expected_rows"]) == 3072
    public_logical_hashes = {
        row["opaque_id"]: row["logical_circuit_sha256"] for row in public["circuits"]
    }
    assert all(
        row["logical_circuit_sha256"]
        == public_logical_hashes[row["opaque_id"]]
        and "logical_qpy_sha256" not in row
        for row in lock["expected_rows"]
    )
    assert all(
        provenance["logical_circuit_sha256"] == row["logical_circuit_sha256"]
        for row in lock["expected_rows"]
        for provenance in row["candidate_provenance"].values()
    )
    assert len(lock["catalog_coverage"]["diagnostic_calibration_opaque_ids"]) == 32
    assert len(lock["label_layout_model"]["components"]) == 256
    table_hashes = {
        analysis.sha256_json(component["row_probabilities"])
        for component in lock["label_layout_model"]["components"]
    }
    assert len(table_hashes) == 256
    assert lock["label_layout_model"]["reference_component_id"] in {
        component["component_id"]
        for component in lock["label_layout_model"]["components"]
    }

    covered = set(lock["catalog_coverage"]["dynamic_analysis_opaque_ids"]) | set(
        lock["catalog_coverage"]["diagnostic_calibration_opaque_ids"]
    )
    assert covered == {row["opaque_id"] for row in public["circuits"]}
    power = lock["ideal_primary_label_power"]
    assert power["reference_component_is_unique_top"] is True
    assert power["pooled_sensitivity_expected_log_likelihood_ratio_lower_bound"] > math.log(
        analysis.POOLED_LR_THRESHOLD
    )
    assert all(
        details["sensitivity_expected_log_likelihood_ratio_lower_bound"]
        > math.log(analysis.PER_BACKEND_LR_THRESHOLD)
        for details in power["per_backend"].values()
    )
    analysis.validate_analysis_lock(lock)
