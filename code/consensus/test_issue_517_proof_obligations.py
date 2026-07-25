"""Regression and adversarial tests for the issue #517 proof audit."""

from __future__ import annotations

import copy
import importlib.util
import json
from pathlib import Path
import sys

import pytest


MODULE_DIR = Path(__file__).resolve().parent
ROOT = MODULE_DIR.parents[1]
sys.path.insert(0, str(MODULE_DIR))

import verify_issue_517_proof_obligations as audit  # noqa: E402


RECEIPT_PATH = MODULE_DIR / "runs/issue_517_proof_obligations.json"


@pytest.fixture(scope="module")
def receipt() -> dict:
    return json.loads(RECEIPT_PATH.read_text(encoding="utf-8"))


def test_stored_receipt_is_exactly_recomputable(receipt: dict) -> None:
    audit.verify_receipt(receipt)


def test_transactional_receipt_checks_every_finite_peak(receipt: dict) -> None:
    row = receipt["transactional_confluence"]
    assert row["receipt_id"] == "TXN-DIAMOND-1"
    assert row["finite_state_space"]["cardinality"] == 16
    assert row["checked_peak_count"] > 0
    assert all(row["checks"].values())
    assert all(peak["joined_in_one_step"] for peak in row["peak_receipts"])


def test_every_load_bearing_confluence_hypothesis_has_a_control(receipt: dict) -> None:
    controls = {
        row["removed_hypothesis"]: row
        for row in receipt["negative_controls"]["transactional_confluence"]
    }
    assert set(controls) == {
        "semantic_dependency_complete_reads_and_revalidation",
        "atomic_conflict_component_commit",
        "coherent_canonical_union_collar_payload",
        "aggregate_support_closure_and_fixed_point_remerge",
        "prepared_source_component_stability_or_frozen_batch_admissibility",
        "payload_determined_by_declared_read_snapshot",
        "strict_well_founded_descent",
        "repair_completeness",
        "boundary_sector_holonomy_preservation",
    }
    assert all(row["counterexample_verified"] for row in controls.values())


def test_prepared_lock_transfers_to_every_new_view_quorum(receipt: dict) -> None:
    row = receipt["prepared_certificate_bft"]
    assert row["receipt_id"] == "BFT-LOCK-1"
    assert all(row["checks"].values())
    assert {
        (item["n"], item["f"], item["q"])
        for item in row["finite_parameter_sweep"]
    } == {(1, 0, 1), (4, 1, 3), (7, 2, 5)}
    assert all(
        item["minimum_nonfaulty_lock_transfer_to_view_change"] >= 1
        for item in row["finite_parameter_sweep"]
    )


def test_bft_controls_cover_same_view_cross_view_and_liveness(receipt: dict) -> None:
    controls = receipt["negative_controls"]["bft"]
    assert len(controls) == 15
    assert all(row["counterexample_verified"] for row in controls)
    violated = {row["violated_conclusion"] for row in controls}
    assert "same_view_safety" in violated
    assert "cross_view_safety" in violated
    assert "liveness" in violated


def test_refinement_and_lp_controls_are_exact(receipt: dict) -> None:
    refinement = receipt["refinement_moduli"]
    assert "for every positive integer" in refinement[
        "uniform_inverse_counterfamily"
    ]["witness_constructor"]
    assert "for every positive integer" in refinement[
        "uniform_residual_counterfamily"
    ]["witness_constructor"]
    assert all(refinement["checks"].values())
    assert len(receipt["negative_controls"]["refinement"]) == 8

    lp = receipt["ell_p_guard"]["negative_control"]
    assert lp["p"] == "1/2"
    assert lp["d_p_xz"] == 4
    assert lp["d_p_xy_plus_d_p_yz"] == 2
    assert lp["counterexample_verified"] is True


def test_selector_is_closed_schema_and_legacy_sieve_is_recomputed(receipt: dict) -> None:
    selector = receipt["selector"]
    assert selector["finite_negative_control_count"] == 19
    assert all(selector["checks"].values())
    assert selector["separate_conditional_variational_sieve"]["status"] == (
        "conditional_finite_selector_theorem"
    )

    validator_path = (
        ROOT
        / "code/particles/hierarchy/validators/validate_screen_sieve_certificate.py"
    )
    spec = importlib.util.spec_from_file_location(
        "test_issue_517_screen_sieve_validator",
        validator_path,
    )
    assert spec is not None and spec.loader is not None
    validator = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = validator
    spec.loader.exec_module(validator)

    legacy_path = (
        ROOT
        / "code/particles/hierarchy/certificates/R_screen_sieve_icosahedral_certificate.json"
    )
    legacy = json.loads(legacy_path.read_text(encoding="utf-8"))
    assert all(validator.validate(legacy).values())

    tampered = copy.deepcopy(legacy)
    tampered["checks"]["seidel_class_is_unique"] = False
    checks = validator.validate(tampered)
    assert checks["certificate_exactly_matches_recomputation"] is False
    assert checks["embedded_check_booleans_match_recomputation"] is False


def test_a5_layers_remain_distinct_and_unpromoted(receipt: dict) -> None:
    row = receipt["a5_layer_separation"]
    assert all(row["checks"].values())
    layers = row["layers"]
    assert len({item["receipt_id"] for item in layers.values()}) == 4
    assert all(item["promoted"] is False for item in layers.values())
    subreceipts = row["conditional_matter_subreceipts"]
    assert set(subreceipts) == {
        "spin_lift",
        "refinement_stable_anomaly_algebra",
    }
    assert all(item["exact_checks_pass"] for item in subreceipts.values())
    assert not any(
        item["physical_source_promoted"] for item in subreceipts.values()
    )
