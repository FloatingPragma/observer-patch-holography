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
    } == {(1, 0, 1), (4, 1, 3), (6, 1, 4), (7, 2, 5)}
    assert all(
        item["minimum_nonfaulty_lock_transfer_to_view_change"] >= 1
        for item in row["finite_parameter_sweep"]
    )
    assert all(
        item["minimum_nonfaulty_decision_commit_lock_holders"]
        >= item["q"] - item["f"]
        for item in row["finite_parameter_sweep"]
    )
    assert row["checks"][
        "prepare_signers_provisional_acceptors_and_committers_are_separate"
    ]
    assert row["checks"][
        "valid_new_view_selection_supersedes_local_predecision_locks"
    ]
    assert row["checks"][
        "omitted_higher_orphan_lock_is_reconciled_in_general_case"
    ]
    for parameter in row["finite_parameter_sweep"]:
        example = parameter["orphan_lock_reconciliation_example"]
        assert example is not None
        assert audit.evaluate_orphan_lock_scenario(example["scenario"]) == (
            example["evaluation"]
        )
        assert example["evaluation"]["recovered"] is True
        assert example["evaluation"]["deadlocked"] is False
        assert example["evaluation"]["decision_certificates"] == []
    general_example = next(
        parameter["orphan_lock_reconciliation_example"]
        for parameter in row["finite_parameter_sweep"]
        if parameter["q"] < parameter["n"] - parameter["f"]
    )
    assert general_example["omitted_higher_lock_holders"]
    selected_view = general_example["evaluation"]["new_view"][
        "selected_report"
    ]["view"]
    assert all(
        general_example["evaluation"]["validator_states"][str(validator)][
            "view"
        ]
        > selected_view
        for validator in general_example["omitted_higher_lock_holders"]
    )


def test_bft_controls_cover_same_view_cross_view_and_liveness(receipt: dict) -> None:
    controls = receipt["negative_controls"]["bft"]
    assert len(controls) == 17
    assert all(row["counterexample_verified"] for row in controls)
    assert {
        row["removed_hypothesis"] for row in controls
    } >= {
        "P4_pre_lock_certificate_acknowledgements_are_provisional",
        "P4_valid_new_view_supersedes_local_predecision_commit_lock",
    }
    violated = {row["violated_conclusion"] for row in controls}
    assert "same_view_safety" in violated
    assert "cross_view_safety" in violated
    assert "liveness" in violated


def test_orphan_lock_traces_execute_every_transition_input(receipt: dict) -> None:
    controls = {
        row["removed_hypothesis"]: row
        for row in receipt["negative_controls"]["bft"]
    }
    pre = controls[
        "P4_pre_lock_certificate_acknowledgements_are_provisional"
    ]["finite_witness"]
    partial = controls[
        "P4_valid_new_view_supersedes_local_predecision_commit_lock"
    ]["finite_witness"]
    pre_scenario = pre["scenario"]
    partial_scenario = partial["scenario"]
    pre_evaluation = pre["evaluation"]
    partial_evaluation = partial["evaluation"]
    assert audit.evaluate_orphan_lock_scenario(pre_scenario) == pre_evaluation
    assert (
        audit.evaluate_orphan_lock_scenario(partial_scenario)
        == partial_evaluation
    )
    assert pre_evaluation["deadlocked"] is True
    assert pre_evaluation["maximum_eligible_votes"] == 2
    assert partial_evaluation["deadlocked"] is True
    assert partial_evaluation["new_view"]["selected_report"]["value"] == "B"
    assert partial_evaluation["maximum_eligible_votes"] == 2

    healed_pre = copy.deepcopy(pre_scenario)
    healed_pre["acknowledgements_are_provisional"] = True
    healed_pre_evaluation = audit.evaluate_orphan_lock_scenario(healed_pre)
    assert healed_pre_evaluation["recovered"] is True
    assert all(
        state["kind"] == "unlocked"
        for state in healed_pre_evaluation["validator_states"].values()
    )

    healed_partial = copy.deepcopy(partial_scenario)
    healed_partial["new_view_supersedes_lc_locks"] = True
    healed_partial_evaluation = audit.evaluate_orphan_lock_scenario(
        healed_partial
    )
    assert healed_partial_evaluation["recovered"] is True
    assert healed_partial_evaluation["voting_outcomes"][0][
        "eligible_voters"
    ] == [1, 2, 3]

    commit_ready = copy.deepcopy(partial_scenario)
    commit_ready["certificates"][1]["lock_certificate_recipients"] = [
        0,
        2,
        3,
    ]
    commit_ready_evaluation = audit.evaluate_orphan_lock_scenario(
        commit_ready
    )
    assert commit_ready_evaluation["decision_certificates"] == []

    delivered_to_lower_holder = copy.deepcopy(partial_scenario)
    delivered_to_lower_holder["certificates"][1][
        "lock_certificate_recipients"
    ] = [1, 2]
    assert audit.evaluate_orphan_lock_scenario(
        delivered_to_lower_holder
    )["recovered"] is True

    formed_decision = copy.deepcopy(commit_ready)
    formed_decision["certificates"][1]["committers"] = [0, 2, 3]
    assert audit.evaluate_orphan_lock_scenario(formed_decision)[
        "decision_certificates"
    ] == [{"view": 1, "value": "B", "committers": [0, 2, 3]}]

    def changed_or_rejected(
        scenario: dict,
        baseline: dict,
        path: tuple[object, ...],
        replacement: object,
    ) -> None:
        mutated = copy.deepcopy(scenario)
        cursor: object = mutated
        for key in path[:-1]:
            cursor = cursor[key]  # type: ignore[index]
        cursor[path[-1]] = replacement  # type: ignore[index]
        try:
            result = audit.evaluate_orphan_lock_scenario(mutated)
        except ValueError:
            return
        assert result != baseline

    mutations = [
        (
            "validators",
            pre_scenario,
            pre_evaluation,
            ("validators",),
            [0, 1, 2, 3, 4],
        ),
        (
            "byzantine",
            pre_scenario,
            pre_evaluation,
            ("byzantine",),
            [0, 3],
        ),
        ("q", pre_scenario, pre_evaluation, ("q",), 2),
        (
            "certificates",
            pre_scenario,
            pre_evaluation,
            ("certificates",),
            [copy.deepcopy(pre_scenario["certificates"][0])],
        ),
        (
            "new_view_senders",
            partial_scenario,
            partial_evaluation,
            ("new_view_senders",),
            [0, 1, 3],
        ),
        (
            "fresh_values",
            pre_scenario,
            pre_evaluation,
            ("fresh_values",),
            ["C"],
        ),
        (
            "acknowledgements_are_provisional",
            pre_scenario,
            pre_evaluation,
            ("acknowledgements_are_provisional",),
            True,
        ),
        (
            "new_view_supersedes_lc_locks",
            partial_scenario,
            partial_evaluation,
            ("new_view_supersedes_lc_locks",),
            True,
        ),
        (
            "certificate.view",
            partial_scenario,
            partial_evaluation,
            ("certificates", 1, "view"),
            2,
        ),
        (
            "certificate.value",
            partial_scenario,
            partial_evaluation,
            ("certificates", 1, "value"),
            "D",
        ),
        (
            "certificate.prepare_signers",
            partial_scenario,
            partial_evaluation,
            ("certificates", 1, "prepare_signers"),
            [0, 2, 2],
        ),
        (
            "certificate.acknowledgers",
            pre_scenario,
            pre_evaluation,
            ("certificates", 0, "acknowledgers"),
            [0],
        ),
        (
            "certificate.lock_certificate_recipients",
            partial_scenario,
            partial_evaluation,
            ("certificates", 1, "lock_certificate_recipients"),
            [1, 2],
        ),
        (
            "certificate.committers",
            commit_ready,
            commit_ready_evaluation,
            ("certificates", 1, "committers"),
            [0, 2, 3],
        ),
    ]
    assert {mutation[0] for mutation in mutations} == {
        "validators",
        "byzantine",
        "q",
        "certificates",
        "new_view_senders",
        "fresh_values",
        "acknowledgements_are_provisional",
        "new_view_supersedes_lc_locks",
        "certificate.view",
        "certificate.value",
        "certificate.prepare_signers",
        "certificate.acknowledgers",
        "certificate.lock_certificate_recipients",
        "certificate.committers",
    }
    for _, scenario, baseline, path, replacement in mutations:
        changed_or_rejected(scenario, baseline, path, replacement)


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

    lp_controls = {
        row["removed_hypothesis"]: row
        for row in receipt["ell_p_guard"]["negative_controls"]
    }
    assert set(lp_controls) == {
        "p>=1",
        "finite_or_pairwise_summable_declared_channel_family",
    }
    exponent = lp_controls["p>=1"]
    assert exponent["p"] == "1/2"
    assert exponent["d_p_xz"] == 4
    assert exponent["d_p_xy_plus_d_p_yz"] == 2
    assert exponent["counterexample_verified"] is True
    summability = lp_controls[
        "finite_or_pairwise_summable_declared_channel_family"
    ]
    assert summability["violated_conclusion"] == (
        "finite_valued_pseudometric_status"
    )
    assert summability["counterexample_verified"] is True


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
