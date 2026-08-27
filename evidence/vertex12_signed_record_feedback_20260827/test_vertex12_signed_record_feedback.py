from __future__ import annotations

import ast
import copy
import hashlib
import json
from pathlib import Path

import pytest

from oph_fpe.dynamics import vertex12_signed_record_feedback as producer
from oph_fpe.dynamics import (
    verify_vertex12_signed_record_feedback_independent as independent,
)


def _canonical_bytes(value: object) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")


def _sha(value: object) -> str:
    return "sha256:" + hashlib.sha256(_canonical_bytes(value)).hexdigest()


def _rehash(report: dict) -> None:
    payload = copy.deepcopy(report)
    payload.pop("receipt_sha256", None)
    report["receipt_sha256"] = _sha(payload)


@pytest.fixture(scope="module")
def receipt() -> dict:
    return producer.produce_receipt()


def _assert_independent_rejects(receipt: dict) -> None:
    _rehash(receipt)
    with pytest.raises(independent.IndependentVerificationError):
        independent.verify_report(receipt)


def test_producer_and_independent_replay_pass(receipt: dict) -> None:
    producer_result = producer.verify_receipt(receipt)
    assert producer_result["receipt"] is True
    assert producer_result["status"] == "PASS"

    independent_result = independent.verify_report(receipt)
    assert independent_result["receipt"] is True
    assert independent_result["status"] == "PASS"
    assert independent_result["producer_imported"] is False
    assert independent_result["literal_integer_transition_independently_replayed"] is True
    assert independent_result["checked_record_count"] == 8
    assert independent_result["checked_feedback_event_count"] == 96
    assert independent_result["checked_carrier_port_pair_count"] == 96
    assert independent_result["checked_A5_covariance_count"] == 720
    assert independent_result["checked_feedback_idempotence_count"] == 96
    assert independent_result["checked_feedback_commutation_count"] == 528
    assert independent_result["parent_endpoint_repair_confluence_verified"] is False
    assert independent_result["physical_attachment_verified"] is False


def test_literal_record_read_causes_a_later_bounded_restore(receipt: dict) -> None:
    records = receipt["literal_record_commit_log"]["records"]
    events = receipt["causal_feedback_log"]["events"]
    records_by_carrier = {row["carrier_index"]: row for row in records}
    assert len(records_by_carrier) == 8
    assert len(events) == 96

    for event in events:
        record = records_by_carrier[event["carrier_index"]]
        port = event["port"]
        assert event["record_id"] == record["record_id"]
        assert event["commit_event_index"] < event["probe_event_index"]
        assert event["probe_event_index"] < event["read_event_index"]
        assert event["read_event_index"] < event["write_event_index"]
        assert event["literal_committed_record_port_value"] == record[
            "signed_port_record"
        ][port]
        assert event["live_state_before_feedback"][port] == (
            event["literal_committed_record_port_value"] + 1
        )
        assert event["feedback_delta"] == -1
        assert event["live_state_after_feedback"] == event[
            "live_state_before_probe"
        ]
        assert event["protected_record_after"] == event["protected_record_before"]
        assert event["bounded_local_port_write"] is True
        assert event["exact_preprobe_state_restored"] is True
        assert event["protected_record_unchanged"] is True
        assert event["action_rule_inputs"] == [
            "live_port_value_after_probe",
            "literal_committed_record_port_value",
        ]
        assert event["hash_value_consumed_by_transition_rule"] is False


def test_ablation_and_record_coordinate_counterfactual_are_transparent(
    receipt: dict,
) -> None:
    for event in receipt["causal_feedback_log"]["events"]:
        port = event["port"]
        ablation = event["ablation"]
        assert ablation["feedback_delta"] == 0
        assert ablation["live_state_before_action"] == event[
            "live_state_before_feedback"
        ]
        assert ablation["live_state_after_action"] == event[
            "live_state_before_feedback"
        ]
        assert ablation["live_state_after_action"] != event[
            "live_state_after_feedback"
        ]
        assert ablation["differs_from_actual_later_state"] is True
        assert ablation["preprobe_state_not_restored"] is True

        counterfactual = event["record_coordinate_counterfactual"]
        before = counterfactual["record_before"]
        after = counterfactual["record_after_intervention"]
        changed = [
            index
            for index, pair in enumerate(zip(before, after, strict=True))
            if pair[0] != pair[1]
        ]
        assert changed == [port]
        assert after[port] == before[port] + 1
        assert counterfactual["changed_record_coordinate_count"] == 1
        assert counterfactual["live_state_before_action"] == event[
            "live_state_before_feedback"
        ]
        assert counterfactual["probe_delta"] == event["probe_delta"] == 1
        assert counterfactual["feedback_delta"] == 0
        assert counterfactual["feedback_delta"] != event["feedback_delta"]
        assert counterfactual["nonrecord_inputs_held_fixed"] is True
        assert counterfactual["later_action_differs"] is True


def test_all_twelve_ports_and_full_A5_action_are_covered(receipt: dict) -> None:
    events = receipt["causal_feedback_log"]["events"]
    assert {(row["carrier_index"], row["port"]) for row in events} == {
        (carrier, port) for carrier in range(8) for port in range(12)
    }

    audit = receipt["a5_covariance_audit"]
    rows = audit["rows"]
    assert audit["group_order"] == 60
    assert audit["port_count"] == 12
    assert audit["check_count"] == 720
    assert {(row["group_element_index"], row["source_port"]) for row in rows} == {
        (group, port) for group in range(60) for port in range(12)
    }
    assert all(
        row["literal_coordinate_feedback_commutes_with_port_action"] is True
        and row["source_feedback_delta"]
        == row["target_feedback_delta"]
        == -1
        and row["source_counterfactual_feedback_delta"]
        == row["target_counterfactual_feedback_delta"]
        == 0
        for row in rows
    )


def test_feedback_transactions_have_a_standalone_local_normal_form(
    receipt: dict,
) -> None:
    audit = receipt["local_feedback_confluence_audit"]
    assert audit["scope"] == "standalone_literal_record_reset_transactions_only"
    assert audit["idempotence_check_count"] == 96
    assert len(audit["idempotence_rows"]) == 96
    assert all(
        row["transaction_idempotent"] is True
        and row["other_coordinates_unchanged"] is True
        for row in audit["idempotence_rows"]
    )
    assert audit["disjoint_commutation_check_count"] == 528
    assert len(audit["disjoint_commutation_rows"]) == 528
    assert all(
        row["disjoint_transactions_commute"] is True
        and row["all_coordinates_equal"] is True
        for row in audit["disjoint_commutation_rows"]
    )
    assert audit["all_port_reset_transactions_idempotent"] is True
    assert audit["all_disjoint_port_reset_transactions_commute"] is True
    assert (
        audit[
            "any_serialized_control_schedule_applying_every_port_reaches_literal_record"
        ]
        is True
    )
    assert audit["parent_endpoint_repair_confluence_established"] is False


def test_receipt_keeps_internal_mechanics_separate_from_physics(receipt: dict) -> None:
    attainment = receipt["attainment"]
    assert attainment["INTERNAL_FINITE_OBSERVER_LIKE_SELF_READING_RECEIPT"] is True
    for key in (
        "CANONICAL_A1_A2_A3_SOURCE_SELECTION_RECEIPT",
        "SOURCE_QUALIFIED_PHYSICAL_OBSERVER_RECEIPT",
        "SPATIAL_TRANSLATION_RECEIPT",
        "LABORATORY_RECORD_REALIZATION_RECEIPT",
        "PHYSICAL_PREDICTION_RECEIPT",
    ):
        assert attainment[key] is False
    assert receipt["comparison_data_read"] is False
    assert "no laboratory" in receipt["claim_boundary"]


def test_independent_verifier_does_not_import_the_producer() -> None:
    tree = ast.parse(Path(independent.__file__).read_text(encoding="utf-8"))
    imported: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.append(node.module)
    assert "oph_fpe.dynamics.vertex12_signed_record_feedback" not in imported


def test_integer_boolean_type_confusion_with_stale_inner_digest_is_rejected(
    receipt: dict,
) -> None:
    changed = copy.deepcopy(receipt)
    changed["causal_feedback_log"]["events"][0]["probe_delta"] = True
    # Keep the old event-log digest: Python container equality would otherwise
    # accept this substitution because ``True == 1``.
    _rehash(changed)
    with pytest.raises(independent.IndependentVerificationError):
        independent.verify_report(changed)


@pytest.mark.parametrize(
    "mutation",
    (
        "literal_record_value",
        "feedback_action",
        "hash_as_action_input",
        "read_before_commit",
        "ablation_false_positive",
        "counterfactual_changes_two_record_coordinates",
        "missing_carrier_port_event",
        "A5_target_port",
        "local_commutation",
        "physical_promotion",
        "parent_pin",
    ),
)
def test_rehashed_semantic_tampering_is_rejected(
    receipt: dict, mutation: str
) -> None:
    changed = copy.deepcopy(receipt)
    event = changed["causal_feedback_log"]["events"][0]

    if mutation == "literal_record_value":
        changed["literal_record_commit_log"]["records"][0][
            "signed_port_record"
        ][0] += 1
    elif mutation == "feedback_action":
        event["feedback_delta"] = 0
    elif mutation == "hash_as_action_input":
        event["action_rule_inputs"].append("record_id")
        event["hash_value_consumed_by_transition_rule"] = True
    elif mutation == "read_before_commit":
        event["read_event_index"] = event["commit_event_index"]
    elif mutation == "ablation_false_positive":
        event["ablation"]["live_state_after_action"] = event[
            "live_state_after_feedback"
        ]
    elif mutation == "counterfactual_changes_two_record_coordinates":
        event["record_coordinate_counterfactual"][
            "record_after_intervention"
        ][1] += 1
        event["record_coordinate_counterfactual"][
            "changed_record_coordinate_count"
        ] = 2
    elif mutation == "missing_carrier_port_event":
        changed["causal_feedback_log"]["events"].pop()
        changed["causal_feedback_log"]["event_count"] -= 1
        changed["causal_feedback_log"]["covered_port_coordinate_count"] -= 1
    elif mutation == "A5_target_port":
        changed["a5_covariance_audit"]["rows"][0]["target_port"] = 1
    elif mutation == "local_commutation":
        changed["local_feedback_confluence_audit"][
            "disjoint_commutation_rows"
        ][0]["disjoint_transactions_commute"] = False
    elif mutation == "physical_promotion":
        changed["attainment"]["LABORATORY_RECORD_REALIZATION_RECEIPT"] = True
    elif mutation == "parent_pin":
        changed["parent_pins"]["constructive_source_law"]["receipt_sha256"] = (
            "sha256:" + "0" * 64
        )
    else:  # pragma: no cover - the parameter list is closed above
        raise AssertionError(mutation)

    # Rehash nested custody fields where doing so makes the forgery stronger.
    changed["literal_record_commit_log"]["records_sha256"] = _sha(
        changed["literal_record_commit_log"]["records"]
    )
    changed["causal_feedback_log"]["events_sha256"] = _sha(
        changed["causal_feedback_log"]["events"]
    )
    changed["a5_covariance_audit"]["rows_sha256"] = _sha(
        changed["a5_covariance_audit"]["rows"]
    )
    changed["local_feedback_confluence_audit"][
        "idempotence_rows_sha256"
    ] = _sha(changed["local_feedback_confluence_audit"]["idempotence_rows"])
    changed["local_feedback_confluence_audit"][
        "disjoint_commutation_rows_sha256"
    ] = _sha(
        changed["local_feedback_confluence_audit"][
            "disjoint_commutation_rows"
        ]
    )
    _assert_independent_rejects(changed)
