"""Independent verifier for the literal vertex12 signed-record feedback.

The producer is intentionally not imported.  This module reconstructs the
bounded integer records from the signed source-event directions, replays every
probe/read/write, checks the ablation and one-coordinate intervention, and
checks all 720 A5 port-action squares directly from the pinned parent data.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

from oph_fpe.core.icosahedral import icosahedral_a5_port_permutations


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RECEIPT = (
    ROOT / "data/repair_closure/vertex12_signed_record_feedback_receipt.json"
)
ATOMIC_PARENT = (
    ROOT / "data/repair_closure/vertex12_atomic_port_transfer_receipt.json"
)
SOURCE_LAW_PARENT = (
    ROOT / "data/repair_closure/vertex12_constructive_source_law_receipt.json"
)
PRODUCER_PATH = ROOT / "oph_fpe/dynamics/vertex12_signed_record_feedback.py"
VERIFIER_PATH = Path(__file__).resolve()
TEST_PATH = ROOT / "tests/test_vertex12_signed_record_feedback.py"

SCHEMA = "oph.vertex12-signed-record-feedback-diagnostic.v1"
STATUS = (
    "INTERNAL_BOUNDED_LITERAL_SIGNED_RECORD_CAUSAL_FEEDBACK_ATTAINED__"
    "CANONICAL_SOURCE_SELECTION_AND_PHYSICAL_ATTACHMENT_OPEN"
)
ATOMIC_SCHEMA = "oph.vertex12-atomic-port-transfer-subpacket.v1"
ATOMIC_STATUS = (
    "INTERNAL_VERTEX12_SEAM_MATCHING_PROJECTORS_AND_IN_PROCESS_SNAPSHOT_REREAD_ATTAINED__"
    "SPATIAL_PHYSICAL_BRIDGE_OPEN"
)
SOURCE_LAW_SCHEMA = "oph.vertex12-constructive-source-law-control.v1"
SOURCE_LAW_STATUS = (
    "CONSTRUCTIVE_SOURCE_LAW_CONTROL_ATTAINED__"
    "CANONICAL_A1_A2_A3_SELECTION_AND_PHYSICAL_ATTACHMENT_OPEN"
)
PORTS = 12
AXES = 6
BOUND = 32
PROBE = 1


class IndependentVerificationError(RuntimeError):
    """Raised when serialized evidence fails an independent exact check."""


def _canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")


def _sha(value: Any) -> str:
    return "sha256:" + hashlib.sha256(_canonical_bytes(value)).hexdigest()


def _load(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise IndependentVerificationError(f"unreadable JSON: {path}") from exc
    if not isinstance(value, dict):
        raise IndependentVerificationError(f"JSON root is not an object: {path}")
    return value


def _fail(condition: bool, message: str) -> None:
    if not condition:
        raise IndependentVerificationError(message)


def _fail_unless_json_exact(actual: Any, expected: Any, message: str) -> None:
    """Compare serialized values without Python's ``True == 1`` coercion."""

    _fail(_canonical_bytes(actual) == _canonical_bytes(expected), message)


def _exact_keys(value: Any, expected: set[str], label: str) -> None:
    _fail(isinstance(value, Mapping), f"{label} is not an object")
    _fail(set(value) == expected, f"{label} key set")


def _raw_pin(path: Path) -> dict[str, Any]:
    raw = path.read_bytes()
    return {
        "path": path.relative_to(ROOT).as_posix(),
        "bytes": len(raw),
        "sha256": "sha256:" + hashlib.sha256(raw).hexdigest(),
    }


def _validate_receipt_digest(report: Mapping[str, Any], label: str) -> None:
    material = dict(report)
    digest = material.pop("receipt_sha256", None)
    _fail(isinstance(digest, str) and digest == _sha(material), f"{label} digest")


def _validate_parent_pin(
    pin: Any,
    *,
    path: Path,
    parent: Mapping[str, Any],
    schema: str,
    label: str,
) -> None:
    _exact_keys(pin, {"schema", "status", "receipt_sha256", "raw_pin"}, label)
    _fail(pin["schema"] == parent.get("schema") == schema, f"{label} schema")
    _fail(pin["status"] == parent.get("status"), f"{label} status")
    _fail(
        pin["receipt_sha256"] == parent.get("receipt_sha256"),
        f"{label} receipt",
    )
    _fail_unless_json_exact(pin["raw_pin"], _raw_pin(path), f"{label} raw pin")


def _parent_contract(
    report: Mapping[str, Any],
) -> tuple[dict[str, Any], dict[str, Any]]:
    atomic = _load(ATOMIC_PARENT)
    source = _load(SOURCE_LAW_PARENT)
    _fail(atomic.get("schema") == ATOMIC_SCHEMA, "atomic parent schema")
    _fail(atomic.get("status") == ATOMIC_STATUS, "atomic parent status")
    _fail(atomic.get("issue") == 655, "atomic parent issue")
    _validate_receipt_digest(atomic, "atomic parent")
    operator = atomic.get("atomic_transfer_operator")
    readback = atomic.get("post_repair_in_process_snapshot_reread")
    boundary = atomic.get("quotient_and_spatial_boundary")
    _fail(isinstance(operator, Mapping), "atomic operator missing")
    _fail(operator.get("carrier_count") == 8, "atomic carrier count")
    _fail(operator.get("port_count") == PORTS, "atomic port count")
    _fail(
        operator.get("source_native_internal_seam_partner_operator_receipt")
        is True,
        "atomic internal transfer",
    )
    _fail(
        operator.get("source_native_spatial_translation_receipt") is False,
        "atomic spatial promotion",
    )
    _fail(isinstance(readback, Mapping), "atomic readback missing")
    _fail(
        readback.get("covered_port_coordinate_count") == 8 * PORTS,
        "atomic readback coverage",
    )
    _fail(
        readback.get("every_carrier_full_port_state_committed") is True
        and readback.get("every_carrier_full_port_state_reread_in_process") is True,
        "atomic readback contract",
    )
    _fail(
        isinstance(boundary, Mapping)
        and boundary.get("internal_seam_transfer_is_spatial_translation") is False,
        "atomic spatial boundary",
    )

    _fail(source.get("schema") == SOURCE_LAW_SCHEMA, "source parent schema")
    _fail(source.get("status") == SOURCE_LAW_STATUS, "source parent status")
    _fail(source.get("issue") == 655, "source parent issue")
    _validate_receipt_digest(source, "source parent")
    source_attainment = source.get("attainment")
    _fail(isinstance(source_attainment, Mapping), "source attainment missing")
    for key in (
        "constructive_source_law_capture_root",
        "twelve_raw_steps",
        "same_Q_A5_covariance",
    ):
        _fail(source_attainment.get(key) is True, f"source attainment {key}")
    for key in (
        "canonical_source_selection",
        "spatial_translation",
        "physical_readout",
        "physical_prediction",
    ):
        _fail(source_attainment.get(key) is False, f"source boundary {key}")

    pins = report.get("parent_pins")
    _exact_keys(
        pins,
        {"atomic_port_transfer", "constructive_source_law"},
        "parent pins",
    )
    _validate_parent_pin(
        pins["atomic_port_transfer"],
        path=ATOMIC_PARENT,
        parent=atomic,
        schema=ATOMIC_SCHEMA,
        label="atomic parent pin",
    )
    _validate_parent_pin(
        pins["constructive_source_law"],
        path=SOURCE_LAW_PARENT,
        parent=source,
        schema=SOURCE_LAW_SCHEMA,
        label="source-law parent pin",
    )
    return atomic, source


def _source_geometry(
    source: Mapping[str, Any],
) -> tuple[
    tuple[tuple[int, ...], ...],
    tuple[tuple[int, ...], ...],
    tuple[str, ...],
]:
    law = source["constructive_source_law"]
    alphabet = law["a1_complete_event_alphabet"]
    events = alphabet["event_rows"]
    _fail(alphabet["event_count"] == len(events) == PORTS, "source event count")
    _fail(alphabet["complete_signed_port_orbit"] is True, "signed orbit")
    directions: list[tuple[int, ...]] = []
    event_ids: list[str] = []
    for port, row in enumerate(events):
        direction = row.get("raw_direction_in_Z_power_6")
        _fail(row.get("port") == port, f"source port {port}")
        _fail(
            isinstance(direction, list)
            and len(direction) == AXES
            and all(type(value) is int and value in (-1, 0, 1) for value in direction)
            and sum(abs(value) for value in direction) == 1,
            f"source direction {port}",
        )
        payload = dict(row)
        event_id = payload.pop("event_id", None)
        _fail(event_id == _sha(payload), f"source event digest {port}")
        directions.append(tuple(direction))
        event_ids.append(str(event_id))

    canonical_actions = tuple(
        tuple(int(value) for value in row)
        for row in icosahedral_a5_port_permutations()
    )
    group_rows = law["same_Q_A5_action"]["group_rows"]
    _fail(len(group_rows) == len(canonical_actions) == 60, "A5 row count")
    actions: list[tuple[int, ...]] = []
    for index, (canonical, row) in enumerate(
        zip(canonical_actions, group_rows, strict=True)
    ):
        action = tuple(int(value) for value in row["port_permutation"])
        _fail(row["group_element_index"] == index, f"A5 index {index}")
        _fail(action == canonical, f"A5 action {index}")
        actions.append(action)
    _fail(len(set(actions)) == 60, "A5 action faithfulness")
    return tuple(directions), tuple(actions), tuple(event_ids)


def _expected_record(
    carrier_id: str,
    carrier_index: int,
    directions: Sequence[Sequence[int]],
    source_event_ids: Sequence[str],
) -> dict[str, Any]:
    axes = [carrier_index + axis + 1 for axis in range(AXES)]
    ports = [
        sum(int(row[axis]) * axes[axis] for axis in range(AXES))
        for row in directions
    ]
    material = {
        "schema": "oph.vertex12-literal-signed-port-record.v1",
        "carrier_id": carrier_id,
        "carrier_index": carrier_index,
        "commit_event_index": carrier_index,
        "axis_record_in_Z_power_6": axes,
        "signed_port_record": ports,
        "source_event_ids_by_port": list(source_event_ids),
        "record_value_type": "literal_bounded_integer",
        "state_bound": BOUND,
    }
    return {**material, "record_id": _sha(material)}


def _expected_event(
    record: Mapping[str, Any],
    *,
    port: int,
    ordinal: int,
    carrier_count: int,
) -> dict[str, Any]:
    protected = [int(value) for value in record["signed_port_record"]]
    before = list(protected)
    probed = list(before)
    probed[port] += PROBE
    actual_delta = protected[port] - probed[port]
    after = list(probed)
    after[port] += actual_delta
    counterfactual_record = list(protected)
    counterfactual_record[port] += PROBE
    counterfactual_delta = counterfactual_record[port] - probed[port]
    counterfactual_after = list(probed)
    counterfactual_after[port] += counterfactual_delta
    event_base = carrier_count + 3 * ordinal
    material = {
        "schema": "oph.vertex12-literal-signed-record-feedback-event.v1",
        "carrier_id": record["carrier_id"],
        "carrier_index": record["carrier_index"],
        "port": port,
        "record_id": record["record_id"],
        "commit_event_index": record["commit_event_index"],
        "probe_event_index": event_base,
        "read_event_index": event_base + 1,
        "write_event_index": event_base + 2,
        "protected_record_before": protected,
        "protected_record_after": list(protected),
        "live_state_before_probe": before,
        "probe_delta": PROBE,
        "live_state_before_feedback": probed,
        "literal_committed_record_port_value": protected[port],
        "feedback_delta": actual_delta,
        "live_state_after_feedback": after,
        "action_rule": "delta=literal_committed_record_port_value-live_port_value",
        "action_rule_inputs": [
            "live_port_value_after_probe",
            "literal_committed_record_port_value",
        ],
        "hash_value_consumed_by_transition_rule": False,
        "bounded_integer_state_minimum": -BOUND,
        "bounded_integer_state_maximum": BOUND,
        "bounded_local_port_write": True,
        "exact_preprobe_state_restored": True,
        "protected_record_unchanged": True,
        "ablation": {
            "mode": "omit_record_read_and_feedback_write",
            "live_state_before_action": probed,
            "feedback_delta": 0,
            "live_state_after_action": list(probed),
            "differs_from_actual_later_state": True,
            "preprobe_state_not_restored": True,
        },
        "record_coordinate_counterfactual": {
            "mode": "do_literal_committed_record_coordinate",
            "intervened_port": port,
            "record_before": protected,
            "record_after_intervention": counterfactual_record,
            "changed_record_coordinate_count": 1,
            "live_state_before_action": probed,
            "probe_delta": PROBE,
            "feedback_delta": counterfactual_delta,
            "live_state_after_action": counterfactual_after,
            "nonrecord_inputs_held_fixed": True,
            "later_action_differs": True,
        },
    }
    return {**material, "event_id": _sha(material)}


def _expected_covariance_rows(
    record: Mapping[str, Any], actions: Sequence[Sequence[int]]
) -> list[dict[str, Any]]:
    reference = [int(value) for value in record["signed_port_record"]]
    rows: list[dict[str, Any]] = []
    for group_index, action in enumerate(actions):
        transformed = [0] * PORTS
        for source, target in enumerate(action):
            transformed[int(target)] = reference[source]
        for source, target in enumerate(action):
            source_value = reference[source]
            target_value = transformed[int(target)]
            rows.append(
                {
                    "group_element_index": group_index,
                    "source_port": source,
                    "target_port": int(target),
                    "source_literal_record_value": source_value,
                    "transformed_literal_record_value": target_value,
                    "source_probe_delta": PROBE,
                    "target_probe_delta": PROBE,
                    "source_feedback_delta": source_value - (source_value + PROBE),
                    "target_feedback_delta": target_value - (target_value + PROBE),
                    "source_counterfactual_feedback_delta": 0,
                    "target_counterfactual_feedback_delta": 0,
                    "literal_coordinate_feedback_commutes_with_port_action": True,
                }
            )
    return rows


def _apply_reset_transaction(
    live_state: Sequence[int], record: Sequence[int], port: int
) -> list[int]:
    after_probe = [int(value) for value in live_state]
    after_probe[port] += PROBE
    feedback_delta = int(record[port]) - after_probe[port]
    after_feedback = list(after_probe)
    after_feedback[port] += feedback_delta
    return after_feedback


def _expected_local_confluence_rows(
    records: Sequence[Mapping[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    idempotence_rows: list[dict[str, Any]] = []
    commutation_rows: list[dict[str, Any]] = []
    for record in records:
        carrier_index = int(record["carrier_index"])
        signed_record = [int(value) for value in record["signed_port_record"]]
        control = [
            value + ((carrier_index + port) % 3) - 1
            for port, value in enumerate(signed_record)
        ]
        _fail(
            all(abs(value) + PROBE <= BOUND for value in control),
            "confluence control bound",
        )
        for port in range(PORTS):
            once = _apply_reset_transaction(control, signed_record, port)
            twice = _apply_reset_transaction(once, signed_record, port)
            idempotence_rows.append(
                {
                    "carrier_index": carrier_index,
                    "port": port,
                    "control_port_value": control[port],
                    "literal_record_port_value": signed_record[port],
                    "after_once_port_value": once[port],
                    "after_twice_port_value": twice[port],
                    "other_coordinates_unchanged": all(
                        once[other] == control[other]
                        for other in range(PORTS)
                        if other != port
                    ),
                    "transaction_idempotent": once == twice,
                }
            )
        for left in range(PORTS):
            for right in range(left + 1, PORTS):
                left_then_right = _apply_reset_transaction(
                    _apply_reset_transaction(control, signed_record, left),
                    signed_record,
                    right,
                )
                right_then_left = _apply_reset_transaction(
                    _apply_reset_transaction(control, signed_record, right),
                    signed_record,
                    left,
                )
                commutation_rows.append(
                    {
                        "carrier_index": carrier_index,
                        "left_port": left,
                        "right_port": right,
                        "control_values_at_pair": [control[left], control[right]],
                        "literal_record_values_at_pair": [
                            signed_record[left],
                            signed_record[right],
                        ],
                        "left_then_right_values_at_pair": [
                            left_then_right[left],
                            left_then_right[right],
                        ],
                        "right_then_left_values_at_pair": [
                            right_then_left[left],
                            right_then_left[right],
                        ],
                        "all_coordinates_equal": (
                            left_then_right == right_then_left
                        ),
                        "disjoint_transactions_commute": (
                            left_then_right == right_then_left
                        ),
                    }
                )
    return idempotence_rows, commutation_rows


def verify_report(report: Mapping[str, Any]) -> dict[str, Any]:
    received = json.loads(_canonical_bytes(dict(report)).decode("utf-8"))
    digest = received.pop("receipt_sha256", None)
    _fail(digest == _sha(received), "receipt digest")
    received["receipt_sha256"] = digest
    _exact_keys(
        received,
        {
            "schema",
            "issue",
            "status",
            "comparison_data_read",
            "parent_pins",
            "finite_contract",
            "literal_record_commit_log",
            "causal_feedback_log",
            "a5_covariance_audit",
            "local_feedback_confluence_audit",
            "attainment",
            "implementation_pins",
            "claim_boundary",
            "receipt_sha256",
        },
        "top-level report",
    )
    _fail(received["schema"] == SCHEMA, "schema")
    _fail(received["status"] == STATUS, "status")
    _fail(received["issue"] == 655, "issue")
    _fail(received["comparison_data_read"] is False, "comparison boundary")
    atomic, source = _parent_contract(received)
    directions, actions, event_ids = _source_geometry(source)
    carrier_ids = [
        str(value) for value in atomic["atomic_transfer_operator"]["carrier_ids"]
    ]
    _fail(len(carrier_ids) == len(set(carrier_ids)) == 8, "carrier census")

    contract = received["finite_contract"]
    _fail_unless_json_exact(
        contract,
        {
            "carrier_count": 8,
            "port_count": PORTS,
            "signed_axis_count": AXES,
            "serialized_working_state_envelope": (
                "integer_coordinates_with_absolute_value_at_most_32"
            ),
            "protected_record_domain": "integer_box_[-32,32]^(8x12)",
            "arbitrary_full_box_transition_closure_claimed": False,
            "probe_event": "append_literal_+1_to_one_working_port_coordinate",
            "feedback_event": "append_literal_record_minus_live_integer_delta",
            "feedback_rule": "z_i,p <- z_i,p + (b_i,p-z_i,p)",
            "record_is_full_literal_twelve_port_integer_row": True,
            "record_hash_used_only_for_identity_and_custody": True,
            "hash_value_consumed_by_transition_rule": False,
            "protected_record_is_not_mutated_by_feedback": True,
            "probe_and_feedback_are_sequential_per_coordinate": True,
        },
        "finite contract",
    )

    commit_log = received["literal_record_commit_log"]
    _exact_keys(
        commit_log, {"record_count", "records", "records_sha256"}, "commit log"
    )
    expected_records = [
        _expected_record(carrier_id, index, directions, event_ids)
        for index, carrier_id in enumerate(carrier_ids)
    ]
    _fail(commit_log["record_count"] == len(expected_records) == 8, "record count")
    _fail_unless_json_exact(
        commit_log["records"], expected_records, "literal record rows"
    )
    _fail(
        commit_log["records_sha256"]
        == _sha(commit_log["records"])
        == _sha(expected_records),
        "record log hash",
    )
    _fail(
        all(
            all(type(value) is int and abs(value) <= BOUND for value in row["signed_port_record"])
            for row in expected_records
        ),
        "record bounds",
    )

    feedback_log = received["causal_feedback_log"]
    _exact_keys(
        feedback_log,
        {
            "event_count",
            "events",
            "events_sha256",
            "covered_carrier_count",
            "covered_port_coordinate_count",
            "all_carrier_port_pairs_covered_exactly_once",
        },
        "feedback log",
    )
    expected_events = [
        _expected_event(
            record,
            port=port,
            ordinal=carrier_index * PORTS + port,
            carrier_count=len(expected_records),
        )
        for carrier_index, record in enumerate(expected_records)
        for port in range(PORTS)
    ]
    _fail(feedback_log["event_count"] == len(expected_events) == 96, "event count")
    _fail_unless_json_exact(
        feedback_log["events"], expected_events, "causal feedback rows"
    )
    _fail(
        feedback_log["events_sha256"]
        == _sha(feedback_log["events"])
        == _sha(expected_events),
        "event log hash",
    )
    _fail(feedback_log["covered_carrier_count"] == 8, "carrier coverage")
    _fail(feedback_log["covered_port_coordinate_count"] == 96, "port coverage")
    _fail(
        feedback_log["all_carrier_port_pairs_covered_exactly_once"] is True,
        "pair coverage verdict",
    )
    pair_set = {
        (row["carrier_index"], row["port"]) for row in feedback_log["events"]
    }
    _fail(pair_set == {(i, p) for i in range(8) for p in range(PORTS)}, "pair set")

    covariance = received["a5_covariance_audit"]
    _exact_keys(
        covariance,
        {
            "group_order",
            "port_count",
            "check_count",
            "rows",
            "rows_sha256",
            "full_twelve_port_orbit_checked_for_every_group_element",
            "literal_coordinate_feedback_rule_A5_equivariant",
        },
        "A5 covariance audit",
    )
    expected_covariance = _expected_covariance_rows(expected_records[0], actions)
    _fail(covariance["group_order"] == 60, "A5 group order")
    _fail(covariance["port_count"] == PORTS, "A5 port count")
    _fail(covariance["check_count"] == len(expected_covariance) == 720, "A5 checks")
    _fail_unless_json_exact(
        covariance["rows"], expected_covariance, "A5 covariance rows"
    )
    _fail(
        covariance["rows_sha256"]
        == _sha(covariance["rows"])
        == _sha(expected_covariance),
        "A5 row hash",
    )
    _fail(
        covariance["full_twelve_port_orbit_checked_for_every_group_element"] is True
        and covariance["literal_coordinate_feedback_rule_A5_equivariant"] is True,
        "A5 covariance verdict",
    )

    confluence = received["local_feedback_confluence_audit"]
    _exact_keys(
        confluence,
        {
            "scope",
            "control_state_rule",
            "idempotence_check_count",
            "idempotence_rows",
            "idempotence_rows_sha256",
            "disjoint_commutation_check_count",
            "disjoint_commutation_rows",
            "disjoint_commutation_rows_sha256",
            "all_port_reset_transactions_idempotent",
            "all_disjoint_port_reset_transactions_commute",
            "any_serialized_control_schedule_applying_every_port_reaches_literal_record",
            "parent_endpoint_repair_confluence_established",
        },
        "local feedback confluence",
    )
    expected_idempotence, expected_commutation = (
        _expected_local_confluence_rows(expected_records)
    )
    _fail(
        confluence["scope"]
        == "standalone_literal_record_reset_transactions_only",
        "confluence scope",
    )
    _fail(
        confluence["control_state_rule"]
        == "z_i,p=b_i,p+((carrier_index+port)_mod_3)-1",
        "confluence control rule",
    )
    _fail(
        confluence["idempotence_check_count"]
        == len(expected_idempotence)
        == 96,
        "idempotence count",
    )
    _fail_unless_json_exact(
        confluence["idempotence_rows"],
        expected_idempotence,
        "idempotence rows",
    )
    _fail(
        confluence["idempotence_rows_sha256"]
        == _sha(confluence["idempotence_rows"])
        == _sha(expected_idempotence),
        "idempotence row hash",
    )
    _fail(
        confluence["disjoint_commutation_check_count"]
        == len(expected_commutation)
        == 528,
        "commutation count",
    )
    _fail_unless_json_exact(
        confluence["disjoint_commutation_rows"],
        expected_commutation,
        "commutation rows",
    )
    _fail(
        confluence["disjoint_commutation_rows_sha256"]
        == _sha(confluence["disjoint_commutation_rows"])
        == _sha(expected_commutation),
        "commutation row hash",
    )
    _fail(
        confluence["all_port_reset_transactions_idempotent"] is True
        and confluence["all_disjoint_port_reset_transactions_commute"] is True
        and confluence[
            "any_serialized_control_schedule_applying_every_port_reaches_literal_record"
        ]
        is True
        and confluence["parent_endpoint_repair_confluence_established"] is False,
        "local confluence verdict",
    )

    attainment = received["attainment"]
    positive = {
        "LITERAL_SIGNED_RECORD_READ_RECEIPT",
        "READ_AFTER_COMMIT_RECEIPT",
        "BOUNDED_LOCAL_PORT_WRITE_RECEIPT",
        "EXACT_RECORD_CONDITIONED_STATE_RESTORATION_RECEIPT",
        "FEEDBACK_ABLATION_CHANGES_LATER_STATE_RECEIPT",
        "RECORD_COORDINATE_COUNTERFACTUAL_RECEIPT",
        "ALL_TWELVE_PORTS_CAUSALLY_COVERED_RECEIPT",
        "A5_EQUIVARIANT_LITERAL_FEEDBACK_RULE_RECEIPT",
        "IDEMPOTENT_LITERAL_FEEDBACK_TRANSACTION_RECEIPT",
        "DISJOINT_PORT_FEEDBACK_COMMUTATION_RECEIPT",
        "SERIALIZED_CONTROL_FEEDBACK_NORMAL_FORM_RECEIPT",
        "INTERNAL_FINITE_OBSERVER_LIKE_SELF_READING_RECEIPT",
    }
    negative = {
        "CANONICAL_A1_A2_A3_SOURCE_SELECTION_RECEIPT",
        "SOURCE_QUALIFIED_PHYSICAL_OBSERVER_RECEIPT",
        "SPATIAL_TRANSLATION_RECEIPT",
        "LABORATORY_RECORD_REALIZATION_RECEIPT",
        "PHYSICAL_PREDICTION_RECEIPT",
    }
    _exact_keys(attainment, positive | negative, "attainment")
    _fail(all(attainment[key] is True for key in positive), "finite attainment")
    _fail(all(attainment[key] is False for key in negative), "physical boundary")

    pins = received["implementation_pins"]
    _fail_unless_json_exact(
        pins,
        [_raw_pin(path) for path in (PRODUCER_PATH, VERIFIER_PATH, TEST_PATH)],
        "implementation pins",
    )
    claim_boundary = received["claim_boundary"]
    _fail(isinstance(claim_boundary, str) and "no laboratory" in claim_boundary, "claim boundary")
    return {
        "schema": "oph.vertex12-signed-record-feedback-independent-verification.v1",
        "receipt": True,
        "status": "PASS",
        "producer_imported": False,
        "literal_integer_transition_independently_replayed": True,
        "checked_record_count": len(expected_records),
        "checked_feedback_event_count": len(expected_events),
        "checked_carrier_port_pair_count": len(pair_set),
        "checked_A5_covariance_count": len(expected_covariance),
        "checked_feedback_idempotence_count": len(expected_idempotence),
        "checked_feedback_commutation_count": len(expected_commutation),
        "parent_endpoint_repair_confluence_verified": False,
        "physical_attachment_verified": False,
        "claim_boundary": (
            "The independent verifier establishes the internal finite literal-record "
            "feedback mechanics only. It does not select the source law or verify a "
            "physical realization."
        ),
    }


def verify_receipt(path: Path = DEFAULT_RECEIPT) -> dict[str, Any]:
    return verify_report(_load(path))


def _write_json(value: Mapping[str, Any], path: Path | None) -> None:
    rendered = json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if path is None:
        print(rendered, end="")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(rendered, encoding="utf-8")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("receipt", nargs="?", type=Path, default=DEFAULT_RECEIPT)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    try:
        result = verify_receipt(args.receipt)
    except IndependentVerificationError as exc:
        result = {
            "schema": "oph.vertex12-signed-record-feedback-independent-verification.v1",
            "receipt": False,
            "status": "FAIL",
            "reason": str(exc),
            "producer_imported": False,
            "physical_attachment_verified": False,
        }
    _write_json(result, args.output)
    return 0 if result["receipt"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
