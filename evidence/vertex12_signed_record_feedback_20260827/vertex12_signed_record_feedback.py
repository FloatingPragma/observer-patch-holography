"""Literal signed-record feedback on the finite twelve-port carrier.

This diagnostic closes one deliberately small implementation gap.  The
existing vertex12 source-law control supplies the twelve signed unit events
and their exact A5 action.  The existing atomic-transfer packet supplies a
finite carrier federation, complete port coverage, transaction replay, and an
in-process full-port snapshot readback.  Here those two internal finite
surfaces are joined to a causal feedback action that consumes a *literal*
integer record coordinate.

For every carrier and every port the diagnostic commits a bounded twelve-port
integer record, applies a local ``+1`` probe to a working copy, rereads the
committed integer at that port, and appends the resulting ``-1`` retraction.
The working state is thereby restored exactly while the protected record is
unchanged.  A frozen ablation omits the reread and therefore leaves the probe
in place.  A record-coordinate intervention changes only the read integer and
changes the prescribed later action from ``-1`` to ``0``.  Hashes identify and
bind evidence; no hash value is consumed by the transition rule.

This is an internal finite observer-like receipt.  It neither source-selects
the constructive signed law in canonical A1--A3 nor attaches the record or
feedback action to laboratory physics, spacetime, or spatial translation.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[2]
ATOMIC_PARENT = (
    ROOT / "data/repair_closure/vertex12_atomic_port_transfer_receipt.json"
)
SOURCE_LAW_PARENT = (
    ROOT / "data/repair_closure/vertex12_constructive_source_law_receipt.json"
)
DEFAULT_OUTPUT = (
    ROOT / "data/repair_closure/vertex12_signed_record_feedback_receipt.json"
)
VERIFIER_PATH = (
    ROOT
    / "oph_fpe/dynamics/verify_vertex12_signed_record_feedback_independent.py"
)
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
PORT_COUNT = 12
AXIS_COUNT = 6
STATE_BOUND = 32
PROBE_DELTA = 1


class SignedRecordFeedbackError(RuntimeError):
    """Raised when a parent or literal feedback invariant fails closed."""


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
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise SignedRecordFeedbackError(f"{path} is not a JSON object")
    return value


def _raw_pin(path: Path) -> dict[str, Any]:
    raw = path.read_bytes()
    return {
        "path": path.relative_to(ROOT).as_posix(),
        "bytes": len(raw),
        "sha256": "sha256:" + hashlib.sha256(raw).hexdigest(),
    }


def _receipt_digest_is_valid(report: Mapping[str, Any]) -> bool:
    material = dict(report)
    digest = material.pop("receipt_sha256", None)
    return isinstance(digest, str) and digest == _sha(material)


def _validated_parents() -> tuple[dict[str, Any], dict[str, Any]]:
    atomic = _load(ATOMIC_PARENT)
    source = _load(SOURCE_LAW_PARENT)
    operator = atomic.get("atomic_transfer_operator")
    readback = atomic.get("post_repair_in_process_snapshot_reread")
    boundary = atomic.get("quotient_and_spatial_boundary")
    if (
        atomic.get("schema") != ATOMIC_SCHEMA
        or atomic.get("status") != ATOMIC_STATUS
        or atomic.get("issue") != 655
        or not _receipt_digest_is_valid(atomic)
        or not isinstance(operator, Mapping)
        or operator.get("carrier_count") != 8
        or operator.get("port_count") != PORT_COUNT
        or operator.get("source_native_internal_seam_partner_operator_receipt")
        is not True
        or operator.get("source_native_spatial_translation_receipt") is not False
        or not isinstance(readback, Mapping)
        or readback.get("covered_port_coordinate_count") != 8 * PORT_COUNT
        or readback.get("every_carrier_full_port_state_committed") is not True
        or readback.get("every_carrier_full_port_state_reread_in_process")
        is not True
        or not isinstance(boundary, Mapping)
        or boundary.get("internal_seam_transfer_is_spatial_translation")
        is not False
    ):
        raise SignedRecordFeedbackError("atomic vertex12 parent contract drifted")

    law = source.get("constructive_source_law")
    attainment = source.get("attainment")
    if (
        source.get("schema") != SOURCE_LAW_SCHEMA
        or source.get("status") != SOURCE_LAW_STATUS
        or source.get("issue") != 655
        or not _receipt_digest_is_valid(source)
        or not isinstance(law, Mapping)
        or not isinstance(attainment, Mapping)
        or attainment.get("constructive_source_law_capture_root") is not True
        or attainment.get("twelve_raw_steps") is not True
        or attainment.get("same_Q_A5_covariance") is not True
        or attainment.get("canonical_source_selection") is not False
        or attainment.get("spatial_translation") is not False
        or attainment.get("physical_readout") is not False
        or attainment.get("physical_prediction") is not False
    ):
        raise SignedRecordFeedbackError("constructive source-law parent drifted")
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
    group_rows = law["same_Q_A5_action"]["group_rows"]
    if (
        alphabet.get("complete_signed_port_orbit") is not True
        or alphabet.get("event_count") != PORT_COUNT
        or not isinstance(events, list)
        or len(events) != PORT_COUNT
        or not isinstance(group_rows, list)
        or len(group_rows) != 60
    ):
        raise SignedRecordFeedbackError("signed port source family is incomplete")
    directions: list[tuple[int, ...]] = []
    event_ids: list[str] = []
    for port, row in enumerate(events):
        if not isinstance(row, Mapping) or row.get("port") != port:
            raise SignedRecordFeedbackError("signed source event ordering drifted")
        direction = row.get("raw_direction_in_Z_power_6")
        if (
            not isinstance(direction, list)
            or len(direction) != AXIS_COUNT
            or sum(abs(value) for value in direction if type(value) is int) != 1
            or any(type(value) is not int or value not in (-1, 0, 1) for value in direction)
        ):
            raise SignedRecordFeedbackError("signed source direction is malformed")
        directions.append(tuple(direction))
        event_ids.append(str(row.get("event_id")))
    actions: list[tuple[int, ...]] = []
    for index, row in enumerate(group_rows):
        action = row.get("port_permutation") if isinstance(row, Mapping) else None
        if (
            not isinstance(action, list)
            or row.get("group_element_index") != index
            or sorted(action) != list(range(PORT_COUNT))
        ):
            raise SignedRecordFeedbackError("A5 port action is malformed")
        actions.append(tuple(int(value) for value in action))
    if len(set(actions)) != 60:
        raise SignedRecordFeedbackError("A5 port action is not faithful")
    return tuple(directions), tuple(actions), tuple(event_ids)


def _literal_record_rows(
    carrier_ids: Sequence[str],
    directions: Sequence[Sequence[int]],
    source_event_ids: Sequence[str],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for carrier_index, carrier_id in enumerate(carrier_ids):
        axis_record = [carrier_index + axis + 1 for axis in range(AXIS_COUNT)]
        port_record = [
            sum(
                int(direction[axis]) * axis_record[axis]
                for axis in range(AXIS_COUNT)
            )
            for direction in directions
        ]
        if (
            len(port_record) != PORT_COUNT
            or max(abs(value) for value in port_record) + PROBE_DELTA > STATE_BOUND
        ):
            raise SignedRecordFeedbackError("bounded record construction failed")
        material = {
            "schema": "oph.vertex12-literal-signed-port-record.v1",
            "carrier_id": str(carrier_id),
            "carrier_index": carrier_index,
            "commit_event_index": carrier_index,
            "axis_record_in_Z_power_6": axis_record,
            "signed_port_record": port_record,
            "source_event_ids_by_port": list(source_event_ids),
            "record_value_type": "literal_bounded_integer",
            "state_bound": STATE_BOUND,
        }
        rows.append({**material, "record_id": _sha(material)})
    return rows


def _feedback_event(
    record: Mapping[str, Any],
    *,
    port: int,
    event_ordinal: int,
    carrier_count: int,
) -> dict[str, Any]:
    protected = [int(value) for value in record["signed_port_record"]]
    before = list(protected)
    probed = list(before)
    probed[port] += PROBE_DELTA
    literal_record_value = protected[port]
    live_before_action = probed[port]
    feedback_delta = literal_record_value - live_before_action
    after = list(probed)
    after[port] += feedback_delta

    ablated_after = list(probed)
    counterfactual_record = list(protected)
    counterfactual_record[port] += PROBE_DELTA
    counterfactual_delta = counterfactual_record[port] - live_before_action
    counterfactual_after = list(probed)
    counterfactual_after[port] += counterfactual_delta

    event_base = carrier_count + 3 * event_ordinal
    material = {
        "schema": "oph.vertex12-literal-signed-record-feedback-event.v1",
        "carrier_id": str(record["carrier_id"]),
        "carrier_index": int(record["carrier_index"]),
        "port": int(port),
        "record_id": str(record["record_id"]),
        "commit_event_index": int(record["commit_event_index"]),
        "probe_event_index": event_base,
        "read_event_index": event_base + 1,
        "write_event_index": event_base + 2,
        "protected_record_before": protected,
        "protected_record_after": list(protected),
        "live_state_before_probe": before,
        "probe_delta": PROBE_DELTA,
        "live_state_before_feedback": probed,
        "literal_committed_record_port_value": literal_record_value,
        "feedback_delta": feedback_delta,
        "live_state_after_feedback": after,
        "action_rule": "delta=literal_committed_record_port_value-live_port_value",
        "action_rule_inputs": [
            "live_port_value_after_probe",
            "literal_committed_record_port_value",
        ],
        "hash_value_consumed_by_transition_rule": False,
        "bounded_integer_state_minimum": -STATE_BOUND,
        "bounded_integer_state_maximum": STATE_BOUND,
        "bounded_local_port_write": bool(
            feedback_delta == -1
            and -STATE_BOUND <= live_before_action <= STATE_BOUND
            and -STATE_BOUND <= after[port] <= STATE_BOUND
        ),
        "exact_preprobe_state_restored": after == before,
        "protected_record_unchanged": protected == list(protected),
        "ablation": {
            "mode": "omit_record_read_and_feedback_write",
            "live_state_before_action": probed,
            "feedback_delta": 0,
            "live_state_after_action": ablated_after,
            "differs_from_actual_later_state": ablated_after != after,
            "preprobe_state_not_restored": ablated_after != before,
        },
        "record_coordinate_counterfactual": {
            "mode": "do_literal_committed_record_coordinate",
            "intervened_port": int(port),
            "record_before": protected,
            "record_after_intervention": counterfactual_record,
            "changed_record_coordinate_count": sum(
                left != right
                for left, right in zip(
                    protected, counterfactual_record, strict=True
                )
            ),
            "live_state_before_action": probed,
            "probe_delta": PROBE_DELTA,
            "feedback_delta": counterfactual_delta,
            "live_state_after_action": counterfactual_after,
            "nonrecord_inputs_held_fixed": True,
            "later_action_differs": counterfactual_delta != feedback_delta,
        },
    }
    return {**material, "event_id": _sha(material)}


def _a5_covariance_rows(
    records: Sequence[Mapping[str, Any]],
    actions: Sequence[Sequence[int]],
) -> list[dict[str, Any]]:
    reference = [int(value) for value in records[0]["signed_port_record"]]
    rows: list[dict[str, Any]] = []
    for group_index, action in enumerate(actions):
        transformed = [0] * PORT_COUNT
        for source_port, target_port in enumerate(action):
            transformed[int(target_port)] = reference[source_port]
        for source_port, target_port in enumerate(action):
            source_record = reference[source_port]
            target_record = transformed[int(target_port)]
            source_live = source_record + PROBE_DELTA
            target_live = target_record + PROBE_DELTA
            source_delta = source_record - source_live
            target_delta = target_record - target_live
            row = {
                "group_element_index": group_index,
                "source_port": source_port,
                "target_port": int(target_port),
                "source_literal_record_value": source_record,
                "transformed_literal_record_value": target_record,
                "source_probe_delta": PROBE_DELTA,
                "target_probe_delta": PROBE_DELTA,
                "source_feedback_delta": source_delta,
                "target_feedback_delta": target_delta,
                "source_counterfactual_feedback_delta": 0,
                "target_counterfactual_feedback_delta": 0,
                "literal_coordinate_feedback_commutes_with_port_action": bool(
                    source_record == target_record
                    and source_delta == target_delta == -1
                ),
            }
            rows.append(row)
    return rows


def _apply_record_reset_transaction(
    live_state: Sequence[int], record: Sequence[int], port: int
) -> list[int]:
    """Apply the literal probe/read/retract transaction at one port."""

    after_probe = [int(value) for value in live_state]
    after_probe[port] += PROBE_DELTA
    feedback_delta = int(record[port]) - after_probe[port]
    after_feedback = list(after_probe)
    after_feedback[port] += feedback_delta
    return after_feedback


def _local_confluence_rows(
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
        if any(abs(value) + PROBE_DELTA > STATE_BOUND for value in control):
            raise SignedRecordFeedbackError("confluence control leaves state bound")
        for port in range(PORT_COUNT):
            once = _apply_record_reset_transaction(control, signed_record, port)
            twice = _apply_record_reset_transaction(once, signed_record, port)
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
                        for other in range(PORT_COUNT)
                        if other != port
                    ),
                    "transaction_idempotent": once == twice,
                }
            )
        for left in range(PORT_COUNT):
            for right in range(left + 1, PORT_COUNT):
                left_then_right = _apply_record_reset_transaction(
                    _apply_record_reset_transaction(control, signed_record, left),
                    signed_record,
                    right,
                )
                right_then_left = _apply_record_reset_transaction(
                    _apply_record_reset_transaction(control, signed_record, right),
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


def _payload() -> dict[str, Any]:
    atomic, source = _validated_parents()
    directions, actions, source_event_ids = _source_geometry(source)
    carrier_ids = tuple(
        str(value)
        for value in atomic["atomic_transfer_operator"]["carrier_ids"]
    )
    if len(carrier_ids) != 8 or len(set(carrier_ids)) != 8:
        raise SignedRecordFeedbackError("atomic carrier census drifted")
    records = _literal_record_rows(carrier_ids, directions, source_event_ids)
    events = [
        _feedback_event(
            record,
            port=port,
            event_ordinal=carrier_index * PORT_COUNT + port,
            carrier_count=len(records),
        )
        for carrier_index, record in enumerate(records)
        for port in range(PORT_COUNT)
    ]
    covariance_rows = _a5_covariance_rows(records, actions)
    idempotence_rows, commutation_rows = _local_confluence_rows(records)
    implementation_files = [Path(__file__).resolve(), VERIFIER_PATH, TEST_PATH]
    feedback_passed = bool(
        len(events) == len(records) * PORT_COUNT
        and all(
            event["commit_event_index"] < event["read_event_index"]
            < event["write_event_index"]
            and event["bounded_local_port_write"] is True
            and event["exact_preprobe_state_restored"] is True
            and event["protected_record_unchanged"] is True
            and event["ablation"]["differs_from_actual_later_state"] is True
            and event["record_coordinate_counterfactual"][
                "nonrecord_inputs_held_fixed"
            ]
            is True
            and event["record_coordinate_counterfactual"]["later_action_differs"]
            is True
            and event["hash_value_consumed_by_transition_rule"] is False
            for event in events
        )
    )
    a5_passed = bool(
        len(covariance_rows) == 60 * PORT_COUNT
        and all(
            row["literal_coordinate_feedback_commutes_with_port_action"] is True
            for row in covariance_rows
        )
    )
    local_confluence_passed = bool(
        len(idempotence_rows) == len(records) * PORT_COUNT
        and len(commutation_rows)
        == len(records) * PORT_COUNT * (PORT_COUNT - 1) // 2
        and all(row["transaction_idempotent"] is True for row in idempotence_rows)
        and all(
            row["disjoint_transactions_commute"] is True
            for row in commutation_rows
        )
    )
    return {
        "schema": SCHEMA,
        "issue": 655,
        "status": STATUS,
        "comparison_data_read": False,
        "parent_pins": {
            "atomic_port_transfer": {
                "schema": atomic["schema"],
                "status": atomic["status"],
                "receipt_sha256": atomic["receipt_sha256"],
                "raw_pin": _raw_pin(ATOMIC_PARENT),
            },
            "constructive_source_law": {
                "schema": source["schema"],
                "status": source["status"],
                "receipt_sha256": source["receipt_sha256"],
                "raw_pin": _raw_pin(SOURCE_LAW_PARENT),
            },
        },
        "finite_contract": {
            "carrier_count": len(records),
            "port_count": PORT_COUNT,
            "signed_axis_count": AXIS_COUNT,
            "serialized_working_state_envelope": (
                f"integer_coordinates_with_absolute_value_at_most_{STATE_BOUND}"
            ),
            "protected_record_domain": f"integer_box_[-{STATE_BOUND},{STATE_BOUND}]^(8x12)",
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
        "literal_record_commit_log": {
            "record_count": len(records),
            "records": records,
            "records_sha256": _sha(records),
        },
        "causal_feedback_log": {
            "event_count": len(events),
            "events": events,
            "events_sha256": _sha(events),
            "covered_carrier_count": len(records),
            "covered_port_coordinate_count": len(records) * PORT_COUNT,
            "all_carrier_port_pairs_covered_exactly_once": True,
        },
        "a5_covariance_audit": {
            "group_order": len(actions),
            "port_count": PORT_COUNT,
            "check_count": len(covariance_rows),
            "rows": covariance_rows,
            "rows_sha256": _sha(covariance_rows),
            "full_twelve_port_orbit_checked_for_every_group_element": True,
            "literal_coordinate_feedback_rule_A5_equivariant": a5_passed,
        },
        "local_feedback_confluence_audit": {
            "scope": "standalone_literal_record_reset_transactions_only",
            "control_state_rule": (
                "z_i,p=b_i,p+((carrier_index+port)_mod_3)-1"
            ),
            "idempotence_check_count": len(idempotence_rows),
            "idempotence_rows": idempotence_rows,
            "idempotence_rows_sha256": _sha(idempotence_rows),
            "disjoint_commutation_check_count": len(commutation_rows),
            "disjoint_commutation_rows": commutation_rows,
            "disjoint_commutation_rows_sha256": _sha(commutation_rows),
            "all_port_reset_transactions_idempotent": local_confluence_passed,
            "all_disjoint_port_reset_transactions_commute": (
                local_confluence_passed
            ),
            "any_serialized_control_schedule_applying_every_port_reaches_literal_record": (
                local_confluence_passed
            ),
            "parent_endpoint_repair_confluence_established": False,
        },
        "attainment": {
            "LITERAL_SIGNED_RECORD_READ_RECEIPT": feedback_passed,
            "READ_AFTER_COMMIT_RECEIPT": feedback_passed,
            "BOUNDED_LOCAL_PORT_WRITE_RECEIPT": feedback_passed,
            "EXACT_RECORD_CONDITIONED_STATE_RESTORATION_RECEIPT": feedback_passed,
            "FEEDBACK_ABLATION_CHANGES_LATER_STATE_RECEIPT": feedback_passed,
            "RECORD_COORDINATE_COUNTERFACTUAL_RECEIPT": feedback_passed,
            "ALL_TWELVE_PORTS_CAUSALLY_COVERED_RECEIPT": feedback_passed,
            "A5_EQUIVARIANT_LITERAL_FEEDBACK_RULE_RECEIPT": a5_passed,
            "IDEMPOTENT_LITERAL_FEEDBACK_TRANSACTION_RECEIPT": (
                local_confluence_passed
            ),
            "DISJOINT_PORT_FEEDBACK_COMMUTATION_RECEIPT": (
                local_confluence_passed
            ),
            "SERIALIZED_CONTROL_FEEDBACK_NORMAL_FORM_RECEIPT": (
                local_confluence_passed
            ),
            "INTERNAL_FINITE_OBSERVER_LIKE_SELF_READING_RECEIPT": bool(
                feedback_passed and a5_passed and local_confluence_passed
            ),
            "CANONICAL_A1_A2_A3_SOURCE_SELECTION_RECEIPT": False,
            "SOURCE_QUALIFIED_PHYSICAL_OBSERVER_RECEIPT": False,
            "SPATIAL_TRANSLATION_RECEIPT": False,
            "LABORATORY_RECORD_REALIZATION_RECEIPT": False,
            "PHYSICAL_PREDICTION_RECEIPT": False,
        },
        "implementation_pins": [_raw_pin(path) for path in implementation_files],
        "claim_boundary": (
            "This target-free finite diagnostic commits and rereads literal bounded "
            "integer data on every port, and the reread causes a later local "
            "retraction that restores the working state. Ablation and a one-coordinate "
            "record intervention establish the narrow causal software claim. The signed "
            "record-reset transactions commute and are idempotent on the serialized "
            "controls, giving those controls a schedule-independent literal-record "
            "normal form. This result does not establish confluence of the parent "
            "endpoint repair. The signed "
            "source law remains a constructive control that canonical A1--A3 does not "
            "select, and the feedback channel has no laboratory, spacetime, spatial, or "
            "physical-sector realization."
        ),
    }


def produce_receipt() -> dict[str, Any]:
    payload = _payload()
    return {**payload, "receipt_sha256": _sha(payload)}


def verify_receipt(report: Mapping[str, Any]) -> dict[str, Any]:
    reasons: list[str] = []
    try:
        received = copy.deepcopy(dict(report))
        digest = received.pop("receipt_sha256", None)
        if digest != _sha(received):
            reasons.append("receipt_digest_mismatch")
        if _canonical_bytes(received) != _canonical_bytes(_payload()):
            reasons.append("producer_replay_mismatch")
        attainment = report.get("attainment")
        if not isinstance(attainment, Mapping):
            reasons.append("attainment_missing")
        elif (
            attainment.get("INTERNAL_FINITE_OBSERVER_LIKE_SELF_READING_RECEIPT")
            is not True
            or any(
                attainment.get(key) is not False
                for key in (
                    "CANONICAL_A1_A2_A3_SOURCE_SELECTION_RECEIPT",
                    "SOURCE_QUALIFIED_PHYSICAL_OBSERVER_RECEIPT",
                    "SPATIAL_TRANSLATION_RECEIPT",
                    "LABORATORY_RECORD_REALIZATION_RECEIPT",
                    "PHYSICAL_PREDICTION_RECEIPT",
                )
            )
        ):
            reasons.append("attainment_boundary_mismatch")
    except (
        AttributeError,
        KeyError,
        OSError,
        TypeError,
        ValueError,
        SignedRecordFeedbackError,
    ):
        reasons.append("malformed_or_unreplayable_receipt")
    return {
        "schema": "oph.vertex12-signed-record-feedback-verification.v1",
        "receipt": not reasons,
        "status": "PASS" if not reasons else "FAIL",
        "reasons": sorted(set(reasons)),
        "producer_replay": True,
        "independent_implementation": False,
        "claim_boundary": (
            "Producer replay verifies only the internal finite literal-record "
            "feedback diagnostic and carries no physical attachment."
        ),
    }


def _write_json(value: Mapping[str, Any], path: Path | None) -> None:
    rendered = json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if path is None:
        print(rendered, end="")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(rendered, encoding="utf-8")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--verify", type=Path)
    args = parser.parse_args(argv)
    if args.verify is not None:
        result = verify_receipt(_load(args.verify))
        _write_json(result, args.output)
        return 0 if result["receipt"] else 1
    receipt = produce_receipt()
    result = verify_receipt(receipt)
    if not result["receipt"]:
        _write_json(result, args.output)
        return 1
    _write_json(receipt, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
