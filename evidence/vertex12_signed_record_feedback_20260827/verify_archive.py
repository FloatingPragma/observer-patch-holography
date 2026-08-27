#!/usr/bin/env python3
"""Fail-closed standalone verifier for the signed-record feedback archive."""

from __future__ import annotations

import hashlib
import json
import sys
from itertools import combinations
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
CONTROL_FILES = {"README.md", "archive_manifest.json", "verify_archive.py"}


class VerificationError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def load_json(name: str) -> Any:
    with (ROOT / name).open("r", encoding="utf-8") as handle:
        return json.load(handle)


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")


def value_hash(value: Any) -> str:
    return "sha256:" + hashlib.sha256(canonical_bytes(value)).hexdigest()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_inventory(archive: dict[str, Any]) -> None:
    expected = {row["path"]: row for row in archive["inventory"]}
    actual = {
        path.name
        for path in ROOT.iterdir()
        if path.is_file() and path.name not in CONTROL_FILES
    }
    require(actual == set(expected), "archive inventory names differ")
    lines: list[str] = []
    total = 0
    for name in sorted(expected):
        path = ROOT / name
        row = expected[name]
        digest = file_hash(path)
        size = path.stat().st_size
        require(digest == row["sha256"], f"SHA-256 mismatch: {name}")
        require(size == row["bytes"], f"byte count mismatch: {name}")
        total += size
        lines.append(f"{digest}  {size}  {name}\n")
    curated = archive["curated_archive"]
    require(len(expected) == curated["file_count"], "file count mismatch")
    require(total == curated["total_bytes"], "total byte count mismatch")
    require(
        hashlib.sha256("".join(lines).encode("utf-8")).hexdigest()
        == curated["inventory_sha256"],
        "inventory digest mismatch",
    )


def verify_pins(receipt: dict[str, Any]) -> None:
    pins = receipt["implementation_pins"]
    by_name = {Path(pin["path"]).name: pin for pin in pins}
    for name in (
        "vertex12_signed_record_feedback.py",
        "verify_vertex12_signed_record_feedback_independent.py",
        "test_vertex12_signed_record_feedback.py",
    ):
        raw = (ROOT / name).read_bytes()
        pin = by_name[name]
        require(len(raw) == pin["bytes"], f"implementation bytes: {name}")
        require("sha256:" + hashlib.sha256(raw).hexdigest() == pin["sha256"], f"implementation hash: {name}")

    parents = receipt["parent_pins"]
    parent_files = {
        "atomic_port_transfer": "vertex12_atomic_port_transfer_receipt.json",
        "constructive_source_law": "vertex12_constructive_source_law_receipt.json",
    }
    for key, name in parent_files.items():
        parent = load_json(name)
        material = dict(parent)
        digest = material.pop("receipt_sha256")
        require(digest == value_hash(material), f"parent receipt digest: {name}")
        pin = parents[key]
        raw = (ROOT / name).read_bytes()
        require(pin["receipt_sha256"] == digest, f"parent receipt pin: {name}")
        require(pin["raw_pin"]["bytes"] == len(raw), f"parent byte pin: {name}")
        require(pin["raw_pin"]["sha256"] == "sha256:" + hashlib.sha256(raw).hexdigest(), f"parent hash pin: {name}")


def verify_feedback(receipt: dict[str, Any]) -> None:
    material = dict(receipt)
    digest = material.pop("receipt_sha256")
    require(digest == value_hash(material), "receipt digest")
    require(receipt["schema"] == "oph.vertex12-signed-record-feedback-diagnostic.v1", "schema")
    contract = receipt["finite_contract"]
    require(contract["carrier_count"] == 8 and contract["port_count"] == 12, "finite census")
    require(contract["record_is_full_literal_twelve_port_integer_row"] is True, "literal record contract")
    require(contract["hash_value_consumed_by_transition_rule"] is False, "hash transition boundary")
    require(contract["protected_record_is_not_mutated_by_feedback"] is True, "record protection")

    commit_log = receipt["literal_record_commit_log"]
    records = commit_log["records"]
    require(commit_log["record_count"] == len(records) == 8, "record count")
    require(commit_log["records_sha256"] == value_hash(records), "record log hash")
    for index, record in enumerate(records):
        values = record["signed_port_record"]
        require(record["carrier_index"] == index, "record carrier order")
        require(len(values) == 12, "record port count")
        require(all(type(value) is int and abs(value) <= 32 for value in values), "record bound")

    feedback = receipt["causal_feedback_log"]
    events = feedback["events"]
    require(feedback["event_count"] == len(events) == 96, "feedback event count")
    require(feedback["events_sha256"] == value_hash(events), "feedback log hash")
    pairs: set[tuple[int, int]] = set()
    event_ids: set[str] = set()
    for event in events:
        carrier = event["carrier_index"]
        port = event["port"]
        record = records[carrier]["signed_port_record"]
        probed = list(record)
        probed[port] += 1
        counterfactual_record = list(record)
        counterfactual_record[port] += 1
        pairs.add((carrier, port))
        event_ids.add(event["event_id"])
        require(event["live_state_before_probe"] == record, "preprobe state")
        require(event["live_state_before_feedback"] == probed, "probe state")
        require(event["literal_committed_record_port_value"] == record[port], "literal read value")
        require(event["feedback_delta"] == -1, "feedback delta")
        require(event["live_state_after_feedback"] == record, "exact restoration")
        require(event["protected_record_before"] == event["protected_record_after"] == record, "record mutation")
        require(event["exact_preprobe_state_restored"] is True, "restoration receipt")
        require(event["ablation"]["live_state_after_action"] == probed, "ablation state")
        require(event["ablation"]["feedback_delta"] == 0, "ablation delta")
        counterfactual = event["record_coordinate_counterfactual"]
        require(counterfactual["record_before"] == record, "counterfactual source record")
        require(counterfactual["record_after_intervention"] == counterfactual_record, "counterfactual record")
        require(counterfactual["feedback_delta"] == 0, "counterfactual feedback delta")
        require(counterfactual["live_state_after_action"] == probed, "counterfactual later state")
        require(event["commit_event_index"] < event["probe_event_index"] < event["read_event_index"] < event["write_event_index"], "event order")
    require(pairs == {(carrier, port) for carrier in range(8) for port in range(12)}, "carrier-port coverage")
    require(len(event_ids) == 96, "event-id uniqueness")

    covariance = receipt["a5_covariance_audit"]
    rows = covariance["rows"]
    require(covariance["group_order"] == 60, "A5 order")
    require(covariance["check_count"] == len(rows) == 720, "A5 check count")
    require(covariance["rows_sha256"] == value_hash(rows), "A5 rows hash")
    require({(row["group_element_index"], row["source_port"]) for row in rows} == {(group, port) for group in range(60) for port in range(12)}, "A5 source coverage")
    require(all(row["literal_coordinate_feedback_commutes_with_port_action"] is True for row in rows), "A5 covariance row")
    require(all(row["source_feedback_delta"] == row["target_feedback_delta"] == -1 for row in rows), "A5 feedback delta")
    require(all(row["source_counterfactual_feedback_delta"] == row["target_counterfactual_feedback_delta"] == 0 for row in rows), "A5 counterfactual delta")

    confluence = receipt["local_feedback_confluence_audit"]
    idempotence = confluence["idempotence_rows"]
    commutation = confluence["disjoint_commutation_rows"]
    require(confluence["idempotence_check_count"] == len(idempotence) == 96, "idempotence count")
    require(confluence["idempotence_rows_sha256"] == value_hash(idempotence), "idempotence hash")
    require({(row["carrier_index"], row["port"]) for row in idempotence} == pairs, "idempotence coverage")
    require(all(row["after_once_port_value"] == row["after_twice_port_value"] == row["literal_record_port_value"] for row in idempotence), "idempotence values")
    require(all(row["transaction_idempotent"] is True and row["other_coordinates_unchanged"] is True for row in idempotence), "idempotence verdict")
    expected_pairs = {(carrier, left, right) for carrier in range(8) for left, right in combinations(range(12), 2)}
    require(confluence["disjoint_commutation_check_count"] == len(commutation) == 528, "commutation count")
    require(confluence["disjoint_commutation_rows_sha256"] == value_hash(commutation), "commutation hash")
    require({(row["carrier_index"], row["left_port"], row["right_port"]) for row in commutation} == expected_pairs, "commutation coverage")
    require(all(row["disjoint_transactions_commute"] is True and row["all_coordinates_equal"] is True for row in commutation), "commutation verdict")

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
    attainment = receipt["attainment"]
    require(set(attainment) == positive | negative, "attainment key set")
    require(all(attainment[key] is True for key in positive), "finite attainment")
    require(all(attainment[key] is False for key in negative), "physical boundary")
    verify_pins(receipt)


def verify_independent_output() -> None:
    report = load_json("independent_verification.json")
    require(report["receipt"] is True and report["status"] == "PASS", "independent verifier output")
    require(report["producer_imported"] is False, "producer imported by independent verifier")
    require(report["literal_integer_transition_independently_replayed"] is True, "independent transition replay")
    require(report["checked_record_count"] == 8, "independent record count")
    require(report["checked_feedback_event_count"] == 96, "independent event count")
    require(report["checked_A5_covariance_count"] == 720, "independent A5 count")
    require(report["checked_feedback_idempotence_count"] == 96, "independent idempotence count")
    require(report["checked_feedback_commutation_count"] == 528, "independent commutation count")
    require(report["parent_endpoint_repair_confluence_verified"] is False, "parent confluence boundary")
    require(report["physical_attachment_verified"] is False, "physical boundary")


def main() -> int:
    try:
        archive = load_json("archive_manifest.json")
        require(archive["schema"] == "oph.signed_record_feedback_archive.v1", "archive schema")
        verify_inventory(archive)
        verify_feedback(load_json("vertex12_signed_record_feedback_receipt.json"))
        verify_independent_output()
    except (OSError, KeyError, TypeError, ValueError, VerificationError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1
    print("PASS: signed records, causal feedback, controls, covariance, confluence, and custody verify")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
