from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path

import pytest

from capacity_indexed_source_family import (
    BRANCH_IDS,
    CERTIFICATE_PATH,
    PROJECTION_PATH,
    SAMPLE_RUNGS,
    SPECTATOR_MULTIPLICITIES,
    build_capacity_packet,
    build_certificate,
    build_projection,
    canonical_json_bytes,
    evaluate_capacity_packet,
    source_signature,
    tagged_sha256,
)


def test_one_generator_emits_every_positive_rung_without_per_d_table():
    for k in (1, 2, 5, 11):
        for branch in BRANCH_IDS:
            spectators = (1, 2) if branch == "hidden_spectator" else (1,)
            for spectator in spectators:
                packet = build_capacity_packet(
                    branch, k, spectator_multiplicity=spectator
                )
                result = evaluate_capacity_packet(packet)
                assert result["status"] == "PASS"
                assert packet["k"] == k
                assert packet["shared_source_signature_sha256"] == source_signature()


def test_exact_all_rung_formulas_on_a_larger_unadvertised_rung():
    k = 7
    expected = {
        "reversible_identity": 24 * k,
        "copy_collapse_erasure": 24,
        "capped_two_class": 48,
        "hidden_spectator": 24 * k,
    }
    for branch, capacity in expected.items():
        spectator = 3 if branch == "hidden_spectator" else 1
        result = evaluate_capacity_packet(
            build_capacity_packet(
                branch, k, spectator_multiplicity=spectator
            )
        )
        assert result["exact_zero_error_capacity"] == capacity


def test_same_antecedent_branches_have_different_exact_zero_sets():
    projection = build_projection()
    assert [row["branch_id"] for row in projection["branches"]] == list(BRANCH_IDS)
    observed = {
        row["branch_id"]: [
            (zero["k"], zero["spectator_multiplicity"])
            for zero in row["claimed_bounded_zero_set"]
        ]
        for row in projection["branches"]
    }
    assert observed["reversible_identity"] == [(k, 1) for k in SAMPLE_RUNGS]
    assert observed["copy_collapse_erasure"] == [(1, 1)]
    assert observed["capped_two_class"] == [(1, 1), (2, 1)]
    assert observed["hidden_spectator"] == [(k, 1) for k in SAMPLE_RUNGS]


def test_hidden_spectator_changes_raw_d_without_changing_public_capacity():
    for k in SAMPLE_RUNGS:
        capacities = []
        dimensions = []
        for spectator in SPECTATOR_MULTIPLICITIES:
            result = evaluate_capacity_packet(
                build_capacity_packet(
                    "hidden_spectator",
                    k,
                    spectator_multiplicity=spectator,
                )
            )
            capacities.append(result["exact_zero_error_capacity"])
            dimensions.append(result["raw_dimension"])
        assert capacities == [24 * k] * len(SPECTATOR_MULTIPLICITIES)
        assert dimensions == [24 * k * s for s in SPECTATOR_MULTIPLICITIES]


def test_continuation_composition_sewing_and_extensions_are_executable():
    certificate = build_certificate()
    assert certificate["status"] == "BOUNDED_COMPLETION_CLASS_NONIDENTIFIABLE"
    assert certificate["zero_sets_differ"] is True
    assert certificate["physical_n_closure_promoted"] is False
    for branch in certificate["branch_receipts"].values():
        for family in (
            branch["continuation_composition"],
            branch["sewing"],
            branch["extension"],
        ):
            assert family
            assert all(receipt["status"] == "PASS" for receipt in family)
        for receipt in branch["sewing"]:
            assert receipt["source_port_count"] == 12
            assert receipt["source_orientation_count"] == 2
            assert receipt["output_bijection_verified"] is True
            assert receipt["fiber_product_count"] == receipt[
                "expected_public_output_count"
            ]
            assert receipt["left_section_count"] > receipt[
                "seam_label_count"
            ]
            assert receipt["right_section_count"] > receipt[
                "seam_label_count"
            ]


def test_target_flags_are_fail_closed():
    packet = build_capacity_packet("reversible_identity", 2)
    tainted = copy.deepcopy(packet)
    tainted["target_cleanliness"]["measured_cosmological_constant_read"] = True
    raw = dict(tainted)
    raw.pop("packet_sha256")
    tainted["packet_sha256"] = tagged_sha256(canonical_json_bytes(raw))
    with pytest.raises(ValueError, match="target-tainted"):
        evaluate_capacity_packet(tainted)


def test_packet_and_projection_mutations_fail():
    packet = build_capacity_packet("copy_collapse_erasure", 3)
    tampered = copy.deepcopy(packet)
    first = tampered["records"][0]
    tampered["deterministic_collapse_outputs"][first] = "invented-output"
    with pytest.raises(ValueError, match="hash mismatch"):
        evaluate_capacity_packet(tampered)

    projection = build_projection()
    projection["branches"][0]["shared_source_signature_sha256"] = "sha256:" + "0" * 64
    certificate = build_certificate(projection)
    assert certificate["projection_sha256"] == tagged_sha256(
        canonical_json_bytes(projection)
    )
    # Independent replay owns semantic rejection. The producer certificate
    # records exact bytes and never silently repairs the mutated projection.


def test_runtime_receipts_are_canonical_and_byte_exact():
    projection = build_projection()
    certificate = build_certificate(projection)
    assert PROJECTION_PATH.read_bytes() == canonical_json_bytes(projection)
    assert CERTIFICATE_PATH.read_bytes() == canonical_json_bytes(certificate)
    assert json.loads(PROJECTION_PATH.read_text(encoding="ascii")) == projection


def test_vendored_independent_replay_is_bound_to_exact_projection_and_commit():
    runtime = PROJECTION_PATH.parent
    receipt_path = (
        runtime / "capacity_indexed_source_family_independent_receipt.json"
    )
    custody_path = (
        runtime / "capacity_indexed_source_family_independent_custody.json"
    )
    receipt_bytes = receipt_path.read_bytes()
    receipt = json.loads(receipt_bytes)
    custody = json.loads(custody_path.read_text(encoding="utf-8"))

    assert receipt["status"] == "PASS"
    assert receipt["scientific_verdict_replayed"] == (
        "BOUNDED_COMPLETION_CLASS_NONIDENTIFIABLE"
    )
    assert receipt["projection_sha256"] == "sha256:" + hashlib.sha256(
        PROJECTION_PATH.read_bytes()
    ).hexdigest()
    assert receipt["scope"]["producer_implementation_independent"] is True
    assert receipt["scope"]["all_positive_integer_rungs_proved"] is False
    assert custody["status"] == "PASS"
    assert len(custody["commit"]) == 40
    assert all(character in "0123456789abcdef" for character in custody["commit"])

    pins = {row["path"]: row["sha256"] for row in custody["artifacts"]}
    assert pins[
        "data/capacity_readback/capacity_indexed_source_family_projection.json"
    ] == hashlib.sha256(PROJECTION_PATH.read_bytes()).hexdigest()
    assert pins[
        "data/capacity_readback/capacity_indexed_source_family_independent_receipt.json"
    ] == hashlib.sha256(receipt_bytes).hexdigest()
