#!/usr/bin/env python3
"""Tests for the bounded issue-566 source-current capability audit."""

from __future__ import annotations

import copy
import json
from pathlib import Path
import sys

import pytest


MODULE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODULE_DIR))

import source_current_capability_certificate as producer  # noqa: E402
import verify_source_current_capability_independent as independent  # noqa: E402


def source_packet() -> tuple[dict, dict]:
    carrier = producer.load_json(producer.CARRIER_PATH)
    response = producer.load_json(producer.RESPONSE_PATH)
    projection = producer.build_projection(carrier, response)
    return projection, producer.certificate_payload(projection)


def test_committed_packet_replays_exactly_and_independently() -> None:
    projection, receipt = source_packet()
    assert projection == producer.load_json(producer.PROJECTION_PATH)
    assert receipt == producer.load_json(producer.RECEIPT_PATH)
    producer.verify_receipt(projection, receipt)
    report = independent.verify(producer.PROJECTION_PATH, producer.RECEIPT_PATH)
    assert report == {
        "verdict": producer.VERDICT,
        "response_algebra_dimension": 4,
        "proper_recharting_count": 60,
        "response_word_recharting_intersection": 1,
    }


def test_bounded_obstruction_preserves_the_abstract_theorem() -> None:
    _, receipt = source_packet()
    assert receipt["physical_current_source_bridge_attained"] is False
    assert receipt["abstract_forced_lie_type_theorem_preserved"] is True
    assert "single registered adjacency recurrence" in receipt["audit_scope"]
    assert "no bound is asserted" in receipt["audit_scope"]
    algebra = receipt["audit"]["response_word_algebra"]
    assert algebra["basis"] == ["I", "A", "A^2", "A^3"]
    assert algebra["exact_dimension"] == 4
    assert algebra["closure_rank_through_degree_six"] == 4
    assert algebra["commutator_nonzero_count"] == 0
    assert algebra["skew_pairing_positive_definite"] is True
    assert algebra["skew_pairing_leading_principal_minors"] == [
        "12",
        "720",
        "172800",
        "207360000",
    ]
    assert receipt["audit"]["registered_unitary_channel_audit"] == {
        "generator": "L = 5*I - A",
        "generator_in_response_word_algebra": True,
        "single_generator_functional_calculus_dimension_upper_bound": 4,
        "order_sensitive_port_indexed_tangent_available": False,
    }


def test_static_rechartings_are_not_forged_into_response_words() -> None:
    _, receipt = source_packet()
    row = receipt["audit"]["recharting_audit"]
    assert row["incidence_automorphism_count"] == 120
    assert row["proper_recharting_count"] == 60
    assert row["proper_rechartings_in_response_word_algebra"] == 1
    assert row["nonidentity_proper_rechartings_in_response_word_algebra"] == 0
    assert row["static_rechartings_are_ordered_response_words"] is False


def test_acceptance_rows_distinguish_partial_missing_and_impossible() -> None:
    _, receipt = source_packet()
    rows = receipt["audit"]["acceptance_classifications"]
    assert rows["ordered_two_sided_port_response_histories"].startswith("PARTIAL_")
    assert rows["twelve_independent_skew_adjoint_generators"].startswith("IMPOSSIBLE_")
    assert rows["exact_nonabelian_commutator_reconstruction"].startswith("IMPOSSIBLE_")
    assert rows["closed_overlap_words_cover_sixty_rechartings"].startswith("MISSING_")
    assert rows["same_word_projective_implementers"] == "MISSING"
    assert rows["target_firewall"] == "SOURCE_NATIVE_ATTAINED"


def test_trace_order_mutation_fails_closed() -> None:
    projection, _ = source_packet()
    mutant = copy.deepcopy(projection)
    mutant["registered_recurrence"]["steps"][1]["step"] = 2
    body = {key: value for key, value in mutant.items() if key != "projection_sha256"}
    mutant["projection_sha256"] = producer.canonical_sha256(body)
    with pytest.raises(producer.CapabilityError, match="TRACE_ORDER"):
        producer.audit_projection(mutant)


def test_non_source_matrix_mutation_fails_closed() -> None:
    projection, _ = source_packet()
    mutant = copy.deepcopy(projection)
    mutant["registered_recurrence"]["steps"][1]["matrix"][0][0] = 1
    body = {key: value for key, value in mutant.items() if key != "projection_sha256"}
    mutant["projection_sha256"] = producer.canonical_sha256(body)
    with pytest.raises(producer.CapabilityError, match="TRACE_REPLAY"):
        producer.audit_projection(mutant)


def test_receipt_boolean_and_hash_mutations_fail_independent_verification(
    tmp_path: Path,
) -> None:
    projection, receipt = source_packet()
    projection_path = tmp_path / "projection.json"
    receipt_path = tmp_path / "receipt.json"
    projection_path.write_text(json.dumps(projection), encoding="utf-8")
    mutant = copy.deepcopy(receipt)
    mutant["physical_current_source_bridge_attained"] = True
    body = {key: value for key, value in mutant.items() if key != "receipt_sha256"}
    mutant["receipt_sha256"] = producer.canonical_sha256(body)
    receipt_path.write_text(json.dumps(mutant), encoding="utf-8")
    with pytest.raises(independent.VerificationError, match="false source gate"):
        independent.verify(projection_path, receipt_path)


@pytest.mark.parametrize(
    ("mutator", "message"),
    [
        (
            lambda receipt: receipt["audit"]["response_word_algebra"].__setitem__(
                "minimal_polynomial", "A^4 = 0"
            ),
            "reported minimal polynomial",
        ),
        (
            lambda receipt: receipt["audit"]["response_word_algebra"][
                "skew_pairing_gram"
            ][0].__setitem__(0, "13"),
            "reported trace-pairing Gram",
        ),
        (
            lambda receipt: receipt["audit"]["refinement_audit"].__setitem__(
                "recurrence_natural_on_registered_maps", False
            ),
            "reported refinement audit",
        ),
        (
            lambda receipt: receipt["audit"]["acceptance_classifications"].__setitem__(
                "same_word_projective_implementers", "SOURCE_NATIVE_ATTAINED"
            ),
            "acceptance classifications",
        ),
    ],
)
def test_independent_verifier_recomputes_reported_algebra_and_status_fields(
    tmp_path: Path,
    mutator,
    message: str,
) -> None:
    projection, receipt = source_packet()
    projection_path = tmp_path / "projection.json"
    receipt_path = tmp_path / "receipt.json"
    projection_path.write_text(json.dumps(projection), encoding="utf-8")
    mutant = copy.deepcopy(receipt)
    mutator(mutant)
    body = {key: value for key, value in mutant.items() if key != "receipt_sha256"}
    mutant["receipt_sha256"] = producer.canonical_sha256(body)
    receipt_path.write_text(json.dumps(mutant), encoding="utf-8")
    with pytest.raises(independent.VerificationError, match=message):
        independent.verify(projection_path, receipt_path)


def test_producer_does_not_import_the_conditional_current_fixture() -> None:
    source = Path(producer.__file__).read_text(encoding="utf-8")
    assert "port_current_inner_certificate" not in source
    assert "port_current_response_reference" not in source
    assert "charged_double_triplet" not in source
