"""Exact positive and mutation controls for the PR-04 phase-operation receipt.

Three mutation families are guarded fail-closed: a wrong committed literal in
the pinned payload, a wrong declared operation (coefficient or phase drift),
and a tampered count in the receipt.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction

import pytest

import produce_phase_operation_receipt as producer
import verify_phase_operation_receipt as verifier


def load_payload() -> dict:
    return json.loads(producer.PAYLOAD_PATH.read_text())


def payload_sha256() -> str:
    return hashlib.sha256(producer.PAYLOAD_PATH.read_bytes()).hexdigest()


def load_receipt() -> dict:
    return json.loads(producer.RECEIPT_PATH.read_text())


def test_producer_is_deterministic_and_matches_committed_receipt() -> None:
    once = producer.produce()
    twice = producer.produce()
    assert once == twice
    assert once == producer.RECEIPT_PATH.read_bytes()


def test_verifier_passes_on_committed_artifacts() -> None:
    report = verifier.verify(load_receipt(), load_payload(), payload_sha256())
    assert report["contexts_verified"] == 8
    assert report["phase_counts"] == [179, 179]
    assert report["phase_window_lhs"] == 0
    assert report["deterministic_semantics_declared"] is True


def test_receipt_masses_are_positive_run_mass_multiples() -> None:
    for row in load_receipt()["contexts"]:
        assert row["mass"] == row["run_mass_multiple"] * 179
        assert row["run_mass_multiple"] >= 1
        assert row["counts"][0] + row["counts"][1] == row["mass"]


def test_wrong_run_literal_fails_producer() -> None:
    payload = load_payload()
    payload["context_web"]["diagonal_context"]["realized_outcome_counts"] = [112, 67]
    with pytest.raises(producer.ProductionError, match="diagonal counts drift"):
        producer.build_receipt(payload, producer.EXPECTED_PAYLOAD_SHA256)


def test_wrong_run_literal_fails_verifier() -> None:
    payload = load_payload()
    payload["context_web"]["diagonal_context"]["realized_outcome_counts"] = [112, 67]
    with pytest.raises(verifier.VerificationError, match="diagonal counts drift"):
        verifier.verify(load_receipt(), payload, payload_sha256())


def test_wrong_projector_literal_fails_verifier() -> None:
    payload = load_payload()
    payload["context_web"]["conjugated_contexts"][2]["projectors"][0][0][0] = [
        "1/3",
        "0",
    ]
    with pytest.raises(verifier.VerificationError, match="conjugated projector"):
        verifier.verify(load_receipt(), payload, payload_sha256())


def test_wrong_operation_coefficient_is_rejected() -> None:
    payload = load_payload()
    reps = [
        producer.decode_payload_matrix(row["matrix"])
        for row in payload["irrep"]["elements"]
    ]
    projector = producer.mmul(
        producer.mmul(reps[3], producer.RECORD), producer.mtranspose(reps[3])
    )
    commutator = producer.msub(
        producer.mmul(projector, producer.RECORD),
        producer.mmul(producer.RECORD, projector),
    )
    wrong = producer.phase_lift(commutator, coefficient=producer.R3.of(Fraction(1)))
    assert producer.mmul(wrong, wrong) != wrong


def test_wrong_operation_matrix_fails_verifier() -> None:
    receipt = load_receipt()
    receipt["declared_operation"]["matrix"][0][1] = ["0", "0", "1/2", "0"]
    with pytest.raises(verifier.VerificationError, match="operation matrix drift"):
        verifier.verify(receipt, load_payload(), payload_sha256())


def test_tampered_phase_count_fails_verifier() -> None:
    receipt = load_receipt()
    receipt["contexts"][-1]["counts"] = [180, 178]
    with pytest.raises(verifier.VerificationError, match="exact scaled Born weight"):
        verifier.verify(receipt, load_payload(), payload_sha256())


def test_tampered_rotated_count_fails_verifier() -> None:
    receipt = load_receipt()
    receipt["contexts"][3]["counts"] = [316, 400]
    with pytest.raises(verifier.VerificationError, match="exact scaled Born weight"):
        verifier.verify(receipt, load_payload(), payload_sha256())


def test_tampered_mass_multiple_fails_minimality() -> None:
    receipt = load_receipt()
    row = receipt["contexts"][0]
    row["run_mass_multiple"] = 2
    row["mass"] = 358
    row["counts"] = [222, 136]
    with pytest.raises(verifier.VerificationError, match="not minimal"):
        verifier.verify(receipt, load_payload(), payload_sha256())


def test_window_check_rejects_out_of_window_counts() -> None:
    with pytest.raises(verifier.VerificationError, match="receipt window"):
        verifier.check_window([179, 1])


def test_window_check_accepts_committed_counts() -> None:
    assert verifier.check_window([179, 179]) == (0, 3869527488)


def test_stripped_semantics_declaration_fails_verifier() -> None:
    receipt = load_receipt()
    receipt["semantics"] = "counts observed in the laboratory"
    with pytest.raises(verifier.VerificationError, match="semantics declaration"):
        verifier.verify(receipt, load_payload(), payload_sha256())


def test_wrong_decision_disposition_fails_verifier() -> None:
    receipt = load_receipt()
    receipt["recorded_decision"]["disposition"] = "remove"
    with pytest.raises(verifier.VerificationError, match="disposition drift"):
        verifier.verify(receipt, load_payload(), payload_sha256())


def test_payload_hash_mismatch_fails_verifier() -> None:
    with pytest.raises(verifier.VerificationError, match="hash mismatch"):
        verifier.verify(load_receipt(), load_payload(), "0" * 64)
