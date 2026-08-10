"""Exact positive and mutation controls for the B13 phase-lift verifier."""

from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

import pytest

from verify_source_phase_lift import (
    CIDENTITY,
    CONE,
    CZERO,
    Q3,
    VerificationError,
    cconj_transpose,
    cmul,
    complexify,
    decode_qmatrix,
    default_paths,
    flattened_pairing_row,
    phase_lift,
    qmul,
    qsub,
    rank,
    verify_payload,
)


PAYLOAD_PATH, SOURCE_ROOT = default_paths()


def load_payload() -> dict:
    return json.loads(PAYLOAD_PATH.read_text())


def source_commutator(payload: dict):
    diagonal = decode_qmatrix(
        payload["context_web"]["diagonal_context"]["projectors"][0]
    )
    rotated = decode_qmatrix(
        payload["context_web"]["conjugated_contexts"][2]["projectors"][0]
    )
    return qsub(qmul(rotated, diagonal), qmul(diagonal, rotated))


def test_committed_payload_has_exact_rank_jump_and_no_instrument_receipt() -> None:
    report = verify_payload(load_payload(), SOURCE_ROOT)
    assert report["source_hashes_verified"] is True
    assert report["native_operator_span_rank"] == 3
    assert report["phase_lifted_operator_span_rank"] == 4
    assert report["phase_lift_equals_plus_y_projector"] is True
    assert report["rotated_outcome_receipt_present"] is False
    assert report["phase_lift_instrument_receipt_present"] is False


def test_zero_phase_coefficient_does_not_close_tomography() -> None:
    payload = load_payload()
    commutator = source_commutator(payload)
    diagonal = decode_qmatrix(
        payload["context_web"]["diagonal_context"]["projectors"][0]
    )
    rotated = decode_qmatrix(
        payload["context_web"]["conjugated_contexts"][2]["projectors"][0]
    )
    zero_lift = phase_lift(commutator, coefficient=Q3())
    assert (
        rank(
            [
                flattened_pairing_row(CIDENTITY),
                flattened_pairing_row(complexify(diagonal)),
                flattened_pairing_row(complexify(rotated)),
                flattened_pairing_row(zero_lift),
            ]
        )
        == 3
    )


def test_omitting_i_makes_the_lift_nonhermitian() -> None:
    real_commutator_lift = phase_lift(
        source_commutator(load_payload()), include_i=False
    )
    assert cconj_transpose(real_commutator_lift) != real_commutator_lift
    assert cmul(real_commutator_lift, real_commutator_lift) != real_commutator_lift


def test_mutated_irrep_entry_fails_homomorphism_or_orthogonality() -> None:
    payload = load_payload()
    payload["irrep"]["elements"][3]["matrix"][0][0] = ["-2/3", "0"]
    with pytest.raises(VerificationError):
        verify_payload(payload, verify_hashes=False)


def test_mutated_projector_entry_fails_exact_orbit_check() -> None:
    payload = load_payload()
    payload["context_web"]["conjugated_contexts"][2]["projectors"][0][0][0] = [
        "1/3",
        "0",
    ]
    with pytest.raises(VerificationError):
        verify_payload(payload, verify_hashes=False)


def test_fabricated_rotated_counts_fail_closed() -> None:
    payload = load_payload()
    payload["context_web"]["conjugated_contexts"][2]["realized_outcome_counts"] = [
        90,
        89,
    ]
    with pytest.raises(VerificationError, match="not source-realized"):
        verify_payload(payload, verify_hashes=False)


def test_mutated_frequency_boundary_fails_closed() -> None:
    payload = load_payload()
    payload["outcome_frequency_boundary"]["contexts_with_realized_frequencies"].append(
        "conjugated_3"
    )
    with pytest.raises(VerificationError, match="realized-frequency boundary"):
        verify_payload(payload, verify_hashes=False)


def test_mutated_coincidence_class_fails_closed() -> None:
    payload = load_payload()
    payload["context_web"]["context_coincidence_classes"][0].append("conjugated_4")
    with pytest.raises(VerificationError, match="coincidence-class"):
        verify_payload(payload, verify_hashes=False)


def test_mutated_source_hash_fails_custody_check() -> None:
    payload = load_payload()
    payload["provenance"]["gauge_state_sha256"] = "0" * 64
    with pytest.raises(VerificationError, match="provenance hash drift"):
        verify_payload(payload, SOURCE_ROOT)


def test_mutated_source_commit_fails_custody_check() -> None:
    payload = load_payload()
    payload["provenance"]["run_git_commit"] = "0" * 40
    with pytest.raises(VerificationError, match="source commit drift"):
        verify_payload(payload, verify_hashes=False)


def test_wrong_phase_normalization_is_not_a_projector() -> None:
    wrong = phase_lift(source_commutator(load_payload()), coefficient=Q3(Fraction(1)))
    assert cmul(wrong, wrong) != wrong
    assert wrong != CIDENTITY
    assert CONE != CZERO


def test_payload_path_is_outside_the_theorem_repository() -> None:
    """Custody control: the verifier reads the sibling source artifact."""

    theorem_repo = Path(__file__).resolve().parents[2]
    assert theorem_repo not in PAYLOAD_PATH.parents
    assert PAYLOAD_PATH.is_file()
