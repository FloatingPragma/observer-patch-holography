from __future__ import annotations

import copy
from fractions import Fraction
import json
from pathlib import Path
import sys

import pytest


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from finite_born_frame_certificate import (  # noqa: E402
    ANTIPODE,
    CertificateError,
    DEFAULT_FRAME,
    DEFAULT_OUTPUT,
    Q5,
    build_certificate,
    canonical_sha256,
    validate_certificate,
)
from verify_finite_born_frame_independent import (  # noqa: E402
    VerificationError,
    verify,
)


@pytest.fixture(scope="module")
def receipt() -> dict:
    return json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8"))


def rehash(receipt: dict) -> None:
    body = {key: value for key, value in receipt.items() if key != "receipt_sha256"}
    receipt["receipt_sha256"] = canonical_sha256(body)


def test_exact_replay_matches_frozen_receipt(receipt: dict) -> None:
    assert build_certificate() == receipt
    validate_certificate(receipt)


def test_independent_verifier_recomputes_rank_gap(receipt: dict) -> None:
    result = verify(receipt)
    assert result == {
        "status": "PASS",
        "context_affine_dimension": 6,
        "born_affine_dimension": 3,
        "central_affine_dimension": 11,
        "conditional_density_representation_unique": True,
        "all_context_weights_density_representable": False,
        "physical_public_effect_attachment": False,
    }


def test_nonborn_and_nonpositive_controls_are_distinct(receipt: dict) -> None:
    controls = receipt["declared_spinor_projective_branch"]["exact_controls"]
    nonborn = controls["context_additive_nonborn_weight"]
    assert nonborn["all_weights_in_unit_interval"] is True
    assert nonborn["has_trace_one_hermitian_representation"] is False
    assert any(value != "0" for value in nonborn["frame_relation_residuals"])

    nonpositive = controls["born_affine_but_nonpositive_weight"]
    assert nonpositive["all_weights_in_unit_interval"] is True
    assert nonpositive["has_unique_trace_one_hermitian_representation"] is True
    assert nonpositive["has_density_representation"] is False
    valid = controls["valid_density_weight"]
    assert valid["has_unique_density_representation"] is True


def test_central_atoms_do_not_fix_an_ambient_density(receipt: dict) -> None:
    central = receipt["declared_central_atom_branch"]
    assert central["state_on_commutative_algebra"]["unique"] is True
    assert central["ambient_M12_density_representation"]["unique"] is False
    witness = central["ambient_M12_density_representation"]["uniform_weight_counterexample"]
    assert witness["same_declared_atom_weights"] is True
    assert witness["distinct_density_matrices"] is True


def test_false_physical_attachment_is_rejected(receipt: dict) -> None:
    mutant = copy.deepcopy(receipt)
    mutant["source_lineage"]["external_projective_adapter"]["physical_promotion_allowed"] = True
    rehash(mutant)
    with pytest.raises(VerificationError, match="physical promotion"):
        verify(mutant)


def test_false_universal_density_claim_is_rejected(receipt: dict) -> None:
    mutant = copy.deepcopy(receipt)
    mutant["declared_spinor_projective_branch"]["decision"][
        "every_context_admissible_weight_has_a_density_representation"
    ] = True
    rehash(mutant)
    with pytest.raises(VerificationError, match="false density universality"):
        verify(mutant)


def test_tampered_affine_dimension_is_rejected(receipt: dict) -> None:
    mutant = copy.deepcopy(receipt)
    mutant["declared_spinor_projective_branch"]["normalization_and_additivity"]["affine_dimension"] = 3
    rehash(mutant)
    with pytest.raises(VerificationError, match="context dimension"):
        verify(mutant)


def test_tampered_nonpositive_witness_is_rejected(receipt: dict) -> None:
    mutant = copy.deepcopy(receipt)
    mutant["declared_spinor_projective_branch"]["exact_controls"][
        "born_affine_but_nonpositive_weight"
    ]["s"] = ["1/5", "0", "0"]
    rehash(mutant)
    with pytest.raises(VerificationError, match="weights do not match frame"):
        verify(mutant)


def test_wrong_antipode_lineage_fails_closed() -> None:
    mutant = list(ANTIPODE)
    mutant[0], mutant[1] = mutant[1], mutant[0]
    with pytest.raises(CertificateError, match="ANTIPODE"):
        build_certificate(antipode=mutant)


def test_rank_collapsed_frame_fails_closed() -> None:
    positive = [DEFAULT_FRAME[index % 2] for index in range(6)]
    collapsed = tuple(positive) + tuple(
        tuple(-entry for entry in positive[5 - index]) for index in range(6)
    )
    with pytest.raises(CertificateError, match="FRAME_RANK"):
        build_certificate(frame=collapsed)


def test_receipt_hash_detects_unresigned_mutation(receipt: dict) -> None:
    mutant = copy.deepcopy(receipt)
    mutant["closure_assessment"]["positive_born_derivation_obtained"] = True
    with pytest.raises(VerificationError, match="receipt hash"):
        verify(mutant)


def test_q5_sign_controls_are_exact() -> None:
    assert Q5.of(Fraction(3, 5), Fraction(3, 10)).sign() > 0
    assert Q5.of(1) - Q5.of(Fraction(3, 10), Fraction(3, 10)) != Q5.of(0)
    assert Q5.of(-1, 1).sign() > 0
