"""Adversarial tests for the local face-Maxwell producer and replay."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

import local_face_maxwell_action as producer  # noqa: E402
import verify_local_face_maxwell_action_independent as verifier  # noqa: E402


def _resign(receipt: dict) -> None:
    body = {k: v for k, v in receipt.items() if k != "receipt_sha256"}
    receipt["receipt_sha256"] = "sha256:" + hashlib.sha256(
        verifier.canonical(body)).hexdigest()


def _set_nested(receipt: dict, path: tuple[str, ...], value: object) -> None:
    target = receipt
    for key in path[:-1]:
        target = target[key]
    target[path[-1]] = value


def test_producer_matches_committed_receipt() -> None:
    assert producer.verify_committed()["issue"] == 733
    assert json.loads(producer.RECEIPT.read_bytes()) == producer.build_receipt()


def test_separate_independent_replay() -> None:
    assert verifier.main([]) == 0


def test_exact_locality_and_rank_receipt() -> None:
    receipt = producer.build_receipt()
    assert receipt["incidence"]["rank_C"] == 19
    assert receipt["incidence"]["kernel_C_dimension"] == 11
    assert receipt["incidence"]["C_support"] == {"nonzero": 60, "total": 600}
    assert receipt["local_hessian"]["support"]["nonzero"] == 150
    assert receipt["local_hessian"]["support"]["nonzero_per_row"] == [5] * 30
    assert receipt["local_hessian"]["diagonal"] == [2] * 30


def test_orientation_erasure_breaks_boundary_of_boundary() -> None:
    texts = producer.load_sources()
    left = producer.parse_vector(texts["SeamCurrentCarrierQuotient.lean"],
                                 "seamLeft", 30)
    right = producer.parse_vector(texts["SeamCurrentCarrierQuotient.lean"],
                                  "seamRight", 30)
    faces = producer.parse_faces(texts["CoreAxioms.lean"])
    b, c = producer.build_incidence(left, right, faces)
    unsigned = [[abs(x) for x in row] for row in c]
    product = producer.mul(unsigned, producer.transpose(b))
    assert product[0][0] == -2
    assert any(x != 0 for row in product for x in row)


def test_exact_spectra_are_recomputed() -> None:
    spectra = producer.build_receipt()["exact_spectra"]
    assert spectra["C_C_transpose"]["quadratic_polynomial"] == "x^2-6x+4"
    assert spectra["C_C_transpose"]["quadratic_sector_dimension"] == 6
    assert spectra["B_B_transpose"]["quadratic_polynomial"] == "x^2-10x+20"
    assert spectra["B_B_transpose"]["quadratic_sector_dimension"] == 6


def test_resigned_semantic_mutation_fails_independent_replay(tmp_path: Path) -> None:
    mutated = json.loads(producer.RECEIPT.read_bytes())
    mutated["local_hessian"]["support"]["nonzero"] = 151
    _resign(mutated)
    path = tmp_path / "mutated.json"
    path.write_bytes(verifier.canonical(mutated))
    with pytest.raises(verifier.VerificationError):
        verifier.verify(path)


@pytest.mark.parametrize(
    ("field_path", "mutated_value"),
    [
        (("issue",), 999),
        (("exact_spectra", "C_C_transpose", "quadratic_polynomial"),
         "x^2+123x+456"),
        (("exact_spectra", "C_C_transpose", "quadratic_sector_dimension"), 99),
        (("exact_spectra", "C_C_transpose", "rational_eigenspace_dimensions", "2"),
         99),
        (("adversarial_controls", "orientation_erasure_rejected"), False),
        (("adversarial_controls", "unsigned_face_zero_port_boundary_value"), 777),
        (("scalar_green_join", "typed_separately_from_seam_current_source"), False),
        (("scalar_green_join", "charge_type"), "Fin 30 -> rational"),
        (("handoff_interface", "schema"), "oph.local_face_maxwell_action.handoff.v999"),
        (("handoff_interface", "source_types"),
         {"charge": "same", "seam_current": "same"}),
        (("handoff_interface", "statement"), "physical comparison armed"),
    ],
)
def test_resigned_previously_unchecked_semantics_fail_independent_replay(
    tmp_path: Path, field_path: tuple[str, ...], mutated_value: object,
) -> None:
    mutated = json.loads(producer.RECEIPT.read_bytes())
    _set_nested(mutated, field_path, mutated_value)
    _resign(mutated)
    path = tmp_path / ("mutated_" + "_".join(field_path) + ".json")
    path.write_bytes(verifier.canonical(mutated))
    with pytest.raises(verifier.VerificationError):
        verifier.verify(path)


def test_source_pin_drift_fails_closed(tmp_path: Path) -> None:
    source = producer.SOURCES[0][0]
    tampered = tmp_path / source.name
    raw = source.read_bytes().replace(b"![0, 0", b"![1, 0", 1)
    tampered.write_bytes(raw)
    with pytest.raises(verifier.VerificationError):
        verifier.pinned_text(str(tampered), len(raw), producer.SOURCES[0][2])


def test_handoff_is_design_only_and_sources_stay_typed() -> None:
    receipt = producer.build_receipt()
    handoff = receipt["handoff_interface"]
    assert handoff["consumer"] == "instrument lane (issue 737)"
    assert handoff["design_only"] is True
    assert handoff["frozen"] is False
    assert handoff["comparison_permitted"] is False
    assert handoff["source_types"] == {
        "charge": "Fin 12 -> real",
        "seam_current": "Fin 30 -> real",
    }
    assert receipt["physical_boundary"]["open_premise_rows"] == ["PR-53", "PR-54"]
    assert not any(
        receipt["physical_boundary"][key]
        for key in (
            "position_attachment_proved",
            "physical_source_identification_proved",
            "laboratory_readout_proved",
            "continuum_limit_proved",
        )
    )
