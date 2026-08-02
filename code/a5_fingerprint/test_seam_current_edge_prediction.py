from __future__ import annotations

import json
from fractions import Fraction

import pytest

import seam_current_edge_prediction as prediction


def test_source_native_edge_ray_is_exact_and_physical_producer_stays_open() -> None:
    receipt = prediction.build_receipt()
    assert receipt["status"] == prediction.STATUS
    assert receipt["producer_scope"]["physical_producer_closed"] is False
    assert receipt["producer_scope"]["frozen_prediction_registered"] is False
    gates = receipt["promotion_gates"]
    assert gates["all_discharged"] is False
    assert gates["physical_producer_closed"] is False
    assert gates["comparison_eligible"] is False
    assert len(gates["gates"]) == 9
    assert {row["status"] for row in gates["gates"]} == {"OPEN"}
    scale_gates = [row for row in gates["gates"] if row["owner"] == "#664"]
    assert {row["gate"] for row in scale_gates} == {
        "finite physical carrier scale",
        "source-derived positive physical scale lower bound",
    }


def test_exact_seam_geometry_binds_to_the_edge_orbit() -> None:
    source = prediction.build_receipt()["exact_source_result"]
    geometry = source["finite_geometry_replay"]
    assert geometry["source_ports"] == 12
    assert geometry["source_seams"] == 30
    assert geometry["port_degree"] == 5
    assert geometry["unoriented_axes"] == 15
    assert geometry["signed_edge_directions"] == 30
    assert geometry["directed_seam_labels"] == 60
    assert geometry["directed_labels_per_signed_direction"] == 2
    assert geometry["even_moments_on_unit_sphere"] == {
        "sum_w_dot_n_squared": "10",
        "sum_w_dot_n_fourth": "6",
        "sum_w_dot_n_sixth": "30/7 - (2/7) I6(n)",
    }
    assert source["seam_current_image"].startswith("D6")
    assert "do not identify" in source["scope_boundary"]


def test_lambda_coefficients_and_relations_are_recomputed_exactly() -> None:
    candidate = prediction.build_receipt()["conditional_physical_candidate"]
    assert "Λ_a" in candidate["operator"]
    assert "omega" not in candidate["operator"].lower()
    coefficients = candidate["coefficients"]
    assert Fraction(coefficients["C4_over_a2"]) == Fraction(-1, 20)
    assert Fraction(coefficients["B0_over_a4"]) == Fraction(1, 840)
    assert Fraction(coefficients["B6_over_a4"]) == Fraction(-1, 12600)
    relations = candidate["scale_free_relations"]
    assert Fraction(relations["B0_over_C4_squared"]) == Fraction(10, 21)
    assert Fraction(relations["B6_over_C4_squared"]) == Fraction(-2, 63)
    assert Fraction(relations["B6_over_B0"]) == Fraction(-1, 15)


def test_target_and_comparison_data_are_unread_and_unarmed() -> None:
    receipt = prediction.build_receipt()
    boundary = receipt["exposure_and_custody_boundary"]
    assert boundary["target_values_read"] is False
    assert boundary["comparison_data_read"] is False
    assert boundary["public_measurement_read"] is False
    assert boundary["comparison_permitted"] is False
    assert boundary["comparison_inputs"] == []
    assert boundary["comparison_state"].startswith("INELIGIBLE_UNARMED")
    assert receipt["prospective_decision_rule"]["scope_of_failure"].endswith(
        "branch forced and exclusive"
    )
    assert "a>0 by itself is insufficient" in receipt["prospective_decision_rule"][
        "no_null_verdict"
    ]


def test_fz11_is_a_distinct_untouched_branch() -> None:
    separation = prediction.build_receipt()["fz11_separation"]
    assert separation["fz11_register_id"] == "FZ-11"
    assert separation["fz11_prediction_receipt_read"] is False
    assert separation["fz11_bytes_modified"] is False
    assert separation["supersedes_fz11"] is False
    assert separation["vertex_B6_over_C4_squared"] == "32/315"
    assert separation["edge_B6_over_C4_squared"] == "-2/63"
    assert separation["opposite_rank_six_sign"] is True


def test_committed_receipt_is_canonical_byte_exact_and_self_hashed() -> None:
    committed = prediction.verify_committed_receipt()
    assert prediction.RECEIPT_PATH.read_bytes() == prediction.base.canonical_json_bytes(
        committed
    )


@pytest.mark.parametrize("field", ["schema", "status", "receipt_sha256"])
def test_typed_geometry_parent_mutation_fails_closed(
    monkeypatch: pytest.MonkeyPatch, tmp_path, field: str
) -> None:
    parent = json.loads(prediction.base.RECEIPT_PATH.read_text())
    parent[field] = "tampered"
    path = tmp_path / "geometry_parent.json"
    path.write_bytes(prediction.base.canonical_json_bytes(parent))
    monkeypatch.setattr(prediction.base, "RECEIPT_PATH", path)
    with pytest.raises(prediction.base.FingerprintError):
        prediction.build_receipt()


def test_claim_bearing_seam_source_mutation_fails_closed(
    monkeypatch: pytest.MonkeyPatch, tmp_path
) -> None:
    source = prediction.SEAM_PROOF_PATH.read_text(encoding="utf-8")
    source = source.replace(
        "theorem exists_seamCurrent_iff_even",
        "theorem mutated_exists_seamCurrent_iff_even",
        1,
    )
    path = tmp_path / "SeamCurrentCarrierQuotient.lean"
    path.write_text(source, encoding="utf-8")
    monkeypatch.setattr(prediction, "SEAM_PROOF_PATH", path)
    with pytest.raises(
        prediction.base.FingerprintError, match="claim-bearing source drift"
    ):
        prediction.build_receipt()


def test_claim_bearing_orbit_source_mutation_fails_closed(
    monkeypatch: pytest.MonkeyPatch, tmp_path
) -> None:
    source = prediction.ORBIT_PROOF_PATH.read_text(encoding="utf-8")
    source = source.replace("| .edge30 => -1 / 12600", "| .edge30 => 1 / 12600", 1)
    path = tmp_path / "A5OrbitRaySeparation.lean"
    path.write_text(source, encoding="utf-8")
    monkeypatch.setattr(prediction, "ORBIT_PROOF_PATH", path)
    with pytest.raises(
        prediction.base.FingerprintError, match="claim-bearing source drift"
    ):
        prediction.build_receipt()


def test_claim_bearing_seam_moment_source_mutation_fails_closed(
    monkeypatch: pytest.MonkeyPatch, tmp_path
) -> None:
    source = prediction.SEAM_MOMENT_PROOF_PATH.read_text(encoding="utf-8")
    source = source.replace(
        "theorem seamMoment6_eq", "theorem mutated_seamMoment6_eq", 1
    )
    path = tmp_path / "SeamCurrentEdge30Moment.lean"
    path.write_text(source, encoding="utf-8")
    monkeypatch.setattr(prediction, "SEAM_MOMENT_PROOF_PATH", path)
    with pytest.raises(
        prediction.base.FingerprintError, match="claim-bearing source drift"
    ):
        prediction.build_receipt()


def test_nonclaim_source_byte_mutation_breaks_committed_pin(
    monkeypatch: pytest.MonkeyPatch, tmp_path
) -> None:
    source = prediction.SEAM_PROOF_PATH.read_text(encoding="utf-8") + "\n/- mutation -/\n"
    path = tmp_path / "SeamCurrentCarrierQuotient.lean"
    path.write_text(source, encoding="utf-8")
    monkeypatch.setattr(prediction, "SEAM_PROOF_PATH", path)
    with pytest.raises(
        prediction.base.FingerprintError, match="edge receipt parent/source drift"
    ):
        prediction.verify_committed_receipt()


@pytest.mark.parametrize("field", ["schema", "status", "receipt_sha256"])
def test_committed_receipt_type_or_digest_mutation_fails_closed(
    monkeypatch: pytest.MonkeyPatch, tmp_path, field: str
) -> None:
    committed = json.loads(prediction.RECEIPT_PATH.read_text())
    committed[field] = "tampered"
    path = tmp_path / "edge_receipt.json"
    path.write_bytes(prediction.base.canonical_json_bytes(committed))
    monkeypatch.setattr(prediction, "RECEIPT_PATH", path)
    with pytest.raises(prediction.base.FingerprintError):
        prediction.verify_committed_receipt()


def test_geometry_mutation_fails_before_a_candidate_is_emitted(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    vertices = prediction.base.cartesian_vertices()
    mutated = list(vertices)
    first = mutated[0]
    mutated[0] = (
        prediction.base.q5_add(first[0], prediction.base.q5(Fraction(1, 100))),
        first[1],
        first[2],
    )
    monkeypatch.setattr(prediction.base, "cartesian_vertices", lambda: mutated)
    with pytest.raises(prediction.base.FingerprintError):
        prediction.build_receipt()
