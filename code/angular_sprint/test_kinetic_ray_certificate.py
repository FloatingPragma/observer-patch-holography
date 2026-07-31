from __future__ import annotations

from fractions import Fraction

import pytest

import kinetic_ray_certificate as kc


def test_receipt_builds_with_expected_status() -> None:
    receipt = kc.build_receipt()
    assert receipt["status"] == (
        "EXACT_PORT_COORDINATE_RAY_ARITHMETIC__"
        "REFERENCE_RAY_DISTINCT_FORM_NOT_REFUTED"
    )
    assert receipt["comparison_boundary"]["comparison_permitted"] is False
    assert receipt["frozen_rg_statistic"]["frozen_before_comparison"] is True
    assert len(receipt["candidate_rays"]) == 3


def test_frozen_beta_column_carries_exact_values() -> None:
    statistic = kc.frozen_rg_statistic({})
    assert statistic["kinetic_column"] == ["10/3", "2", "2"]
    assert statistic["beta_column"] == ["41/6", "-19/6", "-7"]
    assert statistic["exact_cofactors"] == ["-23/3", "37", "-218/9"]
    assert statistic["integer_zero_locus"] == "69 x1 - 333 x2 + 218 x3 = 0"
    assert "computed at scoring" not in str(statistic)


def test_su3_dimension_weighted_average_is_exact() -> None:
    bands = kc.load_pinned_bands()
    decomposition = kc.ideal_decomposition(bands)
    assert decomposition["su3_blocks_disagree"] is True
    assert decomposition["su3_dimension_weighted_average"] == "15/4+1*sqrt5"
    assert "ad-invariant by trace cyclicity" in decomposition["invariance_statement"]
    assert "not an ad-invariant" not in decomposition["invariance_statement"]
    frame = kc.parse_q5(decomposition["su3_block_values"][0])
    quintet = kc.parse_q5(decomposition["su3_block_values"][1])
    average = kc.q5_scale(
        kc.q5_add(kc.q5_scale(frame, Fraction(3)), kc.q5_scale(quintet, Fraction(5))),
        Fraction(1, 8),
    )
    assert kc.q5_str(average) == decomposition["su3_dimension_weighted_average"]


def test_reference_ray_and_commutant_are_excluded() -> None:
    receipt = kc.build_receipt()
    tests = receipt["ray_tests"]
    assert tests["reference_ray_hit"] is False
    assert tests["commutant_relation_hit"] is False
    for row in tests["results"]:
        assert row["proportional_to_reference_5_3_1_1"] is False
        assert row["quadratic_commutant_k1_eq_3k2_minus_2k3"] is False


def test_proportionality_detector_works_on_a_positive_case() -> None:
    scaled = [kc.q5_str(kc.q5(Fraction(5, 3) * 7)), kc.q5_str(kc.q5(7)), kc.q5_str(kc.q5(7))]
    rays = [{"ray_id": "synthetic", "components": scaled}]
    tests = kc.ray_tests(rays)
    assert tests["results"][0]["proportional_to_reference_5_3_1_1"] is True


def test_commutant_detector_works_on_a_positive_case() -> None:
    components = [kc.q5_str(kc.q5(1)), kc.q5_str(kc.q5(1)), kc.q5_str(kc.q5(1))]
    rays = [{"ray_id": "synthetic", "components": components}]
    tests = kc.ray_tests(rays)
    assert tests["results"][0]["quadratic_commutant_k1_eq_3k2_minus_2k3"] is True


def test_committed_receipt_is_byte_exact() -> None:
    committed = kc.RECEIPT_PATH.read_bytes()
    assert committed == kc.canonical_json_bytes(kc.build_receipt())


def test_pinned_band_drift_fails_closed(tmp_path, monkeypatch: pytest.MonkeyPatch) -> None:
    import json

    tampered = json.loads(kc.PORT_RECEIPT_PATH.read_text(encoding="utf-8"))
    tampered["compactness"]["hilbert_schmidt_pullback_band_coefficients"][
        "frame_band"
    ] = "6 + 1*sqrt(5)"
    path = tmp_path / "port_receipt.json"
    path.write_text(json.dumps(tampered), encoding="utf-8")
    monkeypatch.setattr(kc, "PORT_RECEIPT_PATH", path)
    with pytest.raises(kc.KineticError, match="pinned band coefficient drift"):
        kc.build_receipt()
