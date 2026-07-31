from __future__ import annotations

import json
from fractions import Fraction

import pytest

import spin_six_primitive_port_prediction as prediction


def test_prediction_is_frozen_and_comparison_is_unarmed() -> None:
    receipt = prediction.build_receipt()
    assert receipt["status"] == prediction.STATUS
    assert receipt["prediction_scope"]["type"] == (
        "prospective conditional physical-branch prediction"
    )
    boundary = receipt["exposure_and_custody_boundary"]
    assert boundary["new_comparison_data_read"] is False
    assert boundary["comparison_permitted"] is False
    assert boundary["comparison_state"].startswith("UNARMED")
    assert "WMAP" in boundary["prior_related_exposure"]


def test_scale_free_relations_are_exact() -> None:
    receipt = prediction.build_receipt()
    relations = receipt["exact_prediction"]["scale_free_relations"]
    assert Fraction(relations["B6_over_C4_squared"]) == Fraction(32, 315)
    assert Fraction(relations["B0_over_C4_squared"]) == Fraction(10, 21)
    assert Fraction(relations["B6_over_B0"]) == Fraction(16, 75)
    assert receipt["exact_prediction"]["fit_freedom_after_C4"].startswith(
        "one orientation"
    )


def test_failure_scope_does_not_overclaim_oph_wide_falsification() -> None:
    receipt = prediction.build_receipt()
    scope = receipt["prospective_decision_rule"]["scope_of_failure"]
    assert "primitive twelve-port physical propagation branch" in scope
    assert "only if issue #655 proves" in scope


def test_committed_prediction_receipt_is_byte_exact() -> None:
    committed = prediction.RECEIPT_PATH.read_bytes()
    rebuilt = prediction.build_receipt()
    assert committed == prediction.base.canonical_json_bytes(rebuilt)


@pytest.mark.parametrize("field", ["schema", "status", "receipt_sha256"])
def test_typed_parent_mutation_fails_closed(
    monkeypatch: pytest.MonkeyPatch, tmp_path, field: str
) -> None:
    parent = json.loads(prediction.universality.RECEIPT_PATH.read_text())
    parent[field] = "tampered"
    path = tmp_path / "universality.json"
    path.write_bytes(prediction.base.canonical_json_bytes(parent))
    monkeypatch.setattr(prediction.universality, "RECEIPT_PATH", path)
    with pytest.raises(prediction.base.FingerprintError):
        prediction.build_receipt()
