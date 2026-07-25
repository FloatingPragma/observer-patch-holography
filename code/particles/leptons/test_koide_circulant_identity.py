"""Tests for the exact C3-circulant Koide identity receipt."""

from __future__ import annotations

import json
from decimal import Decimal

import derive_koide_circulant_identity as lane


def test_exact_balance_equivalence() -> None:
    balance = Decimal(1) / Decimal(2).sqrt()
    assert (
        abs(lane.koide_from_modulus_ratio(balance) - Decimal(2) / Decimal(3))
        < Decimal("1e-27")
    )
    assert (
        abs(
            lane.modulus_ratio_from_koide(Decimal(2) / Decimal(3))
            - balance
        )
        < Decimal("1e-27")
    )


def test_phase_independence_has_a_positive_chamber_boundary() -> None:
    balance = 2.0**-0.5
    for phase in (0.0, 0.1, 2.0 / 9.0):
        values = lane.roots(1.0, balance, phase)
        assert min(values) >= 0.0
        assert abs(lane.physical_koide(values) - 2.0 / 3.0) < 1.0e-14
    outside = lane.roots(1.0, balance, 0.4)
    assert min(outside) < 0.0
    assert abs(lane.physical_koide(outside) - 2.0 / 3.0) > 1.0e-3


def test_artifact_separates_identity_model_and_comparison() -> None:
    artifact = lane.build_artifact()
    assert artifact["checks_pass"] is True
    assert artifact["source_only_physical_prediction"] is False
    assert artifact["public_physical_promotion_allowed"] is False
    assert (
        artifact["identity"]["equivalence"]
        == "Q = 2/3 iff |b|/a = 1/sqrt(2)"
    )
    assert (
        artifact["compare_only_PDG_2026_central_coordinate"][
            "significance_claimed"
        ]
        is False
    )
    assert artifact["compare_only_PDG_2026_central_coordinate"][
        "masses_MeV"
    ] == lane.EXPECTED_PDG_2026_MASSES_MEV
    assert Decimal(
        artifact["compare_only_PDG_2026_central_coordinate"]["Q"]
    ) == lane.EXPECTED_PDG_2026_Q
    assert artifact["checks"]["PDG_2026_mass_coordinate_is_exact"] is True
    assert "historically target-informed" in artifact[
        "current_MCPR_coordinate"
    ]["provenance"]


def test_committed_artifact_is_byte_exact() -> None:
    artifact = lane.build_artifact()
    assert lane.DEFAULT_OUT.read_bytes() == lane.canonical_bytes(artifact)
    stored = json.loads(lane.DEFAULT_OUT.read_text(encoding="utf-8"))
    assert stored["artifact"] == "oph_koide_circulant_identity"
