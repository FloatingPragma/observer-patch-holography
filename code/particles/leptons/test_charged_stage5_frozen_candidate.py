"""Tests for the frozen reference-free Stage-5 candidate builder."""

from __future__ import annotations

import json
from decimal import Decimal, localcontext

import derive_charged_stage5_frozen_candidate as lane


def _artifact():
    return lane.build_artifact(json.loads(lane.D10.read_text(encoding="utf-8")))


def test_builder_emits_all_three_positive_ordered_masses():
    rows = _artifact()["candidate_mass_rows"]
    values = [Decimal(rows[name]["value_gev"]) for name in ("electron", "muon", "tau")]
    assert Decimal(0) < values[0] < values[1] < values[2]


def test_determinant_law_is_satisfied():
    artifact = _artifact()
    rows = artifact["candidate_mass_rows"]
    product = Decimal(1)
    for name in ("electron", "muon", "tau"):
        product *= Decimal(rows[name]["value_gev"])
    v = Decimal(artifact["inputs"]["v_from_source_transmutation_gev"])
    expected = v**3 / (Decimal(2) * Decimal(6) ** 14)
    assert abs(product / expected - Decimal(1)) < Decimal("1e-70")


def test_candidate_remains_unpromoted_and_frozen():
    artifact = _artifact()
    assert artifact["reference_data_consumed_by_builder"] is False
    assert artifact["public_prediction_promotion_allowed"] is False
    assert artifact["status"] == "FROZEN_ACCURATE_NUMERICAL_CANDIDATE_NOT_THEOREM_GRADE"
    assert artifact["inputs"]["frozen_pixel_P"] == "1.63094"
    assert artifact["branch_and_scheme_boundary"]["current_public_pixel_used"] is False
    assert artifact["branch_and_scheme_boundary"]["pole_mass_comparison_allowed_here"] is False


def test_frozen_mass_values_cannot_drift_silently():
    rows = _artifact()["candidate_mass_rows"]
    expected = {
        "electron": Decimal("0.000510882243294968223284389471842535"),
        "muon": Decimal("0.105635282870957389616826180472195"),
        "tau": Decimal("1.776579124017314908083042601852472"),
    }
    with localcontext() as context:
        context.prec = 80
        for name, reference in expected.items():
            value = Decimal(rows[name]["value_gev"])
            assert abs(value / reference - Decimal(1)) < Decimal("1e-33")
