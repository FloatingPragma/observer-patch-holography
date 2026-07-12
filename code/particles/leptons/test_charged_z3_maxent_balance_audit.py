"""Tests for the Z3 MaxEnt/Koide conditional theorem and countermodel."""

from __future__ import annotations

import math

import derive_charged_z3_maxent_balance_audit as lane


def test_natural_trace_maxent_does_not_emit_koide():
    ratio = lane.effective_ratio(1.0, 2.0, 0.0)
    assert ratio == 2.0
    assert lane.koide_from_power_ratio(ratio) == 1.0


def test_one_bit_compensation_emits_koide_conditionally():
    ratio = lane.effective_ratio(1.0, 2.0, math.log(2.0))
    assert math.isclose(ratio, 1.0, rel_tol=0.0, abs_tol=1.0e-15)
    assert math.isclose(
        lane.koide_from_power_ratio(ratio), 2.0 / 3.0, rel_tol=0.0, abs_tol=1.0e-15
    )


def test_current_oph_input_remains_open():
    artifact = lane.build_artifact()
    assert artifact["public_koide_promotion_allowed"] is False
    assert artifact["current_oph_boundary"]["one_bit_degeneracy_compensation_derived"] is False
    assert artifact["current_oph_boundary"]["maxent_probability_to_yukawa_power_bridge_derived"] is False
