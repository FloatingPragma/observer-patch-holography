"""Tests for the invariant 12/24 charged-family no-go theorem."""

from __future__ import annotations

import derive_charged_12_24_singlet_no_go as lane


def test_centering_kills_invariant_family_image():
    assert lane.centered([7.0, 7.0, 7.0]) == [0.0, 0.0, 0.0]


def test_live_12_24_lock_cannot_fill_shape_scalars():
    artifact = lane.build_artifact(lane._load(lane.HIERARCHY), lane._load(lane.SCALARS))
    assert artifact["hierarchy_inputs"]["screen_ports"] == 12
    assert artifact["hierarchy_inputs"]["repair_register_slots"] == 24
    forced = artifact["theorem"]["forced_shape_on_invariant_branch"]
    assert forced["eta_source_support_extension_log_per_side"] == 0.0
    assert forced["sigma_source_support_extension_total_log_per_side"] == 0.0
    assert forced["mass_ratio_class"] == "1:1:1"
    assert artifact["public_charged_mass_promotion_allowed"] is False


def test_extension_requires_non_singlet_and_no_target_leak():
    artifact = lane.build_artifact(lane._load(lane.HIERARCHY), lane._load(lane.SCALARS))
    required = artifact["minimal_new_theorem"]["required_objects"]
    assert any("non-singlet" in item for item in required)
    assert any("no-target-leak" in item for item in required)
