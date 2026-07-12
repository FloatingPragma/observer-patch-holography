"""Tests for the OPH dynamics obstruction to a nonuniform port record."""

from __future__ import annotations

import json

import pytest

import derive_charged_nonuniform_port_record_dynamics_no_go as lane


def test_uniform_record_is_unique_convex_minimum():
    uniform = lane.uniform_minimum(12.0, 12)
    assert uniform == [1.0] * 12
    assert lane.convex_cost(uniform) == pytest.approx(12.0)
    assert lane.convex_cost([2.0, 0.0] + [1.0] * 10) > 12.0


def test_live_screen_dynamics_force_zero_non_singlet_moments():
    screen = json.loads(lane.SCREEN.read_text(encoding="utf-8"))
    artifact = lane.build_artifact(screen)
    checks = artifact["theorem"]["checks"]
    assert checks["uniform_minimum_matches_certificate"] is True
    assert checks["first_moment_norm"] < 1.0e-12
    assert checks["quadrupole_norm"] < 1.0e-12
    assert artifact["nonuniform_record_derived"] is False
    assert artifact["public_charged_mass_promotion_allowed"] is False


def test_required_extension_is_source_conditioned_and_target_free():
    screen = json.loads(lane.SCREEN.read_text(encoding="utf-8"))
    extension = lane.build_artifact(screen)["smallest_consistent_extension"]
    assert "observer_conditioned" in extension["id"]
    assert any("no-target-leak" in item for item in extension["required_objects"])
    assert "electron" in extension["not_allowed"]
