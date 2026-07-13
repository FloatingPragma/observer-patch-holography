"""Tests for the entropic conditioned-branch no-go."""

from __future__ import annotations

import derive_charged_entropic_branch_no_go as lane


def test_no_go_is_certified_and_fail_closed():
    artifact = lane.build()
    assert artifact["status"] == "ENTROPIC_CONDITIONED_BRANCH_NO_GO"
    assert artifact["checks_pass"] is True
    assert artifact["minimum_gap"] < 1.0e-8
    assert artifact["promotion_allowed"] is False
