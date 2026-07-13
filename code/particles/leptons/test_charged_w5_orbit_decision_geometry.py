"""Tests for the W5 orbit decision geometry (fast components live, locus from disk)."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import derive_charged_w5_orbit_decision_geometry as lane

ARTIFACT = (
    Path(__file__).resolve().parents[2]
    / "particles" / "runs" / "leptons"
    / "charged_w5_orbit_decision_geometry.json"
)


def test_stratum_theorem_certificates():
    strata = lane.stratum_theorem_certificates()
    assert strata["c5_degenerate"] is True
    assert strata["c3_degenerate"] is True
    assert strata["c2_simple"] is True


def test_projector_and_quadrupole_structure():
    assert np.allclose(lane.P5 @ lane.P5, lane.P5, atol=1e-12)
    assert abs(np.trace(lane.P5) - 5.0) < 1e-9
    rng = np.random.default_rng(5)
    w = lane.P5 @ rng.standard_normal(12)
    assert abs(np.trace(lane.quadrupole(w))) < 1e-10


def test_gradient_matches_finite_difference():
    rng = np.random.default_rng(7)
    w = lane.P5 @ rng.standard_normal(12)
    coefs = (-1.0, 0.1, 0.2, 0.25, 0.25)
    g = lane.gradient(w, coefs)
    eps = 1e-7
    for k in (0, 4, 9):
        w2 = w.copy()
        w2[k] += eps
        fd = (lane.potential(lane.P5 @ w2, coefs) - lane.potential(w, coefs)) / eps
        assert abs((lane.P5 @ np.eye(12)[k]) @ g / max(abs(fd), 1e-12) - fd / max(abs(fd), 1e-12)) < 2e-2 or abs(g[k] - fd) < 1e-4


def test_artifact_on_disk_certifies_the_locus():
    artifact = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert artifact["checks_pass"] is True
    assert artifact["target_locus"]["locus_nonempty_certified"] is True
    point = artifact["target_locus"]["locus_point"]
    assert abs(point["ratio"] - artifact["mcpr_target_shape"]["sorted_gap_ratio"]) < 0.05
    assert artifact["source_coefficient_gate"]["status"] == "open"
    assert artifact["promotion_allowed"] is False
