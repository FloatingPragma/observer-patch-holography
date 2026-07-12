"""Tests for the conditional non-singlet twelve-port attachment theorem."""

from __future__ import annotations

import numpy as np
import pytest

import derive_charged_family_non_singlet_port_attachment as lane


def test_equal_load_has_no_non_singlet_moment():
    first, quadrupole = lane.port_moments([1.0] * 12)
    assert np.linalg.norm(first) < 1.0e-12
    assert np.linalg.norm(quadrupole) < 1.0e-12


def test_moments_are_independent_of_port_ids():
    vertices = lane.icosahedron_vertices()
    weights = np.asarray([0.2, 1.1, -0.4, 0.7, 1.9, -1.2, 0.3, 2.1, -0.8, 0.6, 1.4, -0.1])
    permutation = np.asarray([7, 2, 10, 0, 5, 11, 3, 8, 1, 9, 6, 4])
    first, quadrupole = lane.port_moments(weights, vertices)
    first_p, quadrupole_p = lane.port_moments(weights[permutation], vertices[permutation])
    assert first_p == pytest.approx(first)
    assert quadrupole_p == pytest.approx(quadrupole)


def test_generic_nonuniform_record_emits_two_shape_coordinates():
    weights = [0.2, 1.1, -0.4, 0.7, 1.9, -1.2, 0.3, 2.1, -0.8, 0.6, 1.4, -0.1]
    readout = lane.attachment_readout(weights, x2=-0.5175863354681689)
    assert readout["simple_spectrum"] is True
    assert abs(readout["trace_residual"]) < 1.0e-12
    assert readout["sigma_source_support_extension_total_log_per_side"] is not None
    assert readout["eta_source_support_extension_log_per_side"] is not None


def test_live_receipt_fails_closed_without_nonuniform_source_packet():
    hierarchy = lane._load(lane.HIERARCHY) if hasattr(lane, "_load") else None
    if hierarchy is None:
        import json
        hierarchy = json.loads(lane.HIERARCHY.read_text(encoding="utf-8"))
    artifact = lane.build_artifact(hierarchy, None)
    assert artifact["public_promotion_allowed"] is False
    assert artifact["readout"] is None
    assert artifact["next_exact_missing_object"] == "refinement_stable_source_only_nonuniform_12_port_record"
