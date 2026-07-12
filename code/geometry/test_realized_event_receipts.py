#!/usr/bin/env python3
"""Tests for the realized-record event receipts (#503, #525)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from realized_event_receipts import (  # noqa: E402
    cone_structure_receipt,
    e1_population_receipt,
    e2_separation_receipt,
    e4_cocycle_receipt,
    instrument_realized_events,
    one_ring_step,
    pca_dimension_receipt,
    realized_records,
    tower_meshes,
)


def test_tower_meshes_stay_on_the_sphere():
    tower = tower_meshes(3)
    for _, coords in tower:
        norms = np.linalg.norm(coords[np.any(coords != 0, axis=1)], axis=1)
        assert np.allclose(norms, 1.0, atol=1e-12)


def test_e1_e2_on_realized_germs():
    germs = realized_records(stages=3)
    e1 = e1_population_receipt(germs, 3)
    assert e1["radii_shrinking"] and e1["one_germ_per_stage0_cell"]
    e2 = e2_separation_receipt(germs)
    assert e2["all_separated"]
    assert e2["worst_margin"] > 0.0


def test_e4_moebius_cocycle_at_noise_level():
    germs = realized_records(stages=3)
    e4 = e4_cocycle_receipt(germs)
    assert e4["cocycle_at_noise_level"]


def test_cone_is_lorentzian_with_clock_timelike():
    cone = cone_structure_receipt(n_ticks=40, per_tick=3, seed=11)
    assert cone["signature"] == [1, 2]
    assert cone["timelike_eigenvector_clock_alignment"] > 0.9
    assert cone["classification_rate"] > 0.75
    # the cone interior is underfilled by construction: causal labels need
    # actual intermediate commits, so perfect classification would itself
    # be suspicious
    assert cone["classification_rate"] < 0.99


def test_one_ring_step_matches_mesh_scale():
    tower = tower_meshes(3)
    records, coords = tower[2]
    step = one_ring_step(records, coords)
    assert 0.1 < step < 0.5


def test_realized_sheet_dimension_is_three_not_four():
    germs = realized_records(stages=3)
    pca = pca_dimension_receipt(germs)
    assert pca["chart_pca_dimension"] == 3


def test_instrumented_events_report_bulk_negative():
    report = instrument_realized_events()
    w = report["receipts_witnessed"]
    assert w["e1_screen_population"]
    assert w["e2_certified_separation"]
    assert w["e4_moebius_cocycle"]
    assert w["intrinsic_cone_lorentzian_1p2"]
    assert w["realized_sheet_dimension_measured"]
    # negative verdict: no bulk-depth channel exists on the realized
    # tower, so the rank-four clause must NOT be reported witnessed
    assert w["e3_rank_four_bulk_depth"] is False
    assert any("bulk-depth" in p or "bulk" in p
               for p in report["receipts_pending"])
