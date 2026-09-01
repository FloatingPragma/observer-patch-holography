#!/usr/bin/env python3
"""Tests for the realized cyclic cap-net repair run (#503 nonemptiness gate)."""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from cyclic_cap_net_run import build_tower, run_tower  # noqa: E402


def test_tower_shape():
    tower = build_tower(3)
    assert [len(t) for t in tower] == [20, 80, 320]


def test_cyclic_run_witnesses_topology_and_mesh_families():
    report = run_tower(stages=2, seed=7)
    w = report["receipts_witnessed"]
    assert w["d1_repair_clauses"] is True
    assert w["spherical_incidence_all_stages"] is True
    assert w["mesh_modulus_decreasing"] is True
    assert w["refinement_naturality"] is True
    # provenance must state the explicit branch selection, never a derivation
    assert "branch selection" in report["provenance"]
    assert len(report["receipts_pending"]) >= 4


def test_every_stage_is_a_genuine_repair_run():
    report = run_tower(stages=2, seed=11)
    for stage in report["stages"]:
        d1 = stage["d1_verification"]
        assert d1["terminating"] and d1["conflict_free_normal_form"]
        assert d1["schedule_independent"]
        assert d1["n_schedules_per_seed"] >= 20
        inc = stage["incidence"]
        assert inc["surface_classification"] == "S2"
        assert inc["euler_characteristic"] == 2
