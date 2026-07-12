#!/usr/bin/env python3
"""Tests for the bulk-depth record channel (#503, the 3+1 test)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from bulk_depth_receipts import (  # noqa: E402
    build_cells,
    bulk_point,
    cone_4d_receipt,
    e3_rank_four_receipt,
    multiscale_repair_run,
    pca_bulk_dimension_receipt,
    scale_invariant_speed_check,
)


def test_depth_dictionary_orders_scales():
    cells = build_cells(3)
    rho_by_stage = {}
    for c in cells:
        rho_by_stage.setdefault(c["stage"], []).append(c["rho"])
    means = [np.mean(rho_by_stage[s]) for s in range(3)]
    # finer stages sit closer to the ideal boundary (larger rho)
    assert means[0] < means[1] < means[2]
    # bulk points lie on the unit hyperboloid
    eta = np.diag([-1.0, 1.0, 1.0, 1.0])
    for c in cells[::37]:
        x = bulk_point(c)
        assert abs(x @ eta @ x + 1.0) < 1e-9


def test_hyperbolic_metric_equalizes_lattice_speeds():
    speed = scale_invariant_speed_check(build_cells(3))
    assert speed["scale_invariant_within_factor_2"]


def test_e3_rank_four_and_bulk_pca():
    cells = build_cells(3)
    commits, _ = multiscale_repair_run(cells, n_ticks=30, per_tick=3, seed=5)
    e3 = e3_rank_four_receipt(cells, commits)
    assert e3["frame_rank"] == 4
    assert e3["observability_sigma_min"] > 0.5
    pca = pca_bulk_dimension_receipt(cells, commits, step_time=2.5)
    assert pca["chart_pca_dimension"] == 4


def test_cone_is_1_3_on_inverse_system_dynamics():
    cells = build_cells(3, coupling="inverse_system")
    commits, depends = multiscale_repair_run(cells, n_ticks=60, per_tick=4,
                                             seed=1)
    cone = cone_4d_receipt(cells, commits, depends)
    assert cone["signature"] == [1, 3]
    assert cone["timelike_eigenvector_clock_alignment"] > 0.9
    assert cone["classification_rate"] > 0.8


def test_countermodel_global_coupling_degenerates_depth():
    cells = build_cells(3, coupling="global_namespace")
    commits, depends = multiscale_repair_run(cells, n_ticks=60, per_tick=4,
                                             seed=20260712)
    cone = cone_4d_receipt(cells, commits, depends)
    # infinite depth speed must NOT pass as 3+1: the receipt detects it
    assert cone["signature"] != [1, 3]
