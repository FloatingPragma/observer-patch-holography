#!/usr/bin/env python3
"""Tests for the NuFIT 6.1 correlated-profile scorer."""

from __future__ import annotations

import importlib.util
import math
import pathlib
import sys


HERE = pathlib.Path(__file__).resolve().parent
SCRIPT = HERE / "score_neutrino_nufit61.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("oph_nufit61_scorer", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_bilinear_interpolation_and_out_of_grid_failure() -> None:
    module = _load_module()
    grid = module.RectilinearGrid.from_rows(
        [
            (0.0, 0.0, 0.0),
            (0.0, 2.0, 4.0),
            (2.0, 0.0, 2.0),
            (2.0, 2.0, 6.0),
        ],
        "synthetic",
    )
    assert math.isclose(grid.interpolate(1.0, 1.0)["delta_chi2"], 3.0, abs_tol=1.0e-15)
    try:
        grid.interpolate(-0.1, 1.0)
    except ValueError as exc:
        assert "extrapolation is forbidden" in str(exc)
    else:
        raise AssertionError("out-of-grid interpolation must fail")


def test_candidate_coordinate_transform_and_phase_wrap() -> None:
    module = _load_module()
    candidate = {
        "pmns_observables": {
            "theta12_deg": 30.0,
            "theta13_deg": 10.0,
            "theta23_deg": 45.0,
            "delta_deg": 305.0,
        },
        "dimensionless_ratio_dm21_over_dm32": 0.03,
    }
    coordinates = module._candidate_coordinates(candidate)
    assert math.isclose(coordinates["sin2_theta12"], 0.25, abs_tol=1.0e-15)
    assert math.isclose(coordinates["sin2_theta23"], 0.5, abs_tol=1.0e-15)
    assert math.isclose(coordinates["delta_cp_deg_wrapped"], -55.0, abs_tol=1.0e-15)


def test_template_kernel_blocks_prediction_promotion() -> None:
    module = _load_module()
    boundary = module._source_boundary(
        {"status": "template", "proof_status": "conjugacy_riesz_candidate"}
    )
    assert boundary["source_only_prediction_eligible"] is False
    assert boundary["prospective_evidence_eligible"] is False
    assert boundary["historical_target_exposure"] is True


def test_source_emitted_kernel_is_necessary_but_history_stays_retrospective() -> None:
    module = _load_module()
    boundary = module._source_boundary(
        {"status": "source_only_frozen", "proof_status": "closed_source_emitted"}
    )
    assert boundary["source_only_prediction_eligible"] is True
    assert boundary["prospective_evidence_eligible"] is False
