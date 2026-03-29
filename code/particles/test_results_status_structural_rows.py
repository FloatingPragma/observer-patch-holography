#!/usr/bin/env python3
"""Guard the structural massless rows on the public `/particles` surface."""

from __future__ import annotations

import importlib.util
import pathlib


ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "particles" / "scripts" / "build_results_status_table.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("build_results_status_table", SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"unable to load {SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_structural_inventory_contains_graviton() -> None:
    module = _load_module()
    particle_ids = {row["particle_id"] for row in module.INVENTORY}
    assert {"photon", "gluon", "graviton"} <= particle_ids


def test_structural_massless_defaults_include_graviton() -> None:
    module = _load_module()
    updated = module.apply_local_candidate_overrides({})
    assert updated["m_gamma"] == 0.0
    assert updated["m_gluon"] == 0.0
    assert updated["m_graviton"] == 0.0


def test_graviton_uses_structural_surface() -> None:
    module = _load_module()
    surface_state = module.build_surface_state(with_hadrons=False)
    row_spec = next(row for row in module.INVENTORY if row["particle_id"] == "graviton")
    assert module.prediction_surface_for_row(row_spec, surface_state, with_hadrons=False) == "particles_structural_massless"
