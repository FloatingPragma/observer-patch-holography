#!/usr/bin/env python3
"""Smoke tests for the real diagnostic lattice backend runner."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys

import pytest

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCRIPT = HERE / "run_lattice_diagnostic_backend.py"
EXPORT = ROOT / "particles" / "runs" / "hadron" / "lattice_diagnostic_backend_export.json"


@pytest.mark.slow
def test_smoke_run_emits_guarded_export(tmp_path) -> None:
    out = tmp_path / "smoke_export.json"
    subprocess.run(
        [sys.executable, str(SCRIPT), "--smoke", "--output", str(out)],
        check=True, cwd=ROOT, timeout=1200)
    payload = json.loads(out.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_lattice_diagnostic_backend_export"
    assert payload["execution_class"] == "real_lattice_diagnostic_toy_scale"
    guards = payload["guards"]
    assert guards["real_lattice_execution"] is True
    assert guards["target_anchored"] is False
    assert guards["production_execution_class"] is False
    assert guards["public_promotion_allowed"] is False
    assert guards["satisfies_issue_425_closure"] is False
    assert payload["channels"]["exported"] == ["pi_iso", "N_iso_direct", "N_iso_exchange"]
    # correlators are present, finite, and the pion channel is positive
    for cfg in payload["correlators_per_config"]:
        for _kappa, block in cfg.items():
            assert all(v > 0.0 for v in block["pi_iso"])
            assert block["max_cg_relative_residual"] < 1e-5


def _existing_export() -> dict | None:
    if not EXPORT.exists():
        return None
    return json.loads(EXPORT.read_text(encoding="utf-8"))


def test_production_export_when_present_is_real_and_non_promoting() -> None:
    payload = _existing_export()
    if payload is None:
        pytest.skip("production diagnostic export not generated yet")
    assert payload["execution_class"] == "real_lattice_diagnostic_toy_scale"
    assert payload["guards"]["target_anchored"] is False
    assert payload["guards"]["public_promotion_allowed"] is False
    assert payload["ensemble"]["plaquette_agrees_with_literature"] is True
    for block in payload["analysis"]:
        am_pi = block["pi_iso"]["am_plateau"]
        am_n = block["N_iso"]["am_plateau"]
        # real toy-scale masses: finite, heavy, nucleon above pion, and NOT
        # tuned to physical values (the pion here is far from 140 MeV scale)
        assert 0.1 < am_pi < 2.0
        assert am_n > am_pi
    hmc = payload["dynamical_branch_validation"]
    assert hmc["trajectories_measured"] >= 8
    assert hmc["acceptance_rate"] > 0.3
