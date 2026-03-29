#!/usr/bin/env python3
"""Guard quark-row promotion against uncertified spread candidates."""

from __future__ import annotations

import importlib.util
import json
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


def test_public_quark_candidate_policy_uses_explicit_surface_gate() -> None:
    module = _load_module()
    base = {
        "m_u": 1.0,
        "m_d": 2.0,
        "m_s": 3.0,
        "m_c": 4.0,
        "m_b": 5.0,
        "m_t": 6.0,
    }
    updated = module.apply_local_candidate_overrides(base)
    forward = json.loads((ROOT / "particles" / "runs" / "flavor" / "forward_yukawas.json").read_text(encoding="utf-8"))
    if forward.get("public_surface_candidate_allowed", False):
        assert updated["m_u"] != 1.0
        assert updated["m_d"] != 2.0
        assert updated["m_s"] != 3.0
        assert updated["m_c"] != 4.0
        assert updated["m_b"] != 5.0
        assert updated["m_t"] == 6.0
    else:
        assert updated["m_u"] == 1.0
        assert updated["m_d"] == 2.0
        assert updated["m_s"] == 3.0
        assert updated["m_c"] == 4.0
        assert updated["m_b"] == 5.0
