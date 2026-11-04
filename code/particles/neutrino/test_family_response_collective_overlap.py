#!/usr/bin/env python3
"""Guard the collective-mode overlap in the neutrino family-response builder."""

from __future__ import annotations

import importlib.util
import json
import pathlib
import subprocess
import sys

import pytest


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "neutrino" / "derive_family_response_tensor.py"
SECTOR = ROOT / "particles" / "runs" / "flavor" / "sector_transport_pushforward.json"
ANCHOR = ROOT / "particles" / "runs" / "neutrino" / "neutrino_scale_anchor.json"
OUTPUT = ROOT / "particles" / "runs" / "neutrino" / "family_response_tensor.json"


def _load_module():
    spec = importlib.util.spec_from_file_location("derive_family_response_tensor", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_zero_overlap_is_rejected_instead_of_replaced_by_unit_overlap() -> None:
    module = _load_module()
    rank_two = module.np.diag([1.0, 1.0, 0.0])
    with pytest.raises(ValueError, match="degenerate"):
        module._collective_overlap(rank_two, module.np.array([0.0, 0.0, 1.0]))
    # A supported collective mode passes through untouched, not normalized to 1.
    assert module._collective_overlap(rank_two, module.np.array([0.0, 1.0, 0.0])) == pytest.approx(1.0)


def test_degenerate_collective_mode_fails_the_producer(tmp_path) -> None:
    sector = json.loads(SECTOR.read_text(encoding="utf-8"))
    sector["sector_response_object"]["nu"]["K_core_majorana_sym"] = [
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0],
    ]
    sector["sector_response_object"]["nu"]["K_core_directed"] = [
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0],
    ]
    anchor = json.loads(ANCHOR.read_text(encoding="utf-8"))
    anchor["collective_mode"] = {
        "basis": "synthetic degenerate control",
        "u_vector": [0.0, 0.0, 1.0],
        "u_uT": [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
    }
    sector_path = tmp_path / "sector.json"
    anchor_path = tmp_path / "anchor.json"
    sector_path.write_text(json.dumps(sector), encoding="utf-8")
    anchor_path.write_text(json.dumps(anchor), encoding="utf-8")
    completed = subprocess.run(
        [sys.executable, str(SCRIPT),
         "--sector", str(sector_path),
         "--scale-anchor", str(anchor_path),
         "--out", str(tmp_path / "family.json")],
        capture_output=True,
        text=True,
    )
    assert completed.returncode != 0
    assert "collective-mode overlap is degenerate" in completed.stderr
    assert not (tmp_path / "family.json").exists()


def test_committed_inputs_still_reproduce_the_committed_receipt(tmp_path) -> None:
    out_path = tmp_path / "family.json"
    subprocess.run(
        [sys.executable, str(SCRIPT), "--out", str(out_path)],
        check=True,
        capture_output=True,
        text=True,
    )
    committed = json.loads(OUTPUT.read_text(encoding="utf-8"))
    fresh = json.loads(out_path.read_text(encoding="utf-8"))
    assert fresh.pop("generated_utc")
    committed.pop("generated_utc")
    assert fresh == committed
    assert committed["u_normalization"] == pytest.approx(1.061872133120721, rel=1.0e-12)
