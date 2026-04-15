#!/usr/bin/env python3
"""Validate the exact split D11 Higgs/top theorem artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys

import pytest


ROOT = pathlib.Path(__file__).resolve().parents[2]
DECLARED_SURFACE_SCRIPT = ROOT / "particles" / "calibration" / "derive_d11_declared_calibration_surface.py"
HIGGS_SCRIPT = ROOT / "particles" / "calibration" / "derive_d11_live_exact_higgs_promotion.py"
TOP_SCRIPT = ROOT / "particles" / "calibration" / "derive_d11_live_exact_top_promotion.py"
NO_GO_SCRIPT = ROOT / "particles" / "calibration" / "derive_d11_fixed_ray_no_go_theorem.py"
PAIR_SCRIPT = ROOT / "particles" / "calibration" / "derive_d11_live_exact_split_pair_theorem.py"
OUTPUT = ROOT / "particles" / "runs" / "calibration" / "d11_live_exact_split_pair_theorem.json"


def test_d11_live_exact_split_pair_theorem_closes_exact_pair() -> None:
    subprocess.run([sys.executable, str(DECLARED_SURFACE_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(HIGGS_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(TOP_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(NO_GO_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(PAIR_SCRIPT)], check=True, cwd=ROOT)

    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_d11_live_exact_split_pair_theorem"
    assert payload["proof_status"] == "closed_target_anchored_live_exact_split_pair"
    assert payload["closure_logic"]["fixed_ray_blocked"] is True
    assert payload["exact_split_pair"]["mH_gev"] == pytest.approx(125.1995304097179, abs=1.0e-12)
    assert payload["exact_split_pair"]["mt_pole_gev"] == pytest.approx(172.3523553288312, abs=1.0e-12)
    assert payload["exact_split_pair"]["w_HT_exact"] == pytest.approx(-0.0003857630977715052, abs=1.0e-18)

