#!/usr/bin/env python3
"""Validate the smaller D10 fiberwise tree-law primitive."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SOURCE_PAIR_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_source_transport_pair.py"
POPULATION_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_population_evaluator.py"
SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_fiberwise_population_tree_law_beneath_single_tree_identity.py"
OUTPUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_fiberwise_population_tree_law_beneath_single_tree_identity.json"


def test_d10_fiberwise_population_tree_law_reduces_to_tau2_scalar() -> None:
    subprocess.run([sys.executable, str(SOURCE_PAIR_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(POPULATION_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_d10_ew_fiberwise_population_tree_law_beneath_single_tree_identity"
    assert payload["status"] == "closed_smaller_primitive"
    assert payload["strictly_smaller_than"] == "EWSinglePostTransportTreeIdentity_D10"
    assert payload["next_single_residual_object"] == "tau2_tree_exact"
    assert payload["anchor_point"]["tau2_tree_exact"] == 0.0
    assert payload["carrier_basis_scalar"]["v_inherited"] > 0.0
