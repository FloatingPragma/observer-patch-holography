#!/usr/bin/env python3
"""Guard the D10 current-carrier tau2 obstruction artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SOURCE_PAIR_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_source_transport_pair.py"
POPULATION_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_population_evaluator.py"
FIBERWISE_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_fiberwise_population_tree_law_beneath_single_tree_identity.py"
SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_tau2_current_carrier_obstruction.py"
OUTPUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_tau2_current_carrier_obstruction.json"


def test_d10_tau2_current_carrier_obstruction_is_emitted() -> None:
    subprocess.run([sys.executable, str(SOURCE_PAIR_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(POPULATION_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(FIBERWISE_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_d10_ew_tau2_current_carrier_obstruction"
    assert payload["status"] == "closed_smaller_primitive"
    assert payload["proof_status"] == "no_single_tau2_on_closed_current_carrier_can_hit_exact_W_and_exact_Z"
    assert payload["next_single_residual_object"] == "delta_n_tree_exact"
    assert payload["direction_obstruction"]["single_tau2_possible"] is False
