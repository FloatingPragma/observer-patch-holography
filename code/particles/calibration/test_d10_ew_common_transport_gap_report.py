#!/usr/bin/env python3
"""Guard the D10 common-transport gap report."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_common_transport_gap_report.py"
OUTPUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_common_transport_gap_report.json"


def test_d10_common_transport_gap_report_default_p() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_d10_ew_common_transport_gap_report"
    assert payload["status"] == "compare_only_diagnostic"
    assert payload["verdict"]["classification"] == "fails_required_benchmark_transport"
    assert payload["observable_gap_to_benchmark"]["max_mass_gap_mev"] > 10.0


def test_d10_common_transport_gap_report_alex_p() -> None:
    alex_p = "1.6309681897"
    subprocess.run([sys.executable, str(SCRIPT), "--p", alex_p], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_d10_ew_common_transport_gap_report"
    assert payload["p"] == float(alex_p)
    assert payload["verdict"]["classification"] == "fails_required_benchmark_transport"
    assert payload["observable_gap_to_benchmark"]["max_mass_gap_mev"] > 10.0
