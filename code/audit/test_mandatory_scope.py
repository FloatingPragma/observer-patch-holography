"""Scope guard for the mandatory suite (issue #507).

Proves that ordinary ``pytest`` discovery cannot silently treat a cloud or
hardware experiment as mandatory science. The two optional-lane mechanisms
have different contracts, and each is asserted literally:

- the IBM / Qiskit hardware lane is excluded from collection entirely
  (``collect_ignore_glob`` in its conftest), because importing it without the
  IBM extras breaks at import time;
- the legacy arXiv D10 helpers are collected but explicitly skip-marked with
  an opt-in reason, so a direct invocation without the extras fails clearly
  instead of running silently.

Both lanes must exist on disk so the guard cannot pass vacuously.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

OPT_IN_VARIABLES = ["OPH_RUN_IBM", "OPH_RUN_LEGACY_D10", "OPH_LEGACY_PARTICLE_DIR"]

IBM_LANE = "code/ibm_quantum_cloud"
D10_LANE = "code/particles/calibration/test_d10_ew_exactness_audit.py"


def run_pytest(args: list[str]) -> subprocess.CompletedProcess:
    env = {k: v for k, v in os.environ.items() if k not in OPT_IN_VARIABLES}
    return subprocess.run(
        [sys.executable, "-m", "pytest", *args],
        cwd=REPO_ROOT,
        env=env,
        capture_output=True,
        text=True,
    )


def test_ibm_lane_exists_and_stays_out_of_mandatory_collection():
    assert (REPO_ROOT / IBM_LANE).is_dir(), (
        f"optional lane {IBM_LANE} is gone; update the scope guard"
    )
    result = run_pytest(["--collect-only", "-q", "code"])
    assert result.returncode == 0, (
        "mandatory collection failed from a clean environment:\n"
        f"{result.stdout}\n{result.stderr}"
    )
    assert "ibm_quantum_cloud" not in result.stdout, (
        "IBM hardware lane tests were collected into the mandatory suite"
    )


def test_legacy_d10_lane_skips_clearly_without_opt_in():
    assert (REPO_ROOT / D10_LANE).exists(), (
        f"optional lane {D10_LANE} is gone; update the scope guard"
    )
    result = run_pytest(["-rs", D10_LANE])
    assert result.returncode == 0, (
        f"direct D10 invocation errored instead of skipping:\n{result.stdout}"
    )
    assert "OPH_RUN_LEGACY_D10" in result.stdout, (
        "D10 skip output does not name the opt-in variable:\n" + result.stdout
    )
    assert " passed" not in result.stdout and " failed" not in result.stdout, (
        "D10 lane executed without its opt-in extras:\n" + result.stdout
    )
