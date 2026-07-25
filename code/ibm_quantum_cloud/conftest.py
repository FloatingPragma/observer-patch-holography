"""Collection policy for the optional IBM / Qiskit hardware lane.

These tests import qiskit (and, on some hosts, POSIX-only modules such as
``fcntl``) at import time. The lane is an optional hardware surface, not part of
the mandatory clean-clone scientific suite. Broad collection ignores it. A
direct invocation without the explicit opt-in or complete extras exits with an
actionable message instead of reporting the ambiguous pytest result "no tests
ran".

Opt in after installing ``requirements-ibm.txt``::

    OPH_RUN_IBM=1 python -m pytest code/ibm_quantum_cloud

This mirrors the existing opt-in gate for the legacy arXiv D10 helpers
(``OPH_RUN_LEGACY_D10``) in ``code/particles/conftest.py``.
"""

from __future__ import annotations

import importlib.util
import os
from pathlib import Path

import pytest


LANE_ROOT = Path(__file__).resolve().parent
REQUIRED_MODULES = ("qiskit", "qiskit_aer", "qiskit_ibm_runtime")


def _directly_targets_lane(config: pytest.Config) -> bool:
    invocation_dir = Path(config.invocation_params.dir)
    for raw_arg in config.invocation_params.args:
        if raw_arg.startswith("-"):
            continue
        candidate = (invocation_dir / raw_arg).resolve()
        if candidate == LANE_ROOT or LANE_ROOT in candidate.parents:
            return True
    return False


def pytest_configure(config: pytest.Config) -> None:
    if not _directly_targets_lane(config):
        return
    if os.environ.get("OPH_RUN_IBM") != "1":
        pytest.exit(
            "IBM/Qiskit is an optional hardware lane. Set OPH_RUN_IBM=1 "
            "after installing code/ibm_quantum_cloud/requirements-ibm.txt.",
            returncode=4,
        )
    missing = [
        module
        for module in REQUIRED_MODULES
        if importlib.util.find_spec(module) is None
    ]
    if missing:
        pytest.exit(
            "IBM/Qiskit optional extras are missing: "
            + ", ".join(missing)
            + ". Install code/ibm_quantum_cloud/requirements-ibm.txt.",
            returncode=4,
        )

if os.environ.get("OPH_RUN_IBM") != "1":
    collect_ignore_glob = ["*"]
