"""Collection policy for the optional IBM / Qiskit hardware lane.

These tests import qiskit (and, on some hosts, POSIX-only modules such as
``fcntl``) at import time. The lane is an optional hardware surface, not part of
the mandatory clean-clone scientific suite. Collecting it by default makes a
fresh clone without the IBM extras fail at pytest collection before any
mandatory test runs.

Opt in after installing ``requirements-ibm.txt``::

    OPH_RUN_IBM=1 python -m pytest code/ibm_quantum_cloud

This mirrors the existing opt-in gate for the legacy arXiv D10 helpers
(``OPH_RUN_LEGACY_D10``) in ``code/particles/conftest.py``.
"""

from __future__ import annotations

import os

if os.environ.get("OPH_RUN_IBM") != "1":
    collect_ignore_glob = ["*"]
