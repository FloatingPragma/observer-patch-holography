#!/usr/bin/env python3
"""Kernel-check the typed fixed-federation execution audit module."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
LEAN_ROOT = ROOT / "Lean"
TARGET = "Tower.FixedFederationExecutionAudit"


def main() -> None:
    lake = shutil.which("lake")
    if lake is None:
        raise SystemExit("lake executable is unavailable")
    completed = subprocess.run(
        [lake, "build", TARGET],
        cwd=LEAN_ROOT,
        check=False,
        timeout=600,
    )
    if completed.returncode != 0:
        raise SystemExit(completed.returncode)
    print(f"typed execution audit OK: {TARGET}")


if __name__ == "__main__":
    main()
