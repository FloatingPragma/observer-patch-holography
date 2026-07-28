#!/usr/bin/env python3
"""Rebuild, independently verify, and test the #630 frontier."""

from __future__ import annotations

from pathlib import Path
import subprocess
import sys


PACKAGE = Path(__file__).resolve().parent
REPO_ROOT = PACKAGE.parents[3]


def run(*args: str) -> None:
    subprocess.run(args, cwd=REPO_ROOT, check=True)


def main() -> int:
    run(sys.executable, str(PACKAGE / "build_higgs_yukawa_source_frontier.py"))
    run(sys.executable, str(PACKAGE / "check_higgs_yukawa_source_frontier.py"))
    run(
        sys.executable,
        "-m",
        "pytest",
        "-q",
        str(PACKAGE / "tests" / "test_higgs_yukawa_source_frontier.py"),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
