#!/usr/bin/env python3
"""Run the deterministic producer, independent checker, and focused tests."""

from __future__ import annotations

from pathlib import Path
import subprocess
import sys


PACKAGE = Path(__file__).resolve().parent


def run(*args: str) -> None:
    subprocess.run(args, cwd=PACKAGE, check=True)


def main() -> int:
    run(sys.executable, "build_rg_representation_frontier.py", "--check-byte-exact")
    run(sys.executable, "check_rg_representation_frontier.py")
    run(sys.executable, "-m", "pytest", "-q", "tests/test_rg_representation_frontier.py")
    print("issue #32 RG representation frontier: all checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
