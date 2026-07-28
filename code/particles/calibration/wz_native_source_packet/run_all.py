#!/usr/bin/env python3
"""Run the complete issue-#594 source-frontier validation."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent


def run(*argv: str) -> None:
    subprocess.run(
        [sys.executable, *argv],
        cwd=HERE,
        check=True,
    )


def main() -> int:
    run("build_source_parent_inventory.py", "--check-byte-exact")
    run("check_source_parent_inventory.py")
    run("-m", "pytest", "-q", "tests/test_source_parent_inventory.py")
    print("ISSUE_594_SOURCE_FRONTIER_PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
