#!/usr/bin/env python3
"""Fail when a tracked file links a private repository.

Readers of this public repository cannot open a private one, so no tracked file
links it. Naming a private project is allowed; a URL into it is not.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PRIVATE = ("github.com/" + "muellerberndt/openwave",)
LIMIT = 8 * 1024 * 1024


def violations(root: Path = ROOT) -> list[str]:
    listed = subprocess.run(["git", "ls-files", "-z"], cwd=root, capture_output=True, check=True).stdout
    found = []
    for name in filter(None, listed.decode().split("\0")):
        path = root / name
        if not path.is_file() or path.stat().st_size > LIMIT:
            continue
        data = path.read_bytes().lower()
        found += [f"{name}: links {target}" for target in PRIVATE if target.encode() in data]
    return found


def main() -> int:
    found = violations()
    for line in found:
        print(line)
    print(f"private links: {len(found)} file(s)")
    return 1 if found else 0


if __name__ == "__main__":
    sys.exit(main())
