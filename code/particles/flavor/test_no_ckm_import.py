#!/usr/bin/env python3
"""Fail if the flavor sandbox reintroduces hardcoded CKM-angle placeholders."""

from __future__ import annotations

import argparse
import pathlib
import re
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_TARGETS = [
    ROOT / "particles" / "flavor" / "build_forward_yukawas.py",
    ROOT / "particles" / "flavor" / "export_blind_forward_artifact.py",
]
FORBIDDEN_PATTERNS = [
    re.compile(r"\bs12\b"),
    re.compile(r"\bs23\b"),
    re.compile(r"\bs13\b"),
    re.compile(r"Wolfenstein", re.IGNORECASE),
    re.compile(r"delta\s*=\s*0"),
]


def main() -> int:
    parser = argparse.ArgumentParser(description="Check that the flavor sandbox does not hardcode CKM-angle imports.")
    parser.add_argument("targets", nargs="*", help="Files to scan. Defaults to the /particles flavor lane.")
    args = parser.parse_args()

    targets = [pathlib.Path(item) for item in args.targets] if args.targets else DEFAULT_TARGETS
    failures: list[str] = []

    for target in targets:
        text = target.read_text(encoding="utf-8")
        for pattern in FORBIDDEN_PATTERNS:
            if pattern.search(text):
                failures.append(f"{target}: matched forbidden pattern {pattern.pattern!r}")

    if failures:
        for item in failures:
            print(item, file=sys.stderr)
        return 1

    print("no forbidden CKM-import patterns found")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
