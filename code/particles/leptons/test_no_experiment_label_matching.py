#!/usr/bin/env python3
"""Fail if the charged-lepton lane resolves families by experimental labels."""

from __future__ import annotations

import argparse
import pathlib
import re
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_TARGETS = [
    ROOT / "particles" / "leptons" / "derive_lepton_descent_tensor.py",
    ROOT / "particles" / "leptons" / "derive_lepton_channel_norm.py",
    ROOT / "particles" / "leptons" / "build_forward_charged_leptons.py",
    ROOT / "particles" / "leptons" / "run_qed_ew_completion.py",
]
FORBIDDEN_PATTERNS = [
    re.compile(r"\belectron\b", re.IGNORECASE),
    re.compile(r"\bmuon\b", re.IGNORECASE),
    re.compile(r"\btau\b", re.IGNORECASE),
    re.compile(r"\bmeasured\b", re.IGNORECASE),
    re.compile(r"\bPDG\b"),
    re.compile(r"\bbest\s*fit\b", re.IGNORECASE),
    re.compile(r"\bclosest\b", re.IGNORECASE),
]


def main() -> int:
    parser = argparse.ArgumentParser(description="Check that the charged-lepton lane does not match labels to experiment.")
    parser.add_argument("targets", nargs="*", help="Files to scan. Defaults to the /particles charged-lepton lane.")
    args = parser.parse_args()

    targets = [pathlib.Path(item) for item in args.targets] if args.targets else DEFAULT_TARGETS
    failures: list[str] = []
    for target in targets:
        text = target.read_text(encoding="utf-8")
        for pattern in FORBIDDEN_PATTERNS:
            if pattern.search(text):
                failures.append(f"{target}: matched forbidden pattern {pattern.pattern!r}")

    if failures:
        for failure in failures:
            print(failure, file=sys.stderr)
        return 1

    print("no experiment-label-matching patterns found")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
