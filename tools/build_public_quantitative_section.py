#!/usr/bin/env python3
"""Generate the bilingual README quantitative-claim blocks."""

from __future__ import annotations

import argparse
from pathlib import Path

import public_surface_claims


def build(root: Path, *, check: bool) -> list[str]:
    root = root.resolve()
    issues, outputs = public_surface_claims.expected_surface_texts(root)
    if issues:
        return issues
    stale: list[str] = []
    for relative, expected in outputs.items():
        path = root / relative
        actual = path.read_text(encoding="utf-8")
        if actual == expected:
            continue
        if check:
            stale.append(f"{relative}: generated quantitative claim block is stale")
        else:
            path.write_text(expected, encoding="utf-8", newline="\n")
    if not check and not stale:
        stale.extend(public_surface_claims.check_repository(root))
    return stale


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=public_surface_claims.ROOT,
        help="repository root (default: root containing this script)",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="fail if either generated block differs from its sources",
    )
    args = parser.parse_args()
    issues = build(args.root, check=args.check)
    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        raise SystemExit(1)
    action = "in sync" if args.check else "generated"
    print(f"public quantitative README sections {action}")


if __name__ == "__main__":
    main()
