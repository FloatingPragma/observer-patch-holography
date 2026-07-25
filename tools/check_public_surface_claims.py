#!/usr/bin/env python3
"""Fail-closed check for public OPH quantitative comparison surfaces."""

from __future__ import annotations

import argparse
from pathlib import Path

import public_surface_claims


def check(root: Path) -> list[str]:
    return public_surface_claims.check_repository(root.resolve())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=public_surface_claims.ROOT,
        help="repository root (default: root containing this script)",
    )
    args = parser.parse_args()
    issues = check(args.root)
    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        raise SystemExit(1)
    print(
        "public quantitative surfaces OK: registry classes, producers, "
        "comparisons, and generated blocks resolve"
    )


if __name__ == "__main__":
    main()
