#!/usr/bin/env python3
"""Fail closed when committed scientific receipts leak developer-home paths."""

from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RECEIPT_ROOTS = (
    ROOT / "code" / "particles" / "runs",
    ROOT / "claims",
)
DEVELOPER_HOME_PATTERNS = (
    re.compile(r"/Users/[^/]+/"),
    re.compile(r"/home/[^/]+/"),
    re.compile(r"[A-Za-z]:[\\/]Users[\\/][^\\/]+[\\/]", re.IGNORECASE),
)


@dataclass(frozen=True)
class PortabilityViolation:
    path: Path
    json_pointer: str
    value: str


def _escape_pointer(token: object) -> str:
    return str(token).replace("~", "~0").replace("/", "~1")


def _is_developer_home_path(value: str) -> bool:
    return any(pattern.search(value) for pattern in DEVELOPER_HOME_PATTERNS)


def _walk(value: Any, pointer: str = "") -> Iterable[tuple[str, str]]:
    if isinstance(value, dict):
        for key, item in value.items():
            yield from _walk(item, f"{pointer}/{_escape_pointer(key)}")
    elif isinstance(value, list):
        for index, item in enumerate(value):
            yield from _walk(item, f"{pointer}/{index}")
    elif isinstance(value, str):
        yield pointer or "/", value


def receipt_paths(roots: Iterable[Path] = DEFAULT_RECEIPT_ROOTS) -> list[Path]:
    paths: set[Path] = set()
    for root in roots:
        if root.is_file() and root.suffix == ".json":
            paths.add(root)
        elif root.exists():
            paths.update(root.rglob("*.json"))
    return sorted(paths)


def find_violations(paths: Iterable[Path]) -> list[PortabilityViolation]:
    violations: list[PortabilityViolation] = []
    for path in paths:
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            raise ValueError(f"invalid JSON receipt {path}: {exc}") from exc
        for pointer, value in _walk(payload):
            if _is_developer_home_path(value):
                violations.append(PortabilityViolation(path, pointer, value))
    return violations


def check(paths: Iterable[Path]) -> None:
    violations = find_violations(paths)
    if not violations:
        return

    def display_path(path: Path) -> str:
        try:
            return path.relative_to(ROOT).as_posix()
        except ValueError:
            return path.as_posix()

    rendered = "\n".join(
        f"{display_path(item.path)}{item.json_pointer}: {item.value}"
        for item in violations
    )
    raise ValueError(
        "developer-home path leaked into committed scientific receipts:\n"
        f"{rendered}"
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "paths",
        nargs="*",
        type=Path,
        help="optional receipt files/directories (defaults to particle runs and claims)",
    )
    args = parser.parse_args()
    roots = tuple(path.resolve() for path in args.paths) or DEFAULT_RECEIPT_ROOTS
    paths = receipt_paths(roots)
    check(paths)
    print(f"receipt portability OK ({len(paths)} JSON artifacts)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
