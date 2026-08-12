"""Small shared helpers for fail-closed JSON loading.

Python's default JSON decoder silently accepts duplicate object keys and keeps
the last value.  Registry and receipt validators must reject that ambiguity.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


class DuplicateKeyError(ValueError):
    """Raised when a JSON object repeats a key."""


def _object_without_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    value: dict[str, Any] = {}
    for key, item in pairs:
        if key in value:
            raise DuplicateKeyError(f"duplicate JSON key {key!r}")
        value[key] = item
    return value


def loads(text: str) -> Any:
    """Decode JSON while rejecting duplicate object keys at every depth."""

    return json.loads(text, object_pairs_hook=_object_without_duplicates)


def load(path: Path) -> Any:
    """Read and strictly decode one UTF-8 JSON file."""

    return loads(path.read_text(encoding="utf-8"))
