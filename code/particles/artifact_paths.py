"""Stable path references for committed particle receipts.

Scientific artifacts are compared byte-for-byte across clean clones.  They
therefore must not serialize the absolute checkout path of the machine that
generated them.  Builders may still use absolute :class:`~pathlib.Path`
objects internally; this module converts path-shaped strings at the receipt
boundary.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[2]
_WINDOWS_ABSOLUTE = re.compile(r"^[A-Za-z]:[\\/]")
_WORKSPACE_REPOS = ("reverse-engineering-reality", "oph-physics-sim")


def canonical_artifact_ref(value: str | Path) -> str:
    """Return a clone-stable reference for a local artifact path."""

    path = Path(value)
    raw = str(value)
    if not path.is_absolute() and not _WINDOWS_ABSOLUTE.match(raw):
        return path.as_posix()

    if path.is_absolute():
        try:
            return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()
        except (OSError, ValueError):
            pass

    normalized = raw.replace("\\", "/")
    parts = tuple(part for part in normalized.split("/") if part)
    for repo_name in _WORKSPACE_REPOS:
        if repo_name in parts:
            index = parts.index(repo_name)
            suffix = "/".join(parts[index + 1 :])
            return f"oph-workspace://{repo_name}/{suffix}"

    # An explicitly supplied external file is still identifiable without
    # leaking the developer's home directory.
    external_name = normalized.rstrip("/").split("/")[-1]
    return f"external-file://{external_name}"


def canonicalize_artifact_paths(value: Any) -> Any:
    """Recursively replace absolute path strings by stable references."""

    if isinstance(value, dict):
        return {
            key: canonicalize_artifact_paths(item)
            for key, item in value.items()
        }
    if isinstance(value, list):
        return [canonicalize_artifact_paths(item) for item in value]
    if isinstance(value, tuple):
        return [canonicalize_artifact_paths(item) for item in value]
    if isinstance(value, str):
        if Path(value).is_absolute() or _WINDOWS_ABSOLUTE.match(value):
            return canonical_artifact_ref(value)
    return value


def portable_json_dumps(value: Any, **kwargs: Any) -> str:
    """Serialize a receipt only after applying the portability boundary."""

    return json.dumps(canonicalize_artifact_paths(value), **kwargs)
