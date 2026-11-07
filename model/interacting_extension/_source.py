"""Immutable source-byte contracts for the finite interacting controls.

The reference commit identifies provenance. Execution depends on the recorded
bytes, so a descendant checkout is admissible when those inputs are unchanged.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
PIN = "3cf3c14810671a16e61b7915f69079a1773a12a2"


def verify_sources(repo: Path, names) -> dict[str, str]:
    """Validate requested inputs and their recorded project import closure.

    This check deliberately runs on every call, including callers with cached
    mathematical data. Git labels alone do not authenticate working-tree bytes.
    """
    repo = Path(repo).resolve()
    manifest = json.loads((HERE / "source_manifest.json").read_text())
    if manifest["source_commit"] != PIN:
        raise ValueError("unexpected interacting source contract")
    if isinstance(names, (str, bytes)):
        raise TypeError("source names must be an iterable of relative paths")
    paths = set(names)
    if not paths:
        raise ValueError("at least one source input is required")
    expected = manifest["source_sha256"]
    pending = list(paths)
    while pending:
        name = pending.pop()
        if not isinstance(name, str) or name not in expected:
            raise ValueError(f"uncontracted source input: {name}")
        for dependency in manifest["import_closure"].get(name, []):
            if dependency not in paths:
                paths.add(dependency)
                pending.append(dependency)
    result = {}
    for name in sorted(paths):
        path = (repo / name).resolve()
        if not path.is_relative_to(repo) or not path.is_file():
            raise ValueError(f"missing or escaping source input: {name}")
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != expected[name]:
            raise ValueError(f"source bytes differ from r2062: {name}")
        result[name] = digest
    return result
