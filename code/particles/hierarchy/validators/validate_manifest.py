#!/usr/bin/env python3
"""Validate the hierarchy bundle manifest hashes."""

from __future__ import annotations

import hashlib
import json
import pathlib
import sys


def main(manifest_path: str = "manifest.json") -> int:
    manifest_file = pathlib.Path(manifest_path)
    root = manifest_file.resolve().parent
    manifest = json.loads(manifest_file.read_text(encoding="utf-8"))
    checks: list[dict[str, object]] = []
    for item in manifest["files"]:
        path = root / item["path"]
        exists = path.is_file()
        actual_sha = hashlib.sha256(path.read_bytes()).hexdigest() if exists else None
        actual_bytes = path.stat().st_size if exists else None
        checks.append(
            {
                "path": item["path"],
                "exists": exists,
                "sha256_matches": actual_sha == item["sha256"],
                "bytes_match": actual_bytes == item["bytes"],
            }
        )
    payload = {"files": checks, "pass": all(c["exists"] and c["sha256_matches"] and c["bytes_match"] for c in checks)}
    print(json.dumps(payload, indent=2))
    return 0 if payload["pass"] else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "manifest.json"))
