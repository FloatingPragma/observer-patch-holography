#!/usr/bin/env python3
"""Normalize a raw HDF5 backend export bundle into the frozen JSON production dump."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from particles.hadron.backend_export_bundle import load_backend_input_artifact
from particles.hadron.production_execution_support import build_backend_manifest, build_production_dump


def _load_json(path: str | Path) -> dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def main() -> int:
    parser = argparse.ArgumentParser(description="Normalize a raw hadron backend export bundle.")
    parser.add_argument("--receipt", required=True)
    parser.add_argument("--payload", required=True)
    parser.add_argument("--backend-bundle", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--manifest-output", required=True)
    args = parser.parse_args()

    receipt = _load_json(args.receipt)
    payload = _load_json(args.payload)
    backend_input = load_backend_input_artifact(args.backend_bundle)
    dump = build_production_dump(receipt, payload, backend_input)
    manifest = build_backend_manifest(receipt, payload, backend_input, backend_input_path=args.backend_bundle)

    Path(args.output).write_text(json.dumps(dump, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    Path(args.manifest_output).write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {args.output}")
    print(f"wrote {args.manifest_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
