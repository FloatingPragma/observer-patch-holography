#!/usr/bin/env python3
"""Emit continuation-only forward Yukawas from the D12 internal-backread descent."""

from __future__ import annotations

import argparse
import importlib.util
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "flavor" / "quark_d12_internal_backread_descent.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_d12_internal_backread_forward_yukawas.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_forward_builder():
    module_path = ROOT / "particles" / "flavor" / "build_forward_yukawas.py"
    spec = importlib.util.spec_from_file_location("build_forward_yukawas", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"unable to load {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.build_artifact


def build_artifact(payload: dict[str, Any]) -> dict[str, Any]:
    build_forward = _load_forward_builder()
    forward = build_forward(payload)
    forward["artifact"] = "oph_quark_d12_internal_backread_forward_yukawas"
    forward["generated_utc"] = _timestamp()
    forward["scope"] = "D12_continuation_internal_backread_only"
    forward["proof_status"] = "internal_backread_continuation_forward_yukawas_emitted"
    forward["public_promotion_allowed"] = False
    forward["predictive_promotion_allowed"] = False
    forward["source_certified"] = False
    forward["public_surface_candidate_allowed"] = False
    forward["public_surface_candidate_status"] = "not_allowed"
    forward["certification_status"] = (
        "forward_matrix_arithmetic_certified_source_uncertified"
    )
    forward["forward_certification_scope"] = (
        "continuation_sidecar_matrix_arithmetic_only"
    )
    forward["promotion_blockers"] = sorted(
        set(forward.get("promotion_blockers", []))
        | {
            "D12_INTERNAL_BACKREAD_NOT_SOURCE_DERIVED",
            "D12_RAY_VALUE_SOURCE_OPEN",
        }
    )
    forward["source_descent_artifact"] = payload.get("artifact")
    metadata = dict(forward.get("metadata", {}))
    metadata["note"] = (
        "Continuation-only forward Yukawas built from the explicit D12 internal-backread descent. "
        "Its matrix arithmetic is certified on continuation-side inputs, but those internal backreads are not source-derived and cannot promote the public frontier."
    )
    forward["metadata"] = metadata
    return forward


def main() -> int:
    parser = argparse.ArgumentParser(description="Build continuation-only forward Yukawas from the D12 internal-backread descent.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    payload = json.loads(Path(args.input).read_text(encoding="utf-8"))
    artifact = build_artifact(payload)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
