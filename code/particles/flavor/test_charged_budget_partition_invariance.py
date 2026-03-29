#!/usr/bin/env python3
"""Ensure the shared charged-budget artifact keeps sector shares explicit and normalized."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "flavor" / "charged_budget_transport.json"
CHARGED_SECTORS = ("u", "d", "e")


def _stream_map(stream: list[dict[str, object]]) -> dict[str, float]:
    mapped: dict[str, float] = {}
    for item in stream:
        refinement = str(item.get("refinement", "snapshot"))
        value = item.get("value")
        if value is None:
            continue
        mapped[refinement] = float(value)
    return mapped


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate shared charged-budget partition invariance.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT), help="Input charged-budget artifact.")
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    beta_payload = dict(payload.get("beta_by_sector_by_refinement", {}))
    beta_streams = {sector: _stream_map(list(beta_payload.get(sector, []))) for sector in CHARGED_SECTORS}
    if not all(beta_streams.values()):
        print("charged-budget artifact is missing sector share streams", file=sys.stderr)
        return 1

    common_refinements = sorted(set.intersection(*(set(stream.keys()) for stream in beta_streams.values())))
    if not common_refinements:
        print("charged-budget artifact has no common refinements", file=sys.stderr)
        return 1

    for refinement in common_refinements:
        total = sum(beta_streams[sector][refinement] for sector in CHARGED_SECTORS)
        if abs(total - 1.0) > 1.0e-12:
            print(f"sector shares do not sum to one at refinement {refinement}", file=sys.stderr)
            return 1

    print("charged-budget shares are explicit and normalized")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
