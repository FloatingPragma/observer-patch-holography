#!/usr/bin/env python3
"""Check the shared-budget reconstruction identity on all charged-sector streams."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "flavor" / "charged_budget_transport.json"


def _stream_map(stream: list[dict[str, object]]) -> dict[str, float]:
    out: dict[str, float] = {}
    for item in stream:
        refinement = str(item.get("refinement", "snapshot"))
        value = item.get("value")
        if value is None:
            continue
        out[refinement] = float(value)
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate the charged-budget reconstruction identity.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT), help="Input charged-budget artifact.")
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    b_by_sector = {
        sector: _stream_map(list(dict(payload.get("B_by_sector_by_refinement", {})).get(sector, [])))
        for sector in ("u", "d", "e")
    }
    g_by_sector = {
        sector: _stream_map(list(dict(payload.get("g_by_sector_by_refinement", {})).get(sector, [])))
        for sector in ("u", "d", "e")
    }
    beta_by_sector = {
        sector: _stream_map(list(dict(payload.get("beta_by_sector_by_refinement", {})).get(sector, [])))
        for sector in ("u", "d", "e")
    }
    b_total = _stream_map(list(payload.get("B_ch_by_refinement", [])))
    common = sorted(
        set(b_total)
        & set.intersection(*(set(stream) for stream in b_by_sector.values()))
        & set.intersection(*(set(stream) for stream in beta_by_sector.values()))
    )
    if not common:
        print("charged-budget reconstruction identity has no common refinements", file=sys.stderr)
        return 1
    for refinement in common:
        for sector in ("u", "d", "e"):
            lhs = b_by_sector[sector][refinement]
            rhs = beta_by_sector[sector][refinement] * b_total[refinement]
            if abs(lhs - rhs) > 1.0e-12:
                print(
                    f"reconstruction identity failed for sector {sector} at refinement {refinement}",
                    file=sys.stderr,
                )
                return 1
            if g_by_sector[sector] and refinement in g_by_sector[sector] and abs(lhs - g_by_sector[sector][refinement]) > 1.0e-12:
                print(
                    f"g_{sector} stream disagrees with B_{sector} at refinement {refinement}",
                    file=sys.stderr,
                )
                return 1
    print("shared-budget reconstruction identity passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
