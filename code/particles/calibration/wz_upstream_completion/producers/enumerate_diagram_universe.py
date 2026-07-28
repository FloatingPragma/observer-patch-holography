#!/usr/bin/env python3
"""Workstream C: one-loop diagram universe for the W/Z pole blocks.

The enumerator is blind to physics beyond the canonical rule table and
the component quantum numbers: for each external self-energy block it
scans the table for admissible one-loop topologies,

* bubble: two trilinear vertices (ext1, X, Y) and (ext2, Xc, Yc) with
  Xc/Yc the conjugate lines of X/Y,
* seagull: one quartic vertex (ext1, ext2, X, Xc),
* tadpole: a trilinear (ext1, ext2, S) with a neutral scalar S joined
  by the S propagator to a loop vertex (S, X, Xc); tadpole diagrams
  are enumerated explicitly because the direct FJ contract retains
  them,

and emits every match with its vertex hashes.  Conservation holds by
construction because every vertex is a table record; the enumerator
re-checks it anyway and aborts on any violation.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "producers"))

from vertex_format import FIELD_TABLE, conservation_violations  # noqa: E402

TABLE_PATH = ROOT / "outputs" / "rule_table_engine_a.json"
OUT_PATH = ROOT / "outputs" / "diagram_universe_wz.json"

EXTERNAL_BLOCKS = (("Wp", "Wm"), ("W3", "W3"), ("W3", "B"), ("B", "B"))

CONJUGATE_OVERRIDES = {"Wp": "Wm", "Wm": "Wp", "Gp": "Gm", "Gm": "Gp"}


def conjugate(label: str) -> str:
    if label in CONJUGATE_OVERRIDES:
        return CONJUGATE_OVERRIDES[label]
    if label.endswith("_bar"):
        return label[:-4]
    base = FIELD_TABLE[label]
    if base["fermion"] or base["ghost"]:
        return f"{label}_bar"
    return label


def vertex_hash(entry: dict[str, Any]) -> str:
    payload = json.dumps(entry, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def build_universe() -> dict[str, Any]:
    table = json.loads(TABLE_PATH.read_text(encoding="utf-8"))
    by_multiset: dict[tuple, list[dict[str, Any]]] = {}
    for entry in table["entries"]:
        by_multiset.setdefault(tuple(sorted(entry["fields"])), []).append(entry)

    line_fields = sorted(
        set(FIELD_TABLE)
        | {f"{name}_bar" for name, data in FIELD_TABLE.items() if data["fermion"] or data["ghost"]}
    )
    neutral_scalars = [name for name, data in FIELD_TABLE.items()
                      if data["spin"] == "scalar" and data["charge"] == "0"]

    def lookup(fields: list[str]) -> list[dict[str, Any]]:
        return by_multiset.get(tuple(sorted(fields)), [])

    blocks: dict[str, Any] = {}
    for ext1, ext2 in EXTERNAL_BLOCKS:
        block_name = f"{ext1}{ext2}"
        bubbles = []
        seen_bubbles = set()
        for x in line_fields:
            for y in line_fields:
                pair_key = tuple(sorted([x, y]))
                if pair_key in seen_bubbles:
                    continue
                v1_records = lookup([ext1, x, y])
                v2_records = lookup([ext2, conjugate(x), conjugate(y)])
                if v1_records and v2_records:
                    seen_bubbles.add(pair_key)
                    bubbles.append({
                        "internal": list(pair_key),
                        "vertex_1": sorted(vertex_hash(r) for r in v1_records),
                        "vertex_2": sorted(vertex_hash(r) for r in v2_records),
                    })
        seagulls = []
        seen_seagulls = set()
        for x in line_fields:
            if x in seen_seagulls:
                continue
            records = lookup([ext1, ext2, x, conjugate(x)])
            if records:
                seen_seagulls.add(x)
                seen_seagulls.add(conjugate(x))
                seagulls.append({
                    "internal": sorted({x, conjugate(x)}),
                    "vertex": sorted(vertex_hash(r) for r in records),
                })
        tadpoles = []
        for scalar in neutral_scalars:
            head_records = lookup([ext1, ext2, scalar])
            if not head_records:
                continue
            loops = []
            for x in line_fields:
                loop_records = lookup([scalar, x, conjugate(x)])
                if loop_records:
                    loops.append({
                        "internal": sorted({x, conjugate(x)}),
                        "vertex": sorted(vertex_hash(r) for r in loop_records),
                    })
            if loops:
                tadpoles.append({
                    "propagator": scalar,
                    "head_vertex": sorted(vertex_hash(r) for r in head_records),
                    "loops": loops,
                })
        blocks[block_name] = {
            "externals": [ext1, ext2],
            "bubbles": bubbles,
            "seagulls": seagulls,
            "tadpoles": tadpoles,
            "counts": {
                "bubbles": len(bubbles),
                "seagulls": len(seagulls),
                "tadpole_loops": sum(len(t["loops"]) for t in tadpoles),
            },
        }

    for entry in table["entries"]:
        violations = conservation_violations(entry["fields"])
        if violations:
            raise SystemExit(f"enumerator: conservation violation {violations} in {entry['fields']}")

    universe = {
        "schema": "diagram_universe_wz.v1",
        "table_digest": table["table_digest"],
        "action_subject_digest": table["action_subject_digest"],
        "external_blocks": blocks,
    }
    payload = json.dumps(universe, sort_keys=True, separators=(",", ":"))
    universe["universe_digest"] = hashlib.sha256(payload.encode("utf-8")).hexdigest()
    return universe


def main() -> int:
    universe = build_universe()
    OUT_PATH.write_text(json.dumps(universe, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    counts = {name: block["counts"] for name, block in universe["external_blocks"].items()}
    print(json.dumps({"status": "WROTE", "counts": counts, "universe_digest": universe["universe_digest"][:24]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
