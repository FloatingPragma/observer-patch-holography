"""Compile the declared golden read menu and a fixed injective host placement."""
from __future__ import annotations

import argparse
from functools import cmp_to_key
import hashlib
import importlib.util
import json
from math import isqrt
from pathlib import Path
import struct

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def source_definition():
    spec = importlib.util.spec_from_file_location("routing_source_definition", ROOT / "evidence/source_net_causal_poset/build_causal_poset.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def prepare(q, output):
    if q not in (3, 13, 21):
        raise ValueError("declared q values are 3 (control), 13 and 21")
    source = source_definition()
    sites = source.site_coordinates(q, 3)
    orbit = source.orbit(q)
    a, b = source.axis_tables(orbit)
    offsets, indices = source.site_graph(q, 3, orbit, a, b, sites)
    # Rank the exact golden coordinates, then use a three-axis Morton order.
    # Both the placement and the routing policy are supplied M1 data.
    order = sorted(range(q), key=cmp_to_key(lambda i, j: source.phi_sign(
        (orbit[i][0] - orbit[j][0], orbit[i][1] - orbit[j][1]))))
    rank = {v: i for i, v in enumerate(order)}
    def morton(row):
        return sum(((rank[int(row[axis])] >> bit) & 1) << (3 * bit + axis)
                   for bit in range(q.bit_length()) for axis in range(3))
    ordered_sites = sorted(range(q**3), key=lambda i: morton(sites[i]))
    hosts = np.empty(q**3, dtype="<u4")
    hosts[ordered_sites] = np.arange(q**3, dtype="<u4")
    level = 3 if q == 3 else 4 if q == 13 else 5
    support_path = (ROOT / "code/source_routing/support_w12_l3.json") if level == 3 else HERE / f"support_w12_l{level}.json"
    support = json.loads(support_path.read_text(encoding="utf-8"))
    centre_axis = min(range(q), key=cmp_to_key(lambda i, j: source.phi_sign(source.phi_sub(
        source.phi_square((orbit[i][0] * 2 - 1, orbit[i][1] * 2)),
        source.phi_square((orbit[j][0] * 2 - 1, orbit[j][1] * 2))))))
    centre = centre_axis * (q*q + q + 1)
    records, *_ = source.source_records(q, sites)
    rounds = isqrt(q) + (isqrt(q)**2 < q)
    pairs = np.asarray(support["glued_pairs"], dtype="<u4")
    header = struct.pack("<8Q", q, q**3, support["carriers"], rounds, centre,
                         len(pairs), len(indices), 1)
    data = (header + hosts.tobytes() + pairs.tobytes() + offsets.astype("<u8").tobytes()
            + indices.astype("<u4").tobytes() + records.astype("<i8").tobytes())
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes(data)
    info = {"q": q, "sites": q**3, "carriers": support["carriers"], "rounds": rounds,
            "centre_site": centre, "logical_reads": len(indices)*rounds,
            "neighbour_sha256": source.neighbour_digest(offsets, indices),
            "input_sha256": hashlib.sha256(data).hexdigest(),
            "support_path": support_path.relative_to(ROOT).as_posix(),
            "support_sha256": hashlib.sha256(support_path.read_bytes()).hexdigest()}
    print(json.dumps(info, sort_keys=True), flush=True)
    return info


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--q", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.q, args.output)
