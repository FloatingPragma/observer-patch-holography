"""Replay labelled, timed round-trip echoes from observer 0 on archived strict writes."""

import argparse
import hashlib
import json
import math
import struct
import zipfile
import zlib
from collections import deque
from functools import lru_cache
from pathlib import Path

import numpy as np

MANIFEST_SHA = "3c0b240ac4464afd9b82251b46a86c6622a87a8a281d0ad159b641413360f824"
ROW = struct.Struct("<H12H6B3I3Bx")
TIMEOUT = 48
EMISSIONS = tuple(range(0, 113, 4))


class ArchiveRadar:
    def __init__(self, root):
        self.root = root
        raw = (root / "manifest.json").read_bytes()
        if hashlib.sha256(raw).hexdigest() != MANIFEST_SHA:
            raise ValueError("Expected published level-13 history")
        self.manifest = json.loads(raw)
        with zipfile.ZipFile(root / "source.zip") as z:
            seams = z.read("strict_history/seams.json")
        if (
            hashlib.sha256(seams).hexdigest()
            != self.manifest["provenance"]["internal_seam_sha256"]
        ):
            raise ValueError("Internal seam table differs from archive")
        self.internal = [[] for _ in range(12)]
        for a, b in json.loads(seams):
            self.internal[a].append(b)
            self.internal[b].append(a)

    @lru_cache(maxsize=256)
    def frame(self, block, frame):
        entry = self.manifest["shards"][block]
        with (self.root / entry["file"]).open("rb") as f:
            f.seek(4 * frame)
            start, end = struct.unpack("<II", f.read(8))
            f.seek(start)
            return zlib.decompress(f.read(end - start))

    @lru_cache(maxsize=8192)
    def row(self, owner):
        block, i = divmod(owner, self.manifest["observers_per_shard"])
        frame, offset = divmod(i, self.manifest["observers_per_frame"])
        return ROW.unpack_from(self.frame(block, frame), offset * ROW.size)

    def links(self, station):
        owner, p = divmod(station, 12)
        r = self.row(owner)
        mate = (r[13 + p // 2] >> (4 * (p % 2))) & 15
        tick = r[1 + p]
        for q in self.internal[p]:
            yield owner * 12 + q, tick if mate == q else 0
        for peer, pq in zip(r[19:22], r[22:25]):
            if pq & 15 == p:
                yield peer * 12 + (pq >> 4), tick if mate == 12 else 0


def local_graph(links, origin, timeout):
    # Any round trip of at most 2r links stays in the all-transmitting r-ball.
    depth = {origin: 0}
    queue = deque([origin])
    while queue:
        u = queue.popleft()
        if depth[u] == timeout // 2:
            continue
        for v, _ in links(u):
            if v not in depth:
                depth[v] = depth[u] + 1
                queue.append(v)
    ports = list(depth)
    index = {u: i for i, u in enumerate(ports)}
    edges = [{index[v]: t for v, t in links(u) if v in index} for u in ports]
    if any(edges[v].get(u) != t for u, row in enumerate(edges) for v, t in row.items()):
        raise ValueError("Reciprocal link or recorded closure disagreement")
    return ports, [depth[u] for u in ports], edges


def outward(edges, emit, timeout):
    distances = [-1] * len(edges)
    distances[0] = 0
    queue = deque([0])
    while queue:
        u = queue.popleft()
        d = distances[u]
        if d >= timeout:
            continue
        for v, closure in edges[u].items():
            if distances[v] < 0 and (not closure or emit + d < closure):
                distances[v] = d + 1
                queue.append(v)
    return distances


def echo_table(edges, targets, emissions, timeout, settled):
    # Earliest arrival dominates later arrivals because links only close.
    # return_wait[t,v] is the shortest remaining time to the sender at time t.
    miss = timeout + 1
    static = [v if v >= 0 else miss for v in outward(edges, settled, timeout)]
    returns = np.full((settled + 1, len(edges)), miss, np.int32)
    returns[-1] = static
    for t in range(settled - 1, -1, -1):
        returns[t, 0] = 0
        for v in range(1, len(edges)):
            returns[t, v] = min(
                miss,
                1
                + min(
                    (
                        int(returns[t + 1, u])
                        for u, closure in edges[v].items()
                        if not closure or t < closure
                    ),
                    default=miss,
                ),
            )
    table = []
    for t in emissions:
        out = outward(edges, t, timeout)
        row = []
        for v in targets:
            d = out[v]
            echo = d + int(returns[min(settled, t + d), v]) if d >= 0 else miss
            row.append(echo if echo <= timeout else None)
        table.append(row)
    return table


def analyze(root):
    source = ArchiveRadar(root)
    ports, depth, edges = local_graph(source.links, 0, TIMEOUT)
    targets = [v for v, p in enumerate(ports) if p // 12 != 0 and depth[v] <= 12]
    echoes = echo_table(
        edges, targets, EMISSIONS, TIMEOUT, source.manifest["settlement_round"]
    )
    cohorts = []
    for low, high in ((1, 4), (5, 8), (9, 12)):
        enrolled = [i for i, v in enumerate(targets) if low <= depth[v] <= high]
        visible = [[i for i in enrolled if row[i] is not None] for row in echoes]
        fixed = [i for i in enrolled if all(row[i] is not None for row in echoes)]
        if not fixed:
            raise ValueError("Empty fixed cohort")
        total = [sum(row[i] for i in fixed) for row in echoes]
        scale = [s / total[0] for s in total]
        step = EMISSIONS[1] - EMISSIONS[0]

        def estimates(k):
            return dict(
                mid_emission=EMISSIONS[k + 1],
                log_rate=(math.log(scale[k + 2]) - math.log(scale[k])) / (2 * step),
                relative_acceleration=(scale[k + 2] - 2 * scale[k + 1] + scale[k])
                / (step**2 * scale[k + 1]),
            )

        cohorts.append(
            dict(
                blank_one_way_range=[low, high],
                enrolled=len(enrolled),
                always_returned=len(fixed),
                returned_per_emission=[len(x) for x in visible],
                round_trip_tick_sums=total,
                range_scale=scale,
                first=estimates(0),
                last=estimates(len(EMISSIONS) - 3),
            )
        )
    return dict(
        manifest_sha256=MANIFEST_SHA,
        dataset="https://dovybzrnzxfo1.cloudfront.net/manifest.json",
        origin_observer=0,
        origin_port=0,
        emission_ticks=EMISSIONS,
        round_trip_timeout=TIMEOUT,
        instrument="Matched repair histories. Each calibrated material relay takes one supplied clock tick. An edge passes on departure before its recorded write. A labelled reflector immediately answers.",
        surveyed_ports=len(ports),
        cohorts=cohorts,
        reflectors=[[ports[v] // 12, ports[v] % 12, depth[v]] for v in targets],
        echo_waits=echoes,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("archive", type=Path)
    print(json.dumps(analyze(p.parse_args().archive), indent=1, allow_nan=False))
