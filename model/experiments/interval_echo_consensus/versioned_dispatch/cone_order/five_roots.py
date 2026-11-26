"""Read a common five-root prefix and test every subsequent set of open routes."""

import argparse
from collections import Counter
from functools import cache
import gzip
from itertools import product
import json
from math import log, sqrt

import networkx as nx
import numpy as np

from .chart import graph_from_bank
from .measure import HERE, READ_CODES, branch, decode, instrument, write_columns

SHOTS = 2**16
RADIUS = sqrt(log(2 * (6 * 32 + 1) / 0.01) / (2 * SHOTS))


def matchings(graph, ports):
    if not ports:
        yield ()
        return
    a = min(ports)
    for b in sorted(set(graph[a]) & ports):
        for tail in matchings(graph, ports - {a, b}):
            yield ((a, b), *tail)


def apparatus(graph):
    for matching in matchings(graph, set(graph)):
        for hub in matching:
            roots = tuple(e for e in matching if e != hub)
            candidates = [
                sorted(
                    tuple(sorted((a, b)))
                    for a in r
                    for b in hub
                    if graph.has_edge(a, b)
                )
                for r in roots
            ]
            for spokes in product(*candidates):
                bench = instrument([*matching, *spokes], width=4)
                try:
                    for port in range(12):
                        for head in bench.heads:
                            bench.route(port, head)
                except nx.NetworkXNoPath:
                    continue
                return roots, spokes, hub, bench
    raise ValueError("No admissible native contacts")


def past_of_hub(scans):
    word, parents, loads = decode(scans, width=4)
    ancestors = []
    for ps in parents:
        ancestors.append(set(ps).union(*(ancestors[j] for j in ps)))
    return word, sorted(ancestors[-1] & set(range(5))), loads


def detector(bench, hub):
    @cache
    def read(a, b):
        bench.start([0] * 12)
        bench.loads[hub[0]], bench.loads[hub[1]] = a, b
        signals = [[bench.observe(p, j) for p in hub] for j in (0, 1)]
        if signals[0] != signals[1]:
            raise ValueError("Gathering receivers disagree")
        return tuple(signals[0])

    return read


def outcomes(roots, word, hub, root_id, read, muted=None):
    launch = [0] * 12
    if root_id is not None:
        for p in roots[root_id]:
            launch[p] = 2
    distribution = {tuple(launch): 1}
    for edge in word:
        if edge == muted:
            continue
        a, b = edge
        following = Counter()
        for loads, count in distribution.items():
            for x, y in write_columns()[loads[a], loads[b]].values():
                dest = list(loads)
                dest[a], dest[b] = x, y
                following[tuple(dest)] += count
        distribution = following
    raw = Counter()
    for loads, count in distribution.items():
        raw[read(loads[hub[0]], loads[hub[1]])] += count
    assert sum(raw.values()) == 2 ** (len(word) - int(muted is not None))
    # In the one-root intervention every supported terminal total is 0 or 1.
    bins = sorted(raw)
    totals = [sum(READ_CODES.index(v) for v in pair) for pair in bins]
    assert all(t in (0, 1) for t in totals)
    counts = [raw[pair] for pair in bins]
    return bins, totals, counts


def acquire():
    frame = json.loads((HERE / "flights.json").read_text())["sampled"]
    roots, spokes, hub, bench = apparatus(graph_from_bank(frame[0]))
    read = detector(bench, hub)
    rng = np.random.Generator(np.random.PCG64(1652602))
    assert set(write_columns()[1, 1].values()) == {(1, 1)}
    rows = []

    def collect(word, root_id, muted=None):
        bins, totals, counts = outcomes(roots, word, hub, root_id, read, muted)
        weights = np.array(counts) / sum(counts)
        sampled = rng.multinomial(SHOTS, weights)
        if abs(np.dot(totals, weights) - np.dot(totals, sampled) / SHOTS) > RADIUS:
            raise ValueError("Receiver bank differs from its predicted mean")
        return dict(root=root_id, bins=bins, branches=counts, sampled=sampled.tolist())

    for mask in range(32):
        ids = [i for i in range(5) if mask & (1 << i)]
        word = [*roots, *(spokes[i] for i in ids), hub]
        trace = branch(bench, word, [1] * 12, 0, "normal")
        inferred, parents, _ = past_of_hub(trace)
        if (
            inferred != word
            or parents != ids
            or (rows and trace[:6] != rows[0]["trace"][:6])
        ):
            raise ValueError("Common prefix or gathering ancestry failed")
        rows.append(
            dict(trace=trace, banks=[collect(word, i) for i in (*range(5), None)])
        )
    full = [*roots, *spokes, hub]
    muted_trace = branch(bench, full, [1] * 12, 0, "normal", muted=spokes[0])
    assert muted_trace == rows[-1]["trace"]
    # On a uniform host parking the echo preserves the two writer-stamp updates.
    assert write_columns()[1, 1] == {0: (1, 1), 1: (1, 1)}
    muted = collect(full, 0, muted=spokes[0])
    return dict(
        roots=roots,
        spokes=spokes,
        hub=hub,
        rows=rows,
        muted=muted,
        muted_trace=muted_trace,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    raw = (json.dumps(acquire(), separators=(",", ":")) + "\n").encode()
    target = HERE / "five_roots.json.gz"
    if args.write:
        target.write_bytes(gzip.compress(raw, mtime=0))
    elif gzip.decompress(target.read_bytes()) != raw:
        raise ValueError("Five-root acquisition has changed")
