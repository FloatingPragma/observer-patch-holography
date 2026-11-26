"""Measure correctable messages using eight unannounced native tree continuations."""

from argparse import ArgumentParser
import gzip
import json
from functools import reduce
from operator import xor

import networkx as nx
import numpy as np

from .apparatus import calibrate, decode, layout, rail_swaps, raw_word, read_words
from .measure import HERE, edges


def syndrome(bits):
    return reduce(xor, (j + 1 for j in range(7) if bits & (1 << j)), 0)


def correct(bits):
    j = syndrome(bits)
    return bits ^ (1 << (j - 1)) if j else bits


def acquire():
    assert rail_swaps() == {"parked": 0, "third": 1, "half": 2, "full": 4}
    rng = np.random.Generator(np.random.PCG64(1262604))
    c = layout(0, size=8)
    grid = nx.Graph(edges(c))
    forest = list(nx.bfs_edges(grid, 0, sort_neighbors=sorted))
    all_edges = edges(c)
    c["inside"] = [
        c["inside"][next(i for i, e in enumerate(all_edges) if set(e) == {u, v})]
        for u, v in forest
    ]
    c.pop("skip")
    stations = calibrate(c)
    patterns = []
    inverse = {}
    for dep in (0, 1):
        for flips in range(128):
            raw = read_words(raw_word(c, flips, dep), stations)
            w, external = decode(raw, 8)
            assert external == 1 - dep and w not in inverse
            inverse[w] = (dep, flips)
            patterns.append(raw)
    assert len(inverse) == 256
    calendars = list(range(8))

    def move(bits, position):
        return bits if position == 0 else bits ^ (1 << (position - 1))

    # Pointer banks are exact sharp sectors. A code below selects an entire two-head
    # read, including all ninety-six source-port reads and the optional exterior rail.
    outcomes = []
    for dep in (0, 1):
        for initial in range(128):
            word, _ = decode(patterns[128 * dep + initial], 8)
            sector, measured = inverse[word]
            assert sector == dep
            tables = {}
            for setting, source in (
                ("unadjusted", initial),
                ("adjusted", correct(measured)),
            ):
                tables[setting] = [128 * dep + move(source, j) for j in calendars]
            outcomes.append(tables)
    banks = []
    for dep_mode in ("parked", "half"):
        branch = np.full(256, 1 / 256)
        if dep_mode == "parked":
            branch[:128] = 1 / 128
            branch[128:] = 0
        for setting in ("unadjusted", "adjusted"):
            for j in calendars:
                started = rng.multinomial(2**17, branch)
                finished = np.zeros(256, dtype=int)
                for ix, n in enumerate(started):
                    end = outcomes[ix][setting][j]
                    word, _ = decode(patterns[end], 8)
                    finished[word] += n
                banks.append(
                    dict(
                        boundary=dep_mode,
                        source=setting,
                        continuation=j,
                        starts=started.tolist(),
                        terminal_words=finished.tolist(),
                    )
                )
    # Changing to diagonal contacts parks the intended transfer; opening two real
    # faulty contacts tests an out-of-menu alternative with the same sharp receiver.
    return dict(
        config=c,
        stations=stations,
        raw=patterns,
        outcomes=outcomes,
        banks=banks,
        controls=dict(
            ordinary=outcomes[0]["adjusted"][1],
            diagonal=0,
            two_errors=3,
            other_parent=7,
            other_error=3,
        ),
    )


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    target = HERE / "checkpoint.json.gz"
    raw = (json.dumps(acquire(), separators=(",", ":")) + "\n").encode()
    if args.write:
        target.write_bytes(gzip.compress(raw, mtime=0))
    elif gzip.decompress(target.read_bytes()) != raw:
        raise ValueError("Checkpoint receipt changed")
