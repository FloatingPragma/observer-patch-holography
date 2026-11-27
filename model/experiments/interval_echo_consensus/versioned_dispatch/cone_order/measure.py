"""Acquire flight counts, sharp load readings and endpoint writer stamps."""

import argparse
from functools import cache
import gzip
import json
from pathlib import Path

import numpy as np
from ...seam_word_meter.finite_program.machine import (
    Bench,
    LAW,
    READ_CODES,
    decode_word,
    template_matching,
)
from .chart import arrival_counts, graph_from_bank
from .order import census

HERE = Path(__file__).parent


@cache
def write_columns():
    pattern = LAW._write_pattern()
    assert np.array_equal(pattern.T @ pattern, 2 * np.eye(9, dtype=int))
    cols = {}
    for k, pair in enumerate(LAW.SRC):
        outcomes = {}
        for ix in np.flatnonzero(pattern[:, k]):
            mem, payload = divmod(int(ix), len(LAW.SRC))
            coin, contrast = LAW.LABELS[mem]
            assert contrast == pair[0] - pair[1] and abs(pattern[ix, k]) == 1
            outcomes[coin] = LAW.SRC[payload]
        assert set(outcomes) == {0, 1}
        cols[pair] = outcomes
    return cols


def instrument(word, width=3):
    matching = template_matching(set(range(12)))
    storage = [
        (12 * owner + a, 12 * owner + b)
        for owner in range(1, 2 * width + 1)
        for a, b in matching
    ]
    bench = Bench(words=[storage[k : k + width] for k in range(0, len(storage), width)])
    # Readout swaps avoid all event contacts, including occupied memories.
    bench.graph.remove_edges_from(word)
    return bench


def snapshot(bench, at):
    readings = [bench.observe(i, at) for i in range(12)]
    stamps = [bench.survey(w, at) for w in bench.words]
    return {
        "load": readings,
        "stamp": [v for word in stamps for pair in word for v in pair],
    }


def branch(bench, word, launch, coins, mode, muted=None):
    bench.start([0] * 12)
    bench.loads[:12] = launch

    def pair():
        return [snapshot(bench, at) for at in (0, 1)]

    scans = [pair()]
    for i, (a, b) in enumerate(word):
        # Each addressed seam is blank. A sharp read on its endpoints
        # conditions the native integer-column ensemble on the two loads.
        # The hidden equal-weight memory branches are never used as parents.
        if (
            (a, b) != muted
            and mode != "parked"
            and not (mode == "parked_bridge" and i == 2)
        ):
            bench.loads[a], bench.loads[b] = write_columns()[
                bench.loads[a], bench.loads[b]
            ][(coins >> i) & 1]
        bench.put(a, i + 1, 0)
        if not (mode == "shutter" and i == len(word) - 1):
            bench.put(b, i + 1, 1)
        scans.append(pair())
    return scans


def decode(scans, width=3):
    versions = [0] * 12
    word, parents, loads = [], [], []
    for event, (left, right) in enumerate(scans):
        if left != right or len(left["stamp"]) != 24 * width or len(left["load"]) != 12:
            raise ValueError("Unmatched station records")
        bits = left["stamp"]
        tags = [
            decode_word([bits[j : j + 2] for j in range(k, k + 2 * width, 2)])
            for k in range(0, len(bits), 2 * width)
        ]
        if any(t is None for t in tags) or any(
            t not in READ_CODES for t in left["load"]
        ):
            raise ValueError("Unresolved material digit")
        current = [READ_CODES.index(t) for t in left["load"]]
        if event:
            changed = tuple(
                k for k, (x, y) in enumerate(zip(tags, versions, strict=True)) if x != y
            )
            if len(changed) != 2 or any(tags[k] != event for k in changed):
                raise ValueError("Incomplete endpoint commit")
            if any(current[k] != loads[-1][k] for k in range(12) if k not in changed):
                raise ValueError("Unstamped payload change")
            word.append(changed)
            parents.append(
                tuple(sorted({versions[k] - 1 for k in changed if versions[k]}))
            )
        elif any(tags):
            raise ValueError("Nonblank initial versions")
        versions = tags
        loads.append(current)
    return word, tuple(parents), loads


def acquire():
    flights = arrival_counts()
    graph = graph_from_bank(flights["sampled"][0])
    if sorted(graph.edges) != sorted(graph_from_bank(flights["sampled"][1]).edges):
        raise ValueError("Receiver frame disagreement")
    scan = census(graph)
    panels = {}
    for name, example in scan["witnesses"].items():
        word = list(map(tuple, example["edges"]))
        bench = instrument(word)
        rows = []
        for launch in ([1] * 12, [k % 3 for k in range(12)]):
            for coins in range(2 ** len(word)):
                rows.append(branch(bench, word, launch, coins, "normal"))
        control = [
            branch(bench, word, [k % 3 for k in range(12)], 0, mode)
            for mode in ("parked", "shutter")
        ]
        panels[name] = {"ordinary": rows, "parked": control[0], "shutter": control[1]}
    return flights, scan, panels


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    flights, scan, panels = acquire()
    for filename, obj in [
        ("flights.json", flights),
        ("census.json", scan),
        ("records.json.gz", panels),
    ]:
        file = HERE / filename
        raw = (json.dumps(obj, separators=(",", ":")) + "\n").encode()
        if args.write:
            file.write_bytes(
                gzip.compress(raw, mtime=0) if file.suffix == ".gz" else raw
            )
        else:
            existing = (
                gzip.decompress(file.read_bytes())
                if file.suffix == ".gz"
                else file.read_bytes()
            )
            if raw != existing:
                raise ValueError(f"Changed acquisition in {filename}")


if __name__ == "__main__":
    main()
