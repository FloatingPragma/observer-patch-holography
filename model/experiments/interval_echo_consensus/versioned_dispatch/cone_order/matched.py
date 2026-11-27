"""Matched local versions with and without an intervening causal path."""

import argparse
import gzip
import json

from .chart import graph_from_bank
from .measure import HERE, branch, instrument


def matched_words(graph):
    edges = sorted(tuple(sorted(e)) for e in graph.edges)
    for a in edges:
        for c in edges:
            if set(a) & set(c):
                continue
            for b in edges:
                joint = set(b) & set(c)
                if not (set(a) & set(b) and len(joint) == 1):
                    continue
                for x in edges:
                    if set(x) & (set(a) | set(b) | set(c)):
                        continue
                    for e in edges:
                        if set(a) & set(e) or set(e) & set(c) != joint:
                            continue
                        for d in edges:
                            if set(d) & (set(a) | set(c)) or not set(d) & set(e):
                                continue
                            return [a, x, b, c], [a, d, e, c]
    raise ValueError("No admissible comparison")


def acquire():
    flights = json.loads((HERE / "flights.json").read_text())
    linked, separate = matched_words(graph_from_bank(flights["sampled"][0]))
    out = {}
    for name, word, mode in [
        ("linked", linked, "normal"),
        ("separate", separate, "normal"),
        ("parked_bridge", linked, "parked_bridge"),
    ]:
        bench = instrument(word)
        runs = {}
        for tag, level in [("empty", 0), ("loaded", 2)]:
            launch = [0] * 12
            for i in word[0]:
                launch[i] = level
            runs[tag] = [branch(bench, word, launch, coin, mode) for coin in range(16)]
        out[name] = runs
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    raw = (json.dumps(acquire(), separators=(",", ":")) + "\n").encode()
    target = HERE / "matched.json.gz"
    if args.write:
        target.write_bytes(gzip.compress(raw, mtime=0))
    elif gzip.decompress(target.read_bytes()) != raw:
        raise ValueError("Matched acquisition has changed")
