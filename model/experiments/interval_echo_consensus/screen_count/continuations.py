"""Recover all finite-word checkpoint permutations from sharp station transfers."""

from argparse import ArgumentParser
from collections import deque
import gzip
import json
from .apparatus import decode, rail_swaps, raw_word, read_words
from .measure import HERE


def closure(generators):
    identity = tuple(range(len(generators[0])))
    words = {identity: ()}
    queue = deque([identity])
    while queue:
        a = queue.popleft()
        for j, g in enumerate(generators):
            nxt = tuple(g[k] for k in a)
            if nxt not in words:
                words[nxt] = words[a] + (j,)
                queue.append(nxt)
    return words


def acquire():
    with gzip.open(HERE / "checkpoint.json.gz", "rt") as src:
        r = json.load(src)
    assert rail_swaps()["full"] == 4
    labels = [decode(row, 8)[0] for row in r["raw"]]
    assert len(set(labels)) == 256
    ops = []
    for j in range(1, 8):
        transition = [-1] * 256
        for i, origin in enumerate(labels):
            transition[origin] = labels[r["outcomes"][i]["unadjusted"][j]]
        ops.append(transition)
    outside = [None] * 256
    transfer = [-1] * 256
    for dep in (0, 1):
        for flips in range(128):
            origin = labels[dep * 128 + flips]
            raw = read_words(raw_word(r["config"], flips, 1 - dep), r["stations"])
            target, external = decode(raw, 8)
            assert external == dep
            outside[origin] = raw
            transfer[origin] = target
    ops.append(transfer)
    summaries = {}
    for title, gen in [("repeat_inside", ops[:7]), ("repeat_boundary", ops)]:
        maps = closure(gen)
        for g in gen:
            assert all(g[g[i]] == i for i in range(256))
        orbit = {i: frozenset(t[i] for t in maps) for i in range(256)}
        classes = sorted((sorted(c) for c in set(orbit.values())), key=min)
        routes = []
        for f, path in maps.items():
            inv = tuple(range(256))
            for j in reversed(path):
                inv = tuple(gen[j][k] for k in inv)
            assert all(inv[f[i]] == i for i in range(256))
            routes.append(dict(labels=path, terminal=f))
        assert sum(map(len, classes)) == 256
        summaries[title] = dict(classes=classes, routes=routes)
    return dict(boundary_readings=outside, generators=ops, families=summaries)


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    raw = (json.dumps(acquire(), separators=(",", ":")) + "\n").encode()
    f = HERE / "continuations.json.gz"
    if args.write:
        f.write_bytes(gzip.compress(raw, mtime=0))
    elif gzip.decompress(f.read_bytes()) != raw:
        raise ValueError("Continuation records have changed")
