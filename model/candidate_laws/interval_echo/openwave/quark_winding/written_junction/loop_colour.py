"""Material-only distinguishability inside intact joint matter and head components."""

import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np
from scipy import sparse

from ...electron_field.flux_string.search import Patch
from ..belt_walls import belt


def component(patch, start, seed):
    # Fixed local Gauss data make the head word determine all material occupations.
    active = [((1 << a) | (1 << b), 1 << a, h) for (a, b), h in patch.bit.items()]
    fields = [seed]
    matter = [start]
    index = {seed: 0}
    for i, word in enumerate(fields):
        occ = matter[i]
        for ab, a, h in active:
            local = occ & ab
            if local and local != ab and ((local == a) == bool(word & h)):
                new = word ^ h
                if new not in index:
                    index[new] = len(fields)
                    fields.append(new)
                    matter.append(occ ^ ab)
                else:
                    assert matter[index[new]] == occ ^ ab
        for _, new in patch.ring_moves(word):
            if new not in index:
                index[new] = len(fields)
                fields.append(new)
                matter.append(occ)
            else:
                assert matter[index[new]] == occ
    return np.array(fields, dtype=np.int64), np.array(matter, dtype=np.int32)


def advance(patch, occ, head, edge, write):
    a, b = edge
    bits = (1 << a) | (1 << b)
    if occ & bits in (0, bits):
        return occ, head
    start = a if occ & 1 << a else b
    dest, word = (
        patch.move(start, head, edge, 1) if write else (b if start == a else a, head)
    )
    return (occ ^ bits, word) if dest != start else (occ, head)


def witness(patch, fields, matter, previous):
    from collections import deque

    counts = np.bincount(previous)
    cls = int(np.flatnonzero(counts > 1)[0])
    i, j = np.flatnonzero(previous == cls)[:2]
    occ = int(matter[i])
    assert occ == matter[j]
    x, y = int(fields[i]), int(fields[j])
    queue = deque([(occ, x, y, [])])
    seen = {(occ, x, y)}
    while queue:
        m, w, v, path = queue.popleft()
        for edge in patch.edges:
            a, aa = advance(patch, m, w, edge, True)
            b, bb = advance(patch, m, v, edge, True)
            route = path + [edge]
            if a != b:
                result = {}
                for name, write in [("native", True), ("head_independent", False)]:
                    traces = []
                    for initial in (x, y):
                        n, h = occ, initial
                        record = [n]
                        for e in route:
                            n, h = advance(patch, n, h, e, write)
                            record.append(n)
                        traces.append(record)
                    result[name] = traces
                return dict(initial_heads=[x, y], seams=route, material_masks=result)
            if a != m and (a, aa, bb) not in seen:
                seen.add((a, aa, bb))
                queue.append((a, aa, bb, route))
    raise ValueError("Missing distinguishing record")


def refine(patch, fields, matter, *, include_witness=True):
    # Sorting accelerates the exact lookup of each allowed destination.
    order = np.argsort(fields)
    fields, matter = fields[order], matter[order]
    _, classes = np.unique(matter, return_inverse=True)
    sizes = []
    moves = []
    for (a, b), h in patch.bit.items():
        aa = (matter >> a) & 1
        bb = (matter >> b) & 1
        bit = (fields & h) != 0
        eligible = (aa != bb) & (aa == bit)
        ids = np.flatnonzero(eligible)
        dest = np.searchsorted(fields, fields[ids] ^ h)
        assert np.array_equal(fields[dest], fields[ids] ^ h)
        moves.append((ids, dest))
    history_before = None
    while True:
        last = classes
        histogram = Counter(np.bincount(classes).tolist())
        sizes.append(
            dict(
                classes=len(np.unique(classes)),
                largest=max(histogram),
                histogram={str(k): v for k, v in sorted(histogram.items())},
            )
        )
        if sizes[-1]["largest"] == 1:
            break
        # Each split uses the previous round's class at every destination.
        # Encoding one key at a time avoids a rectangular n x 30 history table.
        history_before = last
        for ids, dest in moves:
            observable = np.zeros(len(fields), np.int64)
            observable[ids] = last[dest] + 1
            keys = classes.astype(np.int64) * (int(last.max()) + 2) + observable
            _, classes = np.unique(keys, return_inverse=True)
        if int(classes.max()) == int(last.max()):
            break
    rows, cols = (np.concatenate(seq) for seq in zip(*moves))
    graph = sparse.csr_matrix(
        (np.ones(len(rows), np.int8), (rows, cols)), shape=(len(fields), len(fields))
    )
    count, components = sparse.csgraph.connected_components(graph, directed=False)
    partitions = sorted(np.bincount(components).tolist(), reverse=True)
    native = dict(
        components=count,
        dimension_histogram={str(k): v for k, v in Counter(partitions).items()},
    )
    return (
        sizes,
        sum(len(i) for i, _ in moves),
        witness(patch, fields, matter, history_before) if include_witness else None,
        native,
    )


def wall_source(patch):
    order = next(iter(belt.belts().values()))
    seed = patch.seed
    for i in range(10):
        a, b = order[i], order[(i + 1) % 10]
        h = patch.bit[tuple(sorted((a, b)))]
        seed = (seed | h) if ((a < b) == (i % 3 == 0)) else (seed & ~h)
    return order, seed


def run():
    patch = Patch()
    order, seed = wall_source(patch)
    panels = {}
    for name, s in [("minus", (0, 4, 7)), ("plus", (0, 1, 4, 7))]:
        occ = sum(1 << order[i] for i in s)
        words, material = component(patch, occ, seed)
        sizes, transfers, example, native = refine(patch, words, material)
        panels[name] = dict(
            source_ports=[order[i] for i in s],
            joint_states=len(words),
            directed_material_transfers=transfers,
            history_refinements=sizes,
            witness=example,
            without_face_flip=native,
            head_blind_maximum_indistinguishable=sizes[0]["largest"],
        )
    return dict(common_field=seed, panels=panels)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args()
    encoded = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    with args.output.open("x") as f:
        f.write(encoded)
