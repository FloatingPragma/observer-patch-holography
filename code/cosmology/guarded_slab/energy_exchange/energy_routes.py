"""Addressed energy-exchange and coherent readout paths on the native board."""

from itertools import combinations, product
from pathlib import Path
import sys

import networkx as nx

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from guarded_release import route
from guarded_transfer import ROOT, dimer, federation

sys.path.insert(0, str(ROOT / "code/capacity_readback"))
from correctable_public_record_capacity import (
    compound_confusability_graph,
    maximum_independent_set,
)


def local_layout(site, blocked, paths):
    # The helper rails may be reused in sequential assays; occupied support may not.
    off = 12 * site
    lanes = [tuple(off + j for j in path) for path in paths]
    for cold, warm in product(lanes, repeat=2):
        if set(cold[:2]) & blocked:
            continue
        if set(warm[:2]) & (blocked | set(cold)):
            continue
        if set(warm[2:]) & set(cold[:2]):
            continue
        for select in lanes:
            if set(select[:2]) & (blocked | set(cold) | set(warm)):
                continue
            if select[2] in blocked | set(cold[:2]) | set(warm[:2]):
                continue
            return dict(probe=cold, bath=warm, selector=select[:3])
    raise ValueError("No independent probe, reservoir and change reader")


def transposition(graph, board, left, right):
    paths = [nx.shortest_path(graph, a, b) for a, b in zip(left, right, strict=True)]
    labels = list(range(board.n))
    expected = labels.copy()
    for a, b in zip(left, right, strict=True):
        expected[a], expected[b] = expected[b], expected[a]
    for path in paths:
        edges = list(zip(path, path[1:]))
        for a, b in edges + edges[:-1][::-1]:
            labels[a], labels[b] = labels[b], labels[a]
    assert labels == expected
    return paths


def apparatus(level):
    tag = route(level)
    board = federation(level)
    graph = nx.Graph()
    graph.add_edges_from(board.seams)
    blocked = {
        j
        for leg in tag["deliveries"]
        for j in [*leg["tag_path"], leg["vacant_reference"]]
    }
    readers = []
    paths = dimer.local_paths()
    for leg in tag["deliveries"]:
        owner = leg["receiver"]
        off = 12 * owner
        lab = local_layout(owner, blocked, paths)
        local_graph = graph.subgraph(range(off, off + 12))
        lab["contact"] = transposition(
            local_graph, board, lab["probe"][:2], lab["bath"][:2]
        )
        readers.append(lab)
    ring = [
        transposition(graph, board, readers[0]["probe"][:2], readers[i]["probe"][:2])
        for i in (1, 2)
    ]
    # Move the final tag out of the detector helpers before terminal assays.
    # A fourth neighbour houses the autonomous marker read after each stroke.
    last = tag["deliveries"][-1]
    used = {tag["source"], *(d["receiver"] for d in tag["deliveries"])}
    link = next(
        (a, b)
        for edge in board.external
        for a, b in (edge, edge[::-1])
        if a // 12 == last["receiver"] and b // 12 not in used
    )
    dock = link[1]
    arming = transposition(graph, board, [last["arrival"]], [dock])[0]
    clock_site = dock // 12
    reference = next(j for j in board.neighbours[dock] if j // 12 == clock_site)
    off = 12 * clock_site
    clock = next(
        tuple(off + i for i in p[:3])
        for p in paths
        if not ({off + i for i in p[:3]} & {dock, reference})
    )
    assert all(
        not (set(clock) | {dock, reference})
        & (set(lab["probe"]) | set(lab["bath"]) | set(lab["selector"]))
        for lab in readers
    )
    return dict(
        tag=tag,
        readers=readers,
        ring=ring,
        clock=clock,
        arming_path=arming,
        parked_tag=dock,
        arming_reference=reference,
    )


def rotate(word, n=1):
    shift = n % len(word)
    return word[-shift:] + word[:-shift] if shift else word


def public_capacity(protocol):
    words = ["".join(w) for w in product("01", repeat=3)]
    family = protocol["publicness"]["family"]
    assert family == [[0, 1, 2]]

    def evaluate(groups, hidden=False):
        channels = []
        for group in groups:
            rows = [
                {w: {"".join(rotate(w, n)[i] for i in group): 1} for w in words}
                for n in range(3)
            ]
            if hidden:
                from fractions import Fraction

                rows = [
                    {
                        w: {
                            y: sum(row[w].get(y, 0) for row in rows) * Fraction(1, 3)
                            for y in {out for row in rows for out in row[w]}
                        }
                        for w in words
                    }
                ]
            channels.extend(dict(rows=row) for row in rows)
        graph = compound_confusability_graph(words, channels)
        code = maximum_independent_set(graph)
        return dict(M0=len(code), code=code)

    return dict(
        collective=evaluate(family),
        one=evaluate([[i] for i in range(3)]),
        any_two=evaluate(list(combinations(range(3), 2))),
        count_hidden=evaluate(family, True),
    )
