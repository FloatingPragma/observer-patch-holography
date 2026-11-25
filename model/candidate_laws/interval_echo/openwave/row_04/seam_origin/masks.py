"""Classify all local memories by their integer action on the filed rays."""

from collections import Counter, defaultdict
from fractions import Fraction

import numpy as np

from model.executable_model import federation
from .. import open_pair, topology_census
from ..open_pair_release import source_bound


def contacts(board, path, receiver):
    a, b, c = path
    edges = frozenset(tuple(sorted((p, q))) for p in path for q in board.neighbours[p])
    live = frozenset(tuple(sorted((p, q))) for p in (a, c) for q in (b, receiver))
    if not live <= edges or len(live) != 4 or len(edges) != 15:
        raise ValueError("Expected fifteen incident contacts including four conducting routes")
    return edges - live, live


def eigenmasks(board, edges, live, words, sign):
    """For each mask test 2K and C separately, up to their own scalar multiples."""
    source = (1, sign)
    # Assemble columns for single contacts. A spectator has no exchange or
    # contrast on these words. Keep all destinations, including outside home.
    targets = set(words)
    swaps = []
    for edge in edges:
        ends = set(edge)
        moves = []
        for word, weight in zip(words, source, strict=True):
            if len(ends.intersection(word)) == 1:
                target = tuple(sorted(ends.symmetric_difference(word)))
                targets.add(target)
                moves.append((word, target, weight))
        swaps.append(moves)
    ordered = sorted(targets)
    idx = {word: i for i, word in enumerate(ordered)}
    terms = np.zeros((len(edges), len(idx)), dtype=np.int16)
    cost = np.zeros((len(edges), 2), dtype=np.int16)
    for k, moves in enumerate(swaps):
        for word, target, weight in moves:
            terms[k, idx[word]] += weight
            terms[k, idx[target]] -= weight
            cost[k, words.index(word)] += 1
    kinetic = live @ terms
    potential = live @ cost + [sum(6 - board.degree[p] for p in w) for w in words]
    ray = np.zeros(len(idx), dtype=np.int16)
    ray[[idx[w] for w in words]] = source
    return ((kinetic == kinetic[:, idx[words[0]], None] * ray).all(axis=1)
            & (potential[:, 0] == potential[:, 1]))


def local_classes():
    board = federation(0)
    a, b, c = open_pair.PATH
    closed, openings = contacts(board, (a, b, c), open_pair.BRIDGE)
    edges = sorted(closed | openings)
    if closed != set(open_pair.geometry()[0]):
        raise ValueError("Comparison to the filed material failed")
    masks = np.arange(1 << len(edges))[:, None]
    live = (1 - ((masks >> np.arange(len(edges))) & 1)).astype(np.int16)
    writes = len(edges) - live.sum(axis=1)
    transmitting = live[:, [edges.index(e) for e in openings]].all(axis=1)
    classes = {}
    for name, words in (("one_load", ((a,), (c,))), ("pair", ((b, a), (b, c)))):
        words = tuple(tuple(sorted(w)) for w in words)
        for phase, sign in (("minus", -1), ("plus", 1)):
            yes = eigenmasks(board, edges, live, words, sign)
            hist = np.bincount(writes[yes], minlength=len(edges) + 1)
            # Each exact local pattern with k writes has probability 3^(15-k)/4^15.
            mass = sum(int(n) * Fraction(3**(len(edges)-k), 4**len(edges))
                       for k, n in enumerate(hist))
            hits = np.flatnonzero(yes & transmitting)
            classes[f"{name}_{phase}"] = {
                "by_write_count_0_to_15": hist.tolist(),
                "probability": str(mass),
                "four_routes_open_masks": hits.tolist(),
            }
    filed_mask = sum(1 << i for i, edge in enumerate(edges) if edge in closed)
    if classes["pair_minus"]["four_routes_open_masks"] != [filed_mask]:
        raise ArithmeticError("The four-open pair ray does not select precisely the filed mask")
    event = Fraction(3**4, 4**15)
    floor, exhausted = source_bound()
    return {
        "edges_low_bit_first": [list(e) for e in edges],
        "filed_mask": filed_mask, "classes": classes,
        "filed_local_probability": str(event),
        "entire_filed_board_probability": str(event * Fraction(3, 4)**(len(board.seams)-15)),
        "local_and_source_exhausted_probability": str(event * exhausted),
        "local_and_home_lower_bound": str(event * floor),
    }


def occurrence_bounds(events):
    """First two inclusion-exclusion sums without independence of nearby layouts."""
    if len(events) != len(set(events)):
        raise ValueError("Duplicate material event in the layout census")
    p = Fraction(3**4, 4**15)
    n = len(events)
    # Pairs sharing no seam are independent. Count the others only once.
    incident = defaultdict(list)
    signatures = Counter()
    overlaps = 0
    for i, (closed, opened) in enumerate(events):
        neighbours = {j for edge in closed | opened for j in incident[edge]}
        overlaps += len(neighbours)
        for j in neighbours:
            other_closed, other_opened = events[j]
            if closed.isdisjoint(other_opened) and opened.isdisjoint(other_closed):
                signatures[len(closed | other_closed), len(opened | other_opened)] += 1
        for edge in closed | opened:
            incident[edge].append(i)
    disjoint = n * (n - 1) // 2 - overlaps
    second = disjoint * p**2 + sum(k * Fraction(3**b, 4**(a+b))
                                    for (a, b), k in signatures.items())
    first = n * p
    return {"eligible_layouts": n, "expected_matches": str(first),
            "two_event_sum": str(second), "probability_at_least_one_lower": str(max(0, first-second)),
            "probability_at_least_one_upper": str(min(1, first))}


def census():
    topology_census.native_blocks()
    boards = []
    for level in range(3):
        board = federation(level)
        events = [contacts(board, path, receiver) for home in range(board.cells)
                  for path, receiver in topology_census.layouts(board, home)]
        boards.append({"cells": board.cells, **occurrence_bounds(events)})
    return {"filed_address": local_classes(), "independent_coin_boards": boards}
