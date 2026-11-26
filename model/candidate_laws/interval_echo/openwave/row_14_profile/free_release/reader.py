"""Signed terminal reads made on internal seams of one owner after stopping flight."""

from collections import Counter
from itertools import combinations

import networkx as nx
import numpy as np

from model.candidate_laws.interval_echo import candidate
from model.executable_model import cell

TREE = tuple(nx.bfs_edges(nx.Graph(cell().seams), 0))


def calibrate():
    law = candidate.native_law
    echo = law.echo()
    blank = np.eye(law.DIM)[:, :9]
    result = {}
    for d in range(-2, 3):
        mask = np.repeat([False] + [x == d for _, x in law.LABELS], 9)
        instrument = echo @ (mask[:, None] * (echo @ blank))
        expected = blank * np.array([a - b == d for a, b in law.SRC])
        np.testing.assert_allclose(instrument, expected, rtol=0, atol=3e-14)
    for a, b in law.SRC:
        result[f"{a},{b}"] = a - b
    return result


def decode(clicks):
    if len(clicks) != len(TREE):
        raise ValueError("Incomplete comparison tree")
    relative = {0: 0}
    for (a, b), d in zip(TREE, clicks, strict=True):
        if d not in range(-2, 3):
            raise ValueError("Invalid comparison")
        relative[b] = relative[a] - d
    common, amount = Counter(relative.values()).most_common(1)[0]
    charges = [relative[p] - common for p in range(12)]
    if (
        amount < 8
        or min(charges) < -1
        or max(charges) > 1
        or charges.count(1) > 2
        or charges.count(-1) > 2
    ):
        raise ValueError("Read outside the admitted inventory")
    return charges


def census():
    lookup = calibrate()
    count = 0
    for nplus in range(3):
        for plus in combinations(range(12), nplus):
            for nminus in range(3):
                for minus in combinations(sorted(set(range(12)) - set(plus)), nminus):
                    charges = [int(i in plus) - int(i in minus) for i in range(12)]
                    clicks = [
                        lookup[f"{1 + charges[a]},{1 + charges[b]}"] for a, b in TREE
                    ]
                    if decode(clicks) != charges:
                        raise ValueError("Terminal reader failed")
                    count += 1
    return dict(
        internal_seams=TREE, signed_rail_calibration=lookup, decoded_local_words=count
    )
