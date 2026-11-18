"""Exact labelled one-cell scalar carrier from the pinned OpenWave export.

These labels agree with the delivered CELL.npz. Only its scalar incidence
law is reconstructed here; its separate six-component half-links are not
identified with this scalar channel. Matrices are external simulator data.
"""

from itertools import combinations
import math

import numpy as np

UPSTREAM_SOURCE = {
    "commit": "22c17e489bf1d021ad99e45f4ae7ddd84d5e6cd5",
    "path": "candidates/np-solver/export/README.md",
    "sha256": "33c4413e8e72df5004f779485ceafef28a608e6f586900608db1627986816ed5",
    "bytes": 11257,
}
DELIVERED_CARRIER = {
    "path": "candidates/np-solver/pair_transit/delivered/source/CELL.npz",
    "sha256": "9fdc4f4475baa029c1b95f76f1b0a546617ecaef582fbf9f476952c4b258c42f",
    "bytes": 72883,
}
EDGES = (
    (0, 1), (0, 2), (0, 3), (0, 4), (0, 6), (1, 2), (1, 3), (1, 5),
    (1, 7), (2, 4), (2, 5), (2, 8), (3, 6), (3, 7), (3, 9), (4, 6),
    (4, 8), (4, 10), (5, 7), (5, 8), (5, 11), (6, 9), (6, 10),
    (7, 9), (7, 11), (8, 10), (8, 11), (9, 10), (9, 11), (10, 11),
)


def require(condition, message):
    """A production check, including in optimized Python executions."""
    if not condition:
        raise ValueError(message)


def load_source(which="cell"):
    """Generate an independent copy of the declared 30 by 120 incidence."""
    require(type(which) is str and which == "cell", "source domain must be cell")
    edges = np.array(EDGES, dtype=int)
    neighbors = [sorted({x for edge in EDGES if b in edge for x in edge if x != b})
                 for b in range(12)]
    wedges = np.array([(a, b, c) for b in range(12)
                       for a, c in combinations(neighbors[b], 2)], dtype=int)
    index = {edge: k for k, edge in enumerate(EDGES)}
    incidence = np.zeros((30, 120))
    for w, (a, b, c) in enumerate(wedges):
        incidence[index[tuple(sorted((a, b)))], w] = 1 / math.sqrt(2)
        incidence[index[tuple(sorted((b, c)))], w] = -1 / math.sqrt(2)
    law = np.block([[np.eye(30), incidence], [incidence.T, -np.eye(120)]])
    return {"N": 12, "E": 30, "W": 120, "edges": edges, "wedges": wedges,
            "index": index, "degree": np.full(12, 5, dtype=int), "B": incidence,
            "D": law, "H": law.copy(), "m": 1.0}


def chi(data, center=0):
    """The normalized local three-wedge circulation at level -1."""
    require(type(center) is int and 0 <= center < 12, "center must be an integer in [0,12)")
    adjacent = sorted({int(x) for edge in data["edges"] if center in edge
                       for x in edge if x != center})
    a, c, f = adjacent[:3]
    index = {tuple(map(int, row)): k for k, row in enumerate(data["wedges"])}
    result = np.zeros(data["E"] + data["W"])
    for wedge, sign in (((a, center, c), 1), ((c, center, f), 1), ((a, center, f), -1)):
        result[data["E"] + index[wedge]] = sign / math.sqrt(3)
    require(np.linalg.norm(data["D"] @ result + result) < 1e-13, "local circulation is not an eigenmode")
    return result
