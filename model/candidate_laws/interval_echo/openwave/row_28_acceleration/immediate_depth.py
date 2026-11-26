"""Immediate strict-repair depth on every cell word and exterior-write pattern."""
from itertools import combinations

import numpy as np

from model.experiments.interval_echo_consensus.asynchronous import CELL_SEAMS


def immediate_depth_counts():
    words = np.arange(1 << 12)
    ports = (1 << 12) - 1
    edges = [((1 << a) | (1 << b), ((words >> a) ^ (words >> b)) & 1 != 0)
             for a, b in CELL_SEAMS]
    histogram = np.zeros(13, dtype=np.int64)
    for size in range(4):
        for external in combinations(range(12), size):
            exterior = sum(1 << port for port in external)
            remaining = np.full(len(words), ports)
            depth = np.full(len(words), int(bool(exterior)), dtype=int)
            pending = exterior
            epoch = 0
            while True:
                epoch += 1
                blocked = np.full(len(words), pending)
                remove = blocked.copy()
                for endpoints, unequal in edges:
                    ready = unequal & ((remaining & endpoints) == endpoints)
                    enabled = ready & ((blocked & endpoints) == 0)
                    remove[enabled] |= endpoints
                    # Lower-ranked ready neighbours block this round even when
                    # they are themselves blocked at their other endpoint.
                    blocked[ready] |= endpoints
                if not remove.any():
                    break
                depth[remove != 0] = epoch
                remaining &= ~remove
                pending = 0
            histogram += np.bincount(depth, minlength=len(histogram))
    return np.trim_zeros(histogram, trim="b").tolist()


if __name__ == "__main__":
    print(immediate_depth_counts())
