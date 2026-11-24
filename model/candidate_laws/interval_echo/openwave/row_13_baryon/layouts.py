"""Exhaust the three disjoint three-port fundamental frames on a cell."""

from collections import Counter
from itertools import combinations, permutations

import numpy as np

from model.executable_model import carrier


PERMS = tuple(permutations(range(3)))
SIGNS = np.array([1, -1, -1, 1, 1, -1])


def singlet_support(layout):
    return [tuple(sorted(layout[i][p[i]] for i in range(3))) for p in PERMS]


def enumerate_layouts(cell):
    face = np.zeros((12, 12, 12), dtype=bool)
    for sites in cell.faces():
        for abc in permutations(sites):
            face[abc] = True
    symmetries = carrier.automorphisms(cell.seams)
    remaining = set(combinations(range(12), 3))
    counts, best = Counter(), None
    largest = -1
    while remaining:
        a = min(remaining)
        orbit = {tuple(sorted(g[i] for i in a)) for g in symmetries}
        remaining -= orbit
        rest = sorted(set(range(12)) - set(a))
        for group_b in combinations(rest, 3):
            c = np.array(list(permutations(sorted(set(rest) - set(group_b)), 3)))
            for b in permutations(group_b):
                n = np.zeros(len(c), int)
                for i, j, k in PERMS:
                    n += face[a[i], b[j], c[:, k]]
                for value, count in enumerate(np.bincount(n, minlength=7)):
                    counts[value] += len(orbit) * int(count)
                maximum = int(n.max())
                if maximum > largest:
                    largest = maximum
                    best = (a, b, tuple(map(int, c[np.argmax(n)])))
                elif maximum == largest:
                    best = min(best, (a, b, tuple(map(int, c[np.argmax(n)]))))
    assert all(v % 6 == 0 for v in counts.values())
    return {str(i): count // 6 for i, count in sorted(counts.items())}, best
