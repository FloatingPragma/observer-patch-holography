"""Determine the common clock-scale interval of recorded seam events."""

from functools import cmp_to_key, lru_cache
from itertools import permutations
from math import lcm

from networkx.algorithms.isomorphism import GraphMatcher

from .chart import squared_separations


def sign(a):
    return (a > 0) - (a < 0)


def compare(x, y):
    """Order a+b*sqrt(5) without evaluating a square root."""
    a, b = x[0] - y[0], x[1] - y[1]
    if a * b >= 0:
        return sign(a) or sign(b)
    return sign(a) * sign(a * a - 5 * b * b)


def product(x, y):
    return x[0] * y[0] + 5 * x[1] * y[1], x[0] * y[1] + x[1] * y[0]


def independent(midpoints, ds):
    """Twice the Gram of three differences, in the common distance unit."""
    a, *others = midpoints
    h = [
        [tuple(ds[i][a][k] + ds[j][a][k] - ds[i][j][k] for k in (0, 1)) for i in others]
        for j in others
    ]
    det = [0, 0]
    for i, j, k in permutations(range(3)):
        parity = 1 if (j - i) % 3 == 1 else -1
        x, y = product(h[0][i], product(h[1][j], h[2][k]))
        det[0] += parity * x
        det[1] += parity * y
    return det != [0, 0]


def provenance(edges):
    versions = [None] * 12
    parents = []
    for event, (a, b) in enumerate(edges):
        parents.append(
            tuple(sorted({p for p in (versions[a], versions[b]) if p is not None}))
        )
        versions[a] = versions[b] = event
    return tuple(parents)


@lru_cache(None)
def obligations(parents):
    heights, ancestors, lower, upper, same = [], [], [], [], []
    for i, ps in enumerate(parents):
        past = 0
        for j in ps:
            past |= (1 << j) | ancestors[j]
        h = max((heights[j] + 1 for j in ps), default=0)
        heights.append(h)
        ancestors.append(past)
        lower.extend((j, i, h - heights[j]) for j in ps)
        for j in range(i):
            if not ((1 << j) & past):
                dh = abs(h - heights[j])
                if dh:
                    upper.append((j, i, dh))
                else:
                    same.append((j, i))
    return lower, upper, same


class ScaleTest:
    def __init__(self, graph, horizon=4):
        self.edges, self.squared = squared_separations(graph)
        # All finite bounds are (a+b*sqrt(5))/unit.
        base = lcm(*(h * h for h in range(1, horizon)))
        self.unit = 20 * base
        scaled = [
            [
                [(a * (base // (h * h)), b * (base // (h * h))) for a, b in row]
                for row in self.squared
            ]
            for h in range(1, horizon)
        ]
        self.values = sorted(
            {v for table in scaled for row in table for v in row},
            key=cmp_to_key(compare),
        )
        rank = {v: i for i, v in enumerate(self.values)}
        self.ranks = [
            None,
            *[[[rank[v] for v in row] for row in table] for table in scaled],
        ]
        self.infinity = len(self.values)
        assert self.values[0] == (0, 0)

    def interval(self, word, parents):
        lo, up, same = obligations(parents)
        if any(self.squared[word[j]][word[i]] == (0, 0) for j, i in same):
            raise ValueError("Equal-height incomparable events share a position")
        lower = max((self.ranks[h][word[i]][word[j]] for i, j, h in lo), default=0)
        upper = min(
            (self.ranks[h][word[i]][word[j]] for i, j, h in up), default=self.infinity
        )
        return lower, upper

    def bound(self, rank):
        return None if rank == self.infinity else list(self.values[rank])


def census(graph):
    test = ScaleTest(graph)
    first = test.edges[0]
    orbit = {
        tuple(sorted((p[first[0]], p[first[1]])))
        for p in GraphMatcher(graph, graph).isomorphisms_iter()
    }
    if orbit != set(test.edges):
        raise ValueError("First-contact symmetry failed")
    summaries, examples = [], {}
    for n in range(1, 5):
        counts = {"feasible": 0, "tie": 0, "reversed": 0}
        for suffix in permutations(range(1, 30), n - 1):
            w = (0, *suffix)
            ps = provenance([test.edges[i] for i in w])
            lo, hi = test.interval(w, ps)
            status = "feasible" if lo < hi else "tie" if lo == hi else "reversed"
            counts[status] += len(orbit)
            key = f"{n}-{status}"
            if key in examples or n < 3 or not (lo > 0 and hi < test.infinity):
                continue
            if n == 3 and ps.count(()) != 2:
                continue
            if n == 4 and (status == "tie" or not independent(w, test.squared)):
                continue
            examples[key] = {
                "word": list(w),
                "edges": [test.edges[i] for i in w],
                "lower": test.bound(lo),
                "upper": test.bound(hi),
            }
        summaries.append(counts)
    return {"lengths": summaries, "witnesses": examples, "scale_denominator": test.unit}
