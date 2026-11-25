from collections import Counter, deque
from fractions import Fraction
from itertools import product
import json
from math import sqrt

import numpy as np

from . import measure as m


def test_access_and_independent_ranges():
    r = json.loads((m.HERE / "records.json").read_text())
    board, axes, flag = m.layout()
    assert len(set.union(*(set(a) for a in (*axes, flag)))) == 38
    for path in (*axes, flag):
        assert all(q in board.neighbours[p] for p, q in zip(path, path[1:]))
    for axis in r["axes"]:
        assert all(c == [32, 32] for c in axis["directed_arrivals"])
        assert all(c == [0, 0] for c in axis["parked_arrivals"])
    for corner, row in zip(
        product((5, 6), repeat=3), r["native_range_reads"], strict=True
    ):
        q = deque([corner])
        dist = {corner: 1}
        while q:
            pos = q.popleft()
            for delta in product((-1, 0, 1), repeat=3):
                other = tuple(u + v for u, v in zip(pos, delta, strict=True))
                if all(0 <= x < 12 for x in other) and other not in dist:
                    dist[other] = 1 + dist[pos]
                    q.append(other)
        assert [dist[p] for p in product(range(12), repeat=3)] == row


def test_shell_edges_and_distinct_source_sector():
    r = json.loads((m.HERE / "records.json").read_text())
    q = np.arange(12, dtype=int)
    offset = np.minimum(abs(q - 5), abs(q - 6))
    rad = 1 + np.maximum(
        np.maximum(offset[:, None, None], offset[None, :, None]), offset[None, None, :]
    )
    source = Counter()
    seen = set()
    for p, q in m.edges():
        assert (p, q) not in seen
        seen.add((p, q))
        if p == m.READY:
            source[q] += 1
        else:
            assert sum(abs(x - y) for x, y in zip(p, q, strict=True)) == 1
    assert all(n == 1 for n in source.values()) and len(source) == 8
    assert (5, 5, 5) in source and m.READY not in source
    for bank, cut in zip(r["surveyed"], (None, 5, 6), strict=True):
        assert bank["source_arrivals"] == 8
        for n in range(1, m.DEPTH + 1):
            inside = rad <= n
            doors = []
            for axis in range(3):
                crossing = np.diff(inside.astype(int), axis=axis) != 0
                if axis == 0 and cut is not None:
                    crossing[cut] = False
                doors.append(int(crossing.sum()))
            assert sum(doors) == bank["outward"][n]
            if cut is None:
                assert doors == [8 * n * n] * 3


def test_full_vertex_current_and_original_launch_uncertainty():
    r = json.loads((m.HERE / "records.json").read_text())
    s = r["steady"]
    panel = r["walk"]
    assert s == m.steady() and s["full_interior"] == 1000
    assert [Fraction(x) for x in s["outward"]] == [1] * 5
    centre, corner = [Fraction(s["per_edge"][k]) for k in ("1", "-1")]
    mean = Fraction(s["mean"])
    variance = Fraction(s["variance"])
    assert centre > 2 * corner > 0 and 24 * (centre - corner) == mean
    n = r["shots"]
    assert (
        sum(panel["released"])
        == sum(j for _, j in panel["completions"])
        == sum(j for _, j in panel["contrast"])
        == n
    )
    measured = sum(k * j for k, j in panel["contrast"]) / n
    radius = sqrt(float(variance) / (0.01 * n))
    assert 0 < measured - radius < float(mean) < measured + radius
