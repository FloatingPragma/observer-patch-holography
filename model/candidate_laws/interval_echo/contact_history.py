"""Read-from history of changed loads in the equipped patch transport experiment."""
from __future__ import annotations

from collections import Counter
from math import gamma

import numpy as np


class ContactHistory:
    """Two-party receipt on each load-changing contact, ordered by local versions."""

    def __init__(self, cells):
        self.last = [-1] * cells
        self.hand = np.zeros(cells, dtype=np.int64)
        self.events = []

    def record(self, x, y, changed, positive):
        self.hand += positive
        for a, b in zip(x[changed], y[changed], strict=True):
            event = len(self.events)
            self.events.append((a, b, self.last[a], self.last[b], self.hand[a], self.hand[b]))
            self.last[a] = self.last[b] = event

    def array(self):
        return np.asarray(self.events, dtype=np.int64).reshape(-1, 6)


def interval_summary(events, *, minimum=8, maximum=60):
    """Ordering fractions on comparable inclusive intervals in the stated size window.

    Every comparable pair is eligible. Relabeling events in another topological
    order therefore leaves the result unchanged. No dimensional exponent enters.
    """
    past = []
    heights = []
    future = [1 << event for event in range(len(events))]
    for event, (_, _, p, q, _, _) in enumerate(events):
        parents = (past[p] if p >= 0 else 0) | (past[q] if q >= 0 else 0)
        past.append(parents | 1 << event)
        heights.append(1 + max((heights[i] for i in (p, q) if i >= 0), default=0))
    for event in reversed(range(len(events))):
        for parent in events[event, 2:4]:
            if parent >= 0:
                future[parent] |= future[event]

    total = 0.0
    sizes = Counter()
    for end, ancestors in enumerate(past):
        starts = ancestors ^ (1 << end)
        while starts:
            first = starts & -starts
            starts ^= first
            interval = ancestors & future[first.bit_length() - 1]
            size = interval.bit_count()
            if minimum <= size <= maximum:
                comparable = 0
                members = interval
                while members:
                    member = members & -members
                    members ^= member
                    comparable += (past[member.bit_length() - 1] & interval).bit_count() - 1
                total += 2 * comparable / (size * (size - 1))
                sizes[size] += 1
    intervals = sum(sizes.values())
    # Match inclusive interval sizes. Two distinguished endpoints account for
    # 2n-3 comparable pairs. Interior pairs have the flat-diamond probability.
    controls = {}
    if intervals:
        for dimension in (1, 2, 3, 4):
            interior = gamma(dimension + 1) * gamma(dimension / 2) / (2 * gamma(1.5 * dimension))
            controls[dimension] = sum(
                count * (2 * (2 * n - 3) + (n - 2) * (n - 3) * interior) / (n * (n - 1))
                for n, count in sizes.items()
            ) / intervals
    return {"events": len(events), "source_height": max(heights, default=0),
            "intervals": intervals, "minimum_size": minimum, "maximum_size": maximum,
            "mean_ordering_fraction": total / intervals if intervals else None,
            "flat_diamond_expected_fraction": controls}
