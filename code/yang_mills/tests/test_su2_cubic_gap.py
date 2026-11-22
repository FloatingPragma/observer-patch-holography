from collections import deque
from fractions import Fraction
from itertools import product
import json
from pathlib import Path
import sys

import pytest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
from su2_cubic_gap import cubic_bound, run


def girth_through_origin(side):
    # By translation and cubic symmetry every edge has the same shortest return path.
    start, end = (0, 0, 0), (1, 0, 0)
    queue = deque([(start, 0)])
    seen = {start}
    while queue:
        vertex, distance = queue.popleft()
        for axis, step in product(range(3), (-1, 1)):
            neighbor = list(vertex)
            neighbor[axis] = (neighbor[axis] + step) % side
            neighbor = tuple(neighbor)
            if vertex == start and neighbor == end:
                continue
            if neighbor == end:
                return distance + 2
            if neighbor not in seen:
                seen.add(neighbor)
                queue.append((neighbor, distance + 1))
    raise AssertionError("No fundamental loop")


def test_physical_electric_sector_and_magnetic_bound():
    for side in (3, 4, 5, 6, 8, 9):
        gap = Fraction(3, 4) * girth_through_origin(side)
        row = cubic_bound(side, Fraction(1, 512))
        assert Fraction(row["electric_gap"]) == gap
        assert Fraction(row["spectral_gap_lower"]) == max(
            0, gap - Fraction(3 * side**3, 512)
        )
        assert row["positive_bound"] == (side < 8)
    assert run() == json.loads((HERE / "receipts/su2_cubic_gap.json").read_text())


def test_threshold_exactly_and_unsupported_graphs():
    assert not cubic_bound(4, Fraction(1, 64))["positive_bound"]
    assert cubic_bound(4, Fraction(1, 65))["positive_bound"]
    assert Fraction(cubic_bound(6, Fraction(0))["spectral_gap_lower"]) == 3
    for side, magnetic in ((2, Fraction(0)), (4, -1), (4, 0.5)):
        with pytest.raises(ValueError):
            cubic_bound(side, magnetic)
