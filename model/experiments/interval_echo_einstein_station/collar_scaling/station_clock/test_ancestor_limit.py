"""Check ordered pair admission from station addresses and separate sector sums."""

from collections import Counter
from fractions import Fraction
import json
from itertools import permutations

from model.executable_model import federation
from . import ancestor_limit as source


def test_actual_external_routes_and_counted_ancestor_pairs():
    stored = json.loads((source.HERE / "ancestor_limit.json").read_text())
    for level, record in stored["bundled_native_surveys"].items():
        level = int(level)
        board = federation(level)
        m = 2**level
        adj = [set() for _ in range(board.cells)]
        for a, b in board.external:
            adj[a // 12].add(b // 12)
            adj[b // 12].add(a // 12)
        frontier, seen, admitted = {0}, {0}, set()
        for radius in range(1, record["radius"] + 1):
            frontier = set.union(*(adj[i] for i in frontier)) - seen
            seen |= frontier
            if 2 * radius > record["radius"]:
                admitted.update(frontier)
        counts = Counter(x // (m * m) for x in admitted)
        assert [counts[i] for i in range(20)] == record["counts_by_native_ancestor"]
        assert sum(count * (count - 1) for count in counts.values()) == sum(
            i // (m * m) == j // (m * m) for i, j in permutations(admitted, 2)
        )


def test_five_storeys_and_normalized_limit():
    stored = json.loads((source.HERE / "ancestor_limit.json").read_text())
    for row in stored["midpoint_surveys"]:
        r, m = row["radius"], row["frequency"]
        h = row["eccentricity"]
        # One sector has floor(t/2)+1 triangles at depth t. Count its actual
        # nonnegative depths for each of the five possible sector offsets.
        counts = [
            sum(t // 2 + 1 for t in range(r + 1) if r // 2 < t + offset <= r)
            for offset in (0, 1, 1, 2, 2)
        ]
        assert sorted(x for x in counts if x) == sorted(
            x for x in row["ancestor_populations"] if x
        )
        pairs = sum(n * (n - 1) for n in counts)
        assert pairs == row["ordered_pair_count"]
        total = 20 * m * m
        assert Fraction(pairs, total * (total - 1)) / Fraction(r, h) ** 4 == Fraction(
            row["normalized_pair"]
        )
        if r % 4 == 0:
            assert pairs == source.polynomial(r // 4)
    # Leading term in k=R/4, host cells=5R^4, H/m -> 6, ℓ=2/(bR).
    assert (
        Fraction(45 * 6**4, 256 * 400)
        == Fraction(stored["single_region_pair_limit"])
        == Fraction(729, 1280)
    )
