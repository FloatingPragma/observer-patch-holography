"""Independent checks from the saved geometry, words and seam counts."""

import json
from collections import Counter
from itertools import permutations, product
from math import comb
from pathlib import Path

import numpy as np

DATA = json.loads(Path(__file__).with_name("results.json").read_text())


def test_contact_energy_and_raw_records():
    edges = DATA["edges"]
    shots = DATA["shots_per_seam"]
    assert len(edges) == 62
    for trial in DATA["sectors"]:
        for key in ("separated", "three_walls"):
            data = trial[key]
            word = data["word"]
            occ = [(word >> i) & 1 for i in range(30)]
            both = sum(occ[i] * occ[j] for i, j in edges)
            degree_cut = sum(occ[i] != occ[j] for i, j in edges)
            counts = np.array(data["counts"])
            np.testing.assert_array_equal(
                counts.sum(axis=1), np.full(len(edges), shots)
            )
            np.testing.assert_array_equal(
                counts[:, 0], [shots * occ[i] * occ[j] for i, j in edges]
            )
            for k, record in data["energies"].items():
                assert record["mean"] == degree_cut / 2 - 2 * float(k) * both
                observed = (counts[:, 1] - 2 * float(k) * counts[:, 0]).sum() / shots
                assert observed == record["observed"]
        assert trial["sector"]["count"] == comb(30, sum(trial["loads"]))


def test_wall_block_and_filter_bound_from_all_translations():
    edges = DATA["edges"]
    source = {3: (0, 4, 7), 4: (0, 1, 4, 7)}
    translates = {
        n: [sum(1 << ((i + t) % 10) for i in pattern) for t in range(10)]
        for n, pattern in source.items()
    }
    for sample in DATA["sectors"]:
        walls = {
            x | (y << 10) | (z << 20)
            for labels in set(permutations(sample["loads"]))
            for x, y, z in product(*(translates[n] for n in labels))
        }
        assert len(walls) == sample["walls"]["count"] == 3000
        incoming = Counter()
        max_exits = 0
        max_agreements = 0
        for word in walls:
            occupied = [(word >> i) & 1 for i in range(30)]
            max_agreements = max(
                max_agreements, sum(occupied[a] * occupied[b] for a, b in edges)
            )
            exits = set()
            for a, b in edges:
                if occupied[a] != occupied[b]:
                    changed = word ^ ((1 << a) | (1 << b))
                    if changed not in walls:
                        exits.add(changed)
            max_exits = max(max_exits, len(exits))
            incoming.update(exits)
        gap = -16 * max_agreements - sample["separated"]["energies"]["8.0"]["mean"]
        bound = max_exits * max(incoming.values()) / (4 * gap**2)
        assert gap > 0
        assert np.isclose(
            sample["energy_filter"]["any_energy_filtered_wall_probability_upper"], bound
        )
        assert bound < 0.002
