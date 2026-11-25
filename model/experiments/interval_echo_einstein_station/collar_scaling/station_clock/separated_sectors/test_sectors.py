"""Reconstruct the face distances and entropy and orientation readouts independently."""

from collections import defaultdict, deque
from fractions import Fraction
from itertools import combinations
import json
import math

import numpy as np

from model.executable_model import federation
from .... import shared_variance, volume
from .. import ancestor_limit, measure as clock
from . import measure


def distances(graph, starts):
    seen = dict.fromkeys(starts, 0)
    queue = deque(starts)
    while queue:
        here = queue.popleft()
        for there in graph[here]:
            if there not in seen:
                seen[there] = seen[here] + 1
                queue.append(there)
    return seen


def test_distance_height_and_sector_counts():
    vertices, seams = clock.icosahedron_from_coordinates()
    seams = set(seams)
    faces = [f for f in combinations(range(12), 3) if set(combinations(f, 2)) <= seams]
    turns = ancestor_limit.offsets(faces)
    second = min(i for i, offset in turns.items() if offset == 2)
    corner = faces[0][0]
    height = {v: int(v in faces[0] and v != corner) for v in range(12)}
    for level in range(1, 5):
        vertices, children = clock._subdivide(vertices, faces)
        next_height = {v: 2 * value for v, value in height.items()}
        for (a, b, c), (ab, bc, ac) in zip(faces, children[3::4], strict=True):
            for v, pair in [(ab, (a, b)), (bc, (b, c)), (ac, (a, c))]:
                value = sum(height[i] for i in pair)
                assert next_height.setdefault(v, value) == value
        height, faces = next_height, children
        pot = [max(height[v] for v in f) for f in faces]
        owners = defaultdict(list)
        for i, face in enumerate(faces):
            for edge in combinations(sorted(face), 2):
                owners[edge].append(i)
        graph = defaultdict(set)
        for a, b in owners.values():
            graph[a].add(b)
            graph[b].add(a)
            assert abs(pot[a] - pot[b]) <= 1

        block = 4**level
        first_cells = range(block)
        second_cells = range(second * block, (second + 1) * block)
        assert all(pot[y] == 0 for y in second_cells)
        home = distances(graph, [0])
        gap = distances(graph, second_cells)
        for radius in range(2, 2 * 2**level - 1, 2):
            for inner in (-1, radius // 2):
                first = [i for i in first_cells if inner < home[i] <= radius]
                other = [j for j in second_cells if inner < home[j] <= radius]
                assert len(first) == clock.corner_sector(radius) - clock.corner_sector(
                    inner
                )
                assert len(other) == clock.corner_sector(
                    radius - 2
                ) - clock.corner_sector(inner - 2)
                assert all(home[i] <= 2 * pot[i] and pot[i] <= gap[i] for i in first)
                if inner >= 0:
                    assert min(gap[i] for i in first) > radius / 4


def test_native_modular_and_energy_read():
    out = measure.run()
    assert out == json.loads((measure.HERE / "results.json").read_text())
    words = out["words"]
    assert Fraction(words["conditional_entropy_log_base"]) == Fraction(169, 105)
    assert Fraction(words["conditional_entropy_log_multiplier"]) == Fraction(1, 144)
    assert Fraction(words["conditional_orientation_slope"]) == Fraction(1, 36)
    _, _, _, source = volume.radar_source_and_hand(federation(0))
    (reverse, product, shared), _ = shared_variance.source_words(source)
    delta = product - reverse
    np.testing.assert_allclose(delta.sum(axis=0), 0, rtol=0, atol=1e-16)
    np.testing.assert_allclose(delta.sum(axis=1), 0, rtol=0, atol=1e-16)
    kappa = -np.sum(delta * np.log(reverse))
    np.testing.assert_allclose(kappa, math.log(169 / 105) / 144, rtol=0, atol=1e-15)
    # Read cyclic orientations directly from the three allowed energy pairs.
    score = np.array(
        [
            0 if a == b else 1 if (a, b) in {(0, 32), (32, 33), (33, 0)} else -1
            for a in source["source_word_energies"]
            for b in source["source_word_energies"]
        ]
    )
    observed = np.sum((shared - product) * np.outer(score, score))
    np.testing.assert_allclose(observed, 1 / 36, rtol=0, atol=1e-16)
    assert (
        abs(words["held_shared_minus_fresh"]["mean_difference"] - observed)
        < words["held_shared_minus_fresh"]["simultaneous_99pct_radius"]
    )
    assert out["limits"]["ball"]["eligible_per_fourth"] == "81/200"
    assert out["limits"]["annulus"]["eligible_per_fourth"] == "729/3200"


def test_complete_regional_record_by_placement():
    # Level 1 has four descendants per original face. Use the declared
    # midpoint addresses for the completed station clock's two-hop ball.
    vertices, seams = clock.icosahedron_from_coordinates()
    seams = set(seams)
    faces = [f for f in combinations(range(12), 3) if set(combinations(f, 2)) <= seams]
    _, children = clock._subdivide(vertices, faces)
    home = clock.dual_distances(children)
    region = set(np.flatnonzero(home <= 2))
    _, _, _, src = volume.radar_source_and_hand(federation(0))
    raw, _ = shared_variance.source_words(src)
    reverse, fresh, shared = np.rint(raw * 2592).astype(int)
    entries = {}
    for x in range(80):
        for y in range(80):
            if x == y:
                continue
            face_x, face_y = x // 4, y // 4
            cross = {face_x, face_y} == {0, 3}
            q = shared if face_x == face_y else fresh
            p = reverse if cross else q
            addr = (x if x in region else -1, y if y in region else -1)
            remove = tuple(i for i, site in enumerate(addr) if site < 0)
            weights = np.stack([p.sum(axis=remove), q.sum(axis=remove)])
            if addr in entries:
                entries[addr] += weights
            else:
                entries[addr] = weights
    size = 80 * 79 * 2592
    masses = np.concatenate([v.reshape(2, -1) for v in entries.values()], axis=1)
    np.testing.assert_array_equal(masses.sum(axis=1), [size, size])
    delta = (masses[1] - masses[0]) / size
    slope = -np.sum(delta * np.log(masses[0] / size))
    record = json.loads((measure.HERE / "results.json").read_text())
    pair_prob = Fraction(record["surveys"][0]["eligible"])
    np.testing.assert_allclose(
        slope, float(pair_prob) * math.log(169 / 105) / 144, rtol=1e-12
    )
    for addr, (p, q) in entries.items():
        if -1 in addr or {s // 4 for s in addr} != {0, 3}:
            np.testing.assert_array_equal(p, q)
