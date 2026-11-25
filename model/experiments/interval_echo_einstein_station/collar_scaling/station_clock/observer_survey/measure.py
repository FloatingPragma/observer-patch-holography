"""Station-centered annuli and a test of what the gateway incidence determines."""

import argparse
from collections import Counter, deque
from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.sparse.csgraph import shortest_path

from model.executable_model import cell, federation
from ...measure import ProtectedRanges
from .. import measure, word_entropy

HERE = Path(__file__).parent


def flat_ball(radius):
    return 1 + 3 * radius * (radius + 1) // 2


def atlas(level):
    vertices, edges = measure.icosahedron_from_coordinates()
    edges = set(edges)
    faces = [f for f in combinations(range(12), 3) if set(combinations(f, 2)) <= edges]
    for _ in range(level):
        vertices, faces = measure._subdivide(vertices, faces)
    owners = {}
    adjacency = [[] for _ in faces]
    for i, face in enumerate(faces):
        for side in combinations(sorted(face), 2):
            if side in owners:
                j = owners.pop(side)
                adjacency[i].append(j)
                adjacency[j].append(i)
            else:
                owners[side] = i
    if owners or any(len(row) != 3 for row in adjacency):
        raise ValueError("incomplete dual adjacency")
    return adjacency


def annulus(adjacency, center, radius):
    visited = {center}
    layer = {center}
    enrolled = []
    for step in range(1, radius + 1):
        layer = {v for u in layer for v in adjacency[u]} - visited
        visited.update(layer)
        if 2 * step > radius:
            enrolled.extend(layer)
    return enrolled


def all_centers(adjacency, level, radius, centers):
    per_face = 4**level
    pairs = same_face = 0
    sizes = Counter()
    for source in centers:
        region = annulus(adjacency, source, radius)
        n = len(region)
        sizes[n] += 1
        pairs += n * (n - 1)
        same_face += sum(
            c * (c - 1) for c in Counter(v // per_face for v in region).values()
        )
    return dict(
        level=level,
        radius=radius,
        centers=len(centers),
        annulus_size_histogram={str(k): v for k, v in sorted(sizes.items())},
        ordered_pairs=pairs,
        ordered_pairs_in_one_original_face=same_face,
        flat_annulus_count=flat_ball(radius) - flat_ball(radius // 2),
    )


def gate_triples():
    contacts = cell().adjacency().toarray()
    contacts[4, :] = contacts[:, 4] = 0
    distance = shortest_path(contacts, directed=False, unweighted=True)
    choices = {}
    for ports in combinations([v for v in range(12) if v != 4], 3):
        costs = {int(distance[a, b]) for a, b in combinations(ports, 2)}
        for length in (1, 2):
            if costs == {length}:
                choices.setdefault(length, ports)
    return choices


def patched_ranges(adjacency, dual, patch):
    triples = gate_triples()
    gateways = [
        dict(zip(sorted(row), triples[1 if h <= patch else 2], strict=True))
        for row, h in zip(adjacency, dual, strict=True)
    ]
    edges = [
        (12 * u + gateways[u][v], 12 * v + gateways[v][u])
        for u, row in enumerate(adjacency)
        for v in row
        if u < v
    ]
    board = SimpleNamespace(external=np.asarray(edges), cells=len(adjacency))
    return ProtectedRanges(board).survey(0)


def label_panel(level):
    adjacency = atlas(level)
    dual = np.full(len(adjacency), -1, dtype=int)
    dual[0] = 0
    q = deque([0])
    while q:
        u = q.popleft()
        for v in adjacency[u]:
            if dual[v] < 0:
                dual[v] = dual[u] + 1
                q.append(v)
    patch = 2**level
    slow = patched_ranges(adjacency, dual, -1)
    mixed = patched_ranges(adjacency, dual, patch)
    if np.any(slow < 3 * dual - 2) or np.any(slow > 3 * dual + 4):
        raise ValueError(
            "uniform-transit metric violates the predicted distance envelope"
        )
    if np.any(mixed < 3 * (dual - patch) - 3) or np.any(mixed > 3 * dual + 4):
        raise ValueError("fast-patch metric violates the global distance envelope")
    interior = dual <= patch
    if np.any(mixed[interior] < 2 * dual[interior] - 1) or np.any(
        mixed[interior] > 2 * dual[interior] + 5
    ):
        raise ValueError("fast-patch metric violates its inner distance envelope")
    return dict(
        level=level,
        cells=len(dual),
        patch_radius=patch,
        dual_eccentricity=int(dual.max()),
        slow_eccentricity=int(slow.max()),
        mixed_eccentricity=int(mixed.max()),
        changed_station_ranges=int(np.sum(slow != mixed)),
    )


def run():
    samples = []
    for level, radius in ((2, 2), (3, 4), (4, 4)):
        adjacency = atlas(level)
        # Test every center in the small board, then all positions in one
        # original face, keeping the whole sphere available to the paths.
        centers = range(len(adjacency) if level == 2 else 4**level)
        samples.append(all_centers(adjacency, level, radius, centers))
    native_sizes = []
    for level in range(3):
        board = federation(level)
        adjacency = [[] for _ in range(board.cells)]
        for s, t in board.external:
            a, b = s // 12, t // 12
            adjacency[a].append(b)
            adjacency[b].append(a)
        sizes = Counter(len(annulus(adjacency, i, 2)) for i in range(board.cells))
        native_sizes.append({str(k): v for k, v in sorted(sizes.items())})
    bulk = Fraction(81, 40)
    _, kappa = word_entropy.entropy_and_slope(word_entropy.probabilities(), -1)
    slope = word_entropy.endpoints(kappa * bulk.numerator**2 / bulk.denominator**2)
    return dict(
        native_radius_two_annuli=native_sizes,
        atlas_surveys=samples,
        regular_single_limit=str(bulk),
        regular_pair_limit=str(bulk**2),
        regular_closed_ball_single_limit=str(Fraction(27, 10)),
        regular_closed_ball_pair_limit=str(Fraction(27, 10) ** 2),
        mean_regional_entropy_slope_nats=slope,
        comparison_gateway_triples={str(k): list(v) for k, v in gate_triples().items()},
        comparison_labels=[label_panel(level) for level in (2, 3)],
        original_clock_comparison_single_limits=["27/16", "243/64"],
        original_clock_comparison_closed_ball_single_limits=["9/4", "81/16"],
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (HERE / "results.json").write_text(data)
    else:
        print(data, end="")
