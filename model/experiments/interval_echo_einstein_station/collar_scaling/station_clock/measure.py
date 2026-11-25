"""Count protected collar admissions using a completed station hop per clock cycle."""

import argparse
from collections import Counter
from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import shortest_path

from model.candidate_laws.interval_echo import candidate
from model.executable_model import cell, federation
from model.executable_model.carrier import icosahedron_from_coordinates
from model.experiments.mobile_memory.calibrated.geometry import _subdivide
from .. import clock_balls


HERE = Path(__file__).resolve().parent
CLOCK_TICKS = 7
MEMORY_PORT = 4


def cell_paths():
    """BFS paths from the station, never visiting the written word."""
    from collections import deque

    paths = {0: (0,)}
    queue = deque([0])
    local = cell()
    while queue:
        a = queue.popleft()
        for b in local.neighbours[a]:
            if b != MEMORY_PORT and b not in paths:
                paths[b] = (*paths[a], b)
                queue.append(b)
    return paths


def transfer_tables():
    law = candidate.native_law
    exchange, contrast = law.local_term_parts()
    h = exchange + candidate.KAPPA * contrast
    two = np.array([sum(s) == 2 for s in law.SRC])
    changed = h.copy()
    changed[: law.N_SRC, : law.N_SRC] += exchange[: law.N_SRC, : law.N_SRC] * np.outer(
        two, two
    )
    tables = []
    for generator in (h, changed):
        u = expm(-1j * np.pi * generator)
        columns = [law.SRC.index((q, 0)) for q in (1, 2)]
        arrivals = [law.SRC.index((0, q)) for q in (1, 2)]
        tables.append(
            [float(abs(u[a, c]) ** 2) for a, c in zip(arrivals, columns, strict=True)]
        )
    return tables


def scan(adj):
    return shortest_path(adj, directed=False, unweighted=True, indices=0).astype(int)


def graph(pairs, size):
    a, b = pairs.T
    return coo_matrix(
        (np.ones(2 * len(a), dtype=np.uint8), (np.r_[a, b], np.r_[b, a])),
        shape=(size, size),
    ).tocsr()


def hardware(level):
    board = federation(level)
    routes = cell_paths()
    service = Counter()
    edges = np.asarray(board.external)
    spare = next(
        (e for e in cell().seams if not {0, MEMORY_PORT}.intersection(e)), None
    )
    if spare is None:
        raise ValueError("no spare blank contact")
    for a, b in edges:
        for start, end in ((a, b), (b, a)):
            s, t = start // 12, end // 12
            leg = [s * 12 + i for i in routes[start % 12]]
            leg.extend(t * 12 + i for i in routes[end % 12][::-1])
            count = len(leg) - 1
            if count > CLOCK_TICKS or any(q % 12 == MEMORY_PORT for q in leg):
                raise ValueError("route conflicts with the clock or a protected word")
            if any(v not in board.neighbours[u] for u, v in zip(leg, leg[1:])):
                raise ValueError("an unregistered contact entered a station hop")
            # The busy appointments transport the tag. For each remaining
            # appointment address a blank seam disjoint from the tag and word.
            # It acts on |0,0>, an exact zero-energy eigenstate.
            service[CLOCK_TICKS - count] += 1
    d = scan(graph(edges // 12, board.cells))
    return dict(
        level=level,
        directed_routes=int(2 * len(edges)),
        spare_ticks=dict(sorted(service.items())),
        distance_histogram=np.bincount(d).tolist(),
    )


def dual_distances(faces):
    f = np.asarray(faces, dtype=np.int64)
    edges = np.sort(f[:, ((0, 1), (1, 2), (0, 2))], axis=-1).reshape(-1, 2)
    order = np.lexsort((edges[:, 1], edges[:, 0]))
    if not np.array_equal(edges[order[::2]], edges[order[1::2]]):
        raise ValueError("the triangle edges do not form a closed surface")
    return scan(graph(order.reshape(-1, 2) // 3, len(faces)))


def corner_sector(radius):
    """Number of triangles reachable within one sector from its corner triangle."""
    return (radius + 2) ** 2 // 4 if radius >= 0 else 0


def corner_ball(radius):
    return sum(corner_sector(radius - turn) for turn in (0, 1, 1, 2, 2))


def run():
    native = [hardware(level) for level in range(3)]
    vertices, seams = icosahedron_from_coordinates()
    seams = set(seams)
    faces = [f for f in combinations(range(12), 3) if set(combinations(f, 2)) <= seams]
    panels = []
    for level in range(8):
        m = 2**level
        d = dual_distances(faces)
        height = int(d.max())
        if level < 3 and np.bincount(d).tolist() != native[level]["distance_histogram"]:
            raise ValueError(
                "subdivision survey disagrees with the bundled native stations"
            )
        if not 6 * m - 3 <= height <= 6 * m:
            raise ValueError(
                "calibrating distance violates the geometric height estimate"
            )
        if level % 2:
            radius = 2 ** ((level + 1) // 2)
            n = int(np.sum((d > radius // 2) & (d <= radius)))
            if n != corner_ball(radius) - corner_ball(radius // 2):
                raise ValueError("annulus count disagrees with the five-sector proof")
            size = len(faces)
            ell = Fraction(radius, height)
            single = Fraction(n, size) / ell**2
            pair = Fraction(n * (n - 1), size * (size - 1)) / ell**4
            panels.append(
                dict(
                    level=level,
                    cells=size,
                    radius=radius,
                    eccentricity=height,
                    annulus_stations=n,
                    normalized_single=str(single),
                    normalized_pair=str(pair),
                )
            )
        if level < 7:
            vertices, faces = _subdivide(vertices, faces)
    charge, lower = clock_balls.source_contrast()
    pair_limit = Fraction(729, 256)
    normal, changed = transfer_tables()
    return dict(
        hop_appointments=CLOCK_TICKS,
        native_arrivals=normal,
        changed_exchange_arrivals=changed,
        hardware=native,
        midpoint_surveys=panels,
        single_limit="27/16",
        pair_limit=str(pair_limit),
        conditional_collar_contrast_nats=charge,
        collar_contrast_per_fourth_power_limit_nats=charge * float(pair_limit),
        strict_rational_lower_nats=str(lower * pair_limit),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    out = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (HERE / "results.json").write_text(out)
    else:
        print(out, end="")
