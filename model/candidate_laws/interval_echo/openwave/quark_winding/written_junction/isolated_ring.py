"""Full-board invariant written rings with a prepared outward-pointing boundary."""

from functools import lru_cache
from itertools import product

import networkx as nx
import numpy as np

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.executable_model import carrier


@lru_cache
def geometry(level=0):
    board = carrier.federation(level)
    contacts = nx.Graph(board.seams)
    coarse = nx.Graph((a // 12, b // 12) for a, b in board.external)
    gateway = {
        (a // 12, b // 12): a % 12
        for edge in board.external
        for a, b in (edge, edge[::-1])
    }
    inside = nx.Graph(carrier.cell().seams)
    paths = {}
    tours = set()
    for start in coarse:
        for end in coarse[start]:
            for path in nx.all_simple_paths(coarse, start, end, cutoff=4):
                if len(path) == 5 and coarse.subgraph(path).number_of_edges() == 5:
                    tours.add(
                        min(
                            seq[i:] + seq[:i]
                            for v in (path, path[::-1])
                            for seq in [tuple(v)]
                            for i in range(5)
                        )
                    )
    for t in sorted(tours):
        choices = []
        for j, owner in enumerate(t):
            a = gateway[owner, t[j - 1]]
            b = gateway[owner, t[(j + 1) % 5]]
            by_size = {}
            for path in nx.all_simple_paths(inside, a, b, cutoff=6):
                if inside.subgraph(path).number_of_edges() == len(path) - 1:
                    v = tuple(12 * owner + i for i in path)
                    by_size[len(v)] = min(v, by_size.get(len(v), v))
            choices.append(by_size.values())
        for parts in sorted(product(*choices)):
            p = sum(parts, ())
            if len(p) in (12, 15, 18) and len(p) not in paths:
                assert contacts.subgraph(p).number_of_edges() == len(p)
                paths[len(p)] = p
        if len(paths) == 3:
            break
    return board, paths


def shield(board, path):
    region = set(path)
    boundary = [
        (a, b) if a in region else (b, a)
        for a, b in board.seams
        if (a in region) != (b in region)
    ]
    internal = {e for e in board.seams if set(e) <= region}
    # A head points to the only site from which an extra is allowed to hop.
    # Boundary heads point out, and exterior heads point to the larger address.
    fixed = {
        edge: (v if u in region else u if v in region else max(edge))
        for edge in board.seams
        if edge not in internal
        for u, v in [edge]
    }
    graph = [set(row) for row in board.neighbours]
    triangles = [
        (a, b, c) for a, b in board.seams for c in graph[a] & graph[b] if b < c
    ]
    counts = [0, 0, 0, 0]
    for f in triangles:
        counts[len(region.intersection(f))] += 1
        if set(f) <= region:
            continue
        local = [tuple(sorted(e)) for e in zip(f, f[1:] + f[:1])]
        for tips in product(*[(fixed[e],) if e in fixed else e for e in local]):
            assert set(tips) != set(f)
    h0, h1 = law.local_term_parts()
    h = h0 + 8 * h1
    for q in (0, 1):
        # Interior site is first. A zero coin permits the reverse transfer,
        # so that no occupied interior state can load the empty exterior.
        head = np.zeros(11)
        for d in (-1, 1):
            head[1 + law.LABELS.index((0, d))] = 1 / np.sqrt(2)
        z = np.kron(head, np.eye(9)[law.SRC.index((1 + q, 1))])
        np.testing.assert_allclose(h @ z, 8 * z, atol=1e-12)
    return dict(
        owners=board.cells,
        ports=path,
        reference=(6 - board.degree[list(path)]).tolist(),
        boundary_heads=boundary,
        triangles_by_ring_ports=counts,
        frozen_exterior_seams=len(fixed),
        native_idle_energy_per_seam=8,
    )


def certificate(board, path):
    edges = {tuple(sorted(x)) for x in zip(path, path[1:] + path[:1])}
    region = set(path)
    assert {e for e in board.seams if set(e) <= region} == edges
    return shield(board, path)
