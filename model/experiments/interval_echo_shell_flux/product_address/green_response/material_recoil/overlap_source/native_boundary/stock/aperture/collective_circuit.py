"""Find a four-relay ruler built entirely from contacts of the native carrier."""

from functools import lru_cache
from itertools import pairwise

import networkx as nx
import numpy as np

from . import collective_source as source


def edge(a, b):
    return tuple(sorted((int(a), int(b))))


@lru_cache(maxsize=None)
def cell_paths(a, b):
    from model.executable_model.carrier import CELL_SEAMS

    graph = nx.Graph(CELL_SEAMS)
    examples = {}
    # All simple paths in a twelve-port cell. Prefer lexicographic
    # retrieval after the length is fixed by the circuit design.
    for p in nx.all_simple_paths(graph, a, b):
        length = len(p) - 1
        if length not in examples or p < examples[length]:
            examples[length] = p
    return examples


def owner_cycles(board, layout, rng):
    active = {edges[-1] for _, edges in layout}
    closed = {e for _, edges in layout for e in edges[:-1]}
    ext = {
        edge(a // 12, b // 12): (a, b)
        for a, b in board.external
        if (a, b) not in closed
    }
    links = list(ext)
    seen = set()
    for _ in range(500):
        rng.shuffle(links)
        graph = nx.Graph()
        graph.add_edges_from(links)
        for path in nx.cycle_basis(graph, root=0):
            # Start from the least owner and choose the deterministic orientation.
            i = path.index(min(path))
            path = path[i:] + path[:i]
            if path[-1] < path[1]:
                path = path[:1] + path[:0:-1]
            key = tuple(path)
            if key in seen:
                continue
            seen.add(key)
            seams = [ext[edge(path[j - 1], owner)] for j, owner in enumerate(path)]
            if sum(e in active for e in seams) == 4:
                yield path, seams


def try_cycle(board, path, seams, active, clear):
    ends = []
    for i, owner in enumerate(path):
        enter = next(x for x in seams[i] if x // 12 == owner)
        leave = next(x for x in seams[(i + 1) % len(path)] if x // 12 == owner)
        ends.append(cell_paths(enter % 12, leave % 12))
    marked = [i for i, e in enumerate(seams) if e in active]
    # Enlarge only the interior of the longest source-free arc. Keep the
    # other arcs short so the test concerns dependence on four measured tolls.
    gaps = [
        [(a + t) % len(path) for t in range((b - a) % len(path))]
        for a, b in zip(marked, marked[1:] + marked[:1], strict=True)
    ]
    for gap in sorted(gaps, key=len, reverse=True):
        if len(gap) < 3:
            continue
        lengths = [min(p) for p in ends]
        middle = sorted(gap, key=lambda i: abs(gap.index(i) - (len(gap) - 1) / 2))
        for _ in range(sum(max(ends[i]) - lengths[i] for i in middle) + 1):
            ports = [
                12 * owner + k
                for owner, options, length in zip(path, ends, lengths, strict=True)
                for k in options[length]
            ]
            size = len(ports)
            if size % 4 == 1:
                circ = ports + ports
                L = (size + 9) // 2
                used = set(ports)
                for j, X in enumerate(ports):
                    V, Y = circ[j + L - 1 : j + L + 1]
                    if min(clear[X], clear[V]) < 2:
                        continue
                    if any(
                        edge(a, b) in active for a, b in pairwise(circ[j : j + L + 1])
                    ):
                        continue
                    if any(
                        set(board.neighbours[u]) & used
                        != {ports[(k - 1) % size], ports[(k + 1) % size]}
                        for u, k in ((X, j), (V, (j + L - 1) % size))
                    ):
                        continue
                    if set(board.neighbours[X]) & set(board.neighbours[V]):
                        continue
                    return ports, j, L
            for i in middle:
                if lengths[i] < max(ends[i]):
                    lengths[i] += 1
                    break
            else:
                break
    return None


def find():
    from . import collective_geometry as geom

    board = source.tidal.lab.energy.pair.stock.observations.federation(1)
    rng = np.random.default_rng(619971)
    checked = 0
    for mode in ("lexicographic", "dispersed"):
        layout = source.packing(board, mode, 8)
        s = geom.survey(board, layout, False)
        active = set(s.contacts)
        for owners, seams in owner_cycles(board, layout, rng):
            checked += 1
            result = try_cycle(board, owners, seams, active, s.clearance)
            if result:
                return board, layout, s, mode, checked, result
    raise ValueError(f"No circuit found among {checked} four-relay owner cycles")
