"""Construct the complete vacancy-loop image on the active blank load sector."""

import argparse
from collections import Counter, deque
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.executable_model import cell, federation

HERE = Path(__file__).resolve().parent


def path(graph, start, goal, forbidden=()):
    pending = deque([start])
    previous = {start: None}
    blocked = set(forbidden)
    while pending:
        u = pending.popleft()
        if u == goal:
            result = []
            while u is not None:
                result.append(u)
                u = previous[u]
            return result[::-1]
        for v in graph[u]:
            if v not in previous and v not in blocked:
                pending.append(v)
                previous[v] = u
    raise ValueError("disconnected requested path")


def ears(board):
    """Start with a triangle, then add one open ear with distinct old ends."""
    graph = board.neighbours
    a, b = next((a, b) for a in graph[0] for b in graph[0] if a < b and b in graph[a])
    triangle = [0, a, b]
    present = set(triangle)
    frontier = deque((u, v) for u in triangle for v in graph[u] if v not in present)
    witness = []
    while len(present) < board.n:
        while frontier and frontier[0][1] in present:
            frontier.popleft()
        if not frontier:
            raise ValueError("isolated connected component")
        u, w = frontier.popleft()
        queue = deque([w])
        prev = {w: u}
        terminal = None
        while queue and terminal is None:
            z = queue.popleft()
            for v in graph[z]:
                if v in present:
                    if v != u:
                        terminal = (z, v)
                        break
                elif v not in prev:
                    prev[v] = z
                    queue.append(v)
        if terminal is None:
            raise ValueError("component attaches at just one old vertex")
        z, v = terminal
        ear = [v]
        while z != u:
            ear.append(z)
            z = prev[z]
        ear.append(u)
        ear.reverse()
        assert ear[0] in present and ear[-1] in present and ear[0] != ear[-1]
        assert len(set(ear)) == len(ear) and not (set(ear[1:-1]) & present)
        assert all(v in graph[u] for u, v in zip(ear, ear[1:]))
        witness.append(ear)
        present.update(ear[1:-1])
        frontier.extend((u, v) for u in ear[1:-1] for v in graph[u] if v not in present)
    return triangle, witness


def local_gate():
    # The native exact symbolic proof includes all off-blank rows, and the
    # exchange and squared contrast commute on the invariant nine words.
    checks = law.exact_checks()
    exchange, contrast = law._load_ops()
    np.testing.assert_array_equal(exchange @ contrast, contrast @ exchange)
    assert all(int(c) == c for c in contrast.diagonal())
    h = law.seam_term(8)
    np.testing.assert_allclose(expm(-1j * np.pi * h), exchange, rtol=0, atol=2e-14)
    h0, h1 = law.local_term_parts()
    np.testing.assert_allclose(
        (h0 + 8 * h1)[:, :9], np.vstack((h, np.zeros((90, 9)))), atol=4e-14, rtol=0
    )
    return checks


def move(board, walk):
    """Original label -> final port under swaps with the moving vacancy."""
    assert walk[0] == walk[-1] and all(
        v in board.neighbours[u] for u, v in zip(walk, walk[1:])
    )
    labels = np.arange(board.n)
    for a, b in zip(walk, walk[1:]):
        labels[a], labels[b] = labels[b], labels[a]
    perm = np.argsort(labels)
    assert perm[walk[0]] == walk[0]
    return perm


def witness(board):
    """A simple across-owner cycle with a local triangular wedge in two owners."""
    joins = {}
    adj = [[] for _ in range(board.cells)]
    for u, v in board.external:
        a, b = u // 12, v // 12
        assert (a, b) not in joins
        joins[a, b] = (u, v)
        joins[b, a] = (v, u)
        adj[a].append(b)
        adj[b].append(a)
    nxt = min(adj[0])
    alt = [sorted(ns) for ns in adj]
    alt[0].remove(nxt)
    alt[nxt].remove(0)
    owners = [0, *path(alt, nxt, 0)[:-1]]

    def interior(owner, start, last, need):
        def extend(chain):
            v = chain[-1]
            if v == last:
                if not need:
                    return chain, None
                for i in range(1, len(chain) - 2):
                    a, b, c = chain[i : i + 3]
                    for root in sorted(
                        set(board.neighbours[a])
                        & set(board.neighbours[b])
                        & set(board.neighbours[c])
                    ):
                        if root // 12 == owner and root not in (a, b, c):
                            return chain, (root, [a, b, c])
                return None
            if len(chain) >= 7:
                return None
            for w in board.neighbours[v]:
                if w // 12 == owner and w not in chain:
                    candidate = extend([*chain, w])
                    if candidate:
                        return candidate
            return None

        result = extend([start])
        if not result:
            raise ValueError("no suitable triplet on owner segment")
        return result

    cycle = []
    wedges = []
    for i, owner in enumerate(owners):
        start = joins[owners[i - 1], owner][1]
        end = joins[owner, owners[(i + 1) % len(owners)]][0]
        segment, wedge = interior(owner, start, end, i < 2)
        cycle.extend(segment)
        if wedge:
            wedges.append(wedge)
    assert len(cycle) == len(set(cycle))
    cycle.append(cycle[0])
    generators = []
    arrivals = []
    for root, abc in wedges:
        approach = path(board.neighbours, cycle[0], root, abc)
        g = []
        for a, b in (abc[:2], abc[1:]):
            word = approach + [a, b, root] + approach[-2::-1]
            perm = move(board, word)
            assert np.array_equal(
                np.where(perm != np.arange(board.n))[0], sorted((a, b))
            )
            g.append(perm)
        generators.append(g)
        arrivals.append(approach)
    first, last = (w[1] for w in wedges)
    forward = move(board, cycle)
    back = np.argsort(forward)
    orbit = np.arange(board.n)
    direction = 0
    word = None
    for steps in range(1, len(cycle)):
        orbit = back[orbit]
        if np.array_equal(orbit[first], last):
            word = [cycle[0]] + list(reversed(cycle))[1:] * steps
            direction = -1
            break
    assert word is not None
    transport = move(board, word)
    np.testing.assert_array_equal(transport[first], last)
    for old, new in zip(*generators):
        np.testing.assert_array_equal(transport[old[np.argsort(transport)]], new)
    return dict(
        root=cycle[0],
        owners=owners,
        cycle=cycle,
        cycle_power=steps * direction,
        contact_pulses=len(word) - 1,
        source_triple=first,
        receiving_triple=last,
        wedge_roots=[a for a, _ in wedges],
        approaches=arrivals,
    ), transport


def analyze_bright(board, root):
    """Fold the uniform excess state into one sharp port of the punctured graph."""
    start = next(i for i in range(board.n) if i != root and board.degree[i] == 5)
    tree = [start]
    parent = {start: None}
    for v in tree:
        for w in board.neighbours[v]:
            if w != root and w not in parent:
                tree.append(w)
                parent[w] = v
    assert len(tree) == board.n - 1
    mass = {p: 1 for p in tree}
    # Each two-port real rotation sends (sqrt(a),sqrt(b)) to (sqrt(a+b),0).
    rotations = []
    for w in tree[:0:-1]:
        v = parent[w]
        a, b = mass[v], mass[w]
        rotations.append((v, w, a, b))
        mass[v] += b
        mass[w] = 0
    assert mass[start] == board.n - 1 and sum(mass.values()) == mass[start]
    return start, rotations


def report():
    local_gate()
    boards = {}
    for name, board in (
        ("cell", cell()),
        *((str(n), federation(i)) for i, n in enumerate((20, 80, 320))),
    ):
        triangle, added = ears(board)
        m = board.n - 1
        root, rotations = analyze_bright(board, 0)
        data = dict(
            ports=board.n,
            triangle=triangle,
            ear_count=len(added),
            ear_new_vertex_counts=dict(
                sorted(Counter(len(e) - 2 for e in added).items())
            ),
            ear_certificate_sha256=sha256(
                json.dumps([triangle, added], separators=(",", ":")).encode()
            ).hexdigest(),
            permutation_degree=m,
            invariant_dimensions=[1, m - 1],
            bright_port=root,
            bright_rotations=len(rotations),
        )
        if board.cells > 1:
            data["transport"], _ = witness(board)
        boards[name] = data
    return dict(native_coupling=8, protected_swap_dwell="pi", boards=boards)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    text = json.dumps(report(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (HERE / "holonomy.json").write_text(text)
    else:
        print(text, end="")
