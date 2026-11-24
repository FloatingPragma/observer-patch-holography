"""Commission the physical collars on a whole board, conserving every supplied load."""

import argparse
import hashlib
import json
from collections import deque
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.sparse.csgraph import connected_components

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from .experiment import SEEDS, grow

HERE = Path(__file__).resolve().parent


def native():
    law = candidate.native_law
    kappa = candidate.KAPPA
    exchange, contrast = law.local_term_parts()
    blank = np.eye(law.DIM)[:, : law.N_SRC]
    words = law.SRC
    swap = [words.index((b, a)) for a, b in words]
    term = exchange + kappa * contrast
    np.testing.assert_allclose(
        expm(-1j * np.pi * term) @ blank, blank[:, swap], atol=7e-14
    )
    # Move the second endpoint to an addressed degree-five port, phase the
    # private zero-load reference there, then undo that permutation.
    ref_times = np.pi / (8 * kappa)
    reference = blank[:, [words.index((a, 0)) for a in range(3)]]
    phases = np.exp(-1j * ref_times * kappa * np.arange(3) ** 2)
    np.testing.assert_allclose(
        expm(-1j * ref_times * kappa * contrast) @ reference,
        reference * phases,
        atol=1e-14,
    )
    split = expm(-1j * np.pi * term / 2) @ blank[:, words.index((2, 0))]
    split *= np.tile([phases[b] for _, b in words], law.DIM // len(words))
    written = (law.echo() @ split).reshape(-1, law.N_SRC)
    conditional = written[:, words.index((1, 1))]
    np.testing.assert_allclose(np.linalg.norm(conditional), 1, atol=1e-13)
    conditional = conditional * np.exp(
        -1j * np.angle(conditional[np.argmax(abs(conditional))])
    )
    idle = (law.echo() @ blank[:, words.index((1, 1))]).reshape(-1, law.N_SRC)[
        :, words.index((1, 1))
    ]
    empty = (law.echo() @ blank[:, words.index((0, 0))]).reshape(-1, law.N_SRC)[
        :, words.index((0, 0))
    ]
    np.testing.assert_allclose(idle, empty, atol=1e-14)
    residual = 0.0
    for a, b in ((0, 0), (0, 1), (1, 0), (1, 1)):
        word = np.eye(law.N_SRC)[:, words.index((a, b))]
        for collar, energy in [(idle, 0), (conditional, 4 * kappa * (1 - a * b))]:
            col = np.kron(collar, word)
            residual = max(residual, np.linalg.norm(term @ col - energy * col))
    assert residual < 1e-12
    return round(float(residual), 13)


def edge(a, b):
    return tuple(sorted((int(a), int(b))))


def path_to(start, allowed, goal):
    parent = {start: None}
    queue = deque([start])
    while queue:
        tip = queue.popleft()
        if goal(tip):
            route = [tip]
            while parent[route[-1]] is not None:
                route.append(parent[route[-1]])
            return route[::-1]
        for v in allowed[tip]:
            if v not in parent:
                parent[v] = tip
                queue.append(v)
    raise ValueError("No physical path to the requested port")


def prepare(level, seed):
    board = federation(level)
    components = grow(board, seed)
    enrolled = max((r for closed, r in components if not closed), key=len)
    dimers, connectors = set(), set()
    for closed, strand in components:
        dimers.update(edge(a, b) for a, b in strand)
        links = zip(strand, strand[1:] + (strand[:1] if closed else []))
        connectors.update(edge(left[1], right[0]) for left, right in links)
    marked_ports = {site for pair in connectors for site in pair}
    assert len(marked_ports) == 2 * len(connectors)
    seams = {edge(a, b) for a, b in board.seams}
    assert dimers.isdisjoint(connectors) and (dimers | connectors) <= seams
    graph = board.adjacency().tolil()
    for a, b in connectors:
        graph[a, b] = graph[b, a] = 0
    graph = graph.tocsr()
    graph.eliminate_zeros()
    assert connected_components(graph, directed=False)[0] == 1
    allowed = [
        sorted(graph.indices[graph.indptr[i] : graph.indptr[i + 1]].tolist())
        for i in range(board.n)
    ]
    n = len(enrolled)
    available = sorted(set(range(board.n)) - marked_ports)
    assert len(available) >= n
    # Stocks for commissioning are sharp occupations. None contain a wave.
    stock = np.ones(board.n, dtype=int)
    stock[available[:n]] = 0
    for a, b in connectors:
        stock[a], stock[b] = 2, 0
    conserved = int(sum(stock))
    assert conserved == board.n - n
    counts = dict(
        write_conditional=0,
        phase_round_trip_swaps=0,
        degree_phase=0,
        hole_swaps=0,
        orient_swaps=0,
        write_idle=0,
    )
    addresses = hashlib.sha256()

    def log(op, a, b):
        addresses.update(f"{op} {a} {b}\n".encode("ascii"))

    for a, b in sorted(connectors):
        assert (stock[a], stock[b]) == (2, 0)
        route = path_to(b, allowed, lambda v: board.degree[v] == 5)
        path = list(zip(route, route[1:]))
        log("half-contact", a, b)
        for x, y in path:
            log("swap", x, y)
        log("private-phase", route[-1], route[-1])
        for x, y in reversed(path):
            log("swap", x, y)
        log("echo", a, b)
        # Native column above proves this whole service returns every other
        # register and yields the sharp 11 word times its independent collar.
        stock[a], stock[b] = 1, 1
        counts["write_conditional"] += 1
        counts["degree_phase"] += 1
        counts["phase_round_trip_swaps"] += 2 * len(path)
    target_holes = {a for a, b in enrolled}

    def swap_load(a, b):
        stock[a], stock[b] = stock[b], stock[a]
        log("swap", a, b)

    for dest in sorted(target_holes):
        if stock[dest] == 0:
            continue
        route = path_to(
            dest, allowed, lambda v: v not in target_holes and stock[v] == 0
        )
        path = list(zip(route, route[1:]))
        # A transposition fixes all intermediate occupations.
        for a, b in path + list(reversed(path[:-1])):
            swap_load(a, b)
            counts["hole_swaps"] += 1
    assert set(np.flatnonzero(stock == 0)) == target_holes
    mate = {x: y for p in enrolled for x, y in (p, p[::-1])}
    for a, b in sorted(seams - dimers - connectors):
        for v in (a, b):
            if stock[v] == 0:
                assert stock[mate[v]] == 1 and edge(v, mate[v]) in dimers
                swap_load(v, mate[v])
                counts["orient_swaps"] += 1
        assert (stock[a], stock[b]) == (1, 1)
        log("echo-idle", a, b)
        counts["write_idle"] += 1
    # Set a fixed orientation. Any specified product orientation follows with
    # one more addressed dimer swap per flipped site before free release.
    for a, b in enrolled:
        if stock[a] == 1:
            swap_load(a, b)
            counts["orient_swaps"] += 1
    assert set(np.flatnonzero(stock == 0)) == target_holes
    assert sum(stock) == conserved
    return dict(
        cells=board.cells,
        seed=seed,
        dimers_total=len(dimers),
        dimers_enrolled=n,
        ports=board.n,
        load_total=conserved,
        native_services=counts,
        operation_digest=addresses.hexdigest(),
    )


def run():
    return dict(
        max_native_residual=native(),
        boards=[prepare(i, s) for i in range(3) for s in SEEDS],
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", action="store_true")
    args = parser.parse_args()
    observed = run()
    path = HERE / "formation.json"
    if args.record:
        with path.open("x") as f:
            json.dump(observed, f, indent=1)
            f.write("\n")
    elif json.loads(path.read_text()) != observed:
        raise SystemExit("Full-board commissioning differs from the receipt")
    else:
        print("Verified all nine closed native material preparations and load stocks.")
