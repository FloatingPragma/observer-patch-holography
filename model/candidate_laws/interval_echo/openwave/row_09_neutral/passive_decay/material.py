"""Address and write the passive four-dimer decay material using native operations."""

from itertools import product
from collections import Counter
from hashlib import sha256

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.row_30_collective.formation import (
    edge,
    path_to,
)
from model.executable_model import federation


def design(board, length):
    adjacent = board.neighbours
    dimers = sorted(
        (p, q)
        for p in range(board.n)
        for q in adjacent[p]
        if board.degree[p] == board.degree[q]
    )
    by_first = {p: [] for p in range(board.n)}
    for dimer in dimers:
        by_first[dimer[0]].append(dimer)

    def attached(source, weak):
        for rail in sorted(adjacent[source[0]]):
            for target in by_first[rail]:
                if not weak or target[1] in adjacent[source[0]]:
                    yield target

    def fit(chain):
        if len(set(sum(chain, ()))) != 2 * len(chain):
            return None
        a, b, c = chain[:3]
        heads = {edge(a[0], b[0]): 32}
        pins = [(a[1], 32), (b[1], 33)]
        for x, y in [(a, c)] + ([(c, chain[3])] if len(chain) == 4 else []):
            heads[edge(x[0], y[0])] = 33
            heads[edge(x[0], y[1])] = 32
            pins.append((x[1], 32))
        sites = set(sum(chain, ()))
        for rail, cost in pins:
            choices = sorted(
                q
                for q in adjacent[rail]
                if q not in sites and edge(rail, q) not in heads
            )
            if not choices:
                return None
            heads[edge(rail, choices[0])] = cost
        return chain, heads

    for a in dimers:
        for b, c in product(attached(a, False), attached(a, True)):
            if len(set(a + b + c)) != 6:
                continue
            for tail in [()] if length == 3 else ((d,) for d in attached(c, True)):
                result = fit((a, b, c) + tail)
                if result:
                    return result
    return None


def hardware():
    for length in (4, 3):
        for level in range(3):
            board = federation(level)
            result = design(board, length)
            if result:
                return board, *result
    raise ValueError("No material route")


def columns(comparison=False):
    law = candidate.native_law
    exchange, contrast = law.local_term_parts()
    full = exchange + candidate.KAPPA * contrast
    tested = full
    if comparison:
        reverse, _ = law._label_ops()
        extreme = np.diag([int(abs(d) == 2) for _, d in law.LABELS])
        label_shift = extreme @ (np.eye(law.N_LAB) - reverse) / 2
        shift = np.zeros_like(full)
        shift[law.N_SRC :, law.N_SRC :] = np.kron(label_shift, np.eye(law.N_SRC))
        tested = full - law.echo() @ shift @ law.echo().T
    blank = np.eye(law.DIM)[:, : law.N_SRC]
    words = law.SRC
    unit = blank[:, words.index((1, 1))]
    idle = (law.echo() @ unit).reshape(-1, len(words))[:, words.index((1, 1))]
    # Half swap gives a relative +/-i; a phase on rail at a private reference
    # makes either real sign using positive exposure (the other three quarters).
    split = expm(-0.5j * np.pi * full) @ blank[:, words.index((2, 0))]
    residual = 0.0
    for cost, multiplier in ((32, 1), (33, 3)):
        time = multiplier * np.pi / (8 * candidate.KAPPA)
        phase = np.tile(
            [np.exp(-1j * candidate.KAPPA * time * b * b) for _, b in words],
            law.DIM // len(words),
        )
        repaired = (law.echo() @ (split * phase)).reshape(-1, len(words))
        head = repaired[:, words.index((1, 1))]
        np.testing.assert_allclose(np.linalg.norm(head), 1, atol=1e-12)
        for a, b in product(range(2), repeat=2):
            word = np.eye(len(words))[:, words.index((a, b))]
            for label, energy in (
                (head, (32 if comparison else cost) * (1 - a * b)),
                (idle, 0),
            ):
                col = np.kron(label, word)
                residual = max(residual, np.linalg.norm(tested @ col - energy * col))
    assert residual < 2e-12
    return float(residual)


def commission(board, dimers, heads, sharp=False):
    sites = set(sum(dimers, ()))
    seams = {edge(a, b) for a, b in board.seams}
    blanks = {edge(*x) for x in dimers}
    paths = [[] for _ in range(board.n)]
    for a, b in sorted(seams - set(heads)):
        paths[a].append(b)
        paths[b].append(a)
    # Each write consumes one conserved 2 and one 0, producing 11.
    stock = np.ones(board.n, dtype=int)
    outside = sorted(
        set(range(board.n)) - sites - set(p for pair in heads for p in pair)
    )
    m, n = len(heads), len(dimers)
    stock[outside[:m]] = 2
    stock[outside[m : 2 * m + n]] = 0
    total = sum(stock)
    counts, digest = Counter(), sha256()

    def record(op, a, b):
        counts[op] += 1
        digest.update(f"{op} {a} {b}\n".encode("ascii"))

    def swap(a, b):
        assert b in paths[a]
        stock[a], stock[b] = stock[b], stock[a]
        record("blank_swap", a, b)

    def bring(p, value, protect):
        route = path_to(p, paths, lambda q: q not in protect and stock[q] == value)
        pairs = list(zip(route, route[1:]))
        for a, b in pairs + pairs[-2::-1]:
            swap(a, b)

    for (a, b), cost in sorted(heads.items()):
        bring(a, 2, set())
        bring(b, 0, {a})
        assert (stock[a], stock[b]) == (2, 0)
        if sharp:
            neighbor = min(q for q in paths[a] if q != b)
            route = path_to(
                neighbor, paths, lambda q: q not in (a, b) and stock[q] == 0
            )
            chain = list(zip(route, route[1:]))
            moves = chain + chain[-2::-1]
            for x, y in moves:
                swap(x, y)
            assert stock[neighbor] == 0 and (stock[a], stock[b]) == (2, 0)
            record("sharp_write", a, b)
            stock[a], stock[b] = 1, 1
            record("split_borrowed_hole", a, neighbor)
            record("parity_wait", a, b)
            record("recombine_read_restore", a, neighbor)
            for x, y in reversed(moves):
                swap(x, y)
        else:
            record("half_swap", a, b)
            route = path_to(b, paths, lambda q: board.degree[q] == 5)
            pairs = list(zip(route, route[1:]))
            for x, y in pairs:
                swap(x, y)
            record(
                "quarter_phase" if cost == 32 else "three_quarter_phase",
                route[-1],
                route[-1],
            )
            for x, y in reversed(pairs):
                swap(x, y)
            record("write", a, b)
            stock[a], stock[b] = 1, 1
    holes = {pair[0] for pair in dimers}
    for p in sorted(holes):
        if stock[p]:
            bring(p, 0, holes)
    mates = {p: q for x, y in dimers for p, q in [(x, y), (y, x)]}
    for a, b in sorted(seams - blanks - set(heads)):
        if not {a, b} & sites:
            continue
        for p in (a, b):
            if stock[p] == 0:
                swap(p, mates[p])
        assert (stock[a], stock[b]) == (1, 1)
        record("write_idle", a, b)
    for p in sorted(holes):
        if stock[p]:
            swap(p, mates[p])
    assert sum(stock) == total and set(np.flatnonzero(stock == 0)) == holes
    return dict(
        total_load=int(total), operations=dict(counts), digest=digest.hexdigest()
    )
