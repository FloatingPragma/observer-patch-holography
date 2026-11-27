"""Parity-selected unthreaded charged-vacancy square with native commissioning."""

from collections import deque
from itertools import permutations, product
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.executable_model import federation
from .experiment import LAW, SHOTS, TIMES, dark, probabilities, swap
from .measure import serialize
from .passive import local_certificate

START = (0, 1, 0, 2)
WORDS = tuple(sorted(set(permutations(START))))
RING = ((0, 1), (1, 2), (2, 3), (3, 0))
LEFT, RIGHT = ((0, 1), (2, 3)), ((0, 3), (2, 1))
TARGET = tuple(i for i, w in enumerate(WORDS) if w[1] == w[3] == 0)


def mount(level):
    board = federation(level)
    for a in range(board.n):
        for b in board.neighbours[a]:
            for c in board.neighbours[b]:
                for d in board.neighbours[c]:
                    sites = (a, b, c, d)
                    if (
                        min(sites) == a
                        and b < d
                        and len(set(sites)) == 4
                        and a in board.neighbours[d]
                        and all(board.degree[p] == 5 for p in sites)
                    ):
                        return board, sites
    return board, None


def commission(level):
    board, sites = mount(level)
    if sites is None:
        return dict(cells=board.cells, enrolled=False)
    blanks = {tuple(sorted((sites[i], sites[j]))) for i, j in RING}
    cut = [e for e in board.seams if set(sites).intersection(e) and e not in blanks]
    w = START
    operations = []

    def move(goal):
        nonlocal w
        queue = deque([(w, ())])
        seen = {w}
        while queue:
            x, path = queue.popleft()
            if goal(x):
                for i, j in path:
                    operations.append(("swap", sites[i], sites[j]))
                w = x
                return
            for e in RING:
                y = swap(x, *e)
                if y not in seen:
                    queue.append((y, (*path, e)))
                    seen.add(y)
        raise ValueError("No sharp commission path")

    for a, b in cut:

        def ready(x):
            load_a = x[sites.index(a)] if a in sites else 1
            load_b = x[sites.index(b)] if b in sites else 1
            return load_a == load_b

        move(ready)
        operations.append(("equal_write", a, b))
    move(lambda x: x == START)
    return dict(
        cells=board.cells,
        ports=sites,
        enrolled=True,
        cut=cut,
        operations=operations,
        common_private_cost_units=int(
            sum(board.degree == 5) - 4 + sum(n * n for n in w)
        ),
    )


def hamiltonian(words, bonds, k):
    h = np.zeros((len(words), len(words)))
    g = LAW.seam_term(k)
    for j, w in enumerate(words):
        for a, b in bonds:
            pair = (w[a], w[b])
            u = LAW.SRC.index(pair)
            h[j, j] += g[u, u]
            if pair[0] != pair[1]:
                h[words.index(swap(w, a, b)), j] += g[LAW.SRC.index(pair[::-1]), u]
    return h


def trace(psi, seq, k=8, words=WORDS):
    for e in seq:
        psi = expm(-1j * np.pi * hamiltonian(words, (e,), k)) @ psi
    return psi


def inject(sign, k):
    src = np.eye(len(WORDS))[:, WORDS.index(START)]
    adjacent = trace(src, ((2, 3),), k)
    output = (adjacent + sign * trace(adjacent, ((1, 2),), k)) / 2
    p = float(np.vdot(output, output).real)
    return p, trace(output / np.sqrt(p), ((2, 3),), k)


def compare():
    result = {}
    for sign in (1, -1):
        p, src = inject(sign, 8)
        adj = trace(src, ((2, 3),))
        circuit = (*LEFT, *RIGHT[::-1])
        stationary_host = trace(adj, ((1, 2),))
        result[str(sign)] = dict(
            herald=p,
            dark=dark(trace(src, LEFT), trace(src, RIGHT)),
            end_probability=float(sum(abs(trace(src, LEFT)[list(TARGET)]) ** 2)),
            once_dark=dark(trace(src, circuit), src),
            twice_dark=dark(trace(src, circuit * 2), src),
            host_only_dark=dark(adj, stationary_host),
            corrected_dark=dark(trace(src, circuit), trace(stationary_host, ((2, 3),))),
        )
    single = tuple(sorted(set(permutations((0, 1, 1, 1)))))
    ref = np.eye(4)[:, single.index((0, 1, 1, 1))]
    result["one_vacancy_winding_dark"] = dark(trace(ref, RING, 8, single), ref)
    # Fill the vacancy at port 2, leaving the measured 1,2 host pair and its
    # parity coherent at ports 1,3. The comparator uses the same native route.
    one_start = (0, 1, 1, 2)
    marked = tuple(sorted(set(permutations(one_start))))
    launch = np.eye(len(marked))[:, marked.index(one_start)]
    adj = trace(launch, ((2, 3),), 8, marked)
    for sign in (1, -1):
        filtered = (adj + sign * trace(adj, ((1, 2),), 8, marked)) / np.sqrt(2)
        src = trace(filtered, ((2, 3),), 8, marked)
        result[str(sign)]["marked_one_winding_dark"] = dark(
            src, trace(src, RING, 8, marked)
        )
        result[str(sign)]["marked_one_three_windings_dark"] = dark(
            src, trace(src, RING * 3, 8, marked)
        )
    return result


def run():
    local_certificate()
    # The same zero collar can be written on sharp 00 or 11 for internal chords.
    e = LAW.echo()
    collar = []
    for a in (0, 1):
        i = LAW.SRC.index((a, a))
        out = e[:, i].reshape(-1, 9)
        collar.append(out[:, i])
    np.testing.assert_allclose(*collar)
    boards = [commission(i) for i in range(3)]
    if any(not b["enrolled"] for b in boards):
        raise ValueError("Held mounts missing; record refusals explicitly")
    records = []
    banks = [np.random.default_rng(seed) for seed in (703, 1703)]
    for k, t in product((8, 4), TIMES):
        h = hamiltonian(WORDS, RING, k)
        colour = [WORDS.index(tuple(0 if n == 0 else 3 - n for n in w)) for w in WORDS]
        np.testing.assert_array_equal(h, h[colour][:, colour])
        row = []
        for sign in (1, -1):
            herald, src = inject(sign, k)
            psi = expm(-1j * t * h) @ src
            target = float(sum(abs(psi[list(TARGET)]) ** 2))
            row.extend((herald * target, herald * (1 - target)))
        prob = probabilities(row)
        records.append(
            dict(
                k=k,
                time=t,
                joint=prob.tolist(),
                counts=[rng.multinomial(SHOTS, prob).tolist() for rng in banks],
            )
        )
    return dict(trials=SHOTS, boards=boards, addressed=compare(), free=records)


if __name__ == "__main__":
    Path(__file__).with_name("vacancies.json").write_text(
        json.dumps(serialize(run()), indent=1, allow_nan=False) + "\n"
    )
