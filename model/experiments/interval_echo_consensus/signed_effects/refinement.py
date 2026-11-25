"""Coherent signed parent-address lift and its autonomous leakage on the actual boards."""

import argparse
from collections import deque
from fractions import Fraction
import json
import numpy as np
from scipy import sparse
from scipy.linalg import expm

from model.executable_model import federation
from model.experiments.interval_echo_angular.charged_orientation.measure import (
    generator,
)
from model.experiments.mobile_memory.refinement import _cell_seams, certify_step
from .measure import HERE, SIGNS, native_controls


def block_paths(board, start):
    graph = [
        [v - start for v in board.neighbours[start + i] if start <= v < start + 48]
        for i in range(48)
    ]
    all_paths = []
    for i in range(48):
        paths = {i: [i]}
        todo = deque([i])
        while todo:
            u = todo.popleft()
            for v in graph[u]:
                if v not in paths:
                    paths[v] = [*paths[u], v]
                    todo.append(v)
        assert len(paths) == 48
        all_paths.append(paths)
    return all_paths


def protected_lift(board, parent, splitters):
    base = 48 * parent
    paths = block_paths(board, base)
    degrees = board.degree[base : base + 48]
    states = np.zeros((2, 48, 12), complex)
    states[:, :12] = np.eye(12)
    swap_count = 0

    def swap(word):
        nonlocal swap_count
        for a, b in word:
            assert base + b in board.neighbours[base + a]
            states[:, [a, b]] = states[:, [b, a]]
            swap_count += 1

    for a in range(12):
        for sending, receiving in ((0, 1), (0, 2), (1, 3)):
            i, j = 12 * sending + a, 12 * receiving + a
            walk = paths[i][j]
            routing = list(zip(walk[:-2], walk[1:-1]))
            swap(routing)
            p = walk[-2]
            assert base + j in board.neighbours[base + p]
            for k, u in enumerate(splitters):
                states[k, [p, j]] = u @ states[k, [p, j]]
            swap(reversed(routing))
            # Bring the receiving coordinate to an available degree-five reference.
            route = min((w for t, w in paths[j].items() if degrees[t] == 5), key=len)
            routing = list(zip(route, route[1:]))
            swap(routing)
            for k, sign in enumerate(SIGNS):
                frequency = 8 * ((1 + sign) ** 2 - 1)
                states[k, route[-1]] *= np.exp(-1j * frequency * 3 * np.pi / 16)
            swap(reversed(routing))
    # Omitted background phase is the same on all columns and signs.
    target = np.tile(np.eye(12), (4, 1)) / 2
    for arm in states:
        np.testing.assert_allclose(arm, -1j * target, atol=2e-13, rtol=0)
    lift_cost = swap_count
    encoded = states.copy()
    for mask in range(4):
        states[:] = encoded
        characters = [(-1) ** ((child & mask).bit_count()) for child in range(4)]
        for site in range(48):
            if characters[site // 12] == 1:
                continue
            route = min((w for t, w in paths[site].items() if degrees[t] == 5), key=len)
            routing = list(zip(route, route[1:]))
            swap(routing)
            for k, sign in enumerate(SIGNS):
                frequency = 8 * ((1 + sign) ** 2 - 1)
                states[k, route[-1]] *= np.exp(-1j * frequency * np.pi / 8)
            swap(reversed(routing))
        word = target * np.repeat(characters, 12)[:, None]
        for original, arm in zip(encoded, states):
            np.testing.assert_allclose(arm, -1j * word, rtol=0, atol=3e-13)
            np.testing.assert_allclose(
                original.conj().T @ arm, np.eye(12) * int(mask == 0), rtol=0, atol=3e-13
            )
    return lift_cost


def free_leakage(coarse, fine, sign):
    fine_port = np.arange(fine.n)
    parent_port = (fine_port // 48) * 12 + fine_port % 12
    w = sparse.csc_matrix(
        (np.ones(fine.n, dtype=np.int64), (fine_port, parent_port)),
        shape=(fine.n, coarse.n),
    )
    twice = 2 * generator(fine, sign, 8)
    b = twice.astype(np.int64)
    assert not np.any((twice - b).data)
    bw = b @ w
    residual = 4 * bw - w @ (w.T @ bw)  # 16 (1-V V*) H V, V=W/2
    squares = np.asarray(residual.multiply(residual).sum(axis=0)).ravel()
    # Include only the part strictly outside each sending parent.
    coo = residual.tocoo()
    outside = coo.row // 48 != coo.col // 12
    counts = np.zeros(coarse.n, dtype=np.int64)
    np.add.at(counts, coo.col[outside], coo.data[outside] ** 2)
    traces = counts.reshape(-1, 12).sum(axis=1)
    assert traces.min() >= 48  # denominator 256; 3/16 per owner
    return dict(
        nonzero_columns=int(np.count_nonzero(squares)),
        mean_squared_leakage=str(Fraction(int(squares.sum()), 256 * coarse.n)),
        max_squared_leakage=str(Fraction(int(squares.max()), 256)),
        exterior_trace_per_parent_range=[
            str(Fraction(int(x), 256)) for x in (traces.min(), traces.max())
        ],
    )


def report():
    h, _ = native_controls()
    splitters = [expm(-1j * np.pi / 2 * g) for g in h]
    out = []
    for level in (0, 1):
        coarse, fine = federation(level), federation(level + 1)
        certify_step(coarse.cells, _cell_seams(coarse), fine.cells, _cell_seams(fine))
        costs = [protected_lift(fine, p, splitters) for p in range(coarse.cells)]
        out.append(
            dict(
                coarse_owners=coarse.cells,
                fine_owners=fine.cells,
                protected_child_encodings=len(costs),
                half_exchanges_per_parent=36,
                degree_phases_per_parent=36,
                route_swap_range=[min(costs), max(costs)],
                free_by_sign={str(s): free_leakage(coarse, fine, s) for s in SIGNS},
            )
        )
    return dict(
        transitions=out,
        any_unitary_child_rechart_exterior_trace_floor="3/16",
        any_unitary_child_rechart_singular_norm_floor="1/8",
        child_phase_words=4,
        unbiased_phase_word_recovery="1/4",
        half_coherence_phase_word_schedule=[5, 1, 1, 1],
        half_coherence_recovery="5/8",
        half_coherence_entropy_derivative="-3 log(5)/4",
        coherent_entropy_deficit="log(4)",
        coherent_recovery="1",
        dephased_child_recovery="1/4",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    text = json.dumps(report(), indent=2, allow_nan=False) + "\n"
    if parser.parse_args().write:
        (HERE / "refinement.json").write_text(text)
    else:
        print(text, end="")
