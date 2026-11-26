"""Commission an invariant native Y junction and check the physical core reader."""

from collections import deque
from itertools import permutations

import numpy as np
from scipy.linalg import expm

from model.executable_model import federation
from .experiment import (
    LAW,
    RING,
    START,
    TIMES,
    WORDS,
    effective,
    prepare,
    probabilities,
    positions,
    swap,
)


def local_certificate():
    sharp = np.eye(LAW.DIM)[:, : LAW.N_SRC]
    echo = LAW.echo()
    x, v = LAW.local_term_parts()
    index = LAW.SRC.index((1, 1))
    prepared = (echo @ sharp[:, index]).reshape(-1, 9)
    collar = prepared[:, index]
    idle = np.kron(collar[:, None], np.eye(9))
    np.testing.assert_allclose(prepared, np.outer(collar, np.eye(9)[index]), atol=1e-14)
    np.testing.assert_allclose(idle.T @ idle, np.eye(9), atol=1e-14)
    largest = 0.0
    for k in (4, 8):
        term = x + k * v
        largest = max(largest, float(np.linalg.norm(term @ idle, axis=0).max()))
        np.testing.assert_allclose(term @ idle, 0, atol=1e-12)
        np.testing.assert_allclose(term @ sharp, sharp @ LAW.seam_term(k), atol=1e-12)
        ref = sharp[:, [LAW.SRC.index((a, 0)) for a in range(3)]]
        np.testing.assert_allclose(v @ ref, ref @ np.diag([0, 1, 4]), atol=1e-12)
        physical_plate = expm(-1j * np.pi / (4 * k) * k * v) @ ref
        np.testing.assert_allclose(
            physical_plate, ref @ np.diag([1, np.exp(-1j * np.pi / 4), -1]), atol=1e-12
        )
    for col, (a, b) in enumerate(LAW.SRC):
        written = echo @ sharp[:, col]
        for d in range(-2, 3):
            mask = np.r_[np.zeros(9), np.repeat([q == d for _, q in LAW.LABELS], 9)]
            np.testing.assert_allclose(
                echo.T @ (mask * written), int(a - b == d) * sharp[:, col], atol=1e-12
            )
    return largest


def sharp_route(word, goal):
    queue = deque([(word, ())])
    seen = {word}
    while queue:
        before, steps = queue.popleft()
        if goal(before):
            return steps
        for leaf in range(1, len(word)):
            after = swap(before, 0, leaf)
            if after not in seen:
                seen.add(after)
                queue.append((after, (*steps, leaf)))
    raise ValueError("No blank-star route commissions this contact")


def commission(level, arms):
    board = federation(level)
    guide = None
    for centre, neighbours in enumerate(board.neighbours):
        private = [p for p in neighbours if board.degree[p] == 5]
        if board.degree[centre] == 5 and len(private) >= 4:
            guide = (centre, *private[:arms])
            break
    if guide is None:
        return dict(cells=board.cells, arms=arms, enrolled=False)
    enrolled = set(guide)
    blank = {tuple(sorted((guide[0], p))) for p in guide[1:]}
    cut = [e for e in board.seams if set(e) & enrolled and e not in blank]
    start = START if arms == 3 else (*START, 0)
    word = start
    operations = []

    def do_swaps(route):
        nonlocal word
        for b in route:
            operations.append(("swap", guide[0], guide[b]))
            word = swap(word, 0, b)

    for a, b in cut:
        ports = [guide.index(p) for p in (a, b) if p in enrolled]
        do_swaps(sharp_route(word, lambda w: all(w[p] == 1 for p in ports)))
        # Every exterior is a sharp unit and never moves on this calendar.
        operations.append(("equal_write", a, b))
    do_swaps(sharp_route(word, lambda w: w == start))
    assert word == start
    scalar_units = sum(
        board.degree[p] == 5 for p in range(board.n) if p not in enrolled
    )
    scalar_units += sum(n * n for n in word)
    # Three or four joint QND contrast reads on fresh flight copies. Knowing
    # the fixed stock sum gives the centre and then all leaves.
    for w in set(permutations(start)):
        d = [w[0] - v for v in w[1:]]
        recovered_centre = (sum(start) + sum(d)) // len(w)
        assert (recovered_centre, *(recovered_centre - x for x in d)) == w
    # Pull every physical blank and degree column back to the enrolled stock.
    held = []
    if arms == 3:
        for k in (8, 4):
            h = np.zeros((12, 12))
            native = LAW.seam_term(k)
            for j, w in enumerate(WORDS):
                loads = np.ones(board.n, dtype=int)
                loads[list(guide)] = w
                private = int(np.sum(loads[board.degree == 5] ** 2))
                assert private == scalar_units
                # The scalar is common to all words and both arms.
                for a, b in board.seams:
                    if (a, b) in cut:
                        continue
                    pair = (loads[a], loads[b])
                    idx = LAW.SRC.index(pair)
                    h[j, j] += native[idx, idx]
                    if pair != pair[::-1]:
                        moved = swap(w, guide.index(a), guide.index(b))
                        h[WORDS.index(moved), j] += native[
                            LAW.SRC.index(pair[::-1]), idx
                        ]
            energy, modes = np.linalg.eigh(h)
            for sign in (1, -1):
                _, incoming = prepare(sign, k)
                for t in TIMES:
                    wave = modes @ (np.exp(-1j * t * energy) * (modes.T @ incoming))
                    census = probabilities(
                        [
                            sum(
                                abs(wave[j]) ** 2
                                for j, w in enumerate(WORDS)
                                if positions(w) == q
                            )
                            for q in RING
                        ]
                    )
                    predicted = abs(expm(-1j * t * effective(sign, k))[:, 0]) ** 2
                    np.testing.assert_allclose(census, predicted, atol=2e-11)
                    held.append(dict(k=k, sign=sign, time=t, antipode=float(census[3])))
    return dict(
        cells=board.cells,
        arms=arms,
        ports=guide,
        enrolled=True,
        site_loads=start,
        free_held=held,
        quiet_contacts=cut,
        preparation=operations,
        common_reference_units=int(scalar_units),
    )


def run_passive():
    return dict(
        max_quiet_column_residual=local_certificate(),
        mounts=[commission(level, n) for level in range(3) for n in (3, 4)],
    )
