"""Validate first-tick addresses with a separate 240-by-240 pair evolution."""

import json
import math
from pathlib import Path

import numpy as np
from scipy.special import jv
from scipy.stats import binom

from model.executable_model import federation

HERE = Path(__file__).resolve().parent


def probabilities(board, start, exchange):
    # Psi[a,b] has the double at a and the single at b. Set equal-position
    # words to zero. Each actual seam between the two occupied ports is one
    # exchange move, including its reversed pair, not two hops to forbidden words.
    adjacency = board.adjacency().tocsr()
    touch = adjacency.toarray()
    degree = np.asarray(adjacency.sum(axis=0)).ravel()
    diagonal = exchange * (degree[:, None] + degree[None, :] - touch) / 2 - 32 * touch

    def centered(psi):
        hops = adjacency @ psi + (adjacency @ psi.T).T + touch * psi.T
        out = diagonal * psi - exchange / 2 * hops + 10 * psi
        np.fill_diagonal(out, 0)
        return out / 22

    # Gershgorin: native exchange Laplacian/2 lies in [0,12]; static
    # diagonal in [-32,0]. The shifted spectrum is in [-22,22].
    previous = np.zeros((board.n, board.n))
    previous[2, start] = 1
    current = centered(previous)
    evolved = jv(0, 11) * previous - 2j * jv(1, 11) * current
    for k in range(2, 61):
        following = 2 * centered(current) - previous
        evolved += 2 * (-1j) ** k * jv(k, 11) * following
        previous, current = current, following
    # Unitary scalar exp(+5i) cancels from all probabilities.
    np.testing.assert_allclose(np.vdot(evolved, evolved).real, 1, rtol=0, atol=1e-12)
    tick = abs(evolved[[2, 3]]) ** 2
    tick[:, [2, 3]] = 0
    port = tick.sum(axis=0)
    return np.append(port[:12], port[12:].sum())


def run():
    data = json.loads((HERE / "records.json").read_text())
    board = federation(0)
    n = data["shots_per_bit"]
    cases = ((82, 74, 1), (72, 77, 82), (182, 184, 72))
    cache = {}
    tests = 0
    for zero, one, prev in cases:
        for arm, starts, exchange in (
            ("passive", (zero, one), 1),
            ("relay", (prev, one), 1),
            ("sham", (zero, zero), 1),
            ("half_exchange", (zero, one), 0.5),
        ):
            row = data["stations"][str(zero)][arm]["1"]
            keys = row["keys"]
            for bit, start in enumerate(starts):
                key = (start, exchange)
                if key not in cache:
                    cache[key] = probabilities(board, start, exchange)
                p = cache[key]
                # One fixed challenge bank per arm and input. Bonferroni
                # binomial prediction intervals include both tails.
                al = 0.01 / (2 * 3 * 4 * 2 * 13)
                lower = binom.ppf(al, n, p)
                upper = binom.ppf(1 - al, n, p)
                observed = [
                    row["test"][bit][keys.index(f"unmarked:1:{a}")] for a in range(13)
                ]
                assert np.all(lower <= observed) and np.all(observed <= upper)
                tests += len(observed)
    tail = (
        2 * math.exp(11**2 / (4 * 62)) * (5.5**61 / math.factorial(61)) / (1 - 5.5 / 62)
    )
    print(tests, "first-service bins inside simultaneous 99% prediction bounds")
    print("Degree-60 Chebyshev uniform amplitude truncation bound", tail)


if __name__ == "__main__":
    run()
