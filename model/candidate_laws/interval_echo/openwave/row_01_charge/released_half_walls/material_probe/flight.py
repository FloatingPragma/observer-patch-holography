"""Native two-rail instrument with a contrast contact to both wall endpoints."""

from functools import cache
from itertools import product

import numpy as np
from scipy.linalg import expm

from .. import readers, ring

DWELLS = (1 / 32, 1 / 16)
BOARDS = ((0, 14), (1, 14), (2, 14), (0, 16))
ANALYZERS = ("phase", "rail")


@cache
def mount(level, n):
    ports = ring.geometry(level, n)
    board = ring.federation(level)
    on_ring = set(ports)

    def stems(offsets):
        for j in sorted(offsets):
            u, v = ports[j], ports[(j + 1) % n]
            pair = {u, v}
            for a in sorted(
                set(board.neighbours[u]) & set(board.neighbours[v]) - on_ring
            ):
                if set(board.neighbours[a]) & on_ring != pair:
                    continue
                for b in sorted(set(board.neighbours[a]) - on_ring):
                    if b // 12 == a // 12 and not set(board.neighbours[b]) & on_ring:
                        yield (j, a, b)

    def separate(left, right):
        a, b = left
        u, v = right
        return a // 12 != u // 12 and not (
            set(board.neighbours[a]) | set(board.neighbours[b])
        ) & {u, v}

    for left, right in product(stems(range(n - 5, n - 2)), stems(range(2, 5))):
        if left[0] % 2 != right[0] % 2 and separate(left[1:], right[1:]):
            break
    else:
        raise ValueError("No pair of disjoint two-ended contacts")
    options = sorted(
        (min(a, b), max(a, b))
        for a, b in board.seams
        if a // 12 == b // 12
        and not ({a, b} | set(board.neighbours[a]) | set(board.neighbours[b])) & on_ring
    )
    sham = next((a, b) for a in options for b in options if separate(a, b))
    return (left, right), sham


def sensed_sum(level, n, tip, loads, probe_sign, kappa):
    """Diagonal two-path energy difference from all physical and private contacts.

    Inputs are sharp ring loads. Exchange is parked on other links during sensing.
    The path ending at tip a is (1, 1-sign); the other is the reverse.
    """
    board = ring.federation(level)
    ports = ring.geometry(level, n)
    index = {p: j for j, p in enumerate(ports)}
    a, b = tip
    high, low = 1, 1 - probe_sign
    assert b in board.neighbours[a]
    result = 0
    for site, peer, order in ((a, b, 1), (b, a, -1)):
        exterior = set(board.neighbours[site]) - {peer}
        local_loads = [loads[..., index[v]] if v in index else 0 for v in exterior]
        local_loads += [0] * (6 - board.degree[site])
        assert len(local_loads) == 5
        for s in local_loads:
            result += order * ((high - s) ** 2 - (low - s) ** 2)
    return kappa * result


def probability(detuning, dwell, analyzer):
    """First-rail Born probability starting in the even real path superposition."""
    omega = np.sqrt(1 + np.square(detuning))
    if analyzer == "phase":
        signed = -detuning * np.sin(omega * dwell) / omega
    elif analyzer == "rail":
        signed = -detuning * (1 - np.cos(omega * dwell)) / omega**2
    else:
        raise ValueError("Unknown probe analyzer")
    return (1 + signed) / 2


def verify_plates():
    law = ring.candidate.native_law
    exchange, contrast = ring.columns()
    for sign in (-1, 1):
        src = [law.SRC.index(w) for w in ((1, 1 - sign), (1 - sign, 1))]
        h = (exchange + ring.K * contrast)[np.ix_(src, src)]
        pulse = expm(-1j * np.pi / 2 * h)
        # Half exchange followed by an independently calibrated reference plate.
        prepared = (
            np.diag([readers.phase(n, np.pi / (2 * sign)) for n in (1, 1 - sign)])
            @ pulse[:, 0]
        )
        np.testing.assert_allclose(
            np.outer(prepared, prepared.conj()), np.ones((2, 2)) / 2, atol=1e-13
        )
        for kappa, x, y, dwell in product((0, 4, 8), range(3), range(3), DWELLS):
            d = -2 * kappa * sign * (x + y)
            post = expm(-1j * dwell * (h + np.diag([d, 0]))) @ prepared
            for setting, terminal in (("phase", pulse @ post), ("rail", post)):
                np.testing.assert_allclose(
                    abs(terminal[0]) ** 2, probability(d, dwell, setting), atol=2e-13
                )
