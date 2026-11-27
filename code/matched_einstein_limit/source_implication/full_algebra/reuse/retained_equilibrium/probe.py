"""Native first-hit Kraus blocks on one capped survey and one stored-head preparation."""

from itertools import product

import numpy as np
from scipy.sparse.linalg import expm_multiply

from .. import measure as reuse
from ....material_clocks.measure import density_center

AXES = [
    np.eye(2),
    np.array([[0, 1], [1, 0]]),
    np.array([[0, -1j], [1j, 0]]),
    np.diag([1, -1]),
]


def wave_bank(charge, odd):
    """All unitary dwell kernels from the four actually enrolled stations."""
    formation = reuse.station.iso.pairing.formation
    board = formation.federation(1)
    ports = formation.region(board)
    h = reuse.native_precision.generator(board, ports, charge, bool(odd))
    tick = -0.25j * h
    start = np.zeros((h.shape[0], 16))
    for i, port in enumerate(ports):
        start[4 * port : 4 * port + 4, 4 * i : 4 * i + 4] = np.eye(4)
    raw = expm_multiply(
        tick, start, start=0, stop=12, num=13, traceA=tick.diagonal().sum()
    )
    return ports, raw.reshape(13, board.n, 4, 4, 4).transpose(3, 0, 1, 2, 4)


def histories(bank, outbound=(0, 1)):
    """On each record, the returned bank indexes the unmeasured terminal position."""
    ports, w = bank
    i, j = outbound
    home, relay = ports[i], ports[j]
    u = w[[i, j]]

    def first(waves, dest, other):
        survived = [waves[0].copy()]
        hits = [np.zeros((4, 4), complex)]
        for t in range(1, 13):
            s = waves[t].copy()
            for n in range(1, t):
                s -= other[t - n] @ hits[n]
            hits.append(s[dest].copy())
            s[dest] = 0
            survived.append(s)
        return hits, survived

    f, s = first(u[0], relay, u[1])
    g, z = first(u[1], home, u[0])
    yield ("out_timeout", 0, 0), s[12]
    for a in range(1, 13):
        yield ("ret_timeout", a, 0), z[12 - a] @ f[a]
        for b in range(1, 13 - a):
            yield ("success", a, b), (g[b] @ f[a])[None]


def conditional(k, x):
    return (k @ x @ k.conj().swapaxes(-1, -2)).sum(axis=0)


def kyfan(z):
    return np.linalg.eigvalsh(z)[::-1].cumsum()[:3]


def run():
    x, *_ = density_center()
    for q, par in product((-1, 1), range(2)):
        stock = x[par]
        grouped = {}
        gram = np.zeros((4, 4), complex)
        for label, k in histories(wave_bank(q, par)):
            gram += (k.conj().swapaxes(-1, -2) @ k).sum(axis=0)
            y = conditional(k, stock)
            grouped[label] = y
        merged = sum(grouped.values())
        print(
            q,
            par,
            np.linalg.norm(gram - np.eye(4)),
            "none",
            kyfan(merged) - kyfan(stock),
            "79",
            sum(kyfan(v) for v in grouped.values()) - kyfan(stock),
            flush=True,
        )


if __name__ == "__main__":
    run()
