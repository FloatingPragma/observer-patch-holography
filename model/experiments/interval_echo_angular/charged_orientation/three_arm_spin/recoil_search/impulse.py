"""An addressed longitudinal phase gradient using a dark shutter and private port."""

import networkx as nx
import numpy as np
from scipy.linalg import expm

from ..apparatus import LAW, fixture, winding


def certify(phase):
    f = fixture()
    board = winding.flight.federation(0)
    stock = set(f["donor"]) | set(f["ports"])
    i, r, s = next(
        (i, int(r), int(s))
        for i in f["ports"]
        for r in board.neighbours[i]
        if r not in stock and board.degree[r] == 5
        for s in board.neighbours[r]
        if s not in stock and board.degree[s] == 5
    )
    # A separate balanced apparatus stock lives on two isolated blank dimers.
    # Equal degrees give the two orientations of either stored defect the
    # same private cost. Their symmetric blank-dimer state is stationary and
    # cancels against the identically prepared clock host. A private quarter
    # turn and a blank half contact localize the vacancy at s before the kick.
    u, v = next(
        (int(u), int(v))
        for u, v in board.seams
        if not {u, v} & (stock | {r, s}) and board.degree[u] == board.degree[v] == 5
    )
    for ends in ((r, s), (u, v)):
        for a, b in board.seams:
            if {a, b} & set(ends) and {a, b} != set(ends):
                assert set(ends) - {a, b}  # park a defect off a (1,1) shutter write
    relay = nx.single_source_shortest_path(nx.Graph(f["blank"]), i)
    assert set(relay) == set(f["ports"])
    ex, con = LAW.local_term_parts()
    h = ex + 8 * con
    full = expm(-1j * np.pi * h)
    np.testing.assert_allclose(full, full.conj().T, atol=2e-13)
    blank = full[:9, :9]
    swap, _ = LAW._load_ops()
    np.testing.assert_allclose(blank, swap, atol=2e-13)
    rails = [LAW.SRC.index(pair) for pair in ((0, 1), (1, 0))]
    half = expm(-0.5j * np.pi * h[np.ix_(rails, rails)])
    localized = half @ (np.array([1, 1j]) / np.sqrt(2))
    np.testing.assert_allclose(abs(localized), [0, 1], atol=2e-13)
    echo = LAW.echo()
    dark = echo[:, LAW.SRC.index((1, 1))].reshape(11, 9)[:, LAW.SRC.index((1, 1))]
    # Axes are shutter head, guide load, reference load, staging load, and
    # input column. A guide vacancy and a filled guide site span all uses.
    initial = np.zeros((11, 3, 3, 3, 2), complex)
    for col in (0, 1):
        initial[:, col, 1, 0, col] = dark

    def shutter(op, z):
        return (op @ z.reshape(99, -1)).reshape(initial.shape)

    def stage(z):
        return np.einsum("uv,hivc->hiuc", blank, z.reshape(11, 3, 9, 2)).reshape(
            initial.shape
        )

    z = shutter(full, stage(shutter(echo, initial)))
    rate = 8 * np.array(
        [con[LAW.SRC.index((n, 0)), LAW.SRC.index((n, 0))] for n in range(3)]
    )
    np.testing.assert_allclose(rate, 8 * np.arange(3) ** 2, atol=1e-13)
    z *= np.exp(-1j * (phase % (2 * np.pi)) / 8 * rate)[None, None, :, None, None]
    out = shutter(echo, stage(shutter(full.conj().T, z)))
    expected = initial * np.array([1, np.exp(-1j * phase)])
    error = np.max(abs(out - expected))
    assert error < 1e-12
    return dict(
        guide=i,
        private_port=r,
        staging_port=s,
        compensating_dimer=[u, v],
        relay_reaches=len(relay),
        longest_route=max(len(path) - 1 for path in relay.values()),
        two_column_residual=float(error),
    )
