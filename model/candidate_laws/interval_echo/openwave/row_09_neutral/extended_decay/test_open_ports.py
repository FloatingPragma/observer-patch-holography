"""Check the joined generator, finite cold preparation and frozen timing discriminator."""

import json
from itertools import product
import numpy as np
from scipy.linalg import expm
from scipy import sparse
from . import cold_ports, early, joined_outlet as joined, open_ports, outlet, wing
from .material import route, edge
from .measure import HERE
from .outlet_clock import LABELS


def record(name):
    return json.loads((HERE / (name + ".json")).read_text())


def test_every_physical_word_and_contact():
    b, d, cell, available, sites, seams, cut = open_ports.aperture()
    assert available == [73, 75, 76, 77, 80, 81, 83] and len(cell) == 12
    heads = route((3,), crossing=True)[2]
    pins = {p for pair in heads for p in pair} - set(sum(d, ()))
    dg = wing.installation()[2]
    for name in joined.SETTINGS:
        blank = [e for e in seams if name != "exit_cut" or e not in cut]
        words, tail = outlet.connected_wing(b, d, sites, blank)
        ref = dg[
            {"source_cut": "boundary_cut", "no_reversal": "no_reversal"}.get(
                name, "native"
            )
        ]
        h = joined.generator(ref, tail, words, d, name)
        deposits = dict(heads)
        if name == "source_cut":
            deposits[edge(d[2][1], d[3][0])] = 32
        if name == "no_reversal":
            deposits = dict.fromkeys(deposits, 32)
        states = [
            frozenset(
                {*mobile, *(p[1 - bit] for bit, p in zip(bits, d[:4], strict=True))}
            )
            for bits in product((0, 1), repeat=4)
            for mobile in words
        ]
        ix = {state: i for i, state in enumerate(states)}
        diagonals = []
        coords = []
        cols = []
        for i, state in enumerate(states):
            active = state | pins
            diagonal = sum(8 * (6 - b.degree[p]) for p in state)
            diagonal += sum(
                cost * (1 - int(a in active and c in active))
                for (a, c), cost in deposits.items()
            )
            for a, c in (*d[:4], *blank):
                if (a in state) != (c in state):
                    diagonal += 8.5
                    coords.append(ix[state ^ {a, c}])
                    cols.append(i)
            diagonals.append(diagonal)
        difference = np.asarray(diagonals) - np.diag(h)
        np.testing.assert_allclose(difference, difference[0], atol=2e-13, rtol=0)
        measured = sparse.csr_matrix(h - np.diag(np.diag(h)))
        truth = sparse.csr_matrix(
            (np.full(len(coords), -0.5), (coords, cols)), shape=h.shape
        )
        assert abs((measured - truth).data).max(initial=0) < 3e-14


def test_spectral_preparation_against_clock_doubling():
    # A small nondegenerate Hamiltonian checks the full complex zero-bin filter,
    # the embeddings, and the error on the original-admission normalization.
    rng = np.random.default_rng(8927)
    a = rng.normal(size=(24, 24))
    h = (a + a.T) / 2
    state = rng.normal(size=24)
    state /= np.linalg.norm(state)
    k = np.array([[0.6, 0, -0.2j], [0, 0, 0], [0.2j, 0, 0.2]])
    _, _, weights, r = joined.prepare(h, state, k)
    density = np.outer(state, state)
    u = expm(-2j * np.pi * h / outlet.PERIOD)
    for _ in range(int(np.log2(LABELS))):
        density = (density + u @ density @ u.conj().T) / 2
        u = u @ u
    full = np.kron(np.eye(8), k)
    direct = full @ density @ full.conj().T
    np.testing.assert_allclose(
        weights @ weights.conj().T, direct, atol=r["discarded_weight_upper"] + 2e-12
    )
    np.testing.assert_allclose(np.trace(direct).real, r["admission"], atol=2e-12)


def test_timing_enrollment_and_joined_controls():
    assert not early.verdict(record("early"))
    assert early.verdict(record("ticks"), require_each=True)
    assert not open_ports.verdict(record("open_ports"))
    assert not cold_ports.verdict(record("cold_ports"))
    measured = record("joined_outlet")
    joined.decode(measured)
    rows = measured["settings"]["native"]["arms"]
    pulsed, idle, known = (
        np.array(rows[n]["core_by_outlet"]) for n in ("core_pulse", "idle", "D_phase")
    )
    assert min((known - idle)[1:4, :, 2].sum(1)) > 0.014
    assert max((pulsed - idle)[4:, :, 2].sum(1)) < -0.02
    assert min((pulsed - idle)[4:, :, 4].sum(1)) > 0.0139
    for s in ("source_cut", "exit_cut", "no_reversal"):
        settings = measured["settings"][s]["arms"]
        a, b = (np.array(settings[n]["core_by_outlet"]) for n in ("core_pulse", "idle"))
        np.testing.assert_allclose((a - b).sum(1), 0, atol=4e-11, rtol=0)


def test_remote_soft_read_by_sparse_direct_propagation():
    from scipy.sparse.linalg import expm_multiply

    board, dimers, cell, _, sites, contacts, _ = open_ports.aperture()
    sectors = outlet.receiving_cell(
        board, cell, [e for e in contacts if set(e) <= set(cell)]
    )
    words, tail = outlet.connected_wing(board, dimers, sites, contacts)
    full = wing.installation()[2]["native"]
    h = joined.generator(full, tail, words, dimers, "native")
    _, g = np.linalg.eigh(full)
    stock = outlet.inject(
        g[:, 0].reshape(4, 2, 2, 2), words, dimers, cell, sectors[1][2][:, :1]
    ).T.ravel()
    readers = outlet.read_groups(words, sectors, cell)
    k = cold_ports.zero_bin(readers, sectors)
    _, _, columns, prep = joined.prepare(h, stock, k)
    weights = (abs(columns) ** 2).sum(0)
    # Omitting this tiny portion of the positive density can change any
    # original-launch event probability by at most the omitted total weight.
    start = np.searchsorted(np.cumsum(weights), 0.0005)
    error = weights[:start].sum() + prep["discarded_weight_upper"]
    kept = columns[:, start:]
    core, _, _, _ = outlet.hamiltonians(4)
    energies, v = np.linalg.eigh(core)
    pulsed = joined.core_pulse(v) @ kept.reshape(4, -1)
    packed = np.hstack((kept, pulsed.reshape(kept.shape)))
    released = expm_multiply((-1j * joined.WAITS[5]) * sparse.csr_matrix(h), packed)
    # Independent soft-energy projector weights from Fourier bins on the exact
    # one-load eigenlines. A trace over the other physical ports gives each yield.
    _, one, rv = sectors[1]
    response = outlet.clock(
        one,
        one[0],
        np.ptp(energies) / 10,
        labels=LABELS,
        low_edge=outlet.PERIOD / LABELS,
    )
    soft = (rv * response[:, 1]) @ rv.T
    archive = record("joined_outlet")["settings"]["native"]["arms"]
    for name, state in zip(
        ("idle", "core_pulse"), np.split(released, 2, axis=1), strict=True
    ):
        state = state.reshape(8, 2 * len(words), -1)
        probability = 0.0
        for number, rows, _ in readers:
            if number == 1:
                block = state[:, rows].transpose(1, 0, 2).reshape(12, -1)
                probability += np.vdot(block, soft @ block).real
        target = np.array(archive[name]["core_by_outlet"])[5, :, 2].sum()
        assert -1e-10 <= target - probability <= error + 1e-10
