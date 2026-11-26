"""The open, mobile Ising string on native unequal-degree dimers."""

from functools import cache

import numpy as np
from scipy import sparse
from scipy.linalg import eigh
from scipy.special import jv
from model.executable_model import federation
from ...row_30_collective import experiment as strand

J, G, H = 8.0, 8.0, 0.5
LENGTHS = (2, 3, 4, 5, 6, 7, 8, 10, 12)
HOLDS = 2.0 ** np.arange(13) / 128
RELEASE = np.array([0, np.pi / 20, 0.5, 2, 8, 16])


def route(level):
    board = federation(level)
    mate = {a: b for e in board.external for a, b in (e, e[::-1])}

    def choices(owner, entry=None):
        gates = sorted(p for p in range(12 * owner, 12 * owner + 12) if p in mate)
        out = []
        for a in gates if entry is None else [entry]:
            for b in board.neighbours[a]:
                if board.degree[b] != 5:
                    continue
                for c in board.neighbours[b]:
                    if c == a or board.degree[c] != 5:
                        continue
                    for d in gates:
                        if d != a and d in board.neighbours[c]:
                            out.append((a, b, c, d))
        return sorted(out)

    def extend(blocks, owners):
        if len(blocks) == 8:
            return blocks
        entry = mate[blocks[-1][-1]]
        owner = entry // 12
        if owner in owners:
            return None
        for segment in choices(owner, entry):
            answer = extend(blocks + [segment], owners | {owner})
            if answer:
                return answer

    for owner in range(board.cells):
        for first in choices(owner):
            answer = extend([first], {owner})
            if answer:
                sites = [p for block in answer for p in block]
                assert all(b in board.neighbours[a] for a, b in zip(sites, sites[1:]))
                return sites
    raise ValueError("No prespecified eight-owner material route")


def calibration_energies():
    law = strand.candidate.native_law
    strand.native_checks()
    exchange, contrast = law.local_term_parts()
    h = exchange + strand.K * contrast
    eye = np.eye(law.DIM)

    def state(ab):
        return eye[:, law.SRC.index(ab)]

    uniform = law.echo() @ ((state((0, 2)) + state((2, 0))) / np.sqrt(2))
    collar = uniform.reshape(-1, 9)[:, law.SRC.index((1, 1))]

    def cost(ab):
        psi = np.kron(collar, np.eye(9)[:, law.SRC.index(ab)])
        return float(psi @ h @ psi)

    band = [law.SRC.index(ab) for ab in ((0, 1), (1, 0))]
    gap = np.diff(np.linalg.eigvalsh(h[np.ix_(band, band)]))[0]
    # The degree reference has contrast alone with its other port held at zero.
    load = state((1, 0))
    private = float(strand.K * load @ contrast @ load)
    return np.array([cost((0, 0)) - cost((1, 1)), gap, private])


def bits(n):
    return (np.arange(1 << n)[:, None] >> np.arange(n)) & 1


def generate(diagonal, n, h, sparse_output):
    idx = np.arange(len(diagonal))
    mat = sparse.diags(diagonal, format="csr")
    for i in range(n):
        mat += sparse.csr_matrix(
            (np.full(len(idx), -h), (idx, idx ^ (1 << i))), shape=mat.shape
        )
    return mat if sparse_output else mat.toarray()


def core(length, pin, g, h=H):
    n = length - 2
    inside = bits(n)
    extended = np.column_stack(
        (np.full(len(inside), pin), inside, np.full(len(inside), pin))
    )
    walls = (extended[:, 1:] != extended[:, :-1]).sum(axis=1)
    # The outer two fixed lower spins account for two extra domain boundaries.
    diagonal = 2 * J * (walls + 2 * pin) + g * (inside.sum(axis=1) + 2 * pin)
    return generate(diagonal, n, h, False)


@cache
def static_laws(length, g):
    n = length - 2
    e0, v0 = eigh(core(length, 0, g), driver="evd")
    e1, v1 = eigh(core(length, 1, g), driver="evd")
    w0 = v0[0] ** 2
    a0 = np.exp(-1j * HOLDS[:, None] * e0) @ w0
    panels = []
    for basis in ((1 << n) - 1, 0):
        w1 = v1[basis] ** 2
        a1 = np.exp(-1j * HOLDS[:, None] * e1) @ w1
        signal = a1 * a0.conj()
        weight = w0[e0 < e0[0] + 1e-10].sum() * w1[e1 < e1[0] + 1e-10].sum()
        panels.append((np.stack((signal.real, -signal.imag)), float(weight)))
    return float(e1[0] - e0[0]), panels


def release_generator(n, g, h):
    orientations = bits(n)
    walls = (orientations[:, 1:] != orientations[:, :-1]).sum(axis=1)
    # Both physical terminal connectors join degree-five ports. Their
    # unpaired linear fields favour the upper orientation at the two ends.
    diagonal = (
        2 * J * walls
        + g * orientations.sum(axis=1)
        - 2 * J * (orientations[:, 0] + orientations[:, -1])
    )
    return generate(diagonal, n, h, True)


def free_laws(g, h):
    n = 14
    mask = sum(1 << i for i in range(3, 11))
    launch = [mask, ((1 << n) - 1) ^ mask]
    seed = np.zeros((1 << n, 2))
    seed[launch, [0, 1]] = 1
    if h == 0:
        return np.broadcast_to(seed.T, (len(RELEASE), 2, len(seed)))
    full = release_generator(n, g, h)
    diag = full.diagonal()
    low, high = diag.min() - n * h, diag.max() + n * h
    shift = (high + low) / 2
    radius = (high - low) / 2
    scaled = (full - shift * sparse.eye(len(diag), format="csr")) / radius
    arg = RELEASE * radius
    order = int(np.ceil(arg.max() + 20 * np.cbrt(arg.max() + 1) + 50))
    k = np.arange(order + 1)
    coeff = 2 * (-1j) ** (k % 4) * jv(k[None, :], arg[:, None])
    coeff[:, 0] /= 2
    assert abs(coeff[:, -5:]).max() < 1e-20
    p, q = seed, scaled @ seed
    out = coeff[:, 0, None, None] * p + coeff[:, 1, None, None] * q
    for idx in range(2, len(k)):
        p, q = q, 2 * (scaled @ q) - p
        for j in np.flatnonzero(abs(coeff[:, idx]) > 1e-18):
            out[j] += coeff[j, idx] * q
    probabilities = (abs(out) ** 2).transpose(0, 2, 1)
    np.testing.assert_allclose(probabilities.sum(axis=2), 1, atol=2e-10, rtol=0)
    return probabilities


def event_masks():
    orientations = bits(14)
    # The actual signed endpoint jumps must agree with the orientation
    # of the launch. Additional opposite walls are tested inside the interval.
    masks = []
    for src, initial in enumerate((1, 0)):
        edge = (
            (orientations[:, 2] == 1 - initial)
            & (orientations[:, 3] == initial)
            & (orientations[:, 10] == initial)
            & (orientations[:, 11] == 1 - initial)
        )
        s = orientations[:, 3:11]
        internal = (s[:, 1:] != s[:, :-1]).sum(axis=1)
        island = np.zeros(len(s), bool)
        for j in range(1, 4):
            island |= (
                (s[:, j - 1] == initial)
                & (s[:, j + 4] == initial)
                & (s[:, j : j + 4] == 1 - initial).all(axis=1)
            )
        masks.append((edge, edge & (internal >= 2), edge & island))
    return masks
