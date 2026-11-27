"""Compile paired-dimer level reads and measure ordinary formation-energy fringes."""

from collections import deque
import numpy as np
from scipy.linalg import eigh, expm

from .measure import SHOTS

X = np.array([[0.0, 1.0], [1.0, 0.0]])
Z = np.diag([1.0, -1.0])
IDENTITY = np.eye(2)


def rotate(angle):
    # Native blank exchange supplies X. A private load reference or fixed
    # exterior pin supplies Z. Conjugating X by a quarter Z turn gives Y.
    quarter = expm(-1j * np.pi * Z / 4)
    return np.real_if_close(quarter @ expm(-1j * angle * X) @ quarter.conj().T)


def cnot(cost):
    # One addressed forbidden-word collar implements CZ up to a common phase.
    phase = -np.diag(np.exp(-1j * np.pi / cost * np.array([cost, cost, cost, 0.0])))
    had = rotate(np.pi / 4) @ Z
    return np.real_if_close(np.kron(IDENTITY, had) @ phase @ np.kron(IDENTITY, had))


def conditional_rotation(i, j, angle, cost):
    """Map two computational words onto 10,11 using reversible native gates."""
    cn = cnot(cost)
    had = rotate(np.pi / 4) @ Z
    reverse = np.kron(had, had) @ cn @ np.kron(had, had)
    gates = (np.kron(X, IDENTITY), np.kron(IDENTITY, X), cn, reverse)
    initial = (i, j)
    queue = deque([(initial, np.eye(4))])
    seen = {initial}
    while queue:
        pair, perm = queue.popleft()
        if pair == (2, 3):
            break
        for gate in gates:
            new = tuple(np.argmax(abs(gate[:, k])) for k in pair)
            if new not in seen:
                seen.add(new)
                queue.append((new, gate @ perm))
    else:
        raise ValueError(pair)
    half = np.kron(IDENTITY, rotate(angle / 2))
    return perm.T @ half @ cn @ half.T @ cn @ perm


def circuit(target, cost):
    tri = target.copy()
    pulse = np.eye(4)
    for col in range(3):
        for row in range(3, col, -1):
            angle = np.arctan2(-tri[row, col], tri[row - 1, col])
            gate = conditional_rotation(row - 1, row, angle, cost)
            tri = gate @ tri
            pulse = gate @ pulse
    # Each negative terminal sign is a CZ conjugated by two optional rail swaps.
    signs = np.eye(4)
    cn = cnot(cost)
    had = np.kron(IDENTITY, rotate(np.pi / 4) @ Z)
    cz = had @ cn @ had
    for word in np.flatnonzero(np.diag(tri) < 0):
        perm = np.kron(IDENTITY if word & 2 else X, IDENTITY if word & 1 else X)
        signs = signs @ perm @ cz @ perm
    implemented = pulse.T @ signs
    np.testing.assert_allclose(implemented, target, atol=5e-13)
    return implemented


def clocks(h, cost, rng):
    e, v = eigh(h)
    basis = circuit(v, cost)
    res = []
    radius = np.sqrt(4 * np.log(2 * 100 / 0.01) / SHOTS)
    for k in range(1, 4):
        seq = []
        old = 0.0
        for t in 4.0 ** np.arange(7) / 64:
            initial = (basis[:, 0] + basis[:, k]) / np.sqrt(2)
            out = basis.T @ (v @ (np.exp(-1j * e * t) * (v.T @ initial)))
            bits = []
            for ph in (1, 1j):
                p = abs((out[0] + ph * out[k]) / np.sqrt(2)) ** 2
                bits.append(int(rng.binomial(SHOTS, p)))
            read = 2 * np.array(bits) / SHOTS - 1
            angle = np.arctan2(read[1], read[0])
            angle += 2 * np.pi * round((4 * old - angle) / (2 * np.pi))
            old = angle
            seq.append(dict(hold=float(t), counts=bits))
        gap = old / t
        assert abs(gap - (e[k] - e[0])) <= np.arcsin(radius) / t
        res.append(
            dict(
                line=k,
                gap=round(gap, 9),
                radius=float(np.arcsin(radius) / t),
                records=seq,
            )
        )
    return res


def charge(board, dimers):
    legs = dimers[2:]
    reads = []
    for ports in legs:
        assert all(board.degree[p] == 5 for p in ports)
        blocks = []
        j = np.arange(8)
        for a in range(3):
            for b in range(3):
                # Native private references, each with its zero load kept fixed.
                # In the binary code n squared = n. Subtract a unit standard.
                phases = np.exp(-2j * np.pi * j * (a * a + b * b - 1) / 8)
                blocks.append(np.fft.ifft(phases))
        np.testing.assert_allclose(blocks[1], blocks[3])
        for word, target in ((0, 7), (4, 1), (1, 0), (3, 0)):
            assert abs(blocks[word][target]) > 1 - 1e-12
        reads.append(
            dict(
                ports=ports,
                word_00_to_22=[int(np.argmax(abs(p))) for p in blocks],
                coherent_amplitude_error=round(
                    float(np.linalg.norm(blocks[1] - blocks[3])), 13
                ),
            )
        )
    return reads
