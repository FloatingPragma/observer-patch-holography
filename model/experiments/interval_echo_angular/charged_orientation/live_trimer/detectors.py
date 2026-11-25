"""Native body rotations, local process readback and simultaneous count errors."""

from functools import cache

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate

IDENTITY = np.eye(2)
Z = np.diag([1, -1])
X = np.array([[0, 1], [1, 0]])
Y = np.array([[0, -1j], [1j, 0]])
PAULI = (IDENTITY, Z, X, Y)


def written_source():
    """Native echo and a coin herald leave a sharp ending pair and an unread contrast sign."""
    law = candidate.native_law
    outcomes = {}
    for sign in (-1, 1):
        sharp = np.zeros(law.DIM)
        sharp[law.SRC.index((1, 0))] = 1 / np.sqrt(2)
        sharp[law.SRC.index((0, 1))] = sign / np.sqrt(2)
        after = law.echo() @ sharp
        admitted = np.zeros_like(after)
        offset = 9 + 9 * law.LABELS.index((0, 1)) + law.SRC.index((0, 1))
        reverse = 9 + 9 * law.LABELS.index((0, -1)) + law.SRC.index((0, 1))
        admitted[[offset, reverse]] = after[[offset, reverse]]
        assert abs(admitted[offset] - sign * admitted[reverse]) < 1e-14
        prob = admitted @ admitted
        np.testing.assert_allclose(prob, 0.5)
        outcomes[str(sign)] = float(prob)
    return outcomes


@cache
def spin_readers():
    """Sharp, phase and quadrature bases on the blank |20>, |02> seam."""
    law = candidate.native_law
    idx = [law.SRC.index(pair) for pair in [(2, 0), (0, 2)]]
    ex, cost = law.local_term_parts()
    h = (ex + candidate.KAPPA * cost)[np.ix_(idx, idx)]
    pulse = expm(-1j * np.pi * h / 2)
    private = np.diag(
        [np.exp(-1j * 4 * candidate.KAPPA * 3 * np.pi / (8 * candidate.KAPPA)), 1]
    )
    # With the negative i-sine column of this physical pulse the first
    # projector needs the positive or negative Pauli eigenvalue by axis.
    sharp = [np.diag([1, 0]), np.diag([0, 1])]
    result = []
    for u, axis in [(IDENTITY, Z), (private @ pulse, X), (pulse, Y)]:
        got = [u @ p @ u.conj().T for p in sharp]
        if np.trace(got[0] @ axis).real < 0:
            got.reverse()
        np.testing.assert_allclose(
            got, [(IDENTITY + axis) / 2, (IDENTITY - axis) / 2], atol=1e-13
        )
        result.extend(got)
    return result


def measured_partial(counts, shots):
    """Partial transpose of the normalized-input, unnormalized-output Choi state."""
    outs = []
    for triple in np.array(counts) / shots:
        r0 = triple[0, :2].sum()
        r = triple[:, 0] - triple[:, 1]
        outs.append(
            (r0 * IDENTITY + sum(q * p for q, p in zip(r, PAULI[1:], strict=True))) / 2
        )
    mats = [outs[0] + outs[1], outs[0] - outs[1], outs[2] - outs[3], outs[4] - outs[5]]
    pt = sum(np.kron(p.T, m.T) for p, m in zip(PAULI, mats, strict=True)) / 4
    return np.linalg.eigvalsh(pt)


def simultaneous_bound(shots, appointments):
    # Each input uses one [0,1] success read and three [-1,1] Pauli reads.
    delta = np.sqrt(np.log(2 * appointments * 6 * 4 / 0.01) / (2 * shots))
    return (1 + 2 * np.sqrt(3)) * delta


@cache
def body_turn():
    """Turn the ordered ternary triangle load pattern by two native pi contacts."""
    law = candidate.native_law
    ex, cost = law.local_term_parts()
    u = expm(-1j * np.pi * (ex + candidate.KAPPA * cost))[:9, :9]
    expected, _ = law._load_ops()
    np.testing.assert_allclose(u, expected, atol=2e-13)
    permutation = []
    for i in range(3):
        state = np.eye(3)[:, i]
        for a, b in [(0, 2), (2, 1)]:
            state[[a, b]] = state[[b, a]]
        permutation.append(int(np.argmax(state)))
    turn = np.eye(3)[:, permutation]
    np.testing.assert_array_equal(np.linalg.matrix_power(turn, 3), np.eye(3))
    return turn
