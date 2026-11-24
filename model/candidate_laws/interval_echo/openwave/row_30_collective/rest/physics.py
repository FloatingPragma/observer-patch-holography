"""Periodic wall boundary and physical-dimer check of the spectral meter."""

import numpy as np
from scipy.linalg import eigh, svd


def propagation(n, coupling, flip):
    # Dual parity is fixed by the product of the n measured wall signs.
    # The physical dimer-flip sector fixes the sign of the closing dual bond.
    p = (-1) ** n
    forward = np.eye(n, k=-1)
    forward[0, -1] = -flip * p
    u, frequencies, _ = svd(2 * coupling * np.eye(n) - forward, full_matrices=False)
    return u, frequencies


def full_spin_check(n, flip, coupling):
    words = np.arange(1 << n)
    z = 1 - 2 * ((words[:, None] >> np.arange(n)) & 1)
    walls = -z * np.roll(z, 1, axis=1)
    h = np.diag((-coupling * walls.sum(axis=1)).astype(float))
    for a in range(n):
        h[words, words ^ (1 << a)] = -0.5
    # A physical translation-invariant spin-flip constraint selects a real
    # mixed state. Each orbit has one member with its most significant bit 0.
    half = len(words) // 2
    sector = (
        np.eye(len(words))[:, :half] + flip * np.eye(len(words))[::-1, :half]
    ) / np.sqrt(2)
    e, q = eigh(sector.T @ h @ sector)
    u, sigma = propagation(n, coupling, flip)
    differences = []

    def string(left, right):
        matrix = np.zeros(h.shape, complex)
        mask = (1 << right) - (1 << left)
        matrix[words ^ mask, words] = 1j * np.prod(walls[:, left:right], axis=1)
        np.testing.assert_allclose(matrix, matrix.conj().T)
        return q.T @ sector.T @ matrix @ sector @ q

    for left, right, dest, end in [
        (0, 1, 2, 3),
        (0, 2, 1, n - 1),
        (1, n - 1, 1, n - 1),
    ]:
        first, second = string(left, right), string(dest, end)
        for t in (0.37, 2.15, 5.81):
            g = np.exp(-1j * t * e)
            observed = (
                np.sum(first.T * (g.conj()[:, None] * second * g[None, :])).real / half
            )
            correlations = (u * np.cos(sigma * t)) @ u.T
            predicted = np.linalg.det(correlations[np.ix_([left, right], [dest, end])])
            differences.append(abs(float(observed - predicted)))
    assert max(differences) < 1e-11, (n, flip, max(differences))
    return max(differences)
