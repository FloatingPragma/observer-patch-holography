"""Native open strand, a single-port Lüders arrival and two terminal ports."""

from functools import cache

import numpy as np
from scipy.linalg import block_diag, svd

from ...row_30_collective import experiment, packets


def source(n, gradient=0, mirror=False):
    x = np.arange(n + 1)
    envelope = np.zeros(n + 1)
    envelope[1:16] = np.sin(np.pi * x[1:16] / 16) / np.sqrt(8)
    if mirror:
        envelope = envelope[::-1]
    target = envelope * np.exp(1j * gradient * x)
    # Stop at the last nonzero amplitude. Unused rotations on trailing zero
    # amplitudes of the older all-positive packet were never needed.
    a = np.zeros(n + 1, complex)
    a[0] = 1
    last = np.flatnonzero(envelope)[-1]
    for site in range(last):
        c = abs(target[site]) / np.linalg.norm(target[site:])
        s = np.sqrt(max(0, 1 - c * c))
        packets.verify_conditioned_flip(site, n, np.arctan2(s, c))
        a[site], a[site + 1] = c * a[site], 1j * s * a[site]
    a *= np.exp(1j * (np.angle(target) - np.angle(a)))
    np.testing.assert_allclose(a, target, atol=3e-14)
    return a


def bare_covariance(f):
    coherence = f.conj()[:, None] * f
    aa = -2 * coherence.imag
    ab = 2 * coherence.real - np.eye(len(f))
    return np.block([[aa, ab], [-ab.T, aa]])


class Strand:
    def __init__(self, n, j, h=experiment.H):
        self.n = n
        if h == experiment.H:
            self.u, self.e, self.v = experiment.dual(n, j)
        else:
            m = n + 1
            c = 2 * (j * np.eye(m) - h * np.eye(m, k=-1))
            self.u, self.e, vh = svd(c)
            self.v = vh.T
        self.basis = block_diag(self.u, self.v)

    def rows(self, bonds, times):
        """Rows interleave the two physical Majoranas of each requested bond."""
        angle = np.asarray(times)[:, None, None] * self.e
        c, s = np.cos(angle), np.sin(angle)
        u, v = self.u[bonds], self.v[bonds]
        a = np.concatenate((u * c, u * s), axis=2)
        b = np.concatenate((-v * s, v * c), axis=2)
        return np.stack((a, b), axis=2).reshape(len(times), 2 * len(bonds), -1)

    def enroll(self, f, age=2):
        initial = self.basis.T @ bare_covariance(f) @ self.basis
        rows = self.rows(np.arange(len(f)), [age])[0]
        order = np.r_[np.arange(0, len(rows), 2), np.arange(1, len(rows), 2)]
        rows = rows[order]
        return rows @ initial @ rows.T


def condition(g, sign):
    """Project one endpoint dimer orientation, retaining all unmeasured modes."""
    a, b = 0, len(g) // 2
    weight = (1 + sign * g[a, b]) / 2
    if weight <= 0:
        raise ValueError("This trigger branch has zero probability")
    wedge = np.outer(g[:, b], g[:, a]) - np.outer(g[:, a], g[:, b])
    out = g + sign * wedge / (2 * weight)
    out[[a, b], :] = 0
    out[:, [a, b]] = 0
    out[a, b], out[b, a] = sign, -sign
    return weight, out


def pfaffian(matrices, indices):
    """Wick expansion of interleaved endpoint Majoranas over a time batch."""

    @cache
    def rec(idx):
        if not idx:
            return 1.0
        return sum(
            (-1) ** (j + 1) * matrices[:, idx[0], idx[j]] * rec(idx[1:j] + idx[j + 1 :])
            for j in range(1, len(idx))
        )

    return rec(tuple(indices))


def terminal(strand, g, times, depth=7):
    """All four probabilities of the joint first/depth dimer orientation read."""
    initial = strand.basis.T @ g @ strand.basis
    chunks = []
    for start in range(0, len(times), 128):
        rows = strand.rows(np.arange(depth), times[start : start + 128])
        cov = (rows @ initial) @ rows.transpose(0, 2, 1)
        first = cov[:, 0, 1]
        last = pfaffian(cov, range(2 * depth))
        product = pfaffian(cov, range(2, 2 * depth))
        chunks.append(
            np.column_stack(
                [
                    (1 + a * first + b * last + a * b * product) / 4
                    for a in (-1, 1)
                    for b in (-1, 1)
                ]
            )
        )
    return np.vstack(chunks)


def word_law(enrollment, flight, times, depth=7):
    branches = []
    for sign in (-1, 1):
        weight, g = condition(enrollment, sign)
        branches.append(weight * terminal(flight, g, times, depth))
    probabilities = np.column_stack(branches)
    if probabilities.min() < -1e-12 or abs(probabilities.sum(axis=1) - 1).max() > 1e-12:
        raise ValueError("Incomplete probability law")
    p = np.clip(probabilities, 0, 1)
    return p / p.sum(axis=1, keepdims=True)
