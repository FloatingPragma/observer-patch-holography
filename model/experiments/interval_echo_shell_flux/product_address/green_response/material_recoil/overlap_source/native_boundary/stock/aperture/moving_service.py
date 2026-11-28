"""Stationary record laws of the full monitored native exchange on returned material."""

from functools import lru_cache
from fractions import Fraction as Q

import numpy as np
from scipy.linalg import expm
from scipy.sparse.csgraph import connected_components

from . import closed_energy as energy

V = 240


@lru_cache(maxsize=1)
def gates():
    _, _, h, h0 = energy.pair.stock.kernels()
    eye = np.eye(99, dtype=int)
    # At both dwells the commuting stiffness has rephased. h0/4 is a projector.
    c = h - h0
    np.testing.assert_array_equal(c @ h0, h0 @ c)
    np.testing.assert_array_equal(c @ (c - 32 * eye) @ (c - 128 * eye), 0)
    np.testing.assert_array_equal(h0 @ h0, 4 * h0)
    # The material composition is conserved even before the sharp checks:
    # no off-diagonal term of the native Hamiltonian can change it.
    dest, src = np.nonzero(h)
    assert all(
        sorted(divmod(i % 9, 3)) == sorted(divmod(j % 9, 3))
        for i, j in zip(dest, src, strict=True)
    )
    full = (2 * eye - h0) ** 2
    half = 4 * eye + full
    for t, tickets in [(np.pi, full), (np.pi / 2, half)]:
        np.testing.assert_allclose(
            abs(expm(-1j * t * h / 4)) ** 2, tickets / tickets.sum(0), atol=1e-13
        )
        np.testing.assert_array_equal(tickets, tickets.T)
        dest, src = np.nonzero(tickets)
        for i, j in zip(dest, src, strict=True):
            assert sorted(divmod(i % 9, 3)) == sorted(divmod(j % 9, 3))
    return full, half


@lru_cache(maxsize=None)
def stationary_heads(pop):
    _, half = gates()
    counts = np.array(pop)
    draws = counts[:, None] * (counts[None, :] - np.eye(3, dtype=int))
    feasible = draws.ravel() > 0
    k = half.reshape(11, 9, 11, 9)
    adj = k[:, :, :, feasible].sum(axis=(1, 3))
    np.testing.assert_array_equal(adj > 0, adj.T > 0)
    _, label = connected_components(adj, directed=False)
    return [tuple(map(int, np.flatnonzero(label == x))) for x in label]


def signature(word, extra=0, single=0):
    pop = np.bincount(word[:4], minlength=3)
    pop[0] += V - 4 - extra - single
    pop[2] += extra
    pop[1] += single
    key = tuple(map(int, pop))
    comp = stationary_heads(key)
    return key, tuple(comp[int(m)] for m in word[4:])


@lru_cache(maxsize=None)
def pre(pop, labels):
    # A specified ordered pair is a draw from the material without replacement.
    # The stationary head is uniform in its accessible component at fixed pop.
    dist = np.zeros((11, 9), int)
    n = np.array(pop)
    dist[list(labels)] = (n[:, None] * (n[None, :] - np.eye(3, dtype=int))).ravel()
    assert dist.sum() == len(labels) * V * (V - 1)
    return dist.ravel()


def before(rows):
    ans = np.full((4, 99), Q(0), object)
    for (pop, marks), prob in rows.items():
        for j, labels in enumerate(((0,), *marks)):
            counts = pre(pop, labels)
            ans[j] += (prob / int(counts.sum())) * counts
    assert all(ans.sum(1) == 1)
    return ans


@lru_cache(maxsize=1)
def readers():
    full, half = gates()

    # Columns hold the probability mass (in tickets) of each surprise bin.
    def binned(matrix):
        return np.array([k * (matrix == k).sum(0) for k in range(1, 9)], int)

    all_bins = binned(half)
    material_bins = binned(half.reshape(11, 9, 99).sum(0))
    assert np.all(all_bins.sum(0) == 8) and np.all(material_bins.sum(0) == 8)
    marginal = full.reshape(11, 3, 3, 99).sum(0)
    tv = np.zeros((2, 99), int)
    for code in range(99):
        m, ab = divmod(code, 9)
        a, b = divmod(ab, 3)
        for t in (1, 2):
            # Given the preparatory target load and head, decode a bit from
            # the final target load. Average the two boundary orientations.
            left = marginal[:, :, 9 * m + b].sum(0) - marginal[
                :, :, 9 * m + 3 * t + b
            ].sum(0)
            right = marginal[:, :, 9 * m + 3 * a].sum(1) - marginal[
                :, :, 9 * m + 3 * a + t
            ].sum(1)
            tv[t - 1, code] = abs(left).sum() + abs(right).sum()
    assert np.all(tv[:, :9] == 16)
    # A third receiver keeps the entire recorded boundary word and checks
    # for its return after a full pulse. Each entry uses denominator 16.
    return all_bins, material_bins, np.vstack((tv, 4 * full.diagonal()))


def log_coefficients(bins):
    # For each ticket probability k/8, express -log2(k/8) on logs of primes.
    # Rational combinations of these logarithms vanish exactly when every
    # coefficient vanishes, by clearing denominators and unique factorization.
    return bins @ np.array(
        [
            [3, 0, 0, 0],
            [2, 0, 0, 0],
            [3, -1, 0, 0],
            [1, 0, 0, 0],
            [3, 0, -1, 0],
            [2, -1, 0, 0],
            [3, 0, 0, -1],
            [0, 0, 0, 0],
        ]
    )


def bits(bins):
    return np.asarray(bins, float) @ np.log2(8 / np.arange(1, 9))


def word_sweep(pop, labels):
    # Return the exact ticket-bin law. Bits are converted separately.
    count = pre(pop, labels)
    return readers()[0] @ count, 8 * int(count.sum())


@lru_cache(maxsize=None)
def rate_interval(pop, labels):
    from mpmath import iv

    nums, den = word_sweep(pop, labels)
    return sum(
        iv.mpf(int(v)) / den * iv.log(iv.mpf(8) / k) / iv.log(2)
        for k, v in enumerate(nums, 1)
        if v and k != 8
    )


def support_bounds(classes, blank_edges):
    from mpmath import iv

    iv.dps = 35
    out = {}
    for k, words in enumerate(classes):
        cases = {signature(w) for w in words}
        bands = []
        for pop, marks in cases:
            x = blank_edges * rate_interval(pop, (0,)) + sum(
                rate_interval(pop, m) for m in marks
            )
            x = iv.mpf(x)
            bands.append(
                [
                    float(np.nextafter(float(x.a), -np.inf)),
                    float(np.nextafter(float(x.b), np.inf)),
                ]
            )
        if k == 0:
            assert all(max(map(abs, b)) < 1e-50 for b in bands)
            bands = [[0.0, 0.0]]
        out[str(k)] = dict(
            words=len(words),
            stationary_classes=len(cases),
            minimum_lower=min(b[0] for b in bands),
            maximum_upper=max(b[1] for b in bands),
        )
    return out
