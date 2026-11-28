"""Reconstruct the single-defect material statistics from simulated detector banks."""

from hashlib import sha256

import numpy as np
from scipy.linalg import expm
from scipy.special import entr

N = 2**52
# Fewer than 2e8 two-sided tests over both fixtures, material species and controls.
L = np.log(2 * 2e8 / 0.01)
HOEFFDING = np.sqrt(L / (2 * N))


def verify_receiver(contact):
    """Obtain both analysis phases from a blank splitter and a known written clock."""
    blank = contact.stock.kernels()[2][:9, :9]
    base, extra = contact.blocks()["full"][2]
    head = np.ones(len(base)) / np.sqrt(len(base))
    np.testing.assert_allclose(extra @ head / 4, -32 * head)
    dwell = np.pi / 64
    free = expm(-1j * dwell * base / 4) @ head
    occupied = expm(-1j * dwell * (base + extra) / 4) @ head
    np.testing.assert_allclose(occupied, 1j * free, atol=2e-14)
    for load in (0, 2):
        words = [3 * load + 1, 3 + load]
        split = expm(-1j * np.pi * blank[np.ix_(words, words)] / 8)
        for phase, quadrature in ((1, 1j), (np.vdot(free, occupied), 1)):
            # Standard at the second port. It sees two unit loads when the
            # defect is at the first port. Both paths return the same clock.
            arm = split @ np.diag([phase, 1])
            effect = np.outer(arm[0].conj(), arm[0])
            expected = np.array([[1, quadrature], [np.conj(quadrature), 1]]) / 2
            np.testing.assert_allclose(effect, expected, atol=2e-14)


def frequency_error(observed):
    upper = np.minimum(1, observed + HOEFFDING)
    return np.sqrt(2 * upper * L / N) + 2 * L / (3 * N)


def acquire(rho, rng, homogeneous_outside=False):
    n = len(rho)
    counts = rng.multinomial(N, np.diag(rho))
    diag = counts / N
    rec = np.diag(diag.astype(complex))
    i, j = np.triu_indices(n - int(homogeneous_outside), 1)
    pair = rho[i, i] + rho[j, j]
    digest = sha256(counts.astype("<i8").tobytes())
    for phase, unit in ((rho[i, j].real, 1), (rho[i, j].imag, 1j)):
        # The idle outcome is N - hit. The two output ports split the hits.
        hit = rng.binomial(N, pair)
        plus = rng.binomial(hit, 0.5 + phase / pair)
        rec[i, j] += unit * ((2 * plus - hit) / (2 * N))
        for data in (hit, plus):
            digest.update(data.astype("<i8").tobytes())
    rec[j, i] = rec[i, j].conj()
    # For each quadrature the single-trial outcome is +1/2, -1/2, or zero.
    upper_pair = diag[i] + diag[j] + 2 * HOEFFDING
    pair_error = np.sqrt(upper_pair * L / (2 * N)) + 2 * L / (3 * N)
    diag_error = frequency_error(diag)
    frob = np.sqrt(np.sum(diag_error**2) + 4 * np.sum(pair_error**2))
    distance = np.sqrt(n) * frob / 2
    minimum = np.linalg.eigvalsh(rec)[0]
    assert minimum > 0 and distance < 0.1
    np.testing.assert_allclose(rec.trace(), 1, atol=2e-14, rtol=0)
    # A regional restriction is trace preserving and cannot raise trace distance.
    # Its dimension is at most n, including the single homogeneous outcome.
    entropy_error = distance * np.log(n - 1) + entr(distance) + entr(1 - distance)
    errors = np.diag(diag_error)
    errors[i, j] = errors[j, i] = pair_error
    return (
        rec,
        errors,
        dict(
            entropy_error=float(entropy_error),
            trace_distance_bound=float(distance),
            smallest_reconstructed_eigenvalue=float(minimum),
            sha256_counts=digest.hexdigest(),
        ),
    )


def entropy(block, full):
    outside = 0.0 if full else 1 - block.trace().real
    eig = np.linalg.eigvalsh(block)
    assert eig.min() > 0 and outside >= 0
    return np.array(
        [
            sum(entr(eig)) + entr(outside),
            sum(entr(block.diagonal().real)) + entr(outside),
        ]
    )


def product_interval(x, y):
    products = np.outer(x, y)
    return np.array([products.min(), products.max()])


def difference_interval(x, y):
    return np.array([x[0] - y[1], x[1] - y[0]])


def reference_scores(reference, errors, area, membership, material, uncertainty):
    sums = np.bincount(membership, reference)[membership]
    sum_error = np.bincount(membership, errors)[membership]
    assert np.all(reference > errors) and np.all(sums > sum_error)
    surprise = np.log(sums / reference)
    log_error = -np.log1p(-errors / reference) - np.log1p(-sum_error / sums)
    scores = np.column_stack((area, material + surprise[:, None]))
    # Subtracting a fixed offset does not alter the slope of normalized laws.
    scores -= scores.min(axis=0)
    radii = np.column_stack((np.zeros(len(area)), *(2 * [uncertainty + log_error])))
    return scores, radii


def shift_interval(shift, err, scores, radii):
    middle = shift @ scores
    width = err @ abs(scores) + (abs(shift) + err) @ radii
    return np.stack((middle - width, middle + width), axis=1)
