"""All-open native direction blocks and a local signed terminal read."""

import numpy as np
from scipy import sparse
from scipy.special import jv
from scipy.sparse.linalg import eigsh

from ...row_31_orbitals.occupied_write import core

TICKS = np.array([0, 40, 80, 120, 160, 200, 240, 320, 400])


def preparation():
    law = core.LAW
    state = law.echo()[:, law.SRC.index((1, 0))]
    expected = np.zeros(99)
    for coin, word in ((0, (0, 1)), (1, (1, 0))):
        expected[core.local_index((coin, 1), word)] = 1 / np.sqrt(2)
    np.testing.assert_allclose(state, expected, rtol=0, atol=2e-14)
    minus = law.echo()[:, law.SRC.index((0, 1))]
    reversed_state = np.zeros(99)
    for coin, word in ((0, (0, 1)), (1, (1, 0))):
        reversed_state[core.local_index((coin, -1), word)] = 1 / np.sqrt(2)
    np.testing.assert_allclose(minus, reversed_state, rtol=0, atol=2e-14)
    return state


def block(board, seam, strength, odd):
    """Subtract the shared uniform-load energy only. Both parities retain their offset."""
    n = board.n
    a, b = seam
    cut = (board.laplacian().astype(float) / 2).tolil()
    for x, y in ((a, b), (b, a)):
        cut[x, x] -= strength + 0.5
        cut[x, y] += 0.5
    chi = np.zeros(2 * n)
    chi[[b, n + a]] = 1 / np.sqrt(2)
    mu = strength + int(odd)
    rank = sparse.csr_matrix(chi[None, :])
    h = sparse.block_diag([cut, cut], format="csr") + mu * (
        sparse.eye(2 * n, format="csr") - rank.T @ rank
    )
    return h, chi


def evolve(h, init, times, lo, hi):
    """Chebyshev expansion on a certified common spectral interval."""
    times = np.asarray(times)
    center, radius = (lo + hi) / 2, (hi - lo) / 2
    x = (h - center * sparse.eye(h.shape[0], format="csr")) / radius
    z = times * radius
    peak = float(max(z))
    # Past the turning point the exponentially decaying Bessel tail can be
    # bounded by a shifted-contour integral. Use the best common shift.
    last = int(np.ceil(peak + 18 * max(1, peak ** (1 / 3)) + 35))
    if peak:
        shift = np.arccosh((last + 1) / peak)
        tail = (
            2 * np.exp(peak * np.sinh(shift) - (last + 1) * shift) / (-np.expm1(-shift))
        )
        if tail > 1e-13:
            raise ArithmeticError("Chebyshev tail exceeds declared precision")
    coeff = jv(np.arange(last + 1)[:, None], z[None, :])
    real = coeff[0, :, None] * init[None, :]
    imag = np.zeros_like(real)
    previous, current = init, x @ init
    for k in range(1, last + 1):
        target = imag if k % 2 else real
        sign = 1 if k % 4 in (0, 3) else -1
        for i in np.flatnonzero(abs(coeff[k]) > 1e-19):
            target[i] += (2 * sign * coeff[k, i]) * current
        if k != last:
            previous, current = current, 2 * (x @ current) - previous
    result = (real + 1j * imag) * np.exp(-1j * center * times[:, None])
    np.testing.assert_allclose(np.sum(abs(result) ** 2, axis=1), 1, atol=2e-9, rtol=0)
    return result


def clock_trace(board, seam, strength, times):
    waves, audit = [], []
    for odd in (False, True):
        h, chi = block(board, seam, strength, odd)
        waves.append(evolve(h, chi, times, -strength, 6 + strength + 1))
        if strength == 8:
            val, vectors = eigsh(h, k=1, which="SA", v0=chi, tol=1e-13)
            ground = vectors[:, 0]
            residual = np.linalg.norm(h @ ground - val[0] * ground)
            if residual > 1e-9 or val[0] >= -5.5:
                raise ArithmeticError("Ground isolation or residual check failed")
            audit.append((float(val[0]), ground, float(ground @ chi), float(residual)))
    even, odd = waves
    plus = np.sum(abs((even + odd) / 2) ** 2, axis=1)
    minus = np.sum(abs((even - odd) / 2) ** 2, axis=1)
    np.testing.assert_allclose(plus + minus, 1, atol=2e-9, rtol=0)
    info = None
    if audit:
        e, o = audit
        info = dict(
            absolute_ground_energies=[6 * strength + e[0], 6 * strength + o[0]],
            gap=o[0] - e[0],
            ground_weights=[e[2] ** 2, o[2] ** 2],
            signed_gap_coefficient=float(e[2] * o[2] * (e[1] @ o[1])),
            ground_residuals=[e[3], o[3]],
        )
    # Complement is exactly zero by signed-direction invariance, not rounding.
    return np.clip((1 + plus - minus) / 2, 0, 1), info
