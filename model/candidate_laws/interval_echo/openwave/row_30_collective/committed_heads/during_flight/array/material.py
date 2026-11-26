"""Propagate measured wall effects after the last of two connector reads."""

from itertools import combinations

import numpy as np
from scipy.linalg import eigh

from .. import measure as flight


def after_reads(kappa, n, head, sources):
    (attached, cut), _ = flight.generators(kappa, n, head, head + 1)
    psi = np.zeros((1 << n, len(sources)))
    for col, bond in enumerate(sources):
        psi[flight.inherited.initial_word(n, bond), col] = 1
    first = flight.evolve(attached, psi, 1)
    a = flight.evolve(attached, first, 1)
    b = flight.evolve(cut, first, 1)
    # Columns for records 00/11 and 01/10, on each source. Joint
    # amplitudes for either final coin parity are these columns divided by 4.
    return np.column_stack([a + b, a - b])


def rotation(kappa, length):
    # Work between two fixed end spins Z=+1 and Z=-1. Surviving
    # Majoranas are B_left,A_0,B_0,...,B_last,A_right.
    links = np.tile([2.0 * kappa, -1.0], length + 1)[:-1]
    generator = np.diag(links, 1) - np.diag(links, -1)
    return eigh(1j * generator)


def covariance(psi, start, length):
    """Pull each even Majorana bilinear back to the native dimer bit words."""
    full_n = (len(psi) - 1).bit_length()
    if start < 0 or start + length > full_n:
        raise ValueError("segment outside native strand")
    words = np.arange(len(psi), dtype=np.int32)
    phase_parity = np.array([1 - 2 * (i.bit_count() % 2) for i in range(len(psi))])
    # Each (phase, X-mask, Z-mask) denotes i^phase X_mask Z_mask.
    ops = []
    for j in range(length + 2):
        z = 1 << j
        if j > 0:
            ops.append((0, z - 1, z))
        if j <= length:
            ops.append((1, 2 * z - 1, z))
    dim = len(ops)
    result = np.zeros((dim, dim, psi.shape[1]))
    physical = (1 << length) - 1
    for a, b in combinations(range(dim), 2):
        pa, xa, za = ops[a]
        pb, xb, zb = ops[b]
        x, z = xa ^ xb, za ^ zb
        if x & (1 | (1 << (length + 1))):
            raise ValueError("paired instrument rotates a fixed endpoint")
        phase = 1j ** (1 + pa + pb) * (-1) ** (za & xb).bit_count()
        if z & (1 << (length + 1)):
            phase = -phase
        x = ((x >> 1) & physical) << start
        z = ((z >> 1) & physical) << start
        sign = phase_parity[words & z]
        moment = (psi[words ^ x].conj() * psi).T @ sign
        result[a, b] = (phase * moment).real
        result[b, a] = -result[a, b]
    return result


def walls(kappa, psi, start, length, sites, waits):
    freq, modes = rotation(kappa, length)
    moment = covariance(psi, start, length)
    transformed = np.einsum(
        "ai,abk,bj->ijk", modes.conj(), moment, modes, optimize=True
    )
    coeff = []
    for site in sites:
        a = 2 * (site - start)
        coeff.append(
            (modes[a, :, None] * modes[a + 1, None, :].conj())[:, :, None] * transformed
        )
    coef = np.stack(coeff, axis=2).reshape(len(freq) ** 2, -1)
    gaps = (freq[:, None] - freq).ravel()
    chunks = []
    for batch in np.array_split(waits, max(1, (len(waits) + 255) // 256)):
        phase = np.exp(-1j * batch[:, None] * gaps)
        chunks.append((phase @ coef).real)
    return np.concatenate(chunks).reshape(len(waits), len(sites), -1)


def joint(kappa, n, head, sites, sources, waits):
    if any(not head < site < n for site in sites):
        raise ValueError("right-side wall required")
    psi = after_reads(kappa, n, head, sources)
    masses = np.sum(abs(psi) ** 2, axis=0)
    attached = walls(kappa, psi, 0, n, sites, waits)
    cut = walls(kappa, psi, head, n - head, sites, waits)
    hits = (2 * masses + attached + cut) / 32
    pairs = np.stack((hits, masses / 8 - hits), axis=-1)
    # [wait, site, source, 00/01/10/11, hit/miss]
    pairs = pairs.reshape(len(waits), len(sites), 2, len(sources), 2)
    pairs = np.take(pairs, [0, 1, 1, 0], axis=2).swapaxes(2, 3)
    np.testing.assert_allclose(pairs.sum(axis=(-1, -2)), 1, rtol=0, atol=2e-10)
    if pairs.min() < -2e-11:
        raise ValueError("negative native event probability")
    return np.clip(pairs, 0, 1)
