"""Two native receiver readouts with local middle-dimer shutters."""

from itertools import product

import numpy as np
from scipy import sparse
from scipy.linalg import eigh
from scipy.special import jv

from model.executable_model import federation
from openwave.interval_echo.row_08.continuous_receivers import measure as continuous
from .. import experiment, wall_flight
from ..daughter_spectrum import measure as daughter

ARMS = {"open": (), "left": (0,), "right": (1,), "both": (0, 1)}
TIMES = np.array([0, 8, 16, 32, 48, 64])
# Ring bits, three left bits, three right bits, then two receiver bits.
MIDDLE = (1 << 6, 1 << 3)


def geometry():
    ring, indices, wings = wall_flight.extended_geometry()
    paths = [[ring[j], *w] for j, w in zip(indices, wings, strict=True)]
    board = federation(0)
    assert all(
        b in board.neighbours[a] for path in paths for a, b in zip(path, path[1:])
    )
    assert all(len(path) == 7 for path in paths)
    return paths


def setup():
    boundary, initial, ordered = daughter.reference()
    e = eigh(experiment.wing(8, 3, int(boundary[0])), eigvals_only=True)
    omega = float(e[2] - e[0])
    h = continuous.generator(False, omega, continuous.meter.COUPLING)
    seed = np.zeros((h.shape[0], 2))
    seed[::4] = initial
    return h, seed, ordered, omega


def frozen(h, sides):
    # Remove only the off-diagonal swap of each frozen middle dimer. Its
    # two guarded orientation columns have equal local diagonal energy.
    idx = np.arange(h.shape[0])
    result = h.copy()
    for side in sides:
        flip = idx ^ MIDDLE[side]
        values = np.asarray(h[idx, flip]).ravel()
        np.testing.assert_allclose(values, -0.5, atol=1e-14)
        result -= sparse.csr_matrix((values, (idx, flip)), shape=h.shape)
    return result


def flight(h, seed, sides, times=TIMES):
    """Any fixed middle bits partition the exact Hamiltonian into blocks."""
    output = np.zeros((len(times), *seed.shape), dtype=complex)
    indices = np.arange(h.shape[0])
    for bits in product((0, 1), repeat=len(sides)):
        select = np.ones(len(indices), bool)
        for side, bit in zip(sides, bits, strict=True):
            select &= ((indices & MIDDLE[side]) != 0) == bit
        take = indices[select]
        sub = h[take, :][:, take].toarray(order="F")
        offset = sub.diagonal().mean()
        sub[np.diag_indices(len(sub))] -= offset
        e, v = eigh(sub, driver="evd", overwrite_a=True, check_finite=False)
        overlap = v.T @ seed[take]
        amplitudes = (
            np.exp(-1j * e[:, None, None] * times[None, :, None]) * overlap[:, None, :]
        )
        amplitudes = amplitudes.reshape(len(e), -1)
        propagated = v @ amplitudes.real + 1j * (v @ amplitudes.imag)
        block = propagated.reshape(len(e), len(times), -1).transpose(1, 0, 2)
        output[:, take, :] = np.exp(-1j * offset * times)[:, None, None] * block
    return output


def panels(states):
    """Lexicographic four ring bits, two middle bits, two receiver bits."""
    idx = np.arange(states.shape[1])
    labels = (
        16 * (idx >> 8)
        + 8 * ((idx & MIDDLE[0]) != 0)
        + 4 * ((idx & MIDDLE[1]) != 0)
        + (idx & 3)
    )
    return np.array(
        [
            [
                np.bincount(labels, weights=np.abs(f[:, src]) ** 2, minlength=256)
                for src in range(2)
            ]
            for f in states
        ]
    )


def sparse_endpoint(h, seed, exposure):
    diagonal = h.diagonal()
    non_diag = np.asarray(abs(h).sum(axis=1)).ravel() - abs(diagonal)
    lower, upper = (diagonal - non_diag).min(), (diagonal + non_diag).max()
    shift = (upper + lower) / 2
    radius = (upper - lower) / 2
    scaled = (h - shift * sparse.eye(h.shape[0], format="csr")) / radius
    argument = radius * exposure
    order = int(np.ceil(argument + 20 * np.cbrt(argument + 1) + 50))
    k = np.arange(order + 1)
    coeff = 2 * (-1j) ** (k % 4) * jv(k, argument)
    coeff[0] /= 2
    assert max(abs(coeff[-5:])) < 1e-20
    p, q = seed, scaled @ seed
    total = coeff[0] * p + coeff[1] * q
    for c in coeff[2:]:
        p, q = q, 2 * (scaled @ q) - p
        total += c * q
    return total * np.exp(-1j * shift * exposure)


def native(pinned, original, seed):
    h = frozen(original, pinned)
    states = flight(h, seed, pinned)
    np.testing.assert_allclose((abs(states) ** 2).sum(axis=1), 1, atol=4e-11, rtol=0)
    energies = np.sum(seed * (h @ seed), axis=0)
    for state in states:
        np.testing.assert_allclose(
            np.sum(state.conj() * (h @ state), axis=0).real, energies, atol=3e-9, rtol=0
        )
    laws = panels(states)
    return laws, energies
