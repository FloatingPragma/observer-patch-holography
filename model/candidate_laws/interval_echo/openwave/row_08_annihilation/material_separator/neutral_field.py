"""Supplied neutral XX field on six adjacent seam addresses per receiving arm."""

import numpy as np
from scipy.linalg import eigh
from openwave.interval_echo.row_08.continuous_receivers import measure as continuous
from .. import experiment

TIMES = np.array([0, 8, 16, 32, 48, 64, 96, 128])
LENGTH = 6
G = continuous.meter.COUPLING


def field(hop, cut=False):
    idx = np.arange(1 << LENGTH)
    output = np.zeros((len(idx), len(idx)))
    for j in range(LENGTH - 1):
        if cut and j == 2:
            continue
        a, b = 1 << (LENGTH - 1 - j), 1 << (LENGTH - 2 - j)
        diff = ((idx & a) != 0) != ((idx & b) != 0)
        output[idx[diff], idx[diff] ^ (a | b)] = hop
    return output


def calibration(omega, middle, cut):
    path_hop = omega / (2 * np.cos(2 * np.pi / 7))
    full = field(path_hop)
    half = np.array([i for i in range(len(full)) if i.bit_count() == 3])
    _, vec = eigh(full[np.ix_(half, half)], subset_by_index=(0, 0))
    ground = np.zeros(len(full))
    ground[half] = vec[:, 0]
    # Act only on the source attachment field qubit. The distal material is
    # in the lower eigenstate of its conditional, physical two-state block.
    field_source = ground[np.arange(len(full)) ^ (1 << (LENGTH - 1))]
    wing = experiment.wing(8, 3, 0)
    take = [2 * middle, 2 * middle + 1]
    far = wing[np.ix_(take, take)].copy()
    far -= np.trace(far) * np.eye(2) / 2
    e, v = eigh(far)
    np.testing.assert_allclose(e[1] - e[0], omega, atol=1e-12)
    meter0 = np.array([1.0, 0])
    inputs = np.column_stack(
        [np.kron(f, np.kron(v[:, 0], meter0)) for f in (ground, field_source)]
    )
    x = np.array([[0.0, 1], [1, 0]])
    material_meter = (
        np.kron(far, np.eye(2))
        + np.kron(np.eye(2), np.diag([0, omega]))
        + G * np.kron(x, x)
    )
    h = np.kron(field(path_hop, cut), np.eye(4)) + np.kron(
        np.eye(len(full)), material_meter
    )
    idx = np.arange(len(h))
    # Flip the last added field bit and the actual far dimer orientation.
    h[idx, idx ^ 6] += G
    e, v = eigh(h, driver="evd")
    initial_energy = np.sum(inputs * (h @ inputs), axis=0)
    overlap = v.T @ inputs
    trajectories = []
    for t in TIMES:
        state = v @ (np.exp(-1j * t * e)[:, None] * overlap)
        np.testing.assert_allclose(
            np.sum(state.conj() * (h @ state), axis=0).real,
            initial_energy,
            atol=8e-10,
            rtol=0,
        )
        trajectories.append(np.sum(abs(state[1::2]) ** 2, axis=0))
    return np.asarray(trajectories), initial_energy, path_hop
