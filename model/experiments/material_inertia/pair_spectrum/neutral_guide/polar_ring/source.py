"""One-upper input amplitudes from addressed polar-material contacts."""

import numpy as np
from scipy.linalg import expm


def mixing(angle, private):
    """A conditional connector, local rotations and echo implement an XY mixing gate."""
    eye = np.eye(2)
    x = np.array([[0, 1], [1, 0]])
    y = np.array([[0, -1j], [1j, 0]])
    z = np.diag([1, -1])
    xx = np.kron(x, x)
    load = np.array([0, 1]) if private else np.array([1, 0])
    cost = np.diag(32 * (1 - np.kron(load, load)))
    refocused = (cost + xx @ cost @ xx) / 2
    refocused -= np.trace(refocused) / 4 * np.eye(4)
    np.testing.assert_allclose(refocused, -8 * np.kron(z, z), atol=2e-14)
    # One local flip reverses ZZ. For each axis, rotate its eigenbasis onto Z.
    reverse = np.kron(x, eye)
    zz = reverse @ expm(-1j * angle / 16 * refocused) @ reverse
    ry = expm(-0.25j * np.pi * y)
    rx = expm(0.25j * np.pi * x)
    u = np.eye(4)
    for basis in (ry, rx):
        rotate = np.kron(basis, basis)
        u = (rotate @ zz @ rotate.conj().T) @ u
    return u


def prepared(profile):
    n = len(profile)
    p = np.asarray(profile, dtype=complex)
    np.testing.assert_allclose(np.linalg.norm(p), 1, atol=2e-13)
    magnitudes = abs(p)
    v = np.zeros(n, dtype=complex)
    v[0] = 1
    for j in range(n - 1):
        angle = np.arctan2(np.linalg.norm(magnitudes[j + 1 :]), magnitudes[j])
        gate = mixing(angle, private=(j % 2 == 0))
        np.testing.assert_allclose(gate[[0, 3]][:, [1, 2]], 0, atol=1e-14)
        assert np.max(abs(gate[np.ix_([0, 3], [0, 3])] - np.eye(2))) < 2e-14
        v[[j + 1, j]] = gate[np.ix_([1, 2], [1, 2])] @ v[[j + 1, j]]
        v[j + 1] *= 1j
    v *= np.exp(1j * np.angle(p))
    np.testing.assert_allclose(v, p, atol=4e-14)
    out = np.zeros(1 << n, dtype=complex)
    out[1 << np.arange(n)] = v
    return out
