"""The cubic incidence stencil, with an explicit periodic boundary choice."""
from __future__ import annotations
import numpy as np

PAIRS = ((0, 1), (0, 2), (1, 2))


def forward(value, axis, spacing=1.0):
    return (np.roll(value, -1, axis=axis)-value)/spacing


def adjoint(value, axis, spacing=1.0):
    return (np.roll(value, 1, axis=axis)-value)/spacing


def gradient(value, spacing=1.0):
    return np.stack([forward(value, j, spacing) for j in range(3)], axis=-1)


def divergence_adjoint(edges, spacing=1.0):
    return sum(adjoint(edges[..., j], j, spacing) for j in range(3))


def curl(edges, spacing=1.0):
    return np.stack([forward(edges[..., j], i, spacing)-forward(edges[..., i], j, spacing)
                     for i, j in PAIRS], axis=-1)


def curl_adjoint(faces, spacing=1.0):
    edges = np.zeros_like(faces)
    for f, (i, j) in enumerate(PAIRS):
        edges[..., i] -= adjoint(faces[..., f], j, spacing)
        edges[..., j] += adjoint(faces[..., f], i, spacing)
    return edges


def magnetic_divergence(faces, spacing=1.0):
    return (forward(faces[..., 2], 0, spacing)-forward(faces[..., 1], 1, spacing)
            +forward(faces[..., 0], 2, spacing))


def transverse_plane(n, mode=(1, 1, 0), length=2*np.pi, polarization=0):
    """Return a real seam field in ker D*, without replacing staggered phases."""
    a = length/n
    k = 2*np.pi*np.asarray(mode)/length
    d = np.expm1(1j*a*k)/a
    q2 = float(np.vdot(d, d).real)
    if q2 == 0:
        raise ValueError('nonzero lattice mode required')
    projector = np.eye(3)-np.outer(d, d.conj())/q2
    frame = []
    for column in projector.T:
        v = column.copy()
        for old in frame:
            v -= old*np.vdot(old, v)
        norm = np.linalg.norm(v)
        if norm > 1e-10:
            frame.append(v/norm)
    if len(frame) != 2 or polarization not in (0, 1):
        raise ValueError('two transverse polarizations required')
    coordinates = np.indices((n, n, n))
    phase = np.exp(1j*sum(2*np.pi*mode[j]*coordinates[j]/n for j in range(3)))
    e = np.real(phase[..., None]*frame[polarization])
    e /= np.linalg.norm(e)
    return e, a, q2


def verlet(e0, b0, spacing, duration, steps, speed=1.0):
    if steps <= 0 or spacing <= 0 or speed <= 0:
        raise ValueError('positive grid/time parameters required')
    dt = duration/steps
    if dt*speed*np.sqrt(12)/spacing >= 2:
        raise ValueError('step exceeds complete-grid stability bound')
    e, b = e0.copy(), b0.copy()
    for _ in range(steps):
        half = b-dt*curl(e, spacing)/2
        e = e+dt*speed**2*curl_adjoint(half, spacing)
        b = half-dt*curl(e, spacing)/2
    return e, b


def green(n):
    """Mean-zero inverse of D*D; source is one site minus the uniform mode."""
    angle = 2*np.pi*np.fft.fftfreq(n)
    one_axis = 4*np.sin(angle/2)**2
    symbol = one_axis[:, None, None]+one_axis[None, :, None]+one_axis[None, None, :]
    inverse = np.zeros_like(symbol)
    np.divide(1.0, symbol, out=inverse, where=symbol > 0)
    return np.fft.ifftn(inverse).real
