"""Exact plane-wave evolution of prepared packets of the torus law.

The law's update acts on a plane wave exp(i k.x) as an 8 x 8 unitary.
Its mobile band has the eigenvalue exp(-iE) with cos E = (sum of the four
cosines) / 4. A packet is a finite sum of mobile-band waves, so its exact
state at update n is the same sum with every wave multiplied by exp(-inE).
This is how the two-clock experiment evolves packets on carriers with
512 or 1024 positions per spatial axis and 4096 or 8192 private positions.
``test_simulator.py`` compares every reading here with the literal update.
"""
from __future__ import annotations

import numpy as np

from .law import EQUAL, PARTNER, SLOT_AXIS, SLOT_SIGN, coin

TAU = 2 * np.pi


def mobile_band(angles: np.ndarray):
    """Phases E and unit eigenvectors of the mobile band, exp(-iE), for wave numbers (..., 4)."""
    signed = angles[..., SLOT_AXIS] * SLOT_SIGN
    shift = np.exp(-1j * signed[..., PARTNER])
    flipped = shift * EQUAL
    overlap = np.sum(flipped * EQUAL.conj(), axis=-1).real
    root = np.exp(-1j * np.arccos(np.clip(overlap, -1.0, 1.0)))
    vector = EQUAL - root[..., None] * flipped
    degenerate = np.linalg.norm(vector, axis=-1) < 1e-12
    vector[degenerate] = EQUAL
    vector /= np.linalg.norm(vector, axis=-1)[..., None]
    vector *= np.exp(-1j * np.angle(np.sum(EQUAL.conj() * vector, axis=-1)))[..., None]
    evolved = shift * coin(vector, 1.0)[..., PARTNER]
    poles = np.sum(vector.conj() * evolved, axis=-1)
    if np.max(np.abs(evolved - poles[..., None] * vector)) > 1e-12:
        raise ArithmeticError("a prepared wave is no eigenvector of the update")
    return -np.angle(poles), vector


def _window(center: int, width: float, support: int):
    delta = np.arange(-support, support + 1)
    return center + delta, np.exp(-0.5 * (delta / width) ** 2)


class Packet:
    """One prepared excitation: a clock doublet, an optical pulse or a gapless reference.

    ``kind`` is "clock" (private orders m and m + 1, uniform in z), "optical"
    (private order 0, a pulse travelling around the z ring) or "reference"
    (private order 0, uniform in z). ``route`` gives the signed launch
    direction in x and y in units of ``carrier_order``.
    """

    def __init__(self, kind, route=(0, 0), size=512, private_size=4096, carrier_order=16,
                 private_order=128, width=1.5, support=6, transverse_sign=1, origin=(0, 0, 0)):
        if kind not in ("clock", "optical", "reference"):
            raise ValueError("kind is clock, optical or reference")
        self.kind, self.size, self.private_size = kind, size, private_size
        self.origin = np.asarray(origin, int) % size
        windows = [_window(int(d) * carrier_order, width, support) for d in route]
        if kind == "optical":
            windows.append(_window(transverse_sign * carrier_order, width, support))
        else:
            windows.append((np.array([0]), np.ones(1)))
        private = [private_order, private_order + 1] if kind == "clock" else [0]
        windows.append((np.array(private), np.ones(len(private))))
        self.orders = [orders for orders, _ in windows]
        index = np.stack(np.meshgrid(*self.orders, indexing="ij"), axis=-1)
        angles = index * (TAU / np.array([size, size, size, private_size]))
        weights = np.prod(np.stack(np.meshgrid(*(w for _, w in windows), indexing="ij")), axis=0)
        weights = weights / np.linalg.norm(weights.ravel())
        self.energy, vectors = mobile_band(angles)
        self.spinors = vectors * weights[..., None]

    def state(self, t: float) -> np.ndarray:
        return np.exp(-1j * t * self.energy)[..., None] * self.spinors

    def marginal(self, t: float, axis: int) -> np.ndarray:
        """Probability of each position on one spatial axis, every other index traced."""
        amplitudes = np.moveaxis(self.state(t), axis, 0)
        full = np.zeros((self.size,) + amplitudes.shape[1:], complex)
        full[self.orders[axis] % self.size] = amplitudes
        fields = np.fft.ifft(full, axis=0) * np.sqrt(self.size)
        density = np.sum(np.abs(fields) ** 2, axis=tuple(range(1, fields.ndim)))
        return np.roll(density, self.origin[axis])

    def ramsey(self, t: float, site) -> np.ndarray:
        """Probabilities at a spatial address and the four private positions 0, P/4, P/2, 3P/4."""
        offset = np.asarray(site, int) - self.origin
        phases = [np.exp(1j * TAU * offset[a] * self.orders[a] / self.size) for a in range(3)]
        slice_ = np.einsum("i,j,k,ijkps->ps", *phases, self.state(t))
        private = np.exp(0.5j * np.pi * np.outer(np.arange(4), self.orders[3]))
        amplitude = private @ slice_ / np.sqrt(self.size ** 3 * self.private_size)
        return np.sum(np.abs(amplitude) ** 2, axis=-1)


def circular_position(density: np.ndarray) -> float:
    """Envelope position on a ring from its first circular moment (sites)."""
    size = len(density)
    return float(np.angle(density @ np.exp(1j * TAU * np.arange(size) / size)) * size / TAU)


def track(packet: Packet, axis: int, times) -> np.ndarray:
    """Unwrapped envelope positions along one axis at the given update counts."""
    angles = np.array([circular_position(packet.marginal(t, axis)) for t in times]) * TAU / packet.size
    return np.unwrap(angles) * packet.size / TAU


def beat_phase(probabilities: np.ndarray) -> float:
    """The four-aperture reading of the clock hand: atan2(p1 - p3, p0 - p2)."""
    p = np.asarray(probabilities)
    return float(np.arctan2(p[1] - p[3], p[0] - p[2]))


def pair_band(angles: np.ndarray, contact_phase: float):
    """Mobile band of the bound pair: cos E = a cos(phi/2), a the mean of the four cosines."""
    from .law import contact_coin
    signed = angles[..., SLOT_AXIS] * SLOT_SIGN
    shift = np.exp(-1j * signed[..., PARTNER])
    a = np.mean(shift.real, axis=-1)
    energy = np.arccos(np.clip(a * np.cos(contact_phase / 2), -1.0, 1.0))
    direct = np.exp(-0.5j * contact_phase) * shift * EQUAL
    vector = direct - np.exp(1j * energy)[..., None] * EQUAL
    vector /= np.linalg.norm(vector, axis=-1)[..., None]
    vector *= np.exp(-1j * np.angle(np.sum(vector * EQUAL.conj(), axis=-1)))[..., None]
    evolved = shift * contact_coin(vector, contact_phase)[..., PARTNER]
    if np.max(np.abs(evolved - np.exp(-1j * energy)[..., None] * vector)) > 2e-12:
        raise ArithmeticError("a prepared pair wave is no eigenvector of the contact update")
    return energy, vector


class PairPacket(Packet):
    """A bound pair in its two rails, travelling in x and y on the bound sector.

    The rails take the place of the private orders: ``orders[3]`` = (0, 1) labels them,
    and each rail evolves with its own contact phase and extra phase.
    """

    def __init__(self, law, route=(0, 0), size=512, carrier_order=16, width=1.5, support=6,
                 origin=(0, 0, 0)):
        self.kind, self.size, self.private_size = "pair", size, 8 * size
        self.origin = np.asarray(origin, int) % size
        windows = [_window(int(d) * carrier_order, width, support) for d in route]
        self.orders = [w[0] for w in windows] + [np.array([0]), np.array([0, 1])]
        nx, ny = np.meshgrid(*self.orders[:2], indexing="ij")
        angles = np.stack([nx, ny, 0 * nx, 0 * nx], axis=-1) * (TAU / size)
        profile = np.outer(windows[0][1], windows[1][1])
        profile /= np.linalg.norm(profile)
        self.energy = np.zeros(nx.shape + (1, 2))
        self.spinors = np.zeros(nx.shape + (1, 2, 8), complex)
        for rail, (phi, extra) in enumerate(zip(law.rail_contact_phases, law.rail_extra_phases)):
            energy, vector = pair_band(angles, phi)
            self.energy[:, :, 0, rail] = energy + extra
            self.spinors[:, :, 0, rail] = vector * profile[..., None] / np.sqrt(2)

    def rails(self, t: float, site) -> np.ndarray:
        """Amplitudes of the two rails, eight slots each, at one address."""
        offset = np.asarray(site, int) - self.origin
        px = np.exp(1j * TAU * offset[0] * self.orders[0] / self.size)
        py = np.exp(1j * TAU * offset[1] * self.orders[1] / self.size)
        return np.einsum("i,j,ijrs->rs", px, py, self.state(t)[:, :, 0]) / self.size


def rail_phase(rails: np.ndarray) -> float:
    """The four-outcome rail analyzer: atan2(p_i - p_-i, p_1 - p_-1) with p_a = |r0 + a r1|^2 / 4."""
    apertures = np.array([1, 1j, -1, -1j])
    probes = rails[0][None, :] + apertures[:, None] * rails[1][None, :]
    p = np.sum(np.abs(probes) ** 2, axis=-1) / 4
    return float(np.arctan2(p[1] - p[3], p[0] - p[2]))
