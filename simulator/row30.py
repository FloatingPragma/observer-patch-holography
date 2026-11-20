"""OpenWave Row 30: the matter-wave relation of a pair formed in an encounter.

The mature step of a formed pair of the OPH law is the contact coin
D = exp(i phi/2) [(1 + exp(-i phi)) |S><S| - I] on the formed slot vector S, then the
full hop. On a plane wave exp(i k.x) the step keeps the span of S and its routed image,
and on that span its two eigenphases differ by 2 omega with

    cos(omega) = cos(m) (2 w_private + 2 w_spatial sum_i cos k_i),   m = phi / 2,

w the slot weights of S. This module computes omega from the native 8 x 8 step for the
registered probes and applies the observer's test: the kinetic increment
(omega^2 - m^2) / (c^2 |k|^2) has to lie within 3.5% of one, with c^2 the scalar
propagation standard calibrated on the unbound walk. The changed law moves half of the
spatial weight of S into the private slots; it keeps the rest gap and halves the
curvature.
"""
from __future__ import annotations

import itertools
import time

import numpy as np

from .law import LAWS, PARTNER, SLOT_AXIS, SLOT_SIGN, SLOTS

C2_SCALAR = 0.24935023          # scalar propagation standard, radar units^2 per home tick^2
TOLERANCE = 0.035
PERIODS = (64, 128, 256)

# Recorded by the observers in the Row 30 experiment.
RECORDED = {"items": [
    ["ordinary branch, prepared waves", "72 of 72 held comparisons passed"],
    ["ordinary branch, separate point launch", "40 of 40 nonzero forecasts passed, with no material mass fitted"],
    ["compact collision branch", "14 of 14 signed axial and oblique directions passed both bounds at each of periods 64, 128 and 256"],
    ["agreement with the parent calibration", "terminal probabilities within 6e-15"],
    ["limiting squared kinetic coefficient, rail zero", "0.3215% offset at the fixed gap"],
    ["private circumference halved", "each of the first three rest gaps nearly doubles"],
]}


def probes():
    """The fourteen signed order-three directions: six axial, eight sign reflections of (1, 2, 2)."""
    axial = [tuple(s * 3 if i == a else 0 for i in range(3)) for a in range(3) for s in (1, -1)]
    oblique = [(sx, 2 * sy, 2 * sz) for sx, sy, sz in itertools.product((1, -1), repeat=3)]
    return axial + oblique


def formed_vector(changed: bool = False) -> np.ndarray:
    """Equal slot weights, or half the spatial weight moved into the private pair."""
    weights = np.full(SLOTS, 1 / SLOTS)
    if changed:
        spatial = SLOT_AXIS < 3
        moved = weights[spatial].sum() / 2
        weights[spatial] /= 2
        weights[~spatial] += moved / 2
    return np.sqrt(weights).astype(complex)


def mature_step(angles: np.ndarray, phi: float, vector: np.ndarray) -> np.ndarray:
    """The 8 x 8 operator of one mature update on exp(i k.x): routing after the contact coin."""
    projector = np.outer(vector, vector.conj())
    coin = np.exp(0.5j * phi) * ((1 + np.exp(-1j * phi)) * projector - np.eye(SLOTS))
    signed = np.asarray(angles, float)[SLOT_AXIS] * SLOT_SIGN
    routing = np.zeros((SLOTS, SLOTS), complex)
    for slot in range(SLOTS):
        routing[PARTNER[slot], slot] = np.exp(-1j * signed[slot])
    return routing @ coin


def frequency(angles, phi, vector) -> float:
    """omega from the step on span{S, step S}, where its eigenvalues are exp(+i omega), exp(-i omega)."""
    step = mature_step(angles, phi, vector)
    image = step @ vector
    along = np.vdot(vector, image)
    if np.linalg.norm(image - along * vector) < 1e-12:
        return float(abs(np.angle(along)))
    basis, _ = np.linalg.qr(np.stack([vector, image], axis=1))
    restricted = basis.conj().T @ step @ basis
    if np.max(np.abs(step @ basis - basis @ restricted)) > 1e-12:
        raise ArithmeticError("the span of S and its image is not invariant")
    eigenvalues = np.linalg.eigvals(restricted)
    if abs(eigenvalues[0] - np.conj(eigenvalues[1])) > 1e-10:
        raise ArithmeticError("the two eigenphases on the span are not conjugate")
    return float(abs(np.angle(eigenvalues[0])))


def closed_form(angles, phi, vector) -> float:
    weights = np.abs(vector) ** 2
    w_private = weights[SLOT_AXIS == 3].mean()
    w_spatial = weights[SLOT_AXIS < 3].mean()
    cosines = np.cos(np.asarray(angles, float)[:3]).sum()
    return float(np.arccos(np.cos(phi / 2) * (2 * w_private + 2 * w_spatial * cosines)))


def row30_run():
    """Every probe at every period, for the formed vector of the law and for the changed law."""
    started = time.time()
    law = LAWS["oph"]
    phi = law.rail_contact_phases[0]
    rows, worst_identity = [], 0.0
    for changed in (False, True):
        vector = formed_vector(changed)
        rest = frequency(np.zeros(4), phi, vector)
        for period in PERIODS:
            for probe in probes():
                angles = np.array(list(probe) + [0]) * 2 * np.pi / period
                omega = frequency(angles, phi, vector)
                worst_identity = max(worst_identity, abs(omega - closed_form(angles, phi, vector)))
                k2 = float(np.sum(angles[:3] ** 2))
                kinetic = (omega ** 2 - rest ** 2) / (C2_SCALAR * k2)
                rows.append({"law": "changed" if changed else "oph", "period": period,
                             "probe": list(probe), "k": float(np.sqrt(k2)), "omega": omega,
                             "rest": rest, "kinetic": kinetic,
                             "passes": bool(abs(kinetic - 1) <= TOLERANCE)})
    m = phi / 2
    coefficient = 2 * (1 / SLOTS) * m / np.tan(m)
    return {"law": law, "phi": phi, "rows": rows, "rest": rows[0]["rest"],
            "coefficient": coefficient, "c2": C2_SCALAR, "identity_error": worst_identity,
            "recorded": RECORDED, "seconds": time.time() - started}


def matter_frames(size: int = 64, updates: int = 128, stride: int = 2, order: int = 3):
    """A resting and a moving formed-matter packet: amplitude and phase of the S-reading over (x, y)."""
    from .spectral import PairPacket
    law = LAWS["oph"]
    packets = [PairPacket(law, (0, 0), size=size, carrier_order=order, width=1.2, support=4,
                          origin=(size // 4, size // 2, 0)),
               PairPacket(law, (1, 0), size=size, carrier_order=order, width=1.2, support=4,
                          origin=(size // 4, size // 2, 0))]
    equal = np.ones(SLOTS) / np.sqrt(SLOTS)
    frames = []
    for t in range(0, updates + 1, stride):
        layer = []
        for packet in packets:
            state = packet.state(t)[:, :, 0, 0, :] @ equal.conj()
            grid = np.zeros((size, size), complex)
            grid[np.ix_(packet.orders[0] % size, packet.orders[1] % size)] = state
            field = np.fft.ifft2(grid) * size
            field = np.roll(field, tuple(packet.origin[:2]), axis=(0, 1))
            layer.append(field)
        frames.append(layer)
    return np.array(frames), stride
