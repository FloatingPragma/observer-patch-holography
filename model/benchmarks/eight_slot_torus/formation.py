"""Forming a pair through the capture vertex of the torus law, which exports a marker.

The eight-slot torus law gives every junction a periodic wire of neutral modes. When the
two constituents meet head-on at a junction, with constituent A in slot 5 (arrived
moving +z) and B in slot 4 (arrived moving -z), the vertex rotates that channel by
the angle pi/4 into the matched equal-slot pair sum_s e_s |s, s> plus one marker at
the start of the junction's wire. The marker moves one link per update. After the
vertex, both constituents take the contact update.

A marker that has left its junction never returns within the run (the wire is longer
than the run), so every capture event starts its own branch, orthogonal to all
others. A branch lives in the invariant bound sector, where the pair is one walker.
A branch also carries the pair's age, born at zero, and its hop
opens with the age (``law.step_maturing``). Branches whose birth weight is below
1e-14 are dropped and their total is reported. ``test_simulator.py`` compares the
branch evolution with the literal update that carries the marker wires explicitly.
"""
from __future__ import annotations

import time

import numpy as np

from .law import EQUAL, LAWS, contact_coin, hop_angle, route, step_contact
from .spectral import mobile_band

ANGLE = np.pi / 4
COLLISION = (5, 4)
AGE_PERIOD = 8192


def birth_ages(update, births):
    """Bound register values, given the launch-service indices of the actual writes."""
    born = np.asarray(births)
    if (isinstance(update, bool) or not isinstance(update, (int, np.integer)) or update < 0
            or born.dtype.kind not in "iu" or np.any(born < 0)
            or (born.size and int(born.max()) > int(update))):
        raise ValueError("births must be nonnegative integer services through this update")
    residues = born.astype(np.int64, copy=False) % AGE_PERIOD
    return (int(update) % AGE_PERIOD - residues) % AGE_PERIOD


def _strip_packet(size: int, centre: int, order: int, width: float = 1.2, support: int = 4) -> np.ndarray:
    """One constituent on a strip along z: a mobile-band packet of central order ``order``."""
    orders = order + np.arange(-support, support + 1)
    weights = np.exp(-0.5 * ((orders - order) / width) ** 2)
    angles = np.zeros((len(orders), 4))
    angles[:, 2] = 2 * np.pi * orders / size
    _, spinors = mobile_band(angles)
    z = np.arange(size) - centre
    amplitude = np.exp(2j * np.pi * np.outer(z, orders) / size) @ (weights[:, None] * spinors)
    amplitude /= np.linalg.norm(amplitude)
    return amplitude.reshape(1, 1, size, 1, 8)


def first_writes(free: np.ndarray, angle: float = ANGLE):
    """Apply the native empty-wire rotation on every coincident (5,4) channel.

    The complex array at the output gives the amplitude of each distinct
    first-write wire. Its attached matched material ray is always EQUAL.
    """
    shape = free.shape
    if len(shape) != 10 or shape[:5] != shape[5:] or shape[4] != len(EQUAL):
        raise ValueError("two identical four-direction carriers with eight slots required")
    count = int(np.prod(shape[:4]))
    if not np.isfinite(angle):
        raise ValueError("a finite capture angle is required")
    remaining = np.array(free, dtype=np.result_type(free.dtype, np.float64), order="C", copy=True)
    view = remaining.reshape(count, len(EQUAL), count, len(EQUAL))
    here = np.arange(count)
    entry = (here, COLLISION[0], here, COLLISION[1])
    incoming = view[entry]
    view[entry] = np.cos(angle)*incoming
    return remaining, (np.sin(angle)*incoming).reshape(shape[:4])


def capture(free: np.ndarray, angle: float = ANGLE):
    """Apply the vertex to the unmarked two-constituent state; return it and the new branches.

    Only the unmarked sector carries the collision channel, and no marker sits at the
    start of any wire before the vertex acts, so the rotation acts on |5,4; 0> alone.
    """
    remaining, writes = first_writes(free, angle)
    size = free.shape[2]
    if writes.shape != (1, 1, size, 1):
        raise ValueError("the explicit branch-list receiver uses a z strip")
    born = []
    for z in range(size):
        a = writes[0, 0, z, 0]
        # The strip interface lists each incident channel even when the
        # chosen control gives its write exactly zero weight.
        if free[0, 0, z, 0, COLLISION[0], 0, 0, z, 0, COLLISION[1]] == 0:
            continue
        chi = np.zeros((1, 1, size, 1, 8), complex)
        chi[0, 0, z, 0, :] = a * EQUAL
        born.append((z, chi))
    return remaining, born


def formation_run(size: int = 48, updates: int = 72, orders=(7, -5), contact_phase=None,
                  angle: float = ANGLE, support: int = 4, law_key: str = "control-abrupt",
                  separation: int | None = None):
    """Two constituents launched toward each other; capture events, the formed pair, its markers."""
    started = time.time()
    law = LAWS[law_key]
    phi = law.rail_contact_phases[0] if contact_phase is None else contact_phase
    gap = size // 4 if separation is None else separation
    a = _strip_packet(size, size // 2 - gap, orders[0], support=support)
    b = _strip_packet(size, size // 2 + gap, orders[1], support=support)
    free = np.einsum("abcde,fghij->abcdefghij", a, b)
    # Every capture event is one branch: its junction, its birth update and its bound-sector
    # state. Branches are stacked on a leading axis and stepped together; each has its own age.
    sites, births = np.zeros(0, int), np.zeros(0, int)
    states = np.zeros((0, 1, 1, size, 1, 8), complex)
    dropped = 0.0
    joint, bound, collected, captured, clicks = [], [], [], [], []
    for n in range(updates + 1):
        joint.append(np.sum(np.abs(free[0, 0, :, 0, :, 0, 0, :, 0, :]) ** 2, axis=(1, 3)))
        density = np.sum(np.abs(states[:, 0, 0, :, 0, :]) ** 2, axis=-1)
        pair = density.sum(axis=0) if len(states) else np.zeros(size)
        # A coherent collector joins the marker wires of all junctions; its common port
        # leaves the birthplace unresolved. Conditioned on that port, markers emitted at
        # the same update add coherently with amplitude 1 / sqrt(number of junctions).
        if len(states):
            keys, index = np.unique(births, return_inverse=True)
            summed = np.zeros((len(keys), size, 8), complex)
            np.add.at(summed, index, states[:, 0, 0, :, 0, :])
            common = np.sum(np.abs(summed) ** 2, axis=(0, 2)) / size
        else:
            common = np.zeros(size)
        bound.append(pair)
        collected.append(common)
        captured.append(float(pair.sum()))
        clicks.append(float(common.sum()))
        if n == updates:
            break
        free, born = capture(free, angle)
        kept = [(z, chi) for z, chi in born if np.sum(np.abs(chi) ** 2) > 1e-14]
        dropped += sum(float(np.sum(np.abs(chi) ** 2)) for z, chi in born if np.sum(np.abs(chi) ** 2) <= 1e-14)
        if kept:
            sites = np.concatenate([sites, [z for z, _ in kept]])
            births = np.concatenate([births, [n] * len(kept)])
            states = np.concatenate([states, np.stack([chi for _, chi in kept])])
        free = step_contact(free, phi)
        if len(states):
            ages = birth_ages(n, births)
            theta = np.array([hop_angle(int(age), law.maturation) for age in ages])[:, None, None, None, None, None]
            sign = law.maturation[0] if law.maturation else -1
            coined = contact_coin(states, law.rail_contact_phases[0]) * np.exp(-1j * law.rail_extra_phases[0])
            states = np.exp(1j * sign * theta) * (np.cos(theta) * coined - 1j * sign * np.sin(theta) * route(coined))
    weights = np.sum(np.abs(states) ** 2, axis=(1, 2, 3, 4, 5)) if len(states) else np.zeros(0)
    weight = max(float(weights.sum()), 1e-300)
    birth = float(np.sum(sites * weights) / weight) if len(states) else size / 2
    offset = ((np.arange(size) - birth + size / 2) % size) - size / 2
    ahead = lambda d: float(d[offset > 0].sum() / max(d.sum(), 1e-300))
    markers = [[int(z), int(m), float(w)] for z, m, w in zip(sites, births, weights) if w > 1e-10]
    return {"law": law, "size": size, "updates": updates, "orders": list(orders), "angle": angle,
            "contact_phase": phi, "joint": np.array(joint), "bound": np.array(bound),
            "collected": np.array(collected), "captured": captured, "clicks": clicks,
            "markers": markers, "birth_centre": float(birth),
            "forward_fraction": ahead(bound[-1]), "forward_fraction_collected": ahead(collected[-1]),
            "dropped_weight": dropped, "seconds": time.time() - started}
