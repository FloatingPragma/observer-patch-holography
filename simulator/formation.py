"""Step 5: forming a pair through a capture vertex that exports a marker.

The successor law gives every junction a periodic wire of neutral modes. When the
two constituents meet head-on at a junction, with constituent A in slot 5 (arrived
moving +z) and B in slot 4 (arrived moving -z), the vertex rotates that channel by
the angle pi/4 into the matched equal-slot pair sum_s e_s |s, s> plus one marker at
the start of the junction's wire. The marker moves one link per update. After the
vertex, both constituents take the step-4 contact update.

A marker that has left its junction never returns within the run (the wire is longer
than the run), so every capture event starts its own branch, orthogonal to all
others. A branch lives in the invariant bound sector, where the pair is one walker
(``step_bound``). This makes the evolution exact at the cost of one small array per
capture event. ``test_simulator.py`` compares it with the literal update that carries
the marker wires explicitly.
"""
from __future__ import annotations

import time

import numpy as np

from .law import EQUAL, LAWS, step_bound, step_contact
from .spectral import mobile_band

ANGLE = np.pi / 4
COLLISION = (5, 4)


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


def capture(free: np.ndarray, angle: float = ANGLE):
    """Apply the vertex to the unmarked two-constituent state; return it and the new branches.

    Only the unmarked sector carries the collision channel, and no marker sits at the
    start of any wire before the vertex acts, so the rotation acts on |5,4; 0> alone.
    """
    size = free.shape[2]
    free = free.copy()
    born = []
    c, s = np.cos(angle), np.sin(angle)
    for z in range(size):
        a = free[0, 0, z, 0, COLLISION[0], 0, 0, z, 0, COLLISION[1]]
        if a == 0:
            continue
        free[0, 0, z, 0, COLLISION[0], 0, 0, z, 0, COLLISION[1]] = c * a
        chi = np.zeros((1, 1, size, 1, 8), complex)
        chi[0, 0, z, 0, :] = s * a * EQUAL
        born.append((z, chi))
    return free, born


def formation_run(size: int = 48, updates: int = 72, orders=(7, -5), contact_phase=None,
                  angle: float = ANGLE, support: int = 4):
    """Two constituents launched toward each other; capture events, the formed pair, its markers."""
    started = time.time()
    law = LAWS["compact-pair-capture"]
    phi = law.rail_contact_phases[0] if contact_phase is None else contact_phase
    a = _strip_packet(size, size // 4, orders[0], support=support)
    b = _strip_packet(size, 3 * size // 4, orders[1], support=support)
    free = np.einsum("abcde,fghij->abcdefghij", a, b)
    branches = []          # [site, emission update, state]
    joint, bound, collected, captured, clicks = [], [], [], [], []
    for n in range(updates + 1):
        density = np.sum(np.abs(free[0, 0, :, 0, :, 0, 0, :, 0, :]) ** 2, axis=(1, 3))
        joint.append(density)
        pair = np.zeros(size)
        by_time = {}
        for z, born_at, chi in branches:
            pair += np.sum(np.abs(chi[0, 0, :, 0, :]) ** 2, axis=-1)
            by_time[born_at] = by_time.get(born_at, 0) + chi[0, 0, :, 0, :]
        # A coherent collector joins the marker wires of all junctions; its common port
        # leaves the birthplace unresolved. Conditioned on that port, markers emitted at
        # the same update add coherently with amplitude 1 / sqrt(number of junctions).
        common = sum(np.sum(np.abs(v) ** 2, axis=-1) for v in by_time.values()) / size if by_time else np.zeros(size)
        bound.append(pair)
        collected.append(common)
        captured.append(float(pair.sum()))
        clicks.append(float(np.sum(common)))
        if n == updates:
            break
        free, born = capture(free, angle)
        branches.extend([z, n, chi] for z, chi in born)
        free = step_contact(free, phi)
        for branch in branches:
            branch[2] = step_bound(branch[2], phi)
    weight = max(captured[-1], 1e-300)
    birth = sum(z * float(np.sum(np.abs(chi) ** 2)) for z, _, chi in branches) / weight if branches else size / 2
    zs = np.arange(size)
    ahead = lambda d: float(d[zs > birth].sum() / max(d.sum(), 1e-300))
    markers = [[int(z), int(n), float(np.sum(np.abs(chi) ** 2))] for z, n, chi in branches
               if np.sum(np.abs(chi) ** 2) > 1e-16]
    return {"law": law, "size": size, "updates": updates, "orders": list(orders), "angle": angle,
            "contact_phase": phi, "joint": np.array(joint), "bound": np.array(bound),
            "collected": np.array(collected), "captured": captured, "clicks": clicks,
            "markers": markers, "birth_centre": float(birth),
            "forward_fraction": ahead(bound[-1]), "forward_fraction_collected": ahead(collected[-1]),
            "seconds": time.time() - started}
