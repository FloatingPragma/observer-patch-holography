"""One position-independent local instrument for all natal positions on the torus law."""
from __future__ import annotations

import numpy as np

from . import clock_loader as loader
from .. import formation, law


END = 35
# After the first block every link has refocused; each other block
# implements the exact reciprocal swap of the two ends of every link.
BLOCKS = ((0, 3, 5, 6), (6, 14, 22, 28), (28, 30, 32, END))


def packet(carrier, origin):
    """Encoded spatial ray for one material birthplace."""
    carrier, origin = loader.chart(carrier, origin, private_minimum=1)
    shape = np.zeros((*carrier, law.SLOTS), complex)
    opposite = loader.corners(carrier, origin)[-1]
    shape[origin][[0, 2]] = (0.5, -0.5)
    shape[opposite][[1, 3]] = (0.5, -0.5)
    return shape


def _edge_phase(even, odd):
    return np.diag([even, odd] * 4)


def _edge_pulses(begin, cut, turn, end, swap):
    """Euler plates for a two-port orientation shared by every routed edge."""
    a, b, c = (loader._pulse(i, j) for i, j in ((begin, cut), (cut, turn), (turn, end)))
    phase, near, far = loader._split(a, b, np.cos(c)**2 if swap else np.sin(c)**2)
    align = (near, (1j if swap else -1j) * far)
    def hop(theta):
        return np.array([[np.cos(theta), 1j*np.sin(theta)], [1j*np.sin(theta), np.cos(theta)]])
    total = hop(c) @ np.diag(align) @ hop(b) @ np.diag([1, phase]) @ hop(a)
    amplitudes = (total[0, 1], total[1, 0]) if swap else (total[0, 0], total[1, 1])
    finish = [v.conjugate() / abs(v) for v in amplitudes]
    return {cut: _edge_phase(1, phase), turn: _edge_phase(*align), end: _edge_phase(*finish)}


def uniform_schedule():
    """Rail blocks G[j] at the same site of every cell, indexed by local bound age."""
    unit = loader.IDENTITY
    plates = {j: unit for j in range(1, END+1)}
    for block in BLOCKS:
        for j, gate in _edge_pulses(*block, swap=block[0] != 0).items():
            plates[j] = gate @ plates[j]
    # +x,+y ends at the opposite corner; +z,-z returns to the root.
    first = (unit[0]+unit[4])/np.sqrt(2)
    plates[6] = loader._rotation(law.EQUAL, first) @ plates[6]
    plates[28] = loader._rotation(unit[1], unit[2]) @ plates[28]
    root_ray = (unit[0]-unit[2])/np.sqrt(2)
    opposite_ray = (unit[1]-unit[3])/np.sqrt(2)
    terminal = loader._rotation(unit[4], root_ray) @ loader._rotation(unit[3], opposite_ray)
    plates[END] = terminal @ plates[END]
    undo = np.stack([law.contact_coin(unit, -phi).T for phi in loader.SOURCE.rail_contact_phases])
    return {j: undo @ plate if j < END else np.stack([plate, plate]) for j, plate in plates.items()}


def material_service(states, age, plate=None):
    """Apply the native update on both rails and the optional resulting-age plate."""
    if (not isinstance(age, (int, np.integer)) or isinstance(age, bool)
            or not 0 <= age < formation.AGE_PERIOD):
        raise ValueError("a native age register value is required")
    states = np.asarray(states)
    if (states.ndim < 6 or states.shape[0] != 2 or states.shape[-1] != law.SLOTS
            or 0 in states.shape or not np.isfinite(states).all()):
        raise ValueError("two finite rail fields on a nonempty four-axis carrier required")
    result = np.stack([law.step_maturing(states[r], loader.SOURCE, age, r) for r in range(2)])
    if plate is not None:
        result = np.einsum("rij,r...j->r...i", plate, result)
    return result


def advance(states, age, gates):
    """One full native service followed by the local plate at the resulting age."""
    if not isinstance(age, (int, np.integer)) or isinstance(age, bool) or not 0 <= age < 8191:
        raise ValueError("one pre-wrap service requires an integer age in [0,8190]")
    return material_service(states,age,gates.get(age+1))
