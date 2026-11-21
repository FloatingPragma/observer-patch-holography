"""Load an ordinarily born torus pair using the native hops and on-site controls.

Gate tables act on eight matched slots at one site. A table entry at age j
precedes the complete native step j. The entry at 28 is the terminal on-site
rotation. There is no intervention inside a native service or on a link.
"""
from __future__ import annotations

from numbers import Integral

import numpy as np

from .. import formation, law, pair_persistence


SOURCE = law.LAWS["torus"]
LOAD_AGE = 28
IDENTITY = np.eye(law.SLOTS, dtype=complex)
MODES = ("exact", "common", "passive", "phase-cut")


def chart(carrier, root, *, private_minimum=3):
    if len(carrier) != 4 or len(root) != 4 or any(
        isinstance(v, (bool, np.bool_)) or not isinstance(v, Integral) for v in (*carrier, *root)
    ):
        raise ValueError("four integral periods and four integral root coordinates required")
    if any(n < minimum or not 0 <= x < n for n, x, minimum in zip(carrier, root, (3,3,3,private_minimum))):
        raise ValueError("the root and periods must support the declared routes")
    return tuple(map(int, carrier)), tuple(map(int, root))


def neighbour(site, slot, carrier):
    moved = list(site)
    axis = int(law.SLOT_AXIS[slot])
    moved[axis] = (moved[axis] + int(law.SLOT_SIGN[slot])) % carrier[axis]
    return tuple(moved)


def corners(carrier, root):
    """Root, x, y, xy in the positive x and y directions of the carrier."""
    x = neighbour(root, 0, carrier)
    y = neighbour(root, 2, carrier)
    return root, x, y, neighbour(x, 2, carrier)


def circulation(carrier, root):
    """Same signed, routed edge amplitudes as pair_persistence.circulation on this chart."""
    target = np.zeros((*carrier, law.SLOTS), complex)
    template = pair_persistence.circulation(size=3, corner=(0, 0))
    for index, site in enumerate(corners(carrier, root)):
        target[site] = template[index % 2, index // 2, 0, 0]
    return target


def _rotation(start, end):
    # All supplied slot rays here are real, unit vectors. This reflection
    # exactly maps the first ray to the second, with a unitary complement.
    difference = start - end
    norm = np.vdot(difference, difference).real
    if norm == 0:
        return IDENTITY
    return IDENTITY - (2 / norm) * np.outer(difference, difference.conj())


def _pulse(start, end):
    return sum(law.hop_angle(j, SOURCE.maturation) for j in range(start, end))


def _split(a, b, far_power):
    """Two native rotation groups with a phase plate at their common endpoint."""
    u = np.sin(b) * np.cos(a)
    v = np.cos(b) * np.sin(a)
    cosine = (far_power - u * u - v * v) / (2 * u * v)
    if not -1 <= cosine <= 1:
        raise ValueError("the two native pulse angles cannot reach the requested power")
    phase = np.exp(1j * np.arccos(cosine))
    near = np.cos(a) * np.cos(b) - phase * np.sin(a) * np.sin(b)
    far = 1j * (u + phase * v)
    return phase, near.conjugate() / abs(near), far.conjugate() / abs(far)


def schedule(carrier, root, mode="exact"):
    """Return age -> site -> rail-indexed 8x8 column-vector gate, absent = identity."""
    carrier, root = chart(carrier, root)
    if mode not in MODES:
        raise ValueError(f"unknown clock loader {mode!r}")
    if mode == "passive":
        return {}
    star = [neighbour(root, slot, carrier) for slot in range(law.SLOTS)]
    sites = corners(carrier, root)
    loop = circulation(carrier, root)
    c = _pulse(5, 6)
    refocus = _split(_pulse(0, 3), _pulse(3, 5), np.sin(c) ** 2)
    horizontal = _split(_pulse(6, 17), _pulse(17, 22), 0.5)
    vertical = _split(_pulse(22, 25), _pulse(25, LOAD_AGE), 0.5)
    phases = SOURCE.rail_contact_phases
    if mode == "common":
        phases = (np.mean(phases),) * 2
    inverse = np.stack([law.contact_coin(IDENTITY, -phi).T for phi in phases])
    gates = {}
    for j in range(1, LOAD_AGE + 1):
        support = (root, *star) if j < 6 else (sites[:2] if j < 22 else sites)
        layer = {}
        for site in support:
            plate = IDENTITY
            if j == 3 and site != root:
                plate = refocus[0] * IDENTITY
            elif j == 5:
                phase = refocus[1] if site == root else -1j * refocus[2]
                plate = phase * IDENTITY
            elif j == 6 and site == root:
                plate = _rotation(law.EQUAL, IDENTITY[0])
            elif j == 17 and site == sites[1] and mode != "phase-cut":
                plate = horizontal[0] * IDENTITY
            elif j == 22 and site in sites[:2]:
                i = sites.index(site)
                plate = horizontal[i + 1] * _rotation(IDENTITY[i], IDENTITY[2])
            elif j == 25 and site in sites[2:]:
                plate = vertical[0] * IDENTITY
            elif j == LOAD_AGE:
                side = sites.index(site) // 2
                plate = vertical[side + 1] * _rotation(IDENTITY[2 + side], 2 * loop[site])
            layer[site] = inverse @ plate if j < LOAD_AGE else np.stack([plate, plate])
        gates[j] = layer
    return gates


def apply_sites(states, layer):
    """Apply a layer to two rail columns in place. Each gate has one site of support."""
    for site, unitary in layer.items():
        at = (slice(None), *site)
        states[at] = np.einsum("rij,rj->ri", unitary, states[at])


def prepare_clock(carrier, root, mode="exact", *, local_reference=False):
    """Return the loader's rail columns, or its single charge-one Ramsey state.

    With a reference the two components represent (rail 0, memory 1) and
    (rail 1, memory 0). The memory stays at the root and has unit idle evolution.
    """
    carrier, root = chart(carrier, root)
    gates = schedule(carrier, root, mode)
    states = np.zeros((2, *carrier, law.SLOTS), complex)
    states[(slice(None), *root)] = law.EQUAL
    if local_reference:
        states[1] = 0
    for age in range(LOAD_AGE):
        if local_reference and age == 6:
            at = (slice(None), *root)
            exchange = np.array([[1., 1.], [1., -1.]]) / np.sqrt(2)
            states[at] = exchange @ states[at]
        apply_sites(states, gates.get(age, {}))
        for r in range(2):
            states[r] = law.step_maturing(states[r], SOURCE, age, r)
    apply_sites(states, gates.get(LOAD_AGE, {}))
    return states, gates


def ordinary_feeds(carrier, root, focus_a, focus_b):
    """Two independent off-contact sources aimed at the root and its x neighbour."""
    carrier, root = chart(carrier, root, private_minimum=1)
    if any(not np.isfinite(f) or not 0 <= f <= 1 for f in (focus_a, focus_b)):
        raise ValueError("launch fractions in [0,1] required")
    second = neighbour(root, 0, carrier)
    incoming = []
    for f, launch_slot in ((focus_a, 4), (focus_b, 5)):
        state = np.zeros((*carrier, law.SLOTS), complex)
        for site, probability in ((root, f), (second, 1 - f)):
            inlet = neighbour(site, launch_slot ^ 1, carrier)
            state[inlet] = np.sqrt(probability) * law.coin(IDENTITY[launch_slot])
        incoming.append(state)
    return tuple(incoming)


def ordinary_launch(carrier, root, focus_a, focus_b, angle=None):
    """Apply the two literal incoming steps and distinguish all first-wire refusals."""
    carrier, root = chart(carrier, root)
    angle = formation.ANGLE if angle is None else float(angle)
    if not np.isfinite(angle):
        raise ValueError("a finite capture angle is required")
    second = neighbour(root, 0, carrier)
    incoming = [law.step(s, SOURCE) for s in ordinary_feeds(carrier, root, focus_a, focus_b)]
    collision = incoming[0][..., 5] * incoming[1][..., 4]
    marked = np.sin(angle) ** 2 * np.abs(collision) ** 2
    return {"natal": float(marked[root]), "other": float(marked.sum() - marked[root]),
            "unmarked": float(1 - marked.sum())}
