"""The reflected-routing law family: carrier, local update and its members.

A carrier is a four-axis torus. Axes 0, 1 and 2 are spatial addresses; axis 3
is the private cyclic register of each address. Every position holds eight
slots, paired as (0, 1) on axis 0, (2, 3) on axis 1, (4, 5) on axis 2 and
(6, 7) on axis 3. One update applies the coin at every position and then
routes the amplitude of slot s one position along axis s // 2, forward for
even s and backward for odd s, into the partner slot s ^ 1 of the receiving
position. SPEC.md states the law in full.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

SLOTS = 8
PARTNER = np.arange(SLOTS) ^ 1
SLOT_AXIS = np.arange(SLOTS) // 2
SLOT_SIGN = np.where(np.arange(SLOTS) % 2 == 0, 1, -1)
EQUAL = np.ones(SLOTS, complex) / np.sqrt(SLOTS)


@dataclass(frozen=True)
class Law:
    """One member of the family. ``status`` repeats the journal's standing."""

    key: str
    title: str
    equal_phase: complex = 1.0
    encounter_phase: float = 0.0
    access: str = "field-only"
    status: str = "candidate"
    note: str = ""

    def __post_init__(self):
        if abs(abs(self.equal_phase) - 1.0) > 1e-12:
            raise ValueError("the equal-slot phase must have modulus one for a unitary coin")
        if self.access not in ("field-only", "equipped"):
            raise ValueError("access is 'field-only' or 'equipped'")


LAWS = {
    "reflected-routing": Law(
        "reflected-routing", "Reflected routing, field-only access",
        status="candidate",
        note="Survives steps 1 to 3: long-wave clocks keep proper time, the field has a "
             "3+1D Klein-Gordon limit with c = 1/2, and no inhabitant channel outruns it."),
    "reflected-routing-encounter": Law(
        "reflected-routing-encounter", "Reflected routing with an encounter phase",
        encounter_phase=np.pi / 2, status="undecided",
        note="Identical to reflected routing for one excitation. Two excitations at one "
             "spatial address gain the phase chi per update. Every value of chi passes "
             "steps 1 to 3; a two-excitation experiment decides between them."),
    "equipped-relay": Law(
        "equipped-relay", "Reflected routing with an inhabitant-accessible relay",
        access="equipped", status="excluded",
        note="Excluded in step 2 as a law with one causal speed: the addressed relay "
             "carries a freshly chosen bit at one address per update, twice the speed "
             "of the smooth scalar cone."),
    "phase-deformed": Law(
        "phase-deformed", "Equal-slot phase changed by 0.05 rad (control)",
        equal_phase=np.exp(0.05j), status="disfavoured",
        note="Disfavoured in step 1: a fixed nonzero phase makes the gapless branch "
             "quadratic near zero momentum, and moving clocks read 0.76 where the "
             "candidate reads 0.63 to 0.69 at matched speeds."),
    "phase-i": Law(
        "phase-i", "Equal-slot phase i (control)",
        equal_phase=1j, status="control",
        note="The registered phase control of step 1. Its field records miss the "
             "candidate's forecast by an RMS error of 1.25."),
}


def coin(states: np.ndarray, equal_phase: complex = 1.0, axis: int = -1) -> np.ndarray:
    """(1 + lambda) |e><e| - I on the slot axis: reflection of the equal-slot part for lambda = 1."""
    return (1 + equal_phase) * states.mean(axis=axis, keepdims=True) - states


def route(states: np.ndarray, slot_axis: int = -1) -> np.ndarray:
    """Move slot s one position along axis s // 2 into the partner slot s ^ 1.

    The four carrier axes are the four axes immediately before ``slot_axis``.
    """
    slot_axis %= states.ndim
    first = slot_axis - 4
    if first < 0:
        raise ValueError("a state needs four carrier axes before its slot axis")
    out = np.empty_like(states)
    for slot in range(SLOTS):
        source = np.take(states, slot, axis=slot_axis)
        moved = np.roll(source, int(SLOT_SIGN[slot]), axis=first + int(SLOT_AXIS[slot]))
        index = [slice(None)] * states.ndim
        index[slot_axis] = int(PARTNER[slot])
        out[tuple(index)] = moved
    return out


def step(psi: np.ndarray, law: Law) -> np.ndarray:
    """One update of one excitation, psi with shape (N1, N2, N3, P, 8)."""
    if psi.ndim != 5 or psi.shape[-1] != SLOTS:
        raise ValueError("one excitation has shape (N1, N2, N3, P, 8)")
    return route(coin(psi, law.equal_phase))


def encounter_mask(carrier: tuple[int, int, int, int]) -> np.ndarray:
    """True where two excitations share a spatial address (private positions may differ)."""
    n1, n2, n3, p = carrier
    grids = np.indices((n1, n2, n3))
    same = np.ones((n1, n2, n3, n1, n2, n3), bool)
    for axis in range(3):
        a = grids[axis][:, :, :, None, None, None]
        b = grids[axis][None, None, None, :, :, :]
        same &= a == b
    return same[:, :, :, None, None, :, :, :, None, None]


def step_pair(psi2: np.ndarray, law: Law, mask: np.ndarray | None = None) -> np.ndarray:
    """One update of two excitations with shape (N1, N2, N3, P, 8, N1, N2, N3, P, 8).

    The encounter phase acts first, then the coin and routing act on each
    excitation. With encounter phase zero this is the product of two
    independent one-excitation updates.
    """
    if psi2.ndim != 10 or psi2.shape[4] != SLOTS or psi2.shape[9] != SLOTS:
        raise ValueError("two excitations have shape (N1, N2, N3, P, 8) twice")
    if law.encounter_phase:
        if mask is None:
            mask = encounter_mask(psi2.shape[:4])
        psi2 = np.where(mask, np.exp(1j * law.encounter_phase) * psi2, psi2)
    psi2 = coin(psi2, law.equal_phase, axis=4)
    psi2 = coin(psi2, law.equal_phase, axis=9)
    psi2 = route(psi2, slot_axis=4)
    return route(psi2, slot_axis=9)


def dispersion(angles: np.ndarray) -> np.ndarray:
    """Mobile-band phase per update of the candidate: cos E = (sum of the four cosines) / 4."""
    angles = np.asarray(angles, float)
    return np.arccos(np.clip(np.cos(angles).sum(axis=-1) / 4.0, -1.0, 1.0))


def walk_matrix(angles: np.ndarray, equal_phase: complex = 1.0) -> np.ndarray:
    """The one-update operator on a plane wave exp(i k.x) with wave numbers ``angles`` (4,)."""
    angles = np.asarray(angles, float)
    signed = angles[SLOT_AXIS] * SLOT_SIGN
    coin_matrix = (1 + equal_phase) * np.full((SLOTS, SLOTS), 1 / SLOTS) - np.eye(SLOTS)
    routing = np.zeros((SLOTS, SLOTS), complex)
    for slot in range(SLOTS):
        routing[PARTNER[slot], slot] = np.exp(-1j * signed[slot])
    return routing @ coin_matrix
