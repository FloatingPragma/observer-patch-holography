"""The eight-slot torus law: carrier, local update, and the controls used in its experiments.

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
    """The eight-slot torus law, or a control that changes one of its parts."""

    key: str
    title: str
    equal_phase: complex = 1.0
    encounter_phase: float = 0.0
    access: str = "field-only"
    status: str = "law"
    note: str = ""
    rail_contact_phases: tuple = ()
    rail_extra_phases: tuple = ()
    capture: bool = False
    maturation: tuple = ()

    def __post_init__(self):
        if abs(abs(self.equal_phase) - 1.0) > 1e-12:
            raise ValueError("the equal-slot phase must have modulus one for a unitary coin")
        if self.access not in ("field-only", "equipped"):
            raise ValueError("access is 'field-only' or 'equipped'")


# Contact couplings of the torus law: the two rails of one constituent
# meet the other constituent with contact phases PHI0 and PHI0 + DELTA.
PHI0 = np.pi / 16
DELTA = np.pi / 2048

LAWS = {
    "torus": Law(
        "torus", "The eight-slot torus law",
        status="law",
        rail_contact_phases=(PHI0, PHI0 + DELTA), rail_extra_phases=(0.0, 0.0), capture=True,
        maturation=(-1, "cosine", 128),
        note="Each constituent moves by reflected routing. Two constituents that share an address in "
             "matched slots step together under the contact coin, and the relative phase of the two "
             "rails is the pair's clock. A head-on encounter forms a pair through a reversible capture "
             "that sends a marker away on the junction's wire, and the newborn pair opens its hop on a "
             "cosine profile over 128 updates."),
    # Changed laws used as controls in the experiments. Each changes one part of the torus law.
    "control-encounter": Law(
        "control-encounter", "Control: a scalar phase when two constituents meet",
        encounter_phase=np.pi / 2, status="control", rail_contact_phases=(PHI0, PHI0 + DELTA),
        rail_extra_phases=(0.0, 0.0),
        note="Two constituents at one address gain the phase chi per update in place of the contact coin."),
    "control-relay": Law(
        "control-relay", "Control: an addressed relay available to inhabitants",
        access="equipped", status="control", rail_contact_phases=(PHI0, PHI0 + DELTA),
        rail_extra_phases=(0.0, 0.0),
        note="An addressed relay carries a chosen value one address per update, ahead of the smooth cone."),
    "control-phase": Law(
        "control-phase", "Control: equal-slot phase changed by 0.05 rad",
        equal_phase=np.exp(0.05j), status="control",
        note="A fixed nonzero equal-slot phase makes the gapless branch quadratic near zero momentum."),
    "control-phase-i": Law(
        "control-phase-i", "Control: equal-slot phase i",
        equal_phase=1j, status="control",
        note="The equal-slot reflection replaced by the phase i."),
    "control-spectator": Law(
        "control-spectator", "Control: spectator rails",
        status="control",
        rail_contact_phases=(PHI0 + DELTA / 2, PHI0 + DELTA / 2),
        rail_extra_phases=(-DELTA / 4, DELTA / 4),
        note="Equal contact coins and an extra relative phase per update: a clock that counts services."),
    "control-abrupt": Law(
        "control-abrupt", "Control: a newborn pair mobile at once",
        status="control",
        rail_contact_phases=(PHI0, PHI0 + DELTA), rail_extra_phases=(0.0, 0.0), capture=True,
        note="The capture of the torus law without the opening: the newborn pair takes the full hop at once."),
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
    """Mobile-band phase per update of one constituent: cos E = (sum of the four cosines) / 4."""
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


def contact_coin(states: np.ndarray, contact_phase: float, axis: int = -1) -> np.ndarray:
    """D = exp(i phi/2) [(1 + exp(-i phi)) |e><e| - I] on the eight matched pair states."""
    return np.exp(0.5j * contact_phase) * coin(states, np.exp(-1j * contact_phase), axis=axis)


def step_contact(psi2: np.ndarray, contact_phase: float, extra_phase: float = 0.0) -> np.ndarray:
    """One update of two constituents under the contact law, for one rail.

    Off contact each constituent applies the parent reflection. Where both occupy the
    same address and private position, the eight matched slot pairs (s, s) receive the
    contact coin and the phase exp(-i extra); the 56 unmatched slot pairs there are left
    unchanged. Both constituents then route.
    """
    if psi2.ndim != 10 or psi2.shape[4] != SLOTS or psi2.shape[9] != SLOTS:
        raise ValueError("two constituents have shape (N1, N2, N3, P, 8) twice")
    carrier = psi2.shape[:4]
    volume = int(np.prod(carrier))
    coined = coin(coin(psi2, 1.0, axis=4), 1.0, axis=9).reshape(volume, SLOTS, volume, SLOTS)
    flat = psi2.reshape(volume, SLOTS, volume, SLOTS)
    here = np.arange(volume)
    contact = flat[here, :, here, :].copy()
    slots = np.arange(SLOTS)
    matched = contact[:, slots, slots]
    contact[:, slots, slots] = contact_coin(matched, contact_phase) * np.exp(-1j * extra_phase)
    coined[here, :, here, :] = contact
    coined = coined.reshape(psi2.shape)
    return route(route(coined, slot_axis=4), slot_axis=9)


def step_bound(chi: np.ndarray, contact_phase: float, extra_phase: float = 0.0) -> np.ndarray:
    """The same update on the bound sector, where both constituents share address and slot.

    The matched coincident states form an invariant subspace; on it the pair is one walker
    with the contact coin. ``test_simulator.py`` checks this against ``step_contact``.
    """
    return route(contact_coin(chi, contact_phase)) * np.exp(-1j * extra_phase)


def hop_angle(age: int, maturation: tuple) -> float:
    """Hop angle of a pair of the given age: (pi/2) f((age + 1) / M) while immature, pi/2 after."""
    if not maturation:
        return np.pi / 2
    _, profile, duration = maturation
    if age >= duration:
        return np.pi / 2
    u = (age + 1) / duration
    f = u if profile == "linear" else (1 - np.cos(np.pi * u)) / 2
    return np.pi / 2 * f


def partial_hop(chi: np.ndarray, theta: float, sign: int) -> np.ndarray:
    """exp(i s theta) [cos(theta) I - i s sin(theta) R]; at theta = pi/2 this is the routing R."""
    return np.exp(1j * sign * theta) * (np.cos(theta) * chi - 1j * sign * np.sin(theta) * route(chi))


def step_maturing(chi: np.ndarray, law: "Law", age: int, rail: int = 0) -> np.ndarray:
    """One update of a bound pair of the given age: contact coin, then the age's partial hop."""
    phi = law.rail_contact_phases[rail]
    coined = contact_coin(chi, phi) * np.exp(-1j * law.rail_extra_phases[rail])
    sign = law.maturation[0] if law.maturation else -1
    return partial_hop(coined, hop_angle(age, law.maturation), sign)
