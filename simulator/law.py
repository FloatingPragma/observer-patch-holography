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
    rail_contact_phases: tuple = ()
    rail_extra_phases: tuple = ()
    capture: bool = False
    maturation: tuple = ()

    def __post_init__(self):
        if abs(abs(self.equal_phase) - 1.0) > 1e-12:
            raise ValueError("the equal-slot phase must have modulus one for a unitary coin")
        if self.access not in ("field-only", "equipped"):
            raise ValueError("access is 'field-only' or 'equipped'")


# Contact couplings of the step-4 pair laws: the two rails of one constituent
# meet the other constituent with contact phases PHI0 and PHI0 + DELTA.
PHI0 = np.pi / 16
DELTA = np.pi / 2048

LAWS = {
    "reflected-routing": Law(
        "reflected-routing", "Reflected routing, field-only access",
        status="component",
        note="The one-constituent law of all three candidates. Long-wave clocks keep proper time, "
             "the field has a 3+1D Klein-Gordon limit with c = 1/2 (steps 1 to 3), and a compact "
             "free-wave launch carries no message ahead of the smooth cone at the tested budgets."),
    "reflected-routing-encounter": Law(
        "reflected-routing-encounter", "Reflected routing with an encounter phase",
        encounter_phase=np.pi / 2, status="undecided",
        note="Identical to reflected routing for one excitation. Two excitations at one "
             "spatial address gain the phase chi per update. A two-excitation interferometer "
             "tells the values of chi apart (step 4). As a bound clock, a contact phase i on "
             "an axial guide missed its forecasts and captured no separated packets."),
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
    "compact-pair": Law(
        "compact-pair", "Compact pair with motion-coupled rails",
        status="component",
        rail_contact_phases=(PHI0, PHI0 + DELTA), rail_extra_phases=(0.0, 0.0),
        note="Reflected routing for each constituent, and a contact coin on the eight matched "
             "pair states when both share an address, so a matched pair steps together. Its "
             "moving rail clock keeps proper time (step 4), and the finite law predicts its "
             "deviation from the smooth metric. The bound subspace is invariant, so separated "
             "constituents never form the pair."),
    "compact-pair-capture": Law(
        "compact-pair-capture", "Capture vertex, pair mobile at once (comparator)",
        status="excluded",
        rail_contact_phases=(PHI0, PHI0 + DELTA), rail_extra_phases=(0.0, 0.0), capture=True,
        note="The compact-pair contact law plus a reversible vertex at every junction: a head-on "
             "collision turns into the matched pair and a marker that leaves on the junction's own "
             "wire. It forms persistent pairs (step 5), but the newborn pair, mobile at once, leaves "
             "in both directions, and it failed the direction test of step 6."),
    "maturing-cosine-32": Law(
        "maturing-cosine-32", "Maturing pair, minus cosine over 32 updates",
        status="candidate",
        rail_contact_phases=(PHI0, PHI0 + DELTA), rail_extra_phases=(0.0, 0.0), capture=True,
        maturation=(-1, "cosine", 32),
        note="Capture vertex and contact law; the newborn pair carries an age counter and opens its "
             "hop with a reversible partial swap over 32 updates on a smooth cosine profile. Passed "
             "the directional and clock tests of step 6."),
    "maturing-cosine-128": Law(
        "maturing-cosine-128", "Maturing pair, minus cosine over 128 updates",
        status="candidate",
        rail_contact_phases=(PHI0, PHI0 + DELTA), rail_extra_phases=(0.0, 0.0), capture=True,
        maturation=(-1, "cosine", 128),
        note="The same smooth opening spread over 128 updates. Highest directional fraction and "
             "smallest breakup on the step-6 fixtures."),
    "maturing-linear-128": Law(
        "maturing-linear-128", "Maturing pair, minus linear over 128 updates",
        status="candidate",
        rail_contact_phases=(PHI0, PHI0 + DELTA), rail_extra_phases=(0.0, 0.0), capture=True,
        maturation=(-1, "linear", 128),
        note="The hop angle grows by a constant increment over 128 updates. Passed the directional "
             "and clock tests of step 6, with more breakup when the marker returns."),
    "maturing-plus-cosine-32": Law(
        "maturing-plus-cosine-32", "Maturing pair, plus orientation (comparator)",
        status="excluded",
        rail_contact_phases=(PHI0, PHI0 + DELTA), rail_extra_phases=(0.0, 0.0), capture=True,
        maturation=(+1, "cosine", 32),
        note="The opposite orientation of the partial swap. Among the step-6 comparators that "
             "failed the directional test."),
    "compact-pair-spectator": Law(
        "compact-pair-spectator", "Compact pair with spectator rails (control)",
        status="excluded",
        rail_contact_phases=(PHI0 + DELTA / 2, PHI0 + DELTA / 2),
        rail_extra_phases=(-DELTA / 4, DELTA / 4),
        note="Equal contact coins and an extra relative phase per update: the rail clock "
             "counts services and reads one on every journey of step 4."),
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


def contact_coin(states: np.ndarray, contact_phase: float, axis: int = -1) -> np.ndarray:
    """D = exp(i phi/2) [(1 + exp(-i phi)) |e><e| - I] on the eight matched pair states."""
    return np.exp(0.5j * contact_phase) * coin(states, np.exp(-1j * contact_phase), axis=axis)


def step_contact(psi2: np.ndarray, contact_phase: float, extra_phase: float = 0.0) -> np.ndarray:
    """One update of two constituents under the compact-pair contact law, for one rail.

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
