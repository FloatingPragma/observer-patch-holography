"""Exact tetrahedral address-bit law and its total reversible event extension.

The directed-address rule comes from candidate A-address-bit-v1 at revision
0c738545b10b655c72cf2f877ec55b983592f535.
The tetrahedron and load orientation follow that revision's mass-free link
study. This module declares a finite computational construction; the squared
load functional has no supplied physical energy or particle-mass dictionary.
"""

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations


EDGES = tuple(combinations(range(4), 2))
ADDRESSES = tuple(direction for edge in EDGES for direction in (edge, edge[::-1]))
ADDRESS_COUNT = len(ADDRESSES)
STATE_COUNT = 1 << ADDRESS_COUNT
BACKGROUNDS = ((0, 0, 0, 0), (1, 0, 0, 0))
LAWS = ("native", "hold", "uphill")
PRIMARY_REVISION = "0c738545b10b655c72cf2f877ec55b983592f535"
PRIMARY_PATH = "candidates/A-address-bit-v1/a_address_bit_v1.py"


def _integer(value: int, name: str, lower: int, upper: int) -> int:
    if type(value) is not int or not lower <= value < upper:
        raise ValueError(f"{name} must be an integer in [{lower}, {upper})")
    return value


def _background(background: tuple[int, ...]) -> tuple[int, ...]:
    if (
        type(background) is not tuple
        or len(background) != 4
        or any(type(value) is not int for value in background)
        or background not in BACKGROUNDS
    ):
        raise ValueError("background must be one of the two declared integer tuples")
    return background


def _law(law: str) -> str:
    if type(law) is not str or law not in LAWS:
        raise ValueError(f"law must be one of {LAWS}")
    return law


def loads(state: int, background: tuple[int, ...] = BACKGROUNDS[0]) -> tuple[int, ...]:
    """Return in-bit minus out-bit counts plus the declared background.

    For an edge (a,b), flux is bit(b,a) - bit(a,b); the incidence column
    is +1 at a and -1 at b. Both neutral bit pairs (00 and 11) are retained.
    """
    _integer(state, "state", 0, STATE_COUNT)
    _background(background)
    result = list(background)
    for address, (tail, head) in enumerate(ADDRESSES):
        bit = (state >> address) & 1
        result[tail] -= bit
        result[head] += bit
    return tuple(result)


@dataclass(frozen=True)
class Event:
    target: int
    old_bit: int
    new_bit: int
    current: int
    before_loads: tuple[int, ...]
    after_loads: tuple[int, ...]
    energy_fall: int


def native_event(
    state: int,
    address: int,
    background: tuple[int, ...] = BACKGROUNDS[0],
    law: str = "native",
) -> Event:
    """Evaluate one event, retaining no-operation outcomes.

    ``hold`` suppresses every current. ``uphill`` reverses the strict load
    comparison and uses the same carrier, address and record apparatus.
    ``energy_fall`` is the decrease of the squared-load diagnostic only.
    """
    _integer(state, "state", 0, STATE_COUNT)
    _integer(address, "address", 0, ADDRESS_COUNT)
    _background(background)
    _law(law)
    before = loads(state, background)
    tail, head = ADDRESSES[address]
    bit = (state >> address) & 1
    difference = before[tail] - before[head]
    if law == "uphill":
        difference = -difference
    current = 0
    if law != "hold":
        if bit == 0 and difference > 0:
            current = 1
        elif bit == 1 and difference < 0:
            current = -1
    target = state ^ (1 << address) if current else state
    after = list(before)
    after[tail] -= current
    after[head] += current
    return Event(
        target=target,
        old_bit=bit,
        new_bit=bit + current,
        current=current,
        before_loads=before,
        after_loads=tuple(after),
        energy_fall=sum(q * q for q in before) - sum(q * q for q in after),
    )


def forward(
    state: int,
    address: int,
    record: int,
    background: tuple[int, ...] = BACKGROUNDS[0],
    law: str = "native",
) -> tuple[int, int]:
    """Apply W_k(s,b,r)=(s,r XOR g_k(s,b),b), including dirty r.

    The non-address bits s do not change. With a blank record r=0, the
    source equals the native event output and the record stores its old bit.
    A dirty record is part of the total unitary domain, with explicit action.
    The address is retained by the enclosing address register.
    """
    _integer(record, "record", 0, 2)
    event = native_event(state, address, background, law)
    target = event.target ^ (record << address)
    return target, event.old_bit


def inverse(
    target: int,
    address: int,
    record: int,
    background: tuple[int, ...] = BACKGROUNDS[0],
    law: str = "native",
) -> tuple[int, int]:
    """Invert W on every source and record basis state without a reset."""
    _integer(target, "target", 0, STATE_COUNT)
    _integer(address, "address", 0, ADDRESS_COUNT)
    _integer(record, "record", 0, 2)
    _background(background)
    _law(law)
    target_bit = (target >> address) & 1
    source = (target & ~(1 << address)) | (record << address)
    event = native_event(source, address, background, law)
    old_record = target_bit ^ event.new_bit
    return source, old_record


def decode_history(
    endpoint: int,
    addresses: tuple[int, ...],
    old_bits: tuple[int, ...],
    background: tuple[int, ...] = BACKGROUNDS[0],
    law: str = "native",
) -> tuple[int, tuple[tuple[int, int, int], ...]]:
    """Recover the initial word and primary event records from terminal data.

    Each returned record is (seam, coin, seam-oriented current), including
    zero currents. Reverse bit restoration uses only the endpoint and the
    retained address/old-bit prefix. Every reconstructed event must satisfy
    the declared law; an impossible prefix raises ValueError. The native
    interpretation requires fresh blank event memory during preparation.
    """
    _integer(endpoint, "endpoint", 0, STATE_COUNT)
    _background(background)
    _law(law)
    if type(addresses) is not tuple or type(old_bits) is not tuple:
        raise ValueError("addresses and old_bits must be tuples")
    if len(addresses) != len(old_bits):
        raise ValueError("addresses and old_bits must have the same length")
    for address, old_bit in zip(addresses, old_bits):
        _integer(address, "address", 0, ADDRESS_COUNT)
        _integer(old_bit, "old_bit", 0, 2)
    current = endpoint
    records = []
    for index in range(len(addresses) - 1, -1, -1):
        address, old_bit = addresses[index], old_bits[index]
        previous = (current & ~(1 << address)) | (old_bit << address)
        event = native_event(previous, address, background, law)
        if event.target != current:
            raise ValueError(f"impossible declared event history at index {index}")
        seam, coin = divmod(address, 2)
        records.append((seam, coin, -event.current if coin else event.current))
        current = previous
    return current, tuple(reversed(records))


def dependency_addresses(address: int) -> tuple[int, ...]:
    """The event reads only directed bits in the two endpoint stars.

    This support statement does not supply a bounded-degree gate synthesis
    or a spatially local implementation of the enclosing scheduler.
    """
    _integer(address, "address", 0, ADDRESS_COUNT)
    endpoints = frozenset(ADDRESSES[address])
    return tuple(
        index for index, edge in enumerate(ADDRESSES) if endpoints.intersection(edge)
    )


def omitted_record_overlap(
    left: int,
    right: int,
    background: tuple[int, ...] = BACKGROUNDS[0],
    law: str = "native",
) -> Fraction:
    """Gram entry with uniform retained addresses but omitted old-bit record."""
    _integer(left, "left", 0, STATE_COUNT)
    _integer(right, "right", 0, STATE_COUNT)
    _background(background)
    _law(law)
    matches = sum(
        native_event(left, address, background, law).target
        == native_event(right, address, background, law).target
        for address in range(ADDRESS_COUNT)
    )
    return Fraction(matches, ADDRESS_COUNT)


def no_old_bit_collision() -> dict:
    """Explicit off-diagonal Gram witness for the omitted old-bit record."""
    return {
        "background": BACKGROUNDS[0],
        "left": 0,
        "right": 1,
        "address": 0,
        "target": 0,
        "overlap": omitted_record_overlap(0, 1),
    }


def reused_record_counterexample() -> dict:
    """Two events on one record fail to reproduce two fresh-record events.

    Source bit 0 clears in the first event, leaving old bit 1 in memory.
    Reusing that bit restores the source bit during the second event.
    The total map is invertible; the repeated native-event interpretation fails.
    """
    source, address, background = 1, 0, BACKGROUNDS[0]
    first_source, first_record = forward(source, address, 0, background)
    reused_source, reused_record = forward(first_source, address, first_record, background)
    fresh_source, fresh_record = forward(first_source, address, 0, background)
    return {
        "source": source,
        "address": address,
        "background": background,
        "first": (first_source, first_record),
        "reused": (reused_source, reused_record),
        "fresh": (fresh_source, fresh_record),
    }
