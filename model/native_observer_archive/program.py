"""Four local permutation gates on a complete finite observer carrier.

Native events use the unchanged directed-bit source from native_memory_clock.
Pointer writes copy an orthogonal bit label into finite supplied memory. All
gates also act on dirty memory; only the stated blank preparation interprets
the pointer's second bit as a valid-record flag.
"""

from dataclasses import dataclass

from model.native_memory_clock.native import BACKGROUNDS, forward, inverse


PROGRAM_EDGES = 4
REGISTER_CELLS = {
    "S0": 0, "A0": 0, "b1": 0,
    "S1": 1, "A1": 1, "C1": 1, "PA": 1,
    "S2": 2, "A2": 2, "C2": 2,
    "S3": 3, "A3": 3, "C3": 3, "PB": 3, "b2": 3, "D": 3,
    "S4": 4,
}
_SUPPORTS = (
    ("S0", "S1", "A0", "A1", "b1"),
    ("S1", "S2", "A1", "A2", "C1", "C2", "PA"),
    ("S2", "S3", "A2", "A3", "C2", "C3"),
    ("S3", "S4", "C3", "PB", "D", "b2"),
)


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _integer(value, name, upper):
    _require(type(value) is int and 0 <= value < upper,
             f"{name} must be an integer in [0, {upper})")


@dataclass(frozen=True)
class Program:
    background: tuple = (0, 0, 0, 0)
    law: str = "native"
    communication: bool = True
    feedback: bool = True
    write_a: bool = True

    def __post_init__(self):
        _require(type(self.background) is tuple and len(self.background) == 4
                 and all(type(x) is int for x in self.background)
                 and self.background in BACKGROUNDS, "declared integer background required")
        _require(type(self.law) is str and self.law in ("native", "uphill"),
                 "law must be native or uphill")
        for name in ("communication", "feedback", "write_a"):
            _require(type(getattr(self, name)) is bool, name + " must be boolean")


@dataclass(frozen=True, order=True)
class Data:
    buffers: tuple
    addresses: tuple
    couriers: tuple
    pointer_a: int
    pointer_b: int
    old_bits: tuple
    second_address: int


def validate_data(data):
    _require(type(data) is Data, "Data required")
    for values, size, upper, name in (
        (data.buffers, 5, 4096, "buffers"),
        (data.addresses, 4, 12, "addresses"),
        (data.couriers, 3, 4, "couriers"),
        (data.old_bits, 2, 2, "old_bits"),
    ):
        _require(type(values) is tuple and len(values) == size
                 and all(type(x) is int and 0 <= x < upper for x in values), name)
    _integer(data.pointer_a, "pointer_a", 4)
    _integer(data.pointer_b, "pointer_b", 4)
    _integer(data.second_address, "second_address", 12)


def preparation(source, address):
    """One basis component of the uniform first-address preparation."""
    _integer(source, "source", 4096)
    _integer(address, "address", 12)
    return Data((source, 0, 0, 0, 0), (address, 0, 0, 0), (0, 0, 0),
                0, 0, (0, 0), 0)


def gate_support(j):
    """Explicit data support, contained in neighboring cells j and j+1."""
    _integer(j, "program edge", PROGRAM_EDGES)
    return _SUPPORTS[j]


def resources():
    """The complete data carrier, independent of preparation and controls."""
    dimension = 2 ** 72 * 12 ** 5
    return {
        "cells": 5, "program_edges": PROGRAM_EDGES,
        "source_buffer_qubits": 60, "courier_qubits": 6,
        "observer_pointer_qubits": 4, "old_bit_qubits": 2,
        "address_qudits": 5, "address_dimension": 12,
        "full_data_dimension": str(dimension),
        "register_cells": dict(REGISTER_CELLS),
        "gate_supports": [list(support) for support in _SUPPORTS],
    }


def edge(config, data, j, direction=1):
    """Apply an exact gate or its inverse, on every dirty basis state.

    The two pointer bits are value (bit 0) and valid (bit 1). The feedback
    reads the actual receiver pointer after its courier XOR. The second
    address starts at zero in preparation, so its controlled half-turn is
    a nonuniform record-dependent schedule, not a fresh uniform draw.
    """
    _require(type(config) is Program, "Program required")
    validate_data(data)
    _integer(j, "program edge", PROGRAM_EDGES)
    _require(type(direction) is int and direction in (-1, 1), "direction must be -1 or 1")
    s, a, courier, old = (list(data.buffers), list(data.addresses),
                           list(data.couriers), list(data.old_bits))
    pa, pb, d = data.pointer_a, data.pointer_b, data.second_address
    if direction == 1:
        if j == 0:
            s[0], old[0] = forward(s[0], a[0], old[0], config.background, config.law)
            s[0], s[1] = s[1], s[0]
            a[0], a[1] = a[1], a[0]
        elif j == 1:
            if config.write_a:
                pa ^= 2 | ((s[1] >> a[1]) & 1)
            if config.communication:
                courier[0] ^= pa
            s[1], s[2] = s[2], s[1]
            a[1], a[2] = a[2], a[1]
            courier[0], courier[1] = courier[1], courier[0]
        elif j == 2:
            s[2], s[3] = s[3], s[2]
            a[2], a[3] = a[3], a[2]
            courier[1], courier[2] = courier[2], courier[1]
        else:
            pb ^= courier[2]
            if config.feedback and pb & 2:
                d = (d + 6 * (pb & 1)) % 12
            s[3], old[1] = forward(s[3], d, old[1], config.background, config.law)
            s[3], s[4] = s[4], s[3]
    else:
        if j == 0:
            a[0], a[1] = a[1], a[0]
            s[0], s[1] = s[1], s[0]
            s[0], old[0] = inverse(s[0], a[0], old[0], config.background, config.law)
        elif j == 1:
            courier[0], courier[1] = courier[1], courier[0]
            a[1], a[2] = a[2], a[1]
            s[1], s[2] = s[2], s[1]
            if config.communication:
                courier[0] ^= pa
            if config.write_a:
                pa ^= 2 | ((s[1] >> a[1]) & 1)
        elif j == 2:
            courier[1], courier[2] = courier[2], courier[1]
            a[2], a[3] = a[3], a[2]
            s[2], s[3] = s[3], s[2]
        else:
            s[3], s[4] = s[4], s[3]
            s[3], old[1] = inverse(s[3], d, old[1], config.background, config.law)
            if config.feedback and pb & 2:
                d = (d - 6 * (pb & 1)) % 12
            pb ^= courier[2]
    return Data(tuple(s), tuple(a), tuple(courier), pa, pb, tuple(old), d)


def gate_column(config, data, j, direction=1):
    """The complete column of a permutation gate, with its sole unit entry."""
    return {edge(config, data, j, direction): 1}


def prefixes(config, data):
    """All five exact program prefixes, including the unchanged input."""
    _require(type(config) is Program, "Program required")
    validate_data(data)
    out = [data]
    for j in range(PROGRAM_EDGES):
        out.append(edge(config, out[-1], j))
    return tuple(out)


def public_record(data):
    """Complete basis projection of actual registers; no clock inference.

    Order: A-valid, A-value, B-valid, B-value, second address, all source
    buffers, first-address registers, couriers and native old-bit memories.
    Invalid pointer values are retained. Reading every register is a declared
    complete diagnostic instrument, separate from either observer's access.
    """
    validate_data(data)
    return (data.pointer_a >> 1, data.pointer_a & 1,
            data.pointer_b >> 1, data.pointer_b & 1, data.second_address,
            data.buffers, data.addresses, data.couriers, data.old_bits)
