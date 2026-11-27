"""A serial, sharp, two-rail channel between enrolled observer owners."""

import hashlib
import json
from functools import cache

from model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.face_frame.shared_report import (
    certify,
)
from ..boolean_nodes.material import Machine, balanced_pointer


def encode(packet):
    return json.dumps(packet, separators=(",", ":"), sort_keys=True).encode("ascii")


def identifier(packet):
    return hashlib.sha256(encode(packet)).hexdigest()


@cache
def layout():
    return Machine(0, 6)


ternary_pointer = cache(certify)


def read_rail(machine, loads, rail, standard, pointer):
    a, b = pointer
    route = machine.move(rail, a) + machine.move(rail if standard == a else standard, b)
    for i, j in route:
        loads[i], loads[j] = loads[j], loads[i]
    value, reference = loads[a], loads[b]
    click = int(ternary_pointer()[value, reference, False])
    if value != 1 and balanced_pointer("native")[value, reference] != click:
        raise ValueError("Native pointer instruments disagree")
    for i, j in reversed(route):
        loads[i], loads[j] = loads[j], loads[i]
    return click


@cache
def link(sender, receiver, parked=False):
    machine = layout()
    send = machine.rows[sender][0]
    take, _, standard, pointer = machine.rows[receiver]
    path = tuple(s for a, b in zip(send, take) for s in machine.move(a, b))
    applied = [s for s in path if not parked or s[0] // 12 == s[1] // 12]
    raw = []
    for bit in (0, 1):
        bits = [0] * 6
        bits[sender] = bit
        loads = machine.fresh(bits)
        before = loads[:]
        for a, b in applied:
            loads[a], loads[b] = loads[b], loads[a]
        pair = [read_rail(machine, loads, rail, standard[0], pointer) for rail in take]
        for a, b in reversed(applied):
            loads[a], loads[b] = loads[b], loads[a]
        if loads != before or (not parked and sum(pair) != 1):
            raise ValueError("Unrestored transport or illegal balanced word")
        raw.append(pair)
    return {"source": list(send), "destination": list(take), "path": path, "raw": raw}


class Wire:
    def __init__(self):
        self.packets = {}
        self.events = []
        self.links = {}

    def send(self, sender, receiver, packet, parked=False):
        # Every bit slot starts on a fresh sharp dimer. Exhausting the two basis
        # words therefore gives its entire record law, including both pointers.
        key = f"{sender}-{receiver}-{int(parked)}"
        route = link(sender, receiver, parked)
        self.links[key] = route
        word = encode(packet)
        sent_id = identifier(packet)
        self.packets[sent_id] = packet
        # Received bytes and both full raw pointer streams are losslessly
        # reconstructed by the deposited input packet and this route's table.
        self.events.append([sender, receiver, key, sent_id, len(word) * 8])
        decoded = {(1, 0): 0, (0, 1): 1}
        bits = [decoded.get(tuple(row)) for row in route["raw"]]
        if None in bits:
            return None
        byte_map = bytes(
            sum(bits[(v >> i) & 1] << i for i in range(8)) for v in range(256)
        )
        received = word.translate(byte_map)
        try:
            return json.loads(received)
        except (ValueError, UnicodeError):
            return None


class LockMemory:
    """A persistent local scratch dimer selects a certificate in the owner's tape."""

    def __init__(self, i, loads, machine=None):
        self.i = i
        self.loads = loads
        self.machine = layout() if machine is None else machine
        self.tape = {}
        self.events = []

    def read(self, at):
        m = self.machine
        pair = m.rows[self.i][1]
        _, _, standard, head = m.rows[at]
        raw = [read_rail(m, self.loads, port, standard[0], head) for port in pair]
        value = {(1, 0): 0, (0, 1): 1}[tuple(raw)]
        record = None if not self.tape else self.tape[value]
        self.events.append(
            ["read", at, raw, None if record is None else identifier(record)]
        )
        return value, record

    def install(self, cert):
        old, _ = self.read(self.i)
        value = cert[2]
        if value != old:
            a, b = self.machine.rows[self.i][1]
            self.loads[a], self.loads[b] = self.loads[b], self.loads[a]
        self.tape[value] = cert
        self.events.append(["write", self.i, value, identifier(cert)])
        _, output = self.read(self.i)
        if output != cert:
            raise ValueError("Durable readback failed before commit")
        return output
