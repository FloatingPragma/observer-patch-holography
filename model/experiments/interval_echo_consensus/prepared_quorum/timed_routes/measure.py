"""Measure the finite serial quorum calendar and verify its deposited signed packets."""

import argparse
import gzip
import json
from pathlib import Path

from ..certificates import Observer, Rules, Signers
from ..measure import Trial
from ..transport import LockMemory, encode, identifier
from .routes import Mount, inspect_pair, topology

HERE = Path(__file__).parent
BYTE_CAP = 8192
DECODE = {(1, 0): 0, (0, 1): 1}


class ClockWire:
    def __init__(self, mount, members, table):
        self.mount = mount
        self.members = members
        self.table = table
        self.packets = {}
        self.events = []
        self.phase = ""
        self.elapsed = 0

    def send(self, sender, receiver, packet, parked=False):
        route = self.mount.columns(self.members, sender, receiver, parked)
        raw = self.table["parked" if parked else "native"][f"{sender}/{receiver}"]
        if raw["raw"] != route["raw"]:
            raise ValueError("Clocked code differs from route calibration")
        word = encode(packet)
        h = identifier(packet)
        if len(word) > BYTE_CAP:
            raise ValueError("Outside finite byte enrollment")
        self.packets[h] = packet
        mapping = [DECODE.get(tuple(row)) for row in raw["raw"]]
        indices = [word.count(bytes([v])) for v in range(256)]
        ones = sum(v.bit_count() * n for v, n in enumerate(indices))
        pointer_hist = [0] * 4
        for (a, b), count in zip(raw["raw"], (8 * len(word) - ones, ones), strict=True):
            pointer_hist[2 * a + b] += count
        before = self.elapsed
        self.elapsed += 8 * len(word) * raw["slots"]
        self.events.append(
            {
                "phase": self.phase,
                "from": sender,
                "to": receiver,
                "packet": h,
                "bytes": len(word),
                "parked": parked,
                "raw_hist": pointer_hist,
                "start": before,
                "stop": self.elapsed,
            }
        )
        if None in mapping:
            return None
        byte_table = bytes(
            sum(mapping[(v >> i) & 1] << i for i in range(8)) for v in range(256)
        )
        try:
            return json.loads(word.translate(byte_table))
        except (ValueError, UnicodeError):
            return None


class MountedTrial(Trial):
    def __init__(self, mount, auth, wire, stock=None):
        self.wire = wire
        self.rules = Rules(auth, 3)
        self.faulty = set()
        self.log = []
        self.stock = (
            mount.base.fresh([0] * len(mount.base.rows)) if stock is None else stock
        )
        self.nodes = []
        for i, owner in enumerate(mount.homes):
            memory = LockMemory(owner, self.stock, mount.local[owner])
            self.nodes.append(Observer(i, self.rules, self.stock, memory=memory))


def clock_trace(columns, appointments, initial=0):
    """Lossless periodic read stream, including its short prefix and endpoint."""
    if appointments < 1:
        raise ValueError("Clock exposure must be positive")
    word = []
    digit = initial
    while digit not in [x[0] for x in word]:
        response = columns[digit]
        word.append((digit, response))
        digit = DECODE[tuple(response)]
    if digit != initial:
        raise ValueError("Calibration is not periodic from its starting rail")
    repeats, tail = divmod(appointments, len(word))
    raw = [r for _, r in word]
    return {
        "period": raw,
        "repeats": repeats,
        "suffix": raw[:tail],
        "terminal": raw[(appointments - 1) % len(word)],
    }


def acquire_case(mount, fault, auth):
    members = tuple(i for i in range(4) if i != fault)
    labels = [(s, r) for s in members for r in members]
    table = {
        k: {f"{s}/{r}": mount.columns(members, s, r, park) for s, r in labels}
        for k, park in (("native", False), ("parked", True))
    }
    wire = ClockWire(mount, members, table)
    trial = MountedTrial(mount, auth, wire)
    stock = trial.stock
    idle = {}
    for i in members:
        trial.nodes[i].change(fault)
        owner = mount.homes[i]
        m = mount.local[owner]
        blank = stock[:]
        dimer = m.rows[owner][0]
        for r in dimer:
            blank[r] = 2
        idle[str(i)] = inspect_pair(m, blank, dimer, owner)
    view = fault + 1
    leader = view % 4
    wire.phase = "change"
    reports = trial.reports(view, members)
    offer = trial.propose(view, fault % 2, reports)
    cert = offer
    for stage in ("prepare", "accept", "commit"):
        wire.phase = stage
        cert = trial.collect(stage, cert, members, leader)
        if cert is None:
            raise ValueError("Quorum did not complete a phase")
    wire.phase = "decision"
    decisions = [trial.decision(leader, i, cert) for i in members]
    main_stop = wire.elapsed
    wire.phase = "parked control"
    distant = next(i for i in members if i != leader)
    parked_accept = trial.decision(leader, distant, cert, parked=True)
    by_phase = {}
    for event in wire.events:
        key = event["phase"]
        by_phase[key] = by_phase.get(key, 0) + event["stop"] - event["start"]
    widest = max(r["slots"] for r in table["native"].values())
    theta_wire = 6 * 8 * BYTE_CAP * widest
    clocks = mount.clock()
    for event in wire.events:
        n = event["stop"] - event["start"]
        event["clock"] = clock_trace(clocks["tick"], n, event["start"] % 2)
    # One failed view gets its full communication allowance. Processing is
    # supplied separately. The three local peeks above certify an empty mailbox.
    timeout = 9 * theta_wire
    clock_panel = {
        kind: clock_trace(clocks[kind], timeout) for kind in ("tick", "parked")
    }
    return {
        "regions": mount.regions,
        "silent": fault,
        "view": view,
        "members": members,
        "leader": leader,
        "silent_view_mailboxes": idle,
        "timeout_clock": clock_panel,
        "columns": table,
        "cut": mount.cut_control(members),
        "decisions": [trial.nodes[i].decision for i in members],
        "authenticated": decisions,
        "parked_accepted": parked_accept,
        "signature_bank": auth.archive(),
        "packets": wire.packets,
        "traffic": wire.events,
        "local_records": trial.log,
        "lock_reads": {str(i): trial.nodes[i].memory.events for i in members},
        "phase_slots": by_phase,
        "communication_slots_to_finalization": main_stop,
        "worst_bit_slots": widest,
        "finite_phase_wire_bound": theta_wire,
    }


def decimal_tags():
    # A signature has exactly 64 bytes and fixed-length hex representation.
    # Replacing the signature here affects neither field width nor byte count.
    def byte_length(view):
        return len(encode(["change", 0, view, -1, "", "0" * 128]))

    base = byte_length(0)
    widths = {str(v): byte_length(v) for v in (0, 9, 10, 99, 100, 1000, 10**9, 10**32)}
    if any(byte_length(int(v)) - base != len(v) - 1 for v in widths):
        raise ValueError("Serial view representation differs from decimal grammar")
    return {"bytes_by_view": widths, "constant_bytes_excluding_decimal_view": base - 1}


def acquire(enrollment=None):
    board_records = {}
    for level in (0, 1, 2):
        top = topology(level)
        ring = top["pentagons"][0]
        mount = Mount(level, ring)
        records = {}
        for fault in range(4):
            key = str(fault)
            bank = (
                None
                if enrollment is None
                else enrollment["levels"][str(level)]["cases"][key]["signature_bank"]
            )
            records[key] = acquire_case(mount, fault, Signers(bank))
        board_records[str(level)] = {
            "topology": top,
            "clock": mount.clock(),
            "cases": records,
        }
    return {
        "schema": 1,
        "byte_cap": BYTE_CAP,
        "decimal_tags": decimal_tags(),
        "levels": board_records,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--enroll", type=Path, help="create a signed deposit at a fresh path"
    )
    args = parser.parse_args()
    archive = HERE / "records.json.gz"
    if args.enroll is not None:
        if args.enroll.exists():
            raise SystemExit("Use a fresh enrollment path")
        data = acquire()
        raw = encode(data)
        with args.enroll.open("xb") as stream:
            stream.write(gzip.compress(raw, compresslevel=7, mtime=0))
    else:
        raw = gzip.decompress(archive.read_bytes())
        data = acquire(json.loads(raw))
        if encode(data) != raw:
            raise ValueError("Deposited records differ from replay")
    for level, b in data["levels"].items():
        print(
            level,
            b["topology"]["connected_quadruples"],
            [
                (c["decisions"], c["communication_slots_to_finalization"])
                for c in b["cases"].values()
            ],
        )
