"""Read a receiver-local view from ten physical dimers and verify compact frames."""

import argparse
from copy import copy
from dataclasses import dataclass
from pathlib import Path
import gzip
import json

import networkx as nx

from ..certificates import Signers
from ..transport import encode, identifier
from .duplicate import receive
from .measure import DECODE, HERE, MountedTrial
from .routes import Mount, apply_swaps, inspect_pair, owner_graph, read_cost, topology

SIGNED = {"change", "propose", "prepare", "accept", "commit"}
CERTS = {"PC", "LC", "DC"}
WIDTH = 10
MARK = None


def views(packet, local, *, restore):
    """The signed grammar puts the view at index 2 of votes or 1 of certificates."""
    if not isinstance(packet, list):
        return packet
    tag = packet[0] if packet and isinstance(packet[0], str) else ""
    ix = 2 if tag in SIGNED else 1 if tag in CERTS else None
    if tag and (ix is None or len(packet) != (6 if tag in SIGNED else 5)):
        raise ValueError("Unsupported signed frame")
    result = [views(x, local, restore=restore) for x in packet]
    if ix is not None:
        expected = MARK if restore else local
        if result[ix] != expected or type(result[ix]) is not type(expected):
            raise ValueError("Frame needs a different epoch")
        result[ix] = local if restore else MARK
    return result


@dataclass
class EpochRegister:
    cells: list
    digits: list
    machine: object
    home: int
    stock: list

    def read(self):
        if len(self.digits) < WIDTH:
            return {
                "cells": self.cells,
                "digits": self.digits,
                "available": len(self.digits),
            }
        before = self.stock[:]
        raw = [
            inspect_pair(self.machine, self.stock, d, self.home) for d in self.digits
        ]
        if self.stock != before:
            raise ValueError("Epoch read failed restoration")
        return {
            "cells": self.cells,
            "digits": self.digits,
            "raw": raw,
            "local_slots": sum(
                read_cost(self.machine, d, self.home) for d in self.digits
            ),
            "read": sum(DECODE[tuple(r)] << j for j, r in enumerate(raw)),
        }


def register(mount, members, owner, view, stock):
    if not 0 <= view < 2**WIDTH:
        raise ValueError("Outside the mounted view register")
    homes = [mount.homes[i] for i in members]
    blocked = {mount.clock_owner} | set().union(
        *(set(mount.regions[i]) for i in range(4) if i not in members)
    )
    graph = owner_graph(mount.level)
    graph.remove_nodes_from(blocked)
    distance = {
        i: nx.single_source_shortest_path_length(
            graph.subgraph(set(graph) - (set(homes) - {mount.homes[i]})), mount.homes[i]
        )
        for i in members
    }
    patch = sorted(
        v
        for v in graph
        if min((distance[i].get(v, float("inf")), i) for i in members)[1] == owner
    )
    home = mount.homes[owner]
    machine = copy(mount.base)
    machine.lanes = nx.Graph(
        (a, b) for a, b in machine.car.seams if a // 12 in patch and b // 12 in patch
    )
    occupied = {p for row in [machine.rows[home]] for part in row for p in part}
    # Exclude all other enrolled packet stock, including the second cell of region 3.
    occupied.update(
        p
        for p in machine.lanes
        if p // 12 != home and any(p // 12 in z for z in mount.regions)
    )
    digits = []
    for a, b in sorted(machine.lanes.edges):
        if a // 12 == b // 12 and a not in occupied and b not in occupied:
            digits.append((a, b))
            occupied.update((a, b))
            if len(digits) == WIDTH:
                break
    for j, (a, b) in enumerate(digits):
        val = (view >> j) & 1
        stock[a], stock[b] = 2 * val, 2 * (1 - val)
    return EpochRegister(patch, digits, machine, home, stock)


class EpochWire:
    def __init__(self, mount, members, regs):
        self.mount, self.members, self.regs = mount, members, regs
        self.events = []
        self.packets = {}

    def send(self, s, r, packet, parked=False):
        source_read = self.regs[s].read()
        v = source_read["read"]
        tagged = views(packet, v, restore=False)
        word = encode(tagged)
        col = self.mount.columns(self.members, s, r, parked)
        # Every bit appointment stages the data on the same stock as the epoch
        # dimers. Check both physical transfer columns at its current contents.
        stock = self.regs[s].stock
        initial = stock[:]
        a, b = col["source"]
        raw = []
        applied = [
            edge for edge in col["path"] if not parked or edge[0] // 12 == edge[1] // 12
        ]
        for bit in (0, 1):
            stock[a], stock[b] = 2 * bit, 2 * (1 - bit)
            apply_swaps(stock, applied)
            raw.append(
                inspect_pair(
                    self.mount.local[self.mount.homes[r]],
                    stock,
                    col["destination"],
                    self.mount.homes[r],
                )
            )
            apply_swaps(stock, reversed(applied))
            stock[a], stock[b] = initial[a], initial[b]
            if stock != initial:
                raise ValueError("Transit changed register stock")
        if raw != col["raw"]:
            raise ValueError("Data calibration on stocked epoch differs")
        got = receive(raw, word)
        receiver_read = self.regs[r].read()
        incoming = views(got, receiver_read["read"], restore=True)
        key = identifier(packet)
        self.packets[key] = packet
        self.events.append(
            {
                "source": s,
                "destination": r,
                "packet": key,
                "compact_bytes": len(word),
                "absolute_bytes": len(encode(packet)),
                "slots_per_bit": col["slots"],
                "source_raw": source_read["raw"],
                "receiver_raw": receiver_read["raw"],
                "live_data_raw": raw,
                "source_read_slots": source_read["local_slots"],
                "local_read_slots": receiver_read["local_slots"],
            }
        )
        return incoming


def run(mount, v, bank=None):
    leader = v % 4
    silent = (leader - 1) % 4
    members = tuple(i for i in range(4) if i != silent)
    stock = mount.base.fresh([0] * len(mount.base.rows))
    regs = {i: register(mount, members, i, v, stock) for i in members}
    panels = {i: r.read() for i, r in regs.items()}
    if any("read" not in r for r in panels.values()):
        return {"view": v, "insufficient_stock": panels}
    other_stock = mount.base.fresh([0] * len(mount.base.rows))
    wrong = {i: register(mount, members, i, v + 1, other_stock) for i in members}
    wrong_panels = {i: r.read() for i, r in wrong.items()}
    auth = Signers(bank)
    wire = EpochWire(mount, members, regs)
    trial = MountedTrial(mount, auth, wire, stock=stock)
    cert = trial.issue(v, 1, members)
    if not all(trial.decision(leader, i, cert) for i in members):
        raise ValueError("Compact decision refused")
    frame = views(cert, v, restore=False)
    wrong_results = {}
    for i in members:
        delivered = receive(mount.columns(members, leader, i)["raw"], encode(frame))
        wrong_results[str(i)] = trial.rules.cert(
            views(delivered, wrong[i].read()["read"], restore=True), "DC"
        )
    # The retained prior decision supplies an LC in each next-view report.
    report = trial.nodes[leader].change(v + 1)
    if trial.rules.report_signer(v + 1, report) != leader:
        raise ValueError("Invalid prior-view lock control")
    try:
        views(report, v + 1, restore=False)
    except ValueError:
        prior_refused = True
    else:
        prior_refused = False
    return {
        "view": v,
        "registers": panels,
        "changed_registers": wrong_panels,
        "events": wire.events,
        "packets": wire.packets,
        "decision": identifier(cert),
        "prior_lock_report": report,
        "decisions": [trial.nodes[i].decision for i in members],
        "wrong_view_decisions": wrong_results,
        "prior_lock_refused_by_single_view_codec": prior_refused,
        "signature_bank": auth.archive(),
    }


def acquire(deposit=None):
    by_level = {}
    for level in [0, 1, 2]:
        mount = Mount(level, topology(level)["pentagons"][0])
        cases = {}
        for v in (1, 10, 100, 1000):
            bank = (
                None
                if deposit is None
                else deposit[str(level)][str(v)].get("signature_bank")
            )
            cases[str(v)] = run(mount, v, bank)
        by_level[str(level)] = cases
    return by_level


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)

    cli.add_argument(
        "--enroll", type=Path, help="write an independent enrollment to a fresh file"
    )
    args = cli.parse_args()
    f = HERE / "local_view.json.gz"
    if args.enroll is not None:
        if args.enroll.exists():
            raise FileExistsError("Use a fresh enrollment path")
        raw = encode(acquire())
        with args.enroll.open("xb") as stream:
            stream.write(gzip.compress(raw, mtime=0))
    else:
        raw = gzip.decompress(f.read_bytes())
        if encode(acquire(json.loads(raw))) != raw:
            raise ValueError("Local epoch receiver bank differs")
