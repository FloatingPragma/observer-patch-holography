"""Port-confined native packet transfer and an independent material clock."""

from copy import copy
from functools import cache
from itertools import permutations

import networkx as nx

from model.executable_model import federation
from ...boolean_nodes.material import Machine
from ..transport import read_rail


def owner_graph(level):
    return nx.Graph((a // 12, b // 12) for a, b in federation(level).external)


def topology(level):
    graph = owner_graph(level)
    # Grow connected sets by one vertex. Each connected quadruple has a
    # connected predecessor (remove a spanning-tree leaf).
    sets = {frozenset([v]) for v in graph}
    for _ in range(3):
        sets = {
            group | {w}
            for group in sets
            for v in group
            for w in graph[v]
            if w not in group
        }
    survivors = [
        sorted(g)
        for g in sets
        if all(nx.is_connected(graph.subgraph(g - {v})) for v in g)
    ]
    # Canonicalize all simple 5-cycles by their least vertex and orientation.
    paths = [(i,) for i in graph]
    for _ in range(4):
        paths = [
            p + (v,) for p in paths for v in graph[p[-1]] if v > p[0] and v not in p
        ]
    cycles = sorted(p for p in paths if graph.has_edge(p[-1], p[0]) and p[1] < p[-1])
    return {
        "cells": len(graph),
        "connected_quadruples": len(sets),
        "all_triples_connected_singletons": survivors,
        "pentagons": cycles,
    }


def restrict(machine, allowed):
    result = copy(machine)
    result.lanes = machine.lanes.subgraph(
        [p for p in machine.lanes if p // 12 in allowed]
    ).copy()
    return result


def apply_swaps(loads, route):
    for a, b in route:
        loads[a], loads[b] = loads[b], loads[a]


def inspect_pair(machine, stock, dimer, receiver):
    _, _, standard, head = machine.rows[receiver]
    return [read_rail(machine, stock, rail, standard[0], head) for rail in dimer]


def read_cost(machine, dimer, receiver):
    _, _, standard, (a, b) = machine.rows[receiver]
    # Stage each rail and the standard, write, read, invert, and unstage.
    return sum(
        2
        * (
            len(machine.move(r, a))
            + len(machine.move(r if standard[0] == a else standard[0], b))
        )
        + 3
        for r in dimer
    )


class Mount:
    def __init__(self, level, ring):
        self.level = level
        self.regions = ((ring[0],), (ring[1],), (ring[2],), tuple(ring[3:]))
        self.homes = [r[0] for r in self.regions]
        self.clock_owner = min(set(range(federation(level).cells)) - set(ring))
        self.base = Machine(level, max(*ring, self.clock_owner) + 1)
        self.local = {
            i: restrict(self.base, {i}) for i in [*self.homes, self.clock_owner]
        }

    @cache
    def clock(self):
        owner = self.clock_owner
        m = self.local[owner]
        clock = m.rows[owner][0]
        calibration = {}
        for parked in [False, True]:
            bank = []
            for initial in [0, 1]:
                bits = [0] * len(m.rows)
                bits[owner] = initial
                stock = m.fresh(bits)
                before = stock[:]
                if not parked:
                    apply_swaps(stock, [clock])
                response = inspect_pair(m, stock, clock, owner)
                if not parked:
                    apply_swaps(stock, [clock])
                if stock != before:
                    raise ValueError("Clock read failed restoration")
                bank.append(response)
            calibration["parked" if parked else "tick"] = bank
        return {
            "owner": owner,
            "rails": clock,
            "pointer": m.rows[owner][3],
            "read_appointments": read_cost(m, clock, owner),
            **calibration,
        }

    @cache
    def allowed_machine(self, members):
        cells = {cell for i in members for cell in self.regions[i]}
        return restrict(self.base, cells)

    @cache
    def columns(self, members, sender, receiver, parked=False, disabled=()):
        m = self.allowed_machine(members)
        u, v = self.homes[sender], self.homes[receiver]
        src, dst = m.rows[u][0], m.rows[v][0]
        route = tuple(e for a, b in zip(src, dst, strict=True) for e in m.move(a, b))
        blocked = {cell for i in disabled for cell in self.regions[i]}
        applied = tuple(
            (a, b)
            for a, b in route
            if (not parked or a // 12 == b // 12)
            and a // 12 not in blocked
            and b // 12 not in blocked
        )
        raw = []
        for bit in [0, 1]:
            bits = [0] * len(m.rows)
            bits[u] = bit
            stock = m.fresh(bits)
            initial = stock[:]
            apply_swaps(stock, applied)
            raw.append(inspect_pair(self.local[v], stock, dst, v))
            apply_swaps(stock, reversed(applied))
            if stock != initial:
                raise ValueError("Route disturbed its spectators")
        return {
            "source": src,
            "destination": dst,
            "path": route,
            "raw": raw,
            "slots": 2 * len(route) + read_cost(self.local[v], dst, v),
            "receiver_owner": v,
        }

    def cut_control(self, members):
        cells = {cell for i in members for cell in self.regions[i]}
        graph = nx.Graph()
        graph.add_nodes_from(members)
        board = self.base.car
        by_cell = {v: i for i in members for v in self.regions[i]}
        for a, b in board.external:
            x, y = a // 12, b // 12
            if x in cells and y in cells and by_cell[x] != by_cell[y]:
                graph.add_edge(by_cell[x], by_cell[y])
        (middle,) = (i for i in members if graph.degree[i] == 2)
        left = tuple(i for i in members if i != middle)
        isolated = self.allowed_machine(left)
        disconnected = []
        for a, b in permutations(left, 2):
            src = isolated.rows[self.homes[a]][0][0]
            dst = isolated.rows[self.homes[b]][0][0]
            try:
                isolated.move(src, dst)
            except nx.NetworkXNoPath:
                disconnected.append([a, b])
        return {"parked_validator": middle, "unreachable_directions": disconnected}
