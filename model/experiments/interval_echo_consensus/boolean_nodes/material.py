"""Local sharp-register service executed with routed native seam pulses."""

from functools import lru_cache
from itertools import combinations

import networkx as nx
import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo.candidate import KAPPA, native_law
from model.executable_model import CELL_SEAMS, federation


@lru_cache
def balanced_pointer(mode):
    law = native_law
    sx = np.eye(9)[:, [law.SRC.index((b, a)) for a, b in law.SRC]]
    native_swap = expm(-1j * np.pi * law.seam_term(KAPPA))
    np.testing.assert_allclose(native_swap, sx, atol=5e-14)
    echo = {
        "native": law.echo(),
        "reversed_output": law.echo(),
        "parked": np.eye(law.DIM),
        "swap_query": np.kron(np.eye(law.DIM // 9), sx),
    }[mode]
    table = {}
    for i, ab in enumerate(law.SRC):
        if 1 in ab:
            continue  # Boolean queries address only extremal data rails.
        state = echo[:, i]
        bins = np.tile([int((x, y) == (1, 1)) for x, y in law.SRC], law.DIM // 9)
        probs = np.bincount(bins, weights=abs(state) ** 2, minlength=2)
        q = int(np.argmax(probs))
        np.testing.assert_allclose(probs, np.eye(2)[q], atol=1e-13)
        restored = echo.conj().T @ (state * (bins == q))
        np.testing.assert_allclose(restored, np.eye(law.DIM)[:, i], atol=1e-13)
        table[ab] = q
    return table


@lru_cache
def slots(gateways):
    cell = nx.Graph(CELL_SEAMS)
    private = set(cell) - set(gateways)
    local = [e for e in CELL_SEAMS if set(e) <= private]
    for out, scratch in combinations(local, 2):
        if set(out) & set(scratch):
            continue
        for unit in sorted(private - set(out) - set(scratch)):
            reserved = set(out + scratch + (unit,))
            lanes = cell.subgraph(set(cell) - reserved)
            if nx.is_connected(lanes) and all(
                set(cell[p]) & set(lanes) for p in reserved
            ):
                return out, scratch, unit, tuple(sorted(lanes.edges())[0])
    raise ValueError("No connected node layout")


class Machine:
    def __init__(self, level, size):
        board = federation(level)
        if size > board.cells:
            raise ValueError("Choose a larger federation")
        self.rows = []
        self.lanes = nx.Graph()
        self.lanes.add_edges_from(board.seams)
        for owner in range(size):
            g = tuple(i for i in range(12) if board.degree[12 * owner + i] == 6)
            out, sc, c, e = slots(g)
            row = tuple(
                tuple(12 * owner + i for i in part) for part in (out, sc, (c,), e)
            )
            self.rows.append(row)
            self.lanes.remove_nodes_from(row[0] + row[1] + row[2])
        if not nx.is_connected(self.lanes):
            raise ArithmeticError("Transit network cut")
        self.car = board
        self.rail_swaps = 0

    def fresh(self, bits):
        # Unused transit ports deliberately carry nonuniform occupancies.
        loads = [(7 * i + 2) % 3 for i in range(self.car.n)]
        for value, (out, sc, c, _) in zip(bits, self.rows, strict=True):
            loads[out[0]], loads[out[1]] = 2 * value, 2 * (1 - value)
            loads[sc[0]], loads[sc[1]], loads[c[0]] = 0, 2, 2
        return loads

    @lru_cache
    def move(self, u, v):
        if u == v:
            return ()
        tmp = self.lanes.copy()
        for x in (u, v):
            if x not in self.lanes:
                tmp.add_edges_from(
                    (x, y) for y in self.car.neighbours[x] if y in self.lanes
                )
        path = nx.shortest_path(tmp, u, v)
        edges = tuple(zip(path, path[1:]))
        return edges + edges[:-1][::-1]

    def probe(self, loads, ports, edge, mode):
        u, v = ports
        a, b = edge
        if u == v:
            raise ValueError("Query requires different addresses")
        route = self.move(u, a) + self.move(u if v == a else v, b)
        for i, j in route:
            loads[i], loads[j] = loads[j], loads[i]
        q = balanced_pointer(mode)[loads[a], loads[b]]
        for i, j in reversed(route):
            loads[i], loads[j] = loads[j], loads[i]
        self.rail_swaps += 2 * len(route)
        return q

    def query(self, loads, operands, edge, reference, mode):
        # The registered detector reads only the balanced (1,1) output.
        # A mixed pair is already resolved. After a dark read both sharp
        # inputs agree, so comparing the first to a known double resolves it.
        if self.probe(loads, operands, edge, mode):
            return 1
        first_is_zero = self.probe(loads, (operands[0], reference), edge, mode)
        return 0 if first_is_zero else 2

    def attempt(self, loads, owner, node, mode):
        op, args = node
        out, sc, c, edge = self.rows[owner]
        if op == "nand":
            operands = tuple(self.rows[i][0][0] for i in args)
        elif op == "copy":
            operands = (self.rows[args[0]][0][0], c[0])
        else:
            operands = (sc[1 if op == "true" else 0], c[0])
        before = loads.copy()
        value_read = self.query(loads, operands, edge, c[0], mode)
        opens = (value_read != 2) if op == "nand" else (value_read == 2)
        if opens:
            loads[sc[0]], loads[sc[1]] = loads[sc[1]], loads[sc[0]]
        old_rail = out[1] if mode == "reversed_output" else out[0]
        mismatch_read = self.query(loads, (sc[0], old_rail), edge, c[0], mode)
        if mismatch_read == 1:
            loads[out[0]], loads[out[1]] = loads[out[1]], loads[out[0]]
        if opens:
            loads[sc[0]], loads[sc[1]] = loads[sc[1]], loads[sc[0]]
        if any(x != y for i, (x, y) in enumerate(zip(before, loads)) if i not in out):
            raise ArithmeticError("Write outside output")
        return value_read, mismatch_read
