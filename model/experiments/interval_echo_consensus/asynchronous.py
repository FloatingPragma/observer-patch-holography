"""Test a sharp local arbitration service against all orders on a bounded patch."""

import argparse
from functools import cache
from itertools import combinations, product
import json
from pathlib import Path
import random

import numpy as np

from model.candidate_laws.interval_echo import observations
from model.candidate_laws.interval_echo.candidate import native_law
from model.executable_model import CELL_SEAMS, federation
from .__main__ import local_ops


HERE = Path(__file__).resolve().parent


def verify_native_service():
    _, _, write, probes = local_ops()
    identity = np.eye(99)
    for index, (left, right) in enumerate(native_law.SRC):
        source = identity[:, index]
        ready = abs(left - right) == 2
        # Occupied, blank/quiet and blank/ready are separate physical branches.
        for (_, op), selected in zip(probes, (False, not ready, ready), strict=True):
            expected = source if selected else 0
            np.testing.assert_allclose(op @ source, expected, atol=1e-12)
        expected = np.zeros(99)
        if ready:
            for coin in (0, 1):
                collar = 1 + native_law.LABELS.index((coin, left - right))
                expected[9 * collar + native_law.SRC.index((1, 1))] = 1 / np.sqrt(2)
        np.testing.assert_allclose(write @ source, expected, atol=1e-12)
    np.testing.assert_array_equal(probes[0][1][:, 9:], identity[:, 9:])
    for _, op in (*probes[1:], (True, write)):
        np.testing.assert_allclose(op[:, 9:], 0, atol=1e-12)


def priority(edge):
    a, b = edge
    if a // 12 != b // 12:
        return 0
    return 1 + CELL_SEAMS.index((a % 12, b % 12))


class SharpService:
    """A local guard and commit derived from the sharp native instrument."""

    def __init__(self, ports, seams, reverse=False):
        indices = {p: i for i, p in enumerate(ports)}
        self.edges = [(indices[a], indices[b]) for a, b in seams]
        self.rank = [priority(edge) * (-1 if reverse else 1) for edge in seams]
        incident = [[] for _ in ports]
        for i, (a, b) in enumerate(self.edges):
            incident[a].append(i)
            incident[b].append(i)
        if any(len({self.rank[i] for i in edges}) != len(edges) for edges in incident):
            raise ValueError("the local seam chart contains an incident priority tie")
        self.higher = [
            tuple(j for j in (*incident[a], *incident[b]) if self.rank[j] < self.rank[i])
            for i, (a, b) in enumerate(self.edges)
        ]

    def ready(self, i, loads, collars):
        a, b = self.edges[i]
        return collars[i] == 0 and abs(loads[a] - loads[b]) == 2

    def enabled(self, i, loads, collars):
        return self.ready(i, loads, collars) and all(
            not self.ready(j, loads, collars) for j in self.higher[i]
        )

    def commit(self, i, loads, collars):
        a, b = self.edges[i]
        if not self.ready(i, loads, collars):
            raise ValueError("attempted write to an unready seam")
        collars[i] = loads[a] - loads[b]
        loads[a] = loads[b] = 1


def induced_path(car):
    a, b = min(car.external)
    for first in car.neighbours[a]:
        for last in car.neighbours[b]:
            for end in car.neighbours[last]:
                ports = (first, a, b, last, end)
                if len(set(ports)) == 5 and sum(y in car.neighbours[x] for x, y in combinations(ports, 2)) == 4:
                    return ports
    raise ValueError("missing intercell path")


def exhaustive_patch():
    car = federation(0)
    ports = induced_path(car)
    seams = [e for e in car.seams if set(e) <= set(ports)]
    service = SharpService(ports, seams)

    def after(state, i):
        loads, collars = map(list, state)
        service.commit(i, loads, collars)
        return tuple(loads), tuple(collars)

    @cache
    def endpoint(state):
        moves = [i for i in range(4) if service.enabled(i, *state)]
        if not moves:
            if any(service.ready(i, *state) for i in range(4)):
                raise AssertionError("inconsistent terminal")
            return state
        terminals = {endpoint(after(state, i)) for i in moves}
        if len(terminals) != 1:
            raise AssertionError("different terminal records reached")
        return terminals.pop()

    forks = transitions = terminals = 0
    public_states = set()
    for loads in product(range(3), repeat=5):
        for collars in product((0, -2, 2), repeat=4):
            state = (loads, collars)
            normal = endpoint(state)
            occupied = frozenset(i for i, sign in enumerate(collars) if sign)
            ready = frozenset(i for i in range(4) if service.ready(i, *state))
            public_states.add((occupied, ready))
            moves = [i for i in range(4) if service.enabled(i, *state)]
            transitions += len(moves)
            terminals += not moves
            if sum(normal[0]) != sum(loads) or any(
                sign and sign != final for sign, final in zip(collars, normal[1], strict=True)
            ):
                raise AssertionError("protected boundary changed")
            for i in moves:
                next_state = after(state, i)
                updated, marked = next_state
                if sum(x*x for x in loads) - sum(x*x for x in updated) != 2:
                    raise AssertionError("accepted repair lacks strict descent")
                surviving = {j for j in ready if set(service.edges[j]).isdisjoint(service.edges[i])}
                if (surviving != {j for j in range(4) if service.ready(j, *next_state)}
                        or occupied | {i} != {j for j, sign in enumerate(marked) if sign}):
                    raise AssertionError("public successor depends on hidden data")
            for i, j in combinations(moves, 2):
                left, right = after(state, i), after(state, j)
                if not service.enabled(j, *left) or not service.enabled(i, *right) or after(left, j) != after(right, i):
                    raise AssertionError("one-step peak does not join")
                forks += 1

    initial = ((0, 2, 0, 1, 1), (0,) * 4)
    unguarded = [after(initial, i) for i in range(4) if service.ready(i, *initial)]
    if len(set(unguarded)) != 2 or any(service.ready(i, *s) for s in unguarded for i in range(4)):
        raise AssertionError("unguarded fork does not have two distinct terminals")
    return dict(ports=ports, seams=seams, states=endpoint.cache_info().currsize,
                public_quotient_states=len(public_states), normal_states=terminals,
                transitions=transitions, joined_one_step_peaks=forks, control_start=initial,
                unguarded_endpoints=unguarded, guarded_endpoint=endpoint(initial))


def settle(service, initial, order, seed):
    loads, collars = map(list, initial)
    addresses = list(range(len(collars)))
    rng = random.Random(seed)
    if order == "reference":
        # Independent serial greedy matching of the initially discrepant edges.
        for i in sorted(addresses, key=service.rank.__getitem__):
            if service.ready(i, loads, collars):
                service.commit(i, loads, collars)
    else:
        if order == "reverse":
            addresses.reverse()
        while True:
            if order == "shuffled":
                rng.shuffle(addresses)
            writes = 0
            for i in addresses:
                if service.enabled(i, loads, collars):
                    service.commit(i, loads, collars)
                    writes += 1
            if not writes:
                break
    if any(service.ready(i, loads, collars) for i in addresses):
        raise AssertionError("a terminal blank seam has a two-unit discrepancy")
    return loads, collars


def whole_board(level):
    car = federation(level)
    services = [SharpService(range(car.n), car.seams, reverse) for reverse in (False, True)]
    normal, reversed_rank = services
    edge_index = {e: i for i, e in enumerate(car.seams)}
    _, _, _, sign_read = observations.instruments()
    neighbourhoods = [
        {p // 12 for j in (i, *normal.higher[i]) for p in car.seams[j]}
        for i in range(len(car.seams))
    ]
    rows = []
    for seed in (31, 47, 131):
        rng = random.Random(seed)
        initial = ([rng.randrange(3) for _ in range(car.n)],
                   [rng.choice((-2, 2, 0, 0, 0, 0, 0, 0, 0, 0)) for _ in car.seams])
        final = settle(normal, initial, "reference", seed)
        for order in ("forward", "reverse", "shuffled"):
            if settle(normal, initial, order, seed) != final:
                raise AssertionError("schedule changed physical records")
        before, after = initial[1], final[1]
        writes = sum(x != y for x, y in zip(before, after, strict=True))
        if (any(x and x != y for x, y in zip(before, after, strict=True))
                or sum(initial[0]) != sum(final[0])
                or sum(x*x for x in initial[0]) - sum(x*x for x in final[0]) != 2*writes):
            raise AssertionError("protected data or descent failed")
        opposite = settle(reversed_rank, initial, "reference", seed)
        # The public tag reader needs intact (1,1) endpoints and a vacant neighbour
        # on a blank probe seam. Read each such tag on its own checkpoint copy.
        readable = 0
        for i, (a, b) in enumerate(car.seams):
            if before[i] or not after[i]:
                continue
            if (final[0][a], final[0][b]) != (1, 1):
                raise AssertionError("newly written pair lost its public endpoints")
            for port in (a, b):
                if any(final[0][v] == 0 and not after[edge_index[tuple(sorted((port, v)))]]
                       for v in car.neighbours[port]):
                    local_sign = after[i] if port == a else -after[i]
                    tag_index = 2 if local_sign > 0 else 1  # unwrite yields (2,0) or (0,2)
                    if 2 * int(sign_read[tag_index]) != local_sign:
                        raise AssertionError("native vacancy read did not recover the oriented sign")
                    readable += 1
                    break
        rows.append(dict(seed=seed, initial_discrepancies=sum(normal.ready(i, *initial) for i in range(len(after))),
                         writes=writes, external_writes=sum(before[i] != after[i] for i in range(len(after)) if normal.rank[i] == 0),
                         new_signs_with_vacancy_reader=readable,
                         different_collars_after_priority_reversal=sum(x != y for x, y in zip(opposite[1], after, strict=True))))
    return dict(cells=car.cells, seams=len(car.seams), maximum_guard_cells=max(map(len, neighbourhoods)), runs=rows)


def run():
    verify_native_service()
    return dict(native_dimension=99, sharp_inputs=9, occupied_collar_basis_states=90,
                patch=exhaustive_patch(), boards=[whole_board(level) for level in range(3)])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="write a new receipt")
    path = parser.parse_args().output
    result = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if path is not None:
        with path.open("x") as handle:
            handle.write(result)
    elif result != (HERE / "asynchronous.json").read_text():
        raise SystemExit("sharp consensus differs from its deposited receipt")
    print(result, end="")
