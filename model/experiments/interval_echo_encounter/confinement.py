"""Exact single-cell invariant spaces of the free Interval Echo load sectors."""
import json
import math
from collections import Counter, defaultdict, deque
from fractions import Fraction
from functools import reduce
from itertools import combinations
from pathlib import Path

import numpy as np

from model.executable_model.carrier import CELL_SEAMS, automorphisms, federation

HERE = Path(__file__).resolve().parent
RECEIPT = HERE / "confinement.json"


def blocks(gates, window, strength):
    """Integer blocks A=2PHP and B=2(1-P)HP, after the uniform sector shift."""
    states = {"single": [(i,) for i in range(12)], "home": list(combinations(range(12), 2)),
              "seams": list(CELL_SEAMS)}[window]
    index = {s: i for i, s in enumerate(states)}
    moves = (*CELL_SEAMS, *((p, 12 + j) for j, p in enumerate(gates)))
    inside = np.zeros((len(states), len(states)), dtype=int)
    outside = defaultdict(lambda: [0] * len(states))
    for col, state in enumerate(states):
        occupied = set(state)
        if state in CELL_SEAMS:
            inside[col, col] -= 4 * strength
        for a, b in moves:
            if (a in occupied) == (b in occupied):
                continue
            target = tuple(sorted(occupied ^ {a, b}))
            inside[col, col] += 1
            if target in index:
                inside[index[target], col] -= 1
            else:
                outside[target][col] -= 1
    return states, inside, dict(outside)


def primitive(row):
    divisor = reduce(math.gcd, row, 0)
    return tuple(x // divisor for x in row) if divisor else tuple(row)


def maximal_space(inside, exits):
    """Close the rational exit row space under multiplication by A."""
    a = inside.astype(object)
    queue = deque(map(tuple, exits))
    pivots = {}
    while queue:
        row = queue.popleft()
        for j in sorted(pivots):
            if row[j]:
                basis = pivots[j]
                row = primitive([basis[j] * x - row[j] * y for x, y in zip(row, basis)])
        pivot = next((j for j, x in enumerate(row) if x), None)
        if pivot is not None:
            pivots[pivot] = row
            queue.append(tuple(np.asarray(row, dtype=object) @ a))
    free = [j for j in range(len(a)) if j not in pivots]
    columns = []
    for j in free:
        v = [Fraction(i == j) for i in range(len(a))]
        for k in sorted(pivots, reverse=True):
            row = pivots[k]
            v[k] = -sum(x * y for x, y in zip(row, v)) / row[k]
        columns.append(v)
    return free, np.asarray(columns, dtype=object).reshape(len(free), len(a)).T


def restriction(inside, exits):
    free, vectors = maximal_space(inside, exits.values())
    action = inside.astype(object) @ vectors
    reduced = action[free, :]
    outward = np.array(list(exits.values()), dtype=object).reshape(-1, len(inside))
    if np.any(action != vectors @ reduced) or np.any(outward @ vectors):
        raise ArithmeticError("Invariant witness fails at a retained or exit port")
    # Newton identities for det(xI - 2H) on the entire compact space.
    power = np.eye(len(free), dtype=object)
    polynomial, traces = [Fraction(1)], []
    for order in range(1, len(free) + 1):
        power = power @ reduced
        traces.append(np.trace(power))
        polynomial.append(-sum(x * y for x, y in zip(reversed(polynomial), traces)) / order)
    if any(x.denominator != 1 for x in polynomial):
        raise ArithmeticError("Nonintegral factor of a monic integer characteristic polynomial")
    return vectors, [int(c) for c in polynomial]


def gate_profiles(car):
    internal = [set() for _ in range(car.cells)]
    for a, b in car.seams:
        if a // 12 == b // 12:
            internal[a // 12].add((a % 12, b % 12))
    if any(local != set(CELL_SEAMS) for local in internal):
        raise ValueError("A cell differs from the local exchange chart")
    profiles = [[] for _ in range(car.cells)]
    ends = set()
    for a, b in car.external:
        for p in (a, b):
            if p in ends:
                raise ValueError("External seams must form a matching")
            ends.add(p)
            profiles[p // 12].append(p % 12)
    if any(len(g) != 3 for g in profiles):
        raise ValueError("Each cell needs three distinct gateways")
    return [tuple(sorted(g)) for g in profiles]


def certificate():
    protocol = json.loads((HERE / "confinement_protocol.json").read_text())
    symmetries = automorphisms(CELL_SEAMS)
    if len(symmetries) != 120:
        raise ArithmeticError("Incomplete cell symmetry group")
    populations = defaultdict(dict)
    for level in protocol["carrier_levels"]:
        canonical = Counter(min(tuple(sorted(p[i] for i in g)) for p in symmetries)
                            for g in gate_profiles(federation(level)))
        for gates, count in canonical.items():
            populations[gates][str(level)] = count
    profiles = []
    for gates, counts in sorted(populations.items()):
        cases = {}
        for window in ("single", "home", "seams"):
            cases[window] = {}
            strengths = [0] if window == "single" else protocol["strengths"]
            for strength in strengths:
                states, inside, exits = blocks(gates, window, strength)
                vectors, polynomial = restriction(inside, exits)
                cases[window][str(strength)] = {
                    "dimension": vectors.shape[1], "twice_shifted_energy_polynomial": polynomial}
        profiles.append({"gates": list(gates), "cells_by_level": counts, "windows": cases})
    release = protocol["release"]
    car = federation(release["level"])
    gates = gate_profiles(car)[release["cell"]]
    start = next(e for e in CELL_SEAMS if set(e).intersection(gates))
    escape = {}
    for window in ("home", "seams"):
        states, inside, exits = blocks(gates, window, release["strengths"][0])
        vectors, _ = restriction(inside, exits)
        source = states.index(start)
        if any(vectors[source, :]):
            raise ArithmeticError("Gateway release overlaps an exactly confined preparation")
        escape[window] = str(Fraction(sum(row[source]**2 for row in exits.values()), 4))
    products = {}
    for level in protocol["carrier_levels"]:
        dimensions = [profile["windows"]["single"]["0"]["dimension"]
                      for profile in profiles
                      for _ in range(profile["cells_by_level"].get(str(level), 0))]
        products[str(level)] = sum(a*b for a,b in combinations(dimensions, 2))
    return {"profiles": profiles, "distinct_cell_product_dimensions": products,
            "release_local_ports": list(start),
            "release_compact_weight": 0, "quiet_quadratic_coefficients": escape}


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = certificate()
    if args.write:
        from .transport import write_receipt
        write_receipt(RECEIPT, result)
    elif result != json.loads(RECEIPT.read_text()):
        raise ValueError("Exact confinement receipt differs")
    print(json.dumps(result))
