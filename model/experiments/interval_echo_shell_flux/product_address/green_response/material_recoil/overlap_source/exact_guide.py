"""Exact rational witness on all 1,000 addresses of the twelve-port material guides."""

from collections import Counter
import argparse
from fractions import Fraction as F
from itertools import combinations_with_replacement, product
from pathlib import Path
import json

import sympy as sp
from sympy.polys.matrices import DomainMatrix

N = 10
DIRECTIONS = tuple((a, s) for a in range(3) for s in (-1, 1))


def neighbours(x):
    for axis, sign in DIRECTIONS:
        y = list(x)
        y[axis] += sign
        yield tuple(y)


def orbit_solve(nodes, fold, load):
    """Each orbit equation retains neighbour multiplicities. Absent nodes mean zero."""
    index = {x: i for i, x in enumerate(nodes)}
    rows = {}
    for i, x in enumerate(nodes):
        row = Counter({i: 6})
        for y in neighbours(x):
            k = fold(y)
            if k in index:
                row[index[k]] -= 1
        rows[i] = {j: sp.ZZ(v) for j, v in row.items() if v}
    size = len(nodes)
    a = DomainMatrix(rows, (size, size), sp.ZZ)
    rhs = [load(x) * 8 for x in nodes]
    den = int(sp.ilcm(*(v.denominator for v in rhs)))
    b = DomainMatrix(
        {i: {0: sp.ZZ(v * den)} for i, v in enumerate(rhs) if v}, (size, 1), sp.ZZ
    )
    solution, scale = a.solve_den(b, method="rref")
    vector = solution.to_list()
    values = {x: F(int(vector[i][0]), den * int(scale)) for i, x in enumerate(nodes)}
    for x in nodes:
        assert 6 * values[x] - sum(
            values.get(fold(y), 0) for y in neighbours(x)
        ) == 8 * load(x)
    return values


def run():
    n = N
    period = 2 * (n + 1)

    # Pin the arbitrary additive constant of the neutral periodic fundamental
    # solution at the origin. Odd reflection eliminates this constant.
    def periodic(x):
        return tuple(sorted(min(i % period, (-i) % period) for i in x))

    pnodes = list(combinations_with_replacement(range(n + 2), 3))
    kernel = orbit_solve(pnodes[1:], periodic, lambda _: -F(1, period**3))
    origin = (0, 0, 0)
    kernel[origin] = F(0)
    assert -sum(kernel[periodic(y)] for y in neighbours(origin)) == 8 * (
        1 - F(1, period**3)
    )

    def image(x, y):
        return sum(
            s
            * t
            * u
            * kernel[periodic((x[0] - s * y[0], x[1] - t * y[1], x[2] - u * y[2]))]
            for s, t, u in product((-1, 1), repeat=3)
        )

    def folded(x):
        return tuple(sorted(min(a, n + 1 - a) for a in x))

    nodes = list(combinations_with_replacement(range(1, n // 2 + 1), 3))
    mid = (n // 2,) * 3
    v = orbit_solve(nodes, folded, lambda x: F(int(x == mid), 8))
    hit = {}
    # First-hit boundary-value verification for every inequivalent target.
    # This bypasses any finite-time truncation or stochastic sampling.
    vertices = tuple(product(range(1, n + 1), repeat=3))
    for target in nodes:
        diagonal = image(target, target)
        g = {x: image(x, target) / diagonal for x in vertices}
        assert g[target] == 1 and all(0 < a <= 1 for a in g.values())
        assert all(
            6 * g[x] == sum(g.get(y, 0) for y in neighbours(x))
            for x in vertices
            if x != target
        )
        hit[target] = sum(g[x] for x in vertices if folded(x) == mid) / 8
        assert hit[target] == v[target] / diagonal
    stock = {x: p**2 for x, p in hit.items()}
    field = orbit_solve(nodes, folded, stock.__getitem__)
    shells = []
    for k in range(1, n // 2 + 1):
        region = {x for x in vertices if min(folded(x)) > n // 2 - k}
        flux = sum(
            (field[folded(x)] - field.get(folded(y), 0)) / 8
            for x in region
            for y in neighbours(x)
            if y not in region
        )
        faces = sum(y not in region for x in region for y in neighbours(x))
        mass = sum(stock[folded(x)] for x in region)
        assert faces == 24 * k * k and flux == mass and mass > 0
        shells.append(dict(k=k, edges=faces, enclosed_squared_admission=float(mass)))
    return dict(
        interior=n,
        vertices=len(vertices),
        inequivalent_targets=len(nodes),
        periodic_orbits=len(pnodes),
        highest_hit=float(hit[mid]),
        shells=shells,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    receipt = Path(__file__).with_name("exact_guide.json")
    result = run()
    if args.write:
        receipt.write_text(json.dumps(result, indent=2) + "\n")
    else:
        assert json.loads(receipt.read_text()) == result
