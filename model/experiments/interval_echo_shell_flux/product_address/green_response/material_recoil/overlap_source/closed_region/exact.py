"""Rational Cesaro projection of the complete one-cell native dynamics."""

import argparse
from itertools import combinations
import json

import sympy as sp

from model.executable_model import CELL_SEAMS
from model.experiments.interval_echo_einstein_station.mobile_collar_source import (
    measure as formation,
)
from .measure import HERE


def exact_sector(count, k):
    words = list(combinations(range(12), count))
    idx = {frozenset(word): j for j, word in enumerate(words)}
    h = sp.zeros(len(words))
    terms = []
    # Use the independently deposited integer 4H on all 99 native states.
    raw = formation.emission.four_h(k)
    positions = [
        formation.emission.LAW.SRC.index((a, b)) for a in (0, 1) for b in (0, 1)
    ]
    for a, b in CELL_SEAMS:
        term = sp.zeros(len(words))
        for j, word in enumerate(words):
            rest = set(word) - {a, b}
            before = 2 * (a in word) + (b in word)
            for after, v in enumerate(positions):
                coefficient = int(raw[v, positions[before]])
                if coefficient:
                    out = (
                        rest
                        | ({a} if after // 2 else set())
                        | ({b} if after % 2 else set())
                    )
                    term[idx[frozenset(out)], j] = coefficient
        h += term
        terms.append(term)
    return words, h, terms


def project(h, starts):
    """p has the distinct roots of 4H. Sum all same-energy blocks exactly."""
    p = h.charpoly().as_poly().sqf_part().set_domain(sp.QQ)
    order = p.degree()
    t = p.gen
    inverse = sp.invert(p.diff(), p)
    # Lagrange interpolation: sum f(lambda)/p'(lambda) equals the
    # coefficient of t^(m-1) in f mod p. Here we need 1/p'(lambda)^2.
    shifts = []
    power = inverse
    for _ in range(2 * order - 1):
        shifts.append(power.nth(order - 1))
        power = (power * sp.Poly(t, t, domain=sp.QQ)).rem(p)
    weights = sp.Matrix(order, order, lambda i, j: shifts[i + j])
    coefficients = p.all_coeffs()
    result = sp.zeros(h.rows)
    for start in starts:
        v = sp.eye(h.rows)[:, start]
        lists = [v]
        for c in coefficients[1:-1]:
            lists.append(h * lists[-1] + c * v)
        assert h * lists[-1] + coefficients[-1] * v == sp.zeros(h.rows, 1)
        # The residue of (t-H)^(-1)v at lambda is W[lambda^j]/p'(lambda).
        W = sp.Matrix.hstack(*reversed(lists))
        result += W * weights * W.T
    result /= len(starts)
    assert h * result == result * h and result.trace() == 1
    return result, order


def run():
    adj = {i: set() for i in range(12)}
    for a, b in CELL_SEAMS:
        adj[a].add(b)
        adj[b].add(a)
    pairs = [(a, b) for a, b in combinations(sorted(adj[0]), 2) if b not in adj[a]]
    inside = {0}
    regions = []
    while True:
        regions.append(inside)
        if len(inside) == 12:
            break
        inside = inside.union(*(adj[x] for x in inside))
    rows = []
    for k in (8, 0):
        sectors = []
        degrees = []
        for count in (1, 2):
            words, h, terms = exact_sector(count, k)
            starts = (
                [words.index((x,)) for pair in pairs for x in pair]
                if count == 1
                else [words.index(p) for p in pairs]
            )
            rho, m = project(h, starts)
            degrees.append(m)
            # The two separate one-load configurations carry coefficient +1 each.
            factor = 2 if count == 1 else -1
            energy = [
                factor * sum(t[i, j] * rho[j, i] for (i, j) in t.todok()) / 4
                for t in terms
            ]
            occ = [
                factor * sum(rho[i, i] for i, w in enumerate(words) if x in w)
                for x in range(12)
            ]
            sectors.append((energy, occ))
        energy = [a + b for a, b in zip(sectors[0][0], sectors[1][0], strict=True)]
        occupancy = [a + b for a, b in zip(sectors[0][1], sectors[1][1], strict=True)]
        assert sum(energy) == sum(occupancy) == 0
        shells = []
        for r, s in enumerate(regions):
            val = sum(
                (int(a in s) + int(b in s)) * v / 2
                for (a, b), v in zip(CELL_SEAMS, energy, strict=True)
            )
            shells.append(
                dict(radius=r, energy=str(val), load=str(sum(occupancy[x] for x in s)))
            )
        rows.append(dict(stiffness=k, distinct_sector_energies=degrees, shells=shells))
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = run()
    path = HERE / "exact.json"
    if args.write:
        path.write_text(json.dumps(data, indent=2) + "\n")
    else:
        assert json.loads(path.read_text()) == data
