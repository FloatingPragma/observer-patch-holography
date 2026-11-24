"""Bounds for a three-wall state supported on three whole owners after full release."""

import argparse
from itertools import combinations, permutations, product
import json
from pathlib import Path

import numpy as np
import sympy as sp

from model.candidate_laws.interval_echo import currents
from . import external_field, upper
from ...row_13_walls.energy import is_wall


def masks(k):
    return [sum(1 << i for i in sites) for sites in combinations(range(5), k)]


def dark(k, enabled):
    words = masks(k)
    if not enabled or not k:
        return words, sp.eye(len(words))
    parents = masks(k - 1)
    lower = sp.Matrix([[int((w ^ v).bit_count() == 1) for w in words] for v in parents])
    kernel = lower.nullspace()
    if not kernel:
        return words, sp.zeros(len(words))
    v = sp.Matrix.hstack(*kernel)
    proj = v * (v.T * v).inv() * v.T
    assert lower * proj == sp.zeros(len(parents), len(words)) and proj * proj == proj
    return words, proj


def interleave(a, b):
    return sum(
        ((a >> i) & 1) << (2 * i) | ((b >> i) & 1) << (2 * i + 1) for i in range(5)
    )


def local_ceiling(n, both):
    spectra = []
    maxima = []
    for k in range(n + 1):
        a, proj = dark(k, True)
        b, second = dark(n - k, both)
        indices = [
            (i, j)
            for (i, u), (j, v) in product(enumerate(a), enumerate(b))
            if is_wall(interleave(u, v))
        ]
        if not indices:
            continue
        # The nonzero singular values of the two projectors can be found
        # in the small, sharp-wall basis.
        gram = sp.Matrix(
            [
                [proj[i, other_i] * second[j, m] for other_i, m in indices]
                for i, j in indices
            ]
        )
        polynomial = gram.charpoly().as_poly()
        intervals = polynomial.intervals(eps=sp.Rational(1, 10**12))
        lo, hi = intervals[-1][0]
        maxima.append(hi)
        spectra.append(
            dict(
                visible_loads=k,
                dimension=len(indices),
                polynomial=str(sp.factor(polynomial.as_expr())),
                greatest_root_interval=[str(lo), str(hi)],
            )
        )
    return dict(local_ceiling=float(max(maxima)), exact_blocks=spectra)


def pole_escape(words, amplitudes):
    squared = []
    for cell in range(3):
        for parity in (0, 1):
            remaining, inputs = [], []
            for i in range(10 * cell + parity, 10 * cell + 10, 2):
                src = np.flatnonzero((words >> i) & 1)
                remaining.append(words[src] ^ (1 << i))
                inputs.append(amplitudes[src])
            _, inverse = np.unique(np.concatenate(remaining), return_inverse=True)
            out = np.bincount(inverse, weights=np.concatenate(inputs)) / 2
            squared.append(float(out @ out))
    return squared


def run():
    board, layout = external_field.geometry()
    exterior = [
        (v, w)
        for e in board.external
        for v, w in (e, e[::-1])
        if v // 12 in layout["owners"] and w // 12 not in layout["owners"]
    ]
    source_poles = {p for row in layout["poles"] for p in row}
    assert len(exterior) == len(source_poles) == 5
    assert {v for v, _ in exterior} == source_poles and len(
        {w for _, w in exterior}
    ) == 5
    for p in source_poles:
        sites = [s for s in board.neighbours[p] if s // 12 in layout["owners"]]
        assert len(sites) == 5 and all(s in layout["ports"] for s in sites)
    ceilings = {
        str(n): {
            key: local_ceiling(n, all_poles)
            for key, all_poles in [("two", True), ("one", False)]
        }
        for n in (3, 4)
    }
    cases = []
    for n in (10, 11):
        types = (3, 3, 4) if n == 10 else (3, 4, 4)
        bounds = []
        for order in sorted(set(permutations(types))):
            end = (
                ceilings[str(order[0])]["two"]["local_ceiling"]
                * ceilings[str(order[2])]["two"]["local_ceiling"]
            )
            bounds.append(
                dict(
                    loads=order,
                    all_poles_empty=end
                    * ceilings[str(order[1])]["two"]["local_ceiling"],
                    three_whole_owners=end
                    * ceilings[str(order[1])]["one"]["local_ceiling"],
                )
            )
        (words, _, _, _), (_, state, _) = upper.trial_state(
            upper.bins(), layout["bridges"], n
        )
        coeff = pole_escape(words, state)
        cases.append(
            dict(
                material=n,
                orders=bounds,
                all_poles_empty_upper=max(b["all_poles_empty"] for b in bounds),
                three_whole_owners_upper=max(b["three_whole_owners"] for b in bounds),
                trial_escape_per_pole=coeff,
                trial_escape_quadratic=sum(coeff),
                exchange_only_quadratic=sum(coeff),
                no_hopping_quadratic=0,
            )
        )
    return currents._portable(
        dict(exterior_edges=exterior, one_belt=ceilings, sectors=cases)
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    data = run()
    with args.output.open("x") as f:
        f.write(json.dumps(data, indent=2, allow_nan=False) + "\n")
