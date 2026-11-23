"""Spectral enclosures and Koide screening on native single-charge scalar guides."""

import argparse
from collections import Counter
from fractions import Fraction
from itertools import combinations
import json
from math import ceil, sqrt
from pathlib import Path
import sys

import networkx as nx
import numpy as np
from scipy.linalg import eigvalsh

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[6]))
from model.executable_model.carrier import federation, seam_digest

HOP = 4
PIN = 8


def matrices(weights, pinned, edges, denominator):
    """Clear norm shares and their corresponding onsite energies by the same denominator."""
    mass = np.diag(weights)
    stiffness = np.diag([PIN*w if i in pinned else 0 for i, w in enumerate(weights)])
    for a, b in edges:
        stiffness[a, a] += denominator*HOP
        stiffness[b, b] += denominator*HOP
        stiffness[a, b] -= denominator*HOP
        stiffness[b, a] -= denominator*HOP
    return stiffness, mass


def blocks(limit, kind):
    for k in range(limit + 1):
        if kind == "cycle":
            yield (k,), matrices([1] + [2]*k + [1], {0, k+1}, [(i, i+1) for i in range(k+1)], 2)
        elif kind == "end":
            yield (k,), matrices([1] + [2]*k, {0}, [(i, i+1) for i in range(k)], 2)
        elif k:
            for attach in range(k):
                edges = [(i, i+1) for i in range(1, k)] + [(0, 1+attach)]
                yield (k, attach), matrices([1] + [5]*k, {0}, edges, 5)


def certify_positive(stiffness, mass, trial):
    """Rational symmetric elimination. Positive pivots prove strict positivity."""
    integer_form = stiffness*trial.denominator - mass*trial.numerator
    a = [[Fraction(int(x)) for x in row] for row in integer_form]
    for j in range(len(a)):
        assert a[j][j] > 0
        for p in range(j+1, len(a)):
            if a[p][j]:
                scale = a[p][j]/a[j][j]
                for q in range(p, len(a)):
                    a[p][q] -= scale*a[j][q]
                    a[q][p] = a[p][q]


def enclosure(local, ceiling):
    cases = list(local)
    value, worst = min((float(eigvalsh(a, b)[0]), indices) for indices, (a, b) in cases)
    trial = Fraction(ceil(value*1000)-1, 1000)
    for _, (a, b) in cases:
        certify_positive(a, b, trial)
    return {"block_count": len(cases), "lowest_block": list(worst), "lowest_block_energy": value,
            "certified_strict_lower_bound": str(trial), "upper_bound": ceiling,
            "certified_ratio_ceiling": float(ceiling/trial)}


def koide(energies):
    return float(sum(energies)/sum(sqrt(x) for x in energies)**2)


def one_pin_paths(board, high, degree_six):
    """Count all distinct simple physical paths starting with their sole pinned site."""
    count = Counter()
    examples = {}
    def visit(path):
        k = len(path)-1
        count[k] += 1
        examples.setdefault(k, list(path))
        for v in sorted(high[path[-1]]):
            if v not in path:
                visit((*path, v))
    for p, neighbours in enumerate(board.neighbours):
        if p not in degree_six:
            for q in neighbours:
                if q in degree_six:
                    visit((p, q))
    sampled = {}
    best = None
    for k in sorted(count):
        stiffness, norm = matrices([1]*(k+1), {0}, [(i, i+1) for i in range(k)], 1)
        energies = eigvalsh(stiffness, norm)
        assert abs(eigvalsh(stiffness - PIN*np.diag([1]+[0]*k))[0]) < 1e-12
        triples = []
        for triple in combinations(range(k+1), 3):
            e = energies[list(triple)]
            q = koide(e)
            triples.append((abs(q-2/3), triple, q, e))
        if triples:
            error, indices, q, e = min(triples, key=lambda x: (x[0], x[1]))
            entry = {"active_ports": examples[k], "indices": list(indices), "energies": e.tolist(),
                     "ratios": (e/e[0]).tolist(), "Q": q, "absolute_delta_from_balance": error}
            if best is None or error < best["absolute_delta_from_balance"]:
                best = entry
        sampled[k] = {"physical_paths": count[k], "examples": examples[k], "all_triples": len(triples),
                      "E_min": float(energies[0]), "E_max": float(energies[-1])}
    return {"physical_paths_count": sum(count.values()), "inequivalent_spectra": len(count),
            "all_distinct_template_triples_tested": sum(s["all_triples"] for s in sampled.values()),
            "by_unpinned_length": sampled, "closest_to_balance": best}


def census(level):
    board = federation(level)
    assert set(board.degree) == {5, 6}
    degree_six = set(int(v) for v in np.flatnonzero(board.degree == 6))
    high = nx.Graph()
    high.add_nodes_from(degree_six)
    high.add_edges_from((u, v) for u, v in board.seams if u in degree_six and v in degree_six)
    assert max(dict(high.degree).values()) <= 2
    sizes = Counter(len(c) for c in nx.connected_components(high))
    cycle_sizes = sorted(len(c) for c in nx.connected_components(high) if high.subgraph(c).number_of_edges() == len(c))
    limit = max(sizes)
    middle = list(blocks(limit, "cycle"))
    ends = list(blocks(limit, "end"))
    branches = list(blocks(limit, "branched"))
    return {"cells": board.cells, "physical_ports": board.n, "seam_sha256": seam_digest(board),
            "degree_six_component_sizes": dict(sorted(sizes.items())), "all_degree_six_cycle_sizes": cycle_sizes,
            "largest_consecutive_unpinned_run_bound": limit,
            "maximum_distinct_cells_in_a_degree_six_component":
                max(len({v//12 for v in component}) for component in nx.connected_components(high)),
            "pinned_simple_cycles": enclosure(middle, 24),
            "pinned_open_paths": enclosure(middle+ends, 24),
            "pinned_connected_subgraphs": enclosure(branches, 56),
            "untwisted_unpinned_ground": 0,
            "untwisted_unpinned_positive_ratio_upper": 1/(np.tan(np.pi/(2*limit))**2),
            "one_pin_path_screen": one_pin_paths(board, high, degree_six)}


def main():
    args = argparse.ArgumentParser(description=__doc__)
    args.add_argument("--output", type=Path, default=HERE/"results.json")
    out = args.parse_args().output
    result = {"construction": "one scalar configuration per physical port; native degree reference one; fixed external compensator",
              "minimum_max_to_min_ratio_for_exact_positive_Koide": (4+3*sqrt(2))**2,
              "federations": [census(i) for i in range(3)]}
    out.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
