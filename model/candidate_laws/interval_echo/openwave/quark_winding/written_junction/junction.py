"""Native written seam charges, transport rules and energy of a three-arm source."""

import argparse
from itertools import combinations
import json
from pathlib import Path

import networkx as nx
import numpy as np
import sympy as sp
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.interval_echo import currents
from model.executable_model import carrier
from ..belt_walls import belt


def written_blocks(kappa):
    # Lower right block of U G U*: the native integer table is sqrt(2) W.
    w = sp.Matrix(law._write_pattern())
    q = sp.eye(w.rows) - w * w.T / 2
    r, d = map(sp.Matrix, law._label_ops())
    g = sp.kronecker_product((sp.eye(r.rows) - r) / 2 + kappa * d, sp.eye(9))
    h = q * g * q
    result = {}
    for magnitude in (0, 1, 2):
        tags = [i for i, (_, d) in enumerate(law.LABELS) if abs(d) == magnitude]
        size = len(tags)
        coupling = []
        current_entries = []
        for i, lo in enumerate(law.SRC):
            for j, hi in enumerate(law.SRC):
                block = h.extract([9 * m + i for m in tags], [9 * m + j for m in tags])
                if block != sp.zeros(size):
                    coupling.append(block)
                    if lo[1] != hi[1]:
                        current_entries.extend(block.todok().values())
        constraints = (
            sp.Matrix.vstack(
                *(
                    sp.kronecker_product(b, sp.eye(size))
                    - sp.kronecker_product(sp.eye(size), b.T)
                    for b in coupling
                )
            )
            if coupling
            else sp.zeros(0, size**2)
        )
        basis = [v.reshape(size, size) for v in constraints.nullspace()]
        result[str(magnitude)] = dict(
            memory_dimension=size,
            commutant_dimension=len(basis),
            commuting_algebra_abelian=all(
                a * b == b * a for a, b in combinations(basis, 2)
            ),
            crossing_matrix_element_magnitudes=sorted(
                {float(abs(v)) for v in current_entries}
            ),
        )
        if magnitude == 1:
            indices = [9 * m + i for m in tags for i in range(9)]
            local = h.extract(indices, indices)
            charge = sp.diag(*(law.LABELS[m][0] + b for m in tags for _, b in law.SRC))
            result["1"]["right_material_plus_coin_commutes"] = (
                local * charge == charge * local
            )
    return result


def arms(board, length):
    graph = nx.Graph((a // 12, b // 12) for a, b in board.external)
    # Three disjoint geodesics at the selected hub, with mutually distant endpoints.
    paths = nx.single_source_shortest_path(graph, 0)
    groups = [
        [
            p
            for _, p in sorted(paths.items())
            if len(p) == length + 1 and p[1] == neighbour
        ]
        for neighbour in sorted(graph[0])
    ]
    for a in groups[0]:
        for b in groups[1]:
            if set(a[1:]) & set(b[1:]):
                continue
            for c in groups[2]:
                if set(c[1:]) & (set(a[1:]) | set(b[1:])):
                    continue
                ends = (a[-1], b[-1], c[-1])
                if (
                    min(
                        nx.shortest_path_length(graph, i, j)
                        for i, j in combinations(ends, 2)
                    )
                    >= length
                ):
                    return (a, b, c)
    raise ValueError("No three disjoint rays")


def terminal(board, cell, n):
    axes = belt.belts()
    gateways = {p for p in range(12 * cell, 12 * cell + 12) if board.degree[p] == 6}
    for order in axes.values():
        for sites in combinations(range(10), n):
            count = sum(i in sites and j in sites for i, j in belt.ring(10).seams)
            word = [cell * 12 + order[i] for i in sites]
            if count == (0 if n == 3 else 1) and gateways.isdisjoint(word):
                return word
    raise ValueError("No wall avoiding outgoing tips")


def records(board, rays, kappa, contact):
    graph = nx.Graph((a // 12, b // 12) for a, b in board.external)
    bridge = {tuple(sorted((a // 12, b // 12))): (a, b) for a, b in board.external}
    pathsets = {"geodesic": rays}
    for step in reversed(range(len(rays[0]) - 1)):
        a, b = rays[0][step : step + 2]
        g = graph.copy()
        g.remove_edge(a, b)
        g.remove_nodes_from(set(sum((list(r) for r in rays), [])) - {a, b})
        if nx.has_path(g, a, b):
            bypass = nx.shortest_path(g, a, b)
            pathsets["detour"] = [
                rays[0][:step] + bypass + rays[0][step + 2 :],
                *rays[1:],
            ]
            break
    ends = [p[-1] for p in rays]
    result = {}
    for counts in ((3, 3, 4), (3, 4, 4)):
        shapes = [terminal(board, cell, n) for cell, n in zip(ends, counts)]
        sites = set(sum(shapes, []))
        # Every blank edge and the actual degree references on the whole board.
        blank = sum(
            (kappa + 0.5) for a, b in board.seams if (a in sites) != (b in sites)
        ) + sum(kappa * (6 - board.degree[i]) for i in sites)
        panels = {}
        for label, paths in pathsets.items():
            links = [
                bridge[tuple(sorted((a, b)))]
                for path in paths
                for a, b in zip(path, path[1:])
            ]
            assert len(set(links)) == len(links) and all(
                sites.isdisjoint(e) for e in links
            )
            panels[label] = dict(
                paths=paths,
                written_links=len(links),
                blank_energy=blank,
                written_energy=blank + contact * len(links),
                written_zero_energy=blank,
            )
        result[str(sum(counts))] = dict(terminal_words=shapes, panels=panels)
    return result


def run():
    board = carrier.federation(2)
    x, c = law.local_term_parts()
    families = {}
    for kappa in (8, 0):
        h = x + kappa * c
        idx = 9 + 9 * law.LABELS.index((0, 1)) + law.SRC.index((0, 0))
        zero = 9 + 9 * law.LABELS.index((0, 0)) + law.SRC.index((0, 0))
        assert abs(h[zero, zero]) < 1e-12
        arm = np.eye(99)[:, idx]
        times = (0.0, 1 / 16, 0.25, 1.0)
        # Relative phase against an idle equal arm, measured at two analyser phases.
        readings = []
        for t in times:
            u = expm(-1j * t * h)
            overlap = np.vdot(arm, u @ arm)
            readings.append(
                dict(t=t, x=(1 + overlap.real) / 2, y=(1 + overlap.imag) / 2)
            )
        families[str(kappa)] = dict(
            memory=written_blocks(kappa),
            sign_mark_cost=h[idx, idx],
            local_phase_reads=readings,
            distances={
                str(l): records(board, arms(board, l), kappa, h[idx, idx])
                for l in (1, 2, 4, 8)
            },
        )
    return currents._portable(dict(board_cells=board.cells, cases=families))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    data = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    with args.output.open("x") as f:
        f.write(data)
