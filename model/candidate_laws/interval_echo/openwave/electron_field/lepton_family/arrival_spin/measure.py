"""Boundary-head herald followed by a coherent order test of two physical half-turns."""

import argparse
from itertools import combinations
import json
from pathlib import Path
import sys

import networkx as nx
import numpy as np
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[6]))
# The archive replay launches this script by filename from a temporary directory.
# ruff: noqa: E402
from model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.dirac_read import (
    measure as guide,
)
from model.executable_model import cell
from model.executable_model.carrier import icosahedron_from_coordinates, isomorphism
from model.executable_model.symmetry import rotations
from model.experiments.interval_echo_angular.charged_orientation.measure import (
    PAULI,
    compile_turn,
)

TIMES = (0, 4, 8, 12)
SHOTS = 65536


class WrittenTree:
    def __init__(self, board, route, sign):
        self.route = route
        self.sign = sign
        self.base = route[-1] // 12 * 12
        self.receiver = list(range(self.base, self.base + 12))
        local = nx.Graph([(a + self.base, b + self.base) for a, b in cell().seams])
        if not set(local.edges).issubset(
            set(board.seams) | {(b, a) for a, b in board.seams}
        ):
            raise ValueError("Receiver template fails physical seam incidence")
        branches = list(nx.bfs_edges(local, route[-1], sort_neighbors=sorted))
        self.tree = nx.Graph(branches)
        parents = list(zip(route, route[1:])) + branches
        self.ports = [route[0]] + [b for _, b in parents]
        if len(self.ports) != len(set(self.ports)) or set(route) & set(
            self.receiver
        ) != {route[-1]}:
            raise ValueError(
                "The charge guide and receiver must meet at a single entrance"
            )
        self.edge_bits = {tuple(sorted(edge)): 1 << i for i, edge in enumerate(parents)}
        initial = sum(
            self.edge_bits[tuple(sorted((a, b)))]
            for a, b in parents
            if (a < b) == (sign == 1)
        )
        self.word = {route[0]: initial}
        for a, b in parents:
            self.word[b] = self.word[a] ^ self.edge_bits[tuple(sorted((a, b)))]
        self.entrance = tuple(sorted(route[-2:]))
        # Every allowed sharp hop must change only its own physical head.
        reference = self.inward(initial)
        for p in self.ports:
            charge = self.inward(self.word[p]) - reference
            expected = sign * (
                (np.array(self.ports) == p).astype(int)
                - (np.array(self.ports) == route[0]).astype(int)
            )
            np.testing.assert_array_equal(charge, expected)
            for edge in self.edge_bits:
                target, out = self.contact(p, self.word[p], edge)
                if out != self.word[target] or (target != p) != (p in edge):
                    raise ArithmeticError(
                        "Written sector is not closed under each enabled contact"
                    )
        self.initial_inward = self.boundary_head(initial)

    def inward(self, word):
        heads = [a if word & bit else b for (a, b), bit in self.edge_bits.items()]
        return np.array([heads.count(p) for p in self.ports])

    def boundary_head(self, word):
        a, b = self.entrance
        endpoint = a if word & self.edge_bits[self.entrance] else b
        return int(endpoint in self.receiver)

    def herald(self, addresses):
        delta = np.array(
            [self.boundary_head(self.word[p]) - self.initial_inward for p in addresses]
        )
        np.testing.assert_array_equal(
            delta, self.sign * np.isin(addresses, self.receiver)
        )
        return delta == self.sign

    def contact(self, port, word, edge):
        a, b = edge
        bit = self.edge_bits[edge]
        if port in edge and ((port == a) == bool(word & bit)) == (self.sign == 1):
            return (b if port == a else a), word ^ bit
        return port, word

    def pulse_schedule(self, perm):
        pulses = []
        for i, j in compile_turn(perm):
            path = nx.shortest_path(self.tree, self.base + i, self.base + j)
            along = [tuple(sorted(e)) for e in zip(path, path[1:])]
            pulses.extend(along + along[-2::-1])
        return pulses

    def turn(self, perm):
        pulses = self.pulse_schedule(perm)
        op = np.zeros((12, 12))
        for i, start in enumerate(self.receiver):
            p, word = start, self.word[start]
            for e in pulses:
                p, word = self.contact(p, word, e)
            if p != self.base + perm[i] or word != self.word[p]:
                raise ArithmeticError(
                    "Turn misses its target charge or leaves an unmatched head"
                )
            op[p - self.base, i] = 1
        return op, len(pulses)


def turn_pairs(fixture):
    pts, edges = icosahedron_from_coordinates()
    xyz = pts[isomorphism(cell().seams, edges)]
    half = {}
    for perm in sorted(rotations()):
        if all(perm[k] == k for k in range(12)) or any(
            perm[perm[k]] != k for k in range(12)
        ):
            continue
        rot = np.linalg.lstsq(xyz, xyz[list(perm)], rcond=None)[0]
        axis = np.linalg.eigh((rot + rot.T) / 2)[1][:, -1]
        pos, count = fixture.turn(perm)
        spin = -1j * np.einsum("i,ijk->jk", axis, PAULI)
        half[perm] = (pos, axis, spin, count)
    pairs = []
    for (a, ax, s, _), (b, bx, t, _) in combinations(half.values(), 2):
        if abs(ax @ bx) < 1e-10:
            pairs.append((b @ a, a @ b, t @ s, s @ t))
    if len(half) != 15 or len(pairs) != 15:
        raise ArithmeticError("Incomplete set of icosahedral half-turns or axis pairs")
    return pairs, [{"permutation": list(k), "pulses": v[3]} for k, v in half.items()]


def instrument(local, pairs, rng):
    # The input amplitudes retain their unconditioned norm after the flux pointer clicks.
    arrival = float(np.vdot(local, local).real)
    spin_up = np.array([1.0, 0.0])
    rows = []
    for xy, yx, spin_xy, spin_yx in pairs:
        a, b = xy @ local, yx @ local
        records = {}
        for arm, other, label in (
            (a, b, "charge"),
            (a, a, "same_order"),
            (
                np.outer(a.ravel(), spin_xy @ spin_up),
                np.outer(b.ravel(), spin_yx @ spin_up),
                "supplied_spin",
            ),
        ):
            plus = float(np.linalg.norm((arm + other) / 2) ** 2)
            minus = float(np.linalg.norm((arm - other) / 2) ** 2)
            np.testing.assert_allclose(plus + minus, arrival, atol=2e-13, rtol=0)
            records[label] = np.array([plus, minus, max(0.0, 1 - arrival)])
        rows.append(records)
    panel = {}
    for label in rows[0]:
        probs = np.array([row[label] for row in rows])
        sampled = np.maximum(0.0, np.round(probs[0], 12))
        panel[label] = {
            "first_axis_pair": probs[0].tolist(),
            "min_over_axes": probs.min(axis=0).tolist(),
            "max_over_axes": probs.max(axis=0).tolist(),
            "simulated_counts_first_pair": rng.multinomial(
                SHOTS, sampled / sampled.sum()
            ).tolist(),
        }
    return panel


def acquire(h, init, locations, fixture, pairs, rng):
    # Columns are possible family-register states. No family reader operates at the receiver.
    dim = h.shape[0] // len(locations)
    select = fixture.herald(locations)
    inside = np.asarray(locations)[select]
    if len(inside) != len(set(inside)):
        raise ValueError("Capture needs all arriving grades resolved")
    values, vec = np.linalg.eigh(h)
    coeff = vec.conj().T @ init
    records = []
    for t in TIMES:
        psi = (vec @ (np.exp(-1j * t * values) * coeff)).reshape(len(locations), dim)
        np.testing.assert_allclose(np.vdot(psi, psi), 1.0, atol=3e-13)
        local = np.zeros((12, dim), complex)
        local[inside - fixture.base] = psi[select]
        records.append({"time": t, "fringe": instrument(local, pairs, rng)})
    return records


def laplacian(board, f, cut):
    h = np.diag(
        [float(guide.DEGREE_FIVE_COST * (6 - int(board.degree[p]))) for p in f.ports]
    )
    for e in f.edge_bits:
        if cut and e == f.entrance:
            continue
        i, j = [f.ports.index(p) for p in e]
        h[i, i] += guide.HOP
        h[j, j] += guide.HOP
        h[i, j] -= guide.HOP
        h[j, i] -= guide.HOP
    return h


def family_sources(bonds, f, effect, masses):
    n = len(f.route)
    _, chars = np.linalg.eigh(effect)
    kinetic, modes = np.linalg.eigh(bonds.T @ bonds)
    ground = np.r_[np.ones(n) / np.sqrt(n), np.zeros(len(bonds))]
    h = guide.generator(bonds, effect)
    for j, mass in enumerate(masses):
        energy = np.sqrt(kinetic[1] + mass**2)
        site = modes[:, 1] * np.sqrt((energy + mass) / (2 * energy))
        excited = np.r_[site, bonds @ site / (energy + mass)]
        first = np.kron(excited, chars[:, j])
        rest = np.kron(ground, chars[:, j])
        np.testing.assert_allclose(h @ first, energy * first, atol=2e-13)
        np.testing.assert_allclose(h @ rest, mass * rest, atol=2e-13)
        packet = site[-1] * rest - ground[n - 1] * first
        sharp = np.zeros_like(rest)
        sharp[:3] = chars[:, j]
        yield j, sharp, packet / np.linalg.norm(packet)


def run():
    board, route = guide.longest_unpinned_route()
    effect, masses = guide.register_effect()
    bonds = guide.path_incidence(board, route)
    h_family = guide.generator(bonds, effect)
    # A positive contact pulse of duration pi/(2t) realizes the sharp swap with zero stray phase.
    x = np.array([[0.0, 1.0], [1.0, 0.0]])
    np.testing.assert_allclose(expm(-0.5j * np.pi * (np.eye(2) - x)), x, atol=1e-14)
    rng = np.random.default_rng(191178)
    blocks = []
    for sign in (-1, 1):
        f = WrittenTree(board, route, sign)
        pairs, half = turn_pairs(f)
        source = np.eye(len(f.ports))[:, 0]
        rows = {}
        for cut in (False, True):
            rows["cut_written" if cut else "written"] = acquire(
                laplacian(board, f, cut), source, f.ports, f, pairs, rng
            )
        locations = list(route) + list(route[:-1])
        # Disconnect just the external incidence row; preserve its address and the dimensions.
        cut_bonds = bonds.copy()
        cut_bonds[-1] = 0
        cut_h = guide.generator(cut_bonds, effect)
        for j, sharp, packet in family_sources(bonds, f, effect, masses):
            for name, h, state in [
                ("sharp", h_family, sharp),
                ("positive", h_family, packet),
                ("cut_sharp", cut_h, sharp),
            ]:
                rows[f"family_{j}_{name}"] = acquire(h, state, locations, f, pairs, rng)
        blocks.append(
            {"sign": sign, "initial_inward_head": f.initial_inward, "records": rows}
        )
    return {
        "board_cells": board.cells,
        "route": list(route),
        "entrance": list(f.entrance),
        "receiver_tree": sorted([list(sorted(e)) for e in f.tree.edges]),
        "half_turns": half,
        "perpendicular_pairs": len(pairs),
        "outcomes": ["constructive", "destructive", "no signed flux"],
        "family_rest_effects_native": masses.tolist(),
        "shots": SHOTS,
        "seed": 191178,
        "sectors": blocks,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run()
    with args.output.open("x") as out:
        json.dump(result, out, indent=2, allow_nan=False)
        out.write("\n")
