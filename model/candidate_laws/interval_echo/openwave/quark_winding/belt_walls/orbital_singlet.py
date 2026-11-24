"""Collective singlets, local charge and interaction energy in compact native orbitals."""

import argparse
from itertools import combinations, permutations, product
import json
from pathlib import Path

import networkx as nx
import numpy as np
import sympy as sp
from scipy import sparse

from model.executable_model import carrier, rotations
from model.candidate_laws.interval_echo import candidate, currents
from . import localized


def compact_levels(gates):
    words, exits, exchange, _ = localized.sector(gates, 0, 1)
    h = sp.Matrix(exchange.toarray())
    boundary = sp.zeros(len(exits), len(words))
    for i, row in enumerate(exits):
        for j in row:
            boundary[i, j] = 1
    witness = boundary
    for _ in range(len(words) - 1):
        witness = witness.col_join(witness[-boundary.rows :, :] * h)
    basis = sp.Matrix.hstack(*witness.nullspace())
    active = (basis.T * basis).inv() * basis.T * h * basis
    assert basis * active == h * basis
    space = sp.Matrix.hstack(
        *(boundary.col_join(h - 6 * sp.eye(len(words)))).nullspace()
    )
    port = np.zeros((12, space.cols))
    for row, w in enumerate(words):
        port[w.index(2)] = np.array(space[row, :], dtype=float).ravel()
    frame = np.zeros_like(port)
    free = sorted(set(range(12)) - set(gates))
    frame[free] = np.linalg.qr(port[free])[0]
    return (
        dict(
            dimension=basis.cols,
            twice_kinetic_polynomial=str(sp.factor(active.charpoly().as_expr())),
            energy_three_dimension=space.cols,
        ),
        frame,
        space,
    )


def material_read_algebra(basis):
    g = basis.T * basis
    effects = [g.inv() * basis[i, :].T * basis[i, :] for i in range(basis.rows)]
    size = basis.cols
    # Rational matrices in a nonorthogonal basis represent the exact port effects.
    blocks = [
        sp.kronecker_product(sp.eye(size), e) - sp.kronecker_product(e.T, sp.eye(size))
        for e in effects
    ]
    return size**2 - sp.Matrix.vstack(*blocks).rank()


def placements(board):
    graph = nx.Graph((a // 12, b // 12) for a, b in board.external)
    actions = rotations()
    choices = {}
    for c in range(board.cells):
        gates = [i for i in range(12) if board.degree[12 * c + i] == 6]
        for g in actions:
            if sorted(g[i] for i in gates) == [0, 1, 2]:
                choices[c] = g
                break
    distance = dict(nx.all_pairs_shortest_path_length(graph))
    trips = [
        (tuple(sorted(distance[a][b] for a, b in combinations(cs, 2))), cs)
        for cs in combinations(sorted(choices), 3)
    ]
    near = min(trips)
    far = max(trips, key=lambda item: (item[0], tuple(-c for c in item[1])))
    return choices, (near, far)


def one_port_check(board, c, permutation, frame, q):
    slots = np.arange(12 * c, 12 * c + 12)
    modes = np.zeros((board.n, 3))
    modes[slots] = frame[list(permutation)]
    kinetic = board.laplacian() / 2
    potential = candidate.KAPPA * (6 * q * q + 2 * q * (6 - board.degree))
    e = float(3 + candidate.KAPPA * (6 * q * q + 2 * q))
    residual = (kinetic + sparse.diags(potential)) @ modes - e * modes
    return e, float(np.max(abs(residual)))


def joint_residual(board, owners, choices, frame, charges):
    """Apply the full native empty/occupied seam action to a two- or three-body singlet."""
    state = {}
    for sites in product(range(12), repeat=len(charges)):
        rows = np.stack(
            [frame[choices[c][s]] for c, s in zip(owners, sites, strict=True)]
        )
        amplitude = (
            np.linalg.det(rows) / np.sqrt(6)
            if len(charges) == 3
            else rows[0] @ rows[1] / np.sqrt(3)
        )
        if amplitude:
            key = tuple(12 * c + i for c, i in zip(owners, sites, strict=True))
            state[key] = amplitude
    energy = sum(3 + candidate.KAPPA * (6 * q * q + 2 * q) for q in charges)
    diff = {}
    for sites, amplitude in state.items():
        diagonal = 0.0
        for j, (site, q) in enumerate(zip(sites, charges, strict=True)):
            diagonal += board.degree[site] / 2 + candidate.KAPPA * (
                6 * q * q + 2 * q * (6 - board.degree[site])
            )
            for dest in board.neighbours[site]:
                # All nonzero singlet support has mutually nonadjacent defects.
                assert dest not in sites
                target = sites[:j] + (dest,) + sites[j + 1 :]
                diff[target] = diff.get(target, 0.0) - amplitude / 2
        diff[sites] = diff.get(sites, 0.0) + (diagonal - energy) * amplitude
    norm = sum(v * v for v in state.values())
    assert abs(norm - 1) < 1e-12
    return float(np.linalg.norm(list(diff.values())))


def detector():
    # Relative signs and common local frame belong to the preparation.
    triple = np.zeros((3, 3, 3))
    for i, j, k in permutations(range(3)):
        triple[i, j, k] = (-1) ** ((i > j) + (j > k) + (i > k)) / np.sqrt(6)
    meson = np.eye(3) / np.sqrt(3)
    turn = np.diag(np.exp(1j * np.array([np.pi / 3, 2 * np.pi / 3, -np.pi])))
    common = np.einsum("ai,bj,ck,ijk->abc", turn, turn, turn, triple)
    pair = turn @ meson @ turn.conj().T
    separate = np.einsum("ai,ijk->ajk", turn, triple)
    product = turn[0, 0] ** 3

    def bright(overlap):
        return float((1 + overlap.real) / 2)

    return dict(
        triple_singlet_overlap=float(abs(np.vdot(triple, common)) ** 2),
        pair_singlet_overlap=float(abs(np.vdot(meson, pair)) ** 2),
        triple_common_bright=bright(np.vdot(triple, common)),
        triple_one_end_bright=bright(np.vdot(triple, separate)),
        equal_colour_product_common_bright=bright(product),
        local_probabilities=np.einsum("ijk,ijk->i", triple, triple).tolist(),
    )


def run():
    geometries = {}
    for gates, boards in localized.owner_shapes().items():
        levels, frame, basis = compact_levels(gates)
        geometries[str(gates)] = dict(boards=boards, **levels)
        if basis.cols == 3:
            calibrated = frame
            local_commutant = material_read_algebra(basis)
    sectors = []
    for level in (1, 2):
        board = carrier.federation(level)
        options, layouts = placements(board)
        entries = []
        for distances, owners in layouts:
            units = {}
            for q in (-1, 1):
                reads = [
                    one_port_check(board, c, options[c], calibrated, q) for c in owners
                ]
                units[q] = reads[0][0]
                assert max(r for _, r in reads) < 1e-10
            # Inter-owner terms can differ from the sum of independent terms only
            # when both ends of an inter-owner seam are occupied by deviations.
            # Every gateway amplitude is zero in these compact triplets.
            gate_amplitude = max(
                np.max(abs(calibrated[options[c][i]]))
                for c in owners
                for i in range(12)
                if board.degree[12 * c + i] == 6
            )
            assert gate_amplitude < 1e-12
            entries.append(
                dict(
                    owners=owners,
                    distances=distances,
                    unit_energies=units,
                    triple_positive_energy=3 * units[1],
                    triple_negative_energy=3 * units[-1],
                    neutral_pair_energy=units[1] + units[-1],
                    triple_code_dimension=27,
                    pair_code_dimension=9,
                    energy_relative_to_independent_fragments=0.0,
                    maximum_gateway_amplitude=gate_amplitude,
                    triple_full_residual=joint_residual(
                        board, owners, options, calibrated, (1, 1, 1)
                    ),
                    opposite_pair_full_residual=joint_residual(
                        board, owners[:2], options, calibrated, (1, -1)
                    ),
                )
            )
        sectors.append(
            dict(board=board.cells, eligible_cells=len(options), placements=entries)
        )
    return currents._portable(
        dict(
            gates=geometries,
            local_port_read_commutant=local_commutant,
            detector=detector(),
            comparisons=sectors,
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run()
    with args.output.open("x") as f:
        f.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
