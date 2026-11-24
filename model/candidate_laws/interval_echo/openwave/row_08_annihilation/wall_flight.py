"""Far-end calorimetry and material-wall classification after closed-ring extinction."""

import argparse
from itertools import product
import json
import networkx as nx
import numpy as np
from scipy.linalg import eigh
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import expm_multiply
from model.executable_model import federation
from . import walls, experiment as ex

N = 3
TIMES = (0, 1, 2, 4, 8, 16, 32)


def extended_geometry():
    ring, ap, wings = walls.geometry()
    board = federation(0)
    allsites = set(ring + wings[0] + wings[1])
    graph = nx.Graph(board.seams)
    out = []
    for side in wings:
        opts = []
        for a in board.neighbours[side[-1]]:
            for b in board.neighbours[a]:
                if set((a, b)) & allsites or board.degree[a] != board.degree[b]:
                    continue
                opts.append(
                    (-min(nx.shortest_path_length(graph, p, b) for p in ring), a, b)
                )
        _, a, b = min(opts)
        out.append(side + [a, b])
        allsites.update((a, b))
    return ring, ap, out


def regions(k, start):
    states = list(product((0, 1), repeat=N))
    d = len(states)
    mat = np.zeros((d, d))
    local = ex.term(k)
    for j, row in enumerate(states):
        w = walls.loads(row)
        for i in range(start, N):
            a, b = w[2 * i : 2 * i + 2]
            col = ex.LAW.SRC.index((a, b))
            mat[j, j] += local[col, col]
            idx = states.index(tuple(1 - x if m == i else x for m, x in enumerate(row)))
            mat[idx, j] += local[ex.LAW.SRC.index((b, a)), col]
        for i in range(start, N):
            a, b = w[2 * i - 1 : 2 * i + 1]
            col = ex.LAW.SRC.index((a, b))
            mat[j, j] += local[col, col] * (0.5 if i == start else 1)
    return mat


def decode(psi, ops, rng):
    values = []
    spectra = []
    for i, op in enumerate(ops):
        energy, table = eigh(op)
        p = sum(
            (
                np.sum(abs(table.T @ a) ** 2, axis=1)
                if i == 0
                else np.sum(abs(a @ table) ** 2, axis=0)
            )
            for a in psi
        )
        # Group degenerate levels before drawing. An energy meter never selects their basis.
        levels, idx = np.unique(energy.round(9), return_inverse=True)
        weight = np.bincount(idx, weights=p)
        record = walls.sample(np.r_[weight, 1 - weight.sum()], rng)
        mean = np.dot(levels, record["counts"][:-1]) / sum(record["counts"][:-1])
        spectra.append(
            dict(
                energies=levels.tolist(),
                measurement=record,
                measured_conditional_mean=mean,
            )
        )
        values.append(float(energy @ p / p.sum()))
    return dict(conditional_means=values, spectra=spectra)


def run():
    rng = np.random.default_rng(81129)
    ring, ports, wings = extended_geometry()
    h = walls.operators(8, False, N)
    h = csr_matrix(h - np.eye(len(h)) * np.mean(h.diagonal()))
    ref = walls.reference(0)
    bits = list(product((0, 1), repeat=4))
    order = [bits.index(ref), bits.index(tuple(1 - x for x in ref))]
    q = np.asarray(walls.loads(ref))[ports]
    init = [walls.wing_ref(8, b, False, N) for b in q]
    windows = [regions(8, i) for i in (1, 2)]
    cases = []
    for wall in (False, True):
        seq = tuple(z ^ int(wall and i >= 2) for i, z in enumerate(ref))
        x = np.zeros((16, 8, 8))
        x[bits.index(seq)] = np.outer(init[0][2], init[1][2])
        evolved = expm_multiply(-32j * h, x.ravel()).reshape(16, 8, 8)
        found = evolved[order]
        p = float(np.vdot(found, found).real)
        frequencies = [
            [eigh(ex.wing(8, N, b if i == 0 else 2 - b)) for b in q] for i in (0, 1)
        ]
        flags = []
        for arr, pair in zip(found, frequencies, strict=True):
            bright = [np.eye(8) - np.outer(v[:, 0], v[:, 0]) for _, v in pair]
            flags.append(ex.read_joint(arr, *bright))
        flags = walls.sample(np.r_[np.sum(flags, axis=0), 1 - p], rng)
        states = list(product((0, 1), repeat=3))
        # Each internal 00/22 contact (different adjacent dimer bits) is ordered.
        walls_by_word = np.array(
            [sum(a == b for a, b in zip(s, s[1:])) for s in states]
        )
        times = []
        for t in TIMES:
            psi = []
            for a, pair in zip(found, frequencies, strict=True):
                unitaries = [(v * np.exp(-1j * t * (e - e[0]))) @ v.T for e, v in pair]
                psi.append(unitaries[0] @ a @ unitaries[1].T)
            qp = sum(abs(a) ** 2 for a in psi)
            walls_count = np.bincount(
                (walls_by_word[:, None] + walls_by_word[None, :]).ravel(),
                weights=qp.ravel(),
                minlength=5,
            )
            times.append(
                dict(
                    t=t,
                    distal=[decode(psi, [op, op], rng) for op in windows],
                    wing_kink_count=walls.sample(np.r_[walls_count, 1 - p], rng),
                )
            )
        cases.append(
            dict(
                wall_pair=wall,
                ring_empty_probability=p,
                both_wing_excitation=flags,
                times=times,
            )
        )
    return dict(
        ring=ring,
        attachments=ports,
        wings=wings,
        exchange_time=32,
        seed=81129,
        shots=walls.SHOTS,
        cases=cases,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    d = json.loads(json.dumps(run()))
    p = ex.HERE / "wall_flight.json"
    if args.write:
        p.write_text(json.dumps(d, indent=1, allow_nan=False) + "\n")
    else:
        from .verify import compare

        compare(d, json.loads(p.read_text()), "wall_flight")
        print(
            "Verified distant energy transport after the freely disappearing phase walls."
        )


if __name__ == "__main__":
    main()
