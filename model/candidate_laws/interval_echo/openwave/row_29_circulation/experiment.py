"""Source-sign tomography of native material circulation on connected cell guides."""

import argparse
from itertools import combinations
import json
from pathlib import Path

import networkx as nx
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh, expm_multiply

from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.electron_field.replay import compare
from model.executable_model import carrier

HERE = Path(__file__).parent
TIMES = (0, 1, 2, 4, 8, 16, 24, 32)
SHOTS = 2**24
CURRENT = np.sqrt(3) / 2


def laboratories():
    board = carrier.federation(0)
    graph = nx.Graph((a // 12, b // 12) for a, b in board.external)
    a, b = min(tuple(sorted(e)) for e in graph.edges)
    cycles = set()
    # For each deleted edge, all shortest paths between its ends give its shortest cycles.
    for edge in list(graph.edges):
        graph.remove_edge(*edge)
        for path in nx.all_shortest_paths(graph, *edge):
            cycles.add(
                min(
                    tuple(x[i:] + x[:i])
                    for x in (path, path[::-1])
                    for i in range(len(x))
                )
            )
        graph.add_edge(*edge)
    loop = min(cycles, key=lambda p: (len(p), p))
    return ((a,), (a, b), loop)


def guide(owners):
    ids = [12 * i + j for i in owners for j in range(12)]
    index = {p: i for i, p in enumerate(ids)}
    seams = [
        (index[a], index[b])
        for a, b in carrier.federation(0).seams
        if a in index and b in index
    ]
    return carrier.Carrier("native opened cells", len(ids), seams)


def select_source(board):
    cell = carrier.cell()
    graph = nx.Graph(cell.seams)
    sym = nx.algorithms.isomorphism.GraphMatcher(graph, graph)
    gateways = {a for e in board.external for a in e if a < 12}
    stabilizer = [
        p for p in sym.isomorphisms_iter() if all(p[v] == v for v in gateways)
    ]
    witnesses = {}
    for face in cell.faces():
        a, b, c = face
        reverse = {(c, b, a), (b, a, c), (a, c, b)}
        witnesses[face] = [
            p for p in stabilizer if tuple(p[v] for v in face) in reverse
        ]
    source = next(
        (face for face in cell.faces() if not witnesses[face]), cell.faces()[0]
    )
    mirror = witnesses[source]
    return (
        source,
        sorted(gateways),
        [mirror[0][v] for v in range(12)] if mirror else None,
    )


def circulation(words, index, face):
    row, col, sgn = [], [], []
    for n, w in enumerate(words):
        occ = set(w)
        for a, b in zip(face, (*face[1:], face[0])):
            if (a in occ) != (b in occ):
                row.append(index[tuple(sorted(occ ^ {a, b}))])
                col.append(n)
                sgn.append(1 if a in occ else -1)
    return sparse.csr_matrix(
        (0.5j * np.asarray(sgn), (row, col)), shape=(len(words), len(words))
    )


def acquire(board, load, k, source, dest, disconnected, rng):
    words = list(combinations(range(board.n), load))
    index = {w: i for i, w in enumerate(words)}
    pair = PairSector(board) if load == 2 else None

    def hamiltonian(g):
        if load == 1:
            return g.laplacian() / 2
        return (pair if g is board else PairSector(g)).hamiltonian(k)

    h = hamiltonian(board)
    if k and load == 2:
        ev, vec = eigsh(h, k=2, which="SA", tol=2e-13, v0=np.ones(len(words)))
        ground = vec[:, 0]
        if ground.sum() < 0:
            ground = -ground
        assert np.linalg.norm(h @ ground - ev[0] * ground) < 1e-10
    else:
        ground = np.ones(len(words)) / np.sqrt(len(words))
    energy = float(ground @ (h @ ground))
    src = circulation(words, index, source)
    norm = float(np.linalg.norm(src @ ground) ** 2)
    launch = np.column_stack(
        [expm_multiply((-1j * sign) * src, ground) for sign in (1, -1)]
    )
    if disconnected:
        islands = carrier.Carrier(
            "disconnected",
            board.n,
            [(a, b) for a, b in board.seams if a // 12 == b // 12],
        )
        h = hamiltonian(islands)
    if norm < 1e-24 and not disconnected:
        # An exactly dark source on a stationary state only adds a common phase.
        waves = np.broadcast_to(launch, (len(TIMES), *launch.shape))
    else:
        waves = expm_multiply(-1j * h, launch, start=0, stop=32, num=33)[list(TIMES)]
    hw = waves.transpose(1, 0, 2).reshape(len(words), -1)
    reads = []
    for face in sorted(
        {source} | {face for face in board.faces() if face[0] // 12 == dest}
    ):
        j = circulation(words, index, face)
        jx = j @ hw
        mean = np.sum(hw.conj() * jx, axis=0).real.reshape(-1, 2) / CURRENT
        nonzero = np.sum(abs(jx) ** 2, axis=0).reshape(-1, 2) / CURRENT**2
        probs = np.stack(
            ((nonzero + mean) / 2, 1 - nonzero, (nonzero - mean) / 2), axis=-1
        )
        assert np.max(abs(probs.sum(axis=-1) - 1)) < 1e-12 and probs.min() > -1e-11
        counts = [
            [rng.multinomial(SHOTS, np.clip(p.round(12), 0, 1)).tolist() for p in row]
            for row in probs
        ]
        reads.append(
            dict(face=face, probabilities=probs.round(12).tolist(), counts=counts)
        )
    return dict(
        load=load,
        contrast=k,
        bridge_cut=disconnected,
        relative_ground_energy=energy,
        initial_adjacent_probability=float(pair.adjacent @ (ground**2))
        if pair
        else None,
        source_weight=norm,
        reads=reads,
    )


def run():
    rng = np.random.default_rng(291124)
    labs = []
    for owners in laboratories():
        board = guide(owners)
        source, gates, mirror = select_source(board)
        dest = len(owners) // 2
        cases = [
            (1, candidate.KAPPA, False),
            (2, 0, False),
            (2, candidate.KAPPA, False),
        ]
        if len(owners) > 2:
            cases.append((2, candidate.KAPPA, True))
        labs.append(
            dict(
                owners=owners,
                source=source,
                gateway_ports=gates,
                reflection=mirror,
                receiving_owner=owners[dest],
                experiments=[
                    acquire(board, n, k, source, dest, c, rng) for n, k, c in cases
                ],
            )
        )
    contrasts = sum(
        len(e["reads"]) * len(TIMES) for lab in labs for e in lab["experiments"]
    )
    return dict(
        shots=SHOTS,
        times=TIMES,
        simultaneous_contrasts=contrasts,
        simultaneous_error=float(np.sqrt(4 * np.log(2 * contrasts / 0.01) / SHOTS)),
        labs=labs,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    data = run()
    if args.output:
        with args.output.open("x") as f:
            json.dump(data, f, separators=(",", ":"), allow_nan=False)
            f.write("\n")
    else:
        compare(
            json.loads(json.dumps(data)),
            json.loads((HERE / "results.json").read_text()),
            "circulation",
        )
        print("Verified all signed circulation readings and the matched controls.")


if __name__ == "__main__":
    main()
