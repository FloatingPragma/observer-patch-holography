"""Static centre-constrained force between two neutral native defects."""

import argparse
from collections import deque
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from model.executable_model import federation
from model.experiments.interval_echo_encounter.conjugation import free_seam
from model.candidate_laws.contrast_memory_spherical_v1 import law
from .. import neutral_defect as nd
from . import measure


def cell_paths(board):
    edges = {}
    for a, b in board.external:
        edges[a // 12, b // 12] = (a, b)
        edges[b // 12, a // 12] = (b, a)
    paths = {0: [0]}
    queue = deque([0])
    while queue:
        cell = queue.popleft()
        for neighbour in sorted(b for a, b in edges if a == cell and b not in paths):
            paths[neighbour] = paths[cell] + [neighbour]
            queue.append(neighbour)
    return edges, {
        len(p) - 1: p for _, p in sorted(paths.items(), reverse=True) if len(p) > 1
    }


def confined(board, gate, x, y, h, reverse):
    host = gate // 12
    selected = np.flatnonzero((x // 12 == host) & (y // 12 == host))
    xx, yy = x[selected], y[selected]
    aa, bb = (yy, xx) if reverse else (xx, yy)
    close = np.isin(aa, [n for n in board.neighbours[gate] if board.degree[n] == 5])
    source = ((bb == gate) & close).astype(float)
    source /= np.linalg.norm(source)
    delta = (xx == gate).astype(int) - (yy == gate).astype(int)
    return h[selected][:, selected].tocsr(), source, delta


def read(h1, s1, h2, s2, dipole):
    independent = sparse.kronsum(h2, h1, format="csr")
    joint = independent + sparse.diags(dipole.ravel())
    source = np.kron(s1, s2)
    for t in (1 / 128, 1 / 256):
        baseline = np.kron(
            *(
                expm_multiply((-1j * t) * h, s, traceA=-1j * t * h.diagonal().sum())
                for h, s in ((h1, s1), (h2, s2))
            )
        )
        flight = expm_multiply(
            (-1j * t) * joint, source, traceA=-1j * t * joint.diagonal().sum()
        )
        amplitude = baseline.conj() @ flight
        yield t, (1 - amplitude.imag) / 2


def run():
    board = federation(0)
    x, y, _, _, h, _ = nd.generators(board)
    edges, paths = cell_paths(board)
    seam = free_seam(nd.KAPPA)
    local = np.array(
        [
            seam[law.SRC.index((l, r)), law.SRC.index((l, r))]
            for l in range(3)
            for r in range(3)
        ]
    ).reshape(3, 3)
    connected = local - local[:, 1:2] - local[1:2, :] + local[1, 1]
    records = []
    for distance, path in sorted(paths.items()):
        a, _ = edges[path[0], path[1]]
        _, b = edges[path[-2], path[-1]]
        h1, s1, origin = confined(board, a, x, y, h, False)
        for reverse in (False, True):
            h2, s2, target = confined(board, b, x, y, h, reverse)
            dipole = (
                connected[1 + origin[:, None], 1 + target[None, :]]
                if distance == 1
                else np.zeros((len(s1), len(s2)))
            )
            exact = float(s1**2 @ dipole @ s2**2)
            outcomes = []
            for j, (t, prob) in enumerate(read(h1, s1, h2, s2, dipole)):
                click = measure.read_ramsey(
                    prob, 2**36, [measure.SEED, 6, distance, int(reverse), j]
                )
                outcomes.append(
                    dict(
                        time=t,
                        **click,
                        sampled_connected_energy=(2 * click["click"] / 2**36 - 1) / t,
                        standard_error=2 * math.sqrt(prob * (1 - prob) / 2**36) / t,
                    )
                )
            records.append(
                dict(
                    distance=distance,
                    path=path,
                    gateways=[a, b],
                    reverse=reverse,
                    operator_norm=float(abs(dipole).max()),
                    expectation=exact,
                    counts=outcomes,
                )
            )
    return dict(
        protocol_sha256=hashlib.sha256(
            (measure.HERE / "FORCE_PROTOCOL.md").read_bytes()
        ).hexdigest(),
        connected_cross_seam=connected.tolist(),
        configurations_per_pair=132,
        records=records,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run()
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    for a in result["records"]:
        print(
            a["path"],
            a["reverse"],
            a["expectation"],
            a["counts"][1]["sampled_connected_energy"],
        )
