"""Independent sharp-word propagation from the pinned board file, with no producer imports."""
import argparse
import hashlib
import json
from itertools import combinations
from math import sqrt
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

HERE = Path(__file__).resolve().parent
BOARD = HERE.parents[2]/"openwave/interval_echo/row_04/FEDERATION_LEVEL0.npz"
BOARD_SHA = "99870398e3d516f899c600e9bf25ba53c7f513142ba8facbe1cedd94b23f915b"


def check(times):
    if hashlib.sha256(BOARD.read_bytes()).hexdigest() != BOARD_SHA:
        raise ValueError("Unexpected carrier bytes")
    with np.load(BOARD, allow_pickle=False) as data:
        edges = sorted(tuple(sorted(map(int, edge))) for edge in data["edges"])
    n = max(b for _, b in edges)+1
    neighbors = [set() for _ in range(n)]
    for a, b in edges:
        neighbors[a].add(b)
        neighbors[b].add(a)
    width = min(len(neighbors[a])+len(neighbors[b])-2 for a, b in edges)
    positions = list(combinations(range(n), 2))
    index = {pair: j for j, pair in enumerate(positions)}
    rows, cols = [], []
    for j, (a, b) in enumerate(positions):
        for moved, fixed in ((a, b), (b, a)):
            for dest in neighbors[moved]-{fixed}:
                rows.append(j)
                cols.append(index[tuple(sorted((dest, fixed)))])
    adjacency = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(index),)*2)
    exchange = (sparse.diags(np.asarray(adjacency.sum(axis=1)).ravel())-adjacency)/2
    spec = json.loads((HERE/"protocol.json").read_text(encoding="utf-8"))
    edge = tuple(spec["launch"])
    masked = np.array([len(set(pair).intersection(edge)) == 1 for pair in positions], dtype=float)
    adjacent = np.array([b in neighbors[a] for a, b in positions], dtype=float)
    # Every active spectator is odd. Each tagged sign instead has sharp-basis
    # diagonal 1/2. A word of weight r has w-r upward and r downward flips.
    material = exchange-sparse.diags(2*width*adjacent + width*masked/2)
    off = [-sqrt((r+1)*(width-r))/2 for r in range(width)]
    sign_hops = sparse.diags([off, off], [-1, 1], shape=(width+1,)*2)
    h = (sparse.kron(sparse.eye(width+1), material, format="csr")
         + sparse.kron(sign_hops, sparse.diags(masked), format="csr"))
    state = np.zeros(h.shape[0], complex)
    state[index[edge]] = 1
    receipt = json.loads((HERE/"results.json").read_text(encoding="utf-8"))
    seams = [e for e in edges if e[0]//12 == e[1]//12 == edge[0]//12]
    if seams != list(map(tuple, receipt["seams"])) or width != receipt["tagged_signs"]:
        raise ValueError("Detector or tagged register mismatch")
    home = [index[e] for e in seams]
    outside = np.ones(len(index), dtype=bool)
    outside[home] = False
    previous = 0
    errors = {}
    for t in times:
        step = -1j*(t-previous)
        state = expm_multiply(step*h, state, traceA=step*h.diagonal().sum())
        weights = abs(state.reshape(width+1, len(index)))**2
        observed = np.column_stack((weights[:, home], weights[:, outside].sum(axis=1)))
        target = np.array(receipt["appointments"][str(t)]["weight_and_home_or_quiet"])
        if target.shape != observed.shape or not np.isfinite(target).all() or np.any(target < 0):
            raise ValueError("Invalid stored joint distribution")
        errors[t] = float(abs(target-observed).sum())
        if abs(observed.sum()-1) > 1e-9 or errors[t] > 1e-9:
            raise ArithmeticError(f"Independent sharp-word distribution differs at {t}: {errors[t]}")
        previous = t
    return errors


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quick", action="store_true", help="check the complete first appointment")
    args = parser.parse_args()
    times = json.loads((HERE/"protocol.json").read_text(encoding="utf-8"))["exposures"]
    print(json.dumps(check(times[:1] if args.quick else times)))
