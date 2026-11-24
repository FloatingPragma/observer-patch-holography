"""Complete one-owner compact matter sectors on the three native federations."""

import argparse
from collections import Counter, defaultdict
from itertools import combinations
import json
from pathlib import Path

import numpy as np
import sympy as sp
from scipy import sparse

from model.executable_model import carrier, rotations
from model.candidate_laws.interval_echo import candidate

PRIME = 1009


def owner_shapes():
    group = rotations()
    out = defaultdict(list)
    for level in range(3):
        board = carrier.federation(level)
        counts = Counter()
        for cell in range(board.cells):
            gates = [p for p in range(12) if board.degree[12 * cell + p] == 6]
            assert len(gates) == 3
            counts[min(tuple(sorted(g[p] for p in gates)) for g in group)] += 1
        for shape, count in counts.items():
            out[shape].append([board.cells, count])
    return out


def sector(gates, a, b):
    """Symbols 0,1,2 denote host and the two deviations, independently of their actual values."""
    free = [i for i in range(12) if i not in gates]
    words = []
    for first in combinations(free, a):
        for second in combinations(sorted(set(free) - set(first)), b):
            w = bytearray(12)
            for i in first:
                w[i] = 1
            for i in second:
                w[i] = 2
            words.append(bytes(w))
    index = {w: i for i, w in enumerate(words)}
    exits = defaultdict(list)
    rows, cols = [], []
    diff = np.zeros((3, len(words)), dtype=np.int64)
    pairs = ((0, 1), (0, 2), (1, 2))
    edges = carrier.cell().seams
    for j, w in enumerate(words):
        for u, v in edges:
            if w[u] == w[v]:
                continue
            pair = tuple(sorted((w[u], w[v])))
            diff[pairs.index(pair), j] += 1
            target = bytearray(w)
            target[u], target[v] = target[v], target[u]
            key = bytes(target)
            if u in gates or v in gates:
                exits[key].append(j)
            else:
                rows.append(index[key])
                cols.append(j)
    width = len(words)
    exchange = sparse.diags(
        diff.sum(axis=0), dtype=np.int64, format="csr"
    ) - sparse.csr_matrix(
        (np.ones(len(rows), dtype=np.int64), (rows, cols)), shape=(width, width)
    )
    return words, list(exits.values()), exchange, diff


def boundary_kernel(rows, width):
    """Sparse row reduction modulo the prime, followed by back substitution."""
    pivots = {}
    for columns in rows:
        row = {j: 1 for j in columns}
        while row:
            i = min(row)
            if i not in pivots:
                inv = pow(row[i], -1, PRIME)
                pivots[i] = {k: v * inv % PRIME for k, v in row.items()}
                break
            coeff = row[i]
            for k, v in pivots[i].items():
                x = (row.get(k, 0) - coeff * v) % PRIME
                if x:
                    row[k] = x
                else:
                    row.pop(k, None)
    free = sorted(set(range(width)) - set(pivots))
    basis = np.zeros((width, len(free)), dtype=np.int64)
    basis[free] = np.eye(len(free), dtype=np.int64)
    for i in sorted(pivots, reverse=True):
        ks = [k for k in pivots[i] if k != i]
        if ks:
            basis[i] = (-np.array([pivots[i][k] for k in ks]) @ basis[ks]) % PRIME
    return basis


def modular_null(rows):
    rows = np.array(rows % PRIME, copy=True)
    width = rows.shape[1]
    pivots = []
    for i in range(width):
        r = len(pivots)
        hits = np.flatnonzero(rows[r:, i])
        if not len(hits):
            continue
        j = r + hits[0]
        rows[[r, j]] = rows[[j, r]]
        rows[r] = rows[r] * pow(int(rows[r, i]), -1, PRIME) % PRIME
        coeff = rows[:, i].copy()
        coeff[r] = 0
        rows = (rows - coeff[:, None] * rows[r]) % PRIME
        pivots.append(i)
    free = sorted(set(range(width)) - set(pivots))
    out = np.zeros((width, len(free)), dtype=np.int64)
    out[free] = np.eye(len(free), dtype=np.int64)
    out[pivots] = -rows[: len(pivots), free] % PRIME
    return out


def invariant(exits, kernel, h):
    width, dim = kernel.shape
    if not dim:
        return [], 0
    row, col = (
        zip(*[(i, j) for i, values in enumerate(exits) for j in values])
        if exits
        else ([], [])
    )
    cap = sparse.csr_matrix(
        (np.ones(len(row), dtype=np.int64), (row, col)), shape=(len(exits), width)
    )
    flight = kernel
    layers = [dim]
    while dim:
        flight = (h @ flight) % PRIME
        sol = modular_null(cap @ flight)
        reduced = sol.shape[1]
        if reduced == dim:
            # The modular full-rank exclusions hold over the rationals. For a
            # surviving subspace, reconstruct it and check invariance exactly.
            block = sp.Matrix(cap.toarray())
            H = sp.SparseMatrix(
                width,
                width,
                {(int(i), int(j)): int(h[i, j]) for i, j in zip(*h.nonzero())},
            )
            constraints = block
            for _ in range(len(layers)):
                block = block * H
                constraints = constraints.col_join(block)
            while True:
                exact = constraints.nullspace()
                if not exact:
                    return layers, 0
                basis = sp.Matrix.hstack(*exact)
                if basis.row_join(H * basis).rank() == basis.cols:
                    return layers, basis.cols
                block = block * H
                constraints = constraints.col_join(block)
        flight = flight @ sol % PRIME
        dim = reduced
        layers.append(dim)
    return layers, 0


def run():
    result = []
    for gates, boards in owner_shapes().items():
        entries = []
        for a in range(10):
            for b in range(a, 10 - a):
                if not a + b:
                    continue
                words, exits, exchange, diff = sector(gates, a, b)
                kernel = boundary_kernel(exits, len(words))
                # Swapping both nonhost names and reversing all physical loads
                # only permutes histograms. Test every distinct contrast cost.
                choices = {}
                for values in ((0, 1, 2), (0, 2, 1), (1, 0, 2)):
                    squared = tuple(
                        (values[i] - values[j]) ** 2
                        for i, j in ((0, 1), (0, 2), (1, 2))
                    )
                    profile = tuple(np.array(squared) @ diff)
                    choices.setdefault(profile, []).append(values)
                readings = []
                for diag, labels in choices.items():
                    h = exchange + sparse.diags(
                        2 * int(candidate.KAPPA) * np.array(diag),
                        dtype=np.int64,
                        format="csr",
                    )
                    layers, dimension = invariant(exits, kernel, h)
                    readings.append(
                        dict(
                            symbol_values=labels,
                            successive_dimensions=layers,
                            compact_dimension=dimension,
                        )
                    )
                layers, dimension = invariant(exits, kernel, exchange)
                entries.append(
                    dict(
                        counts=[a, b],
                        words=len(words),
                        boundary_kernel_upper=kernel.shape[1],
                        native=readings,
                        exchange_only=dict(
                            successive_dimensions=layers, compact_dimension=dimension
                        ),
                    )
                )
        result.append(dict(gateways=gates, boards=boards, histograms=entries))
    return dict(prime=PRIME, gate_classes=result)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run()
    with args.output.open("x") as f:
        f.write(json.dumps(result, indent=2) + "\n")
