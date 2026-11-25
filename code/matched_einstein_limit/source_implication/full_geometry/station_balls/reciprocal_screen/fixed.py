"""Exact finite-field lower bound on response constraints for a source-independent analyzer."""

from math import comb
import json

import numpy as np

from ....beam_response import measure as camera

PRIME = 65521
ORDER = 10


def powers(h, rings, limit=ORDER):
    stations = camera.indices(rings[0])
    dest = camera.indices(rings[1])
    start = np.zeros((h.shape[0], len(stations)), dtype="int64")
    start[stations] = np.eye(len(stations), dtype="int64")
    seq = []
    for _ in range(limit + 1):
        seq.append(
            start[dest].reshape(5, 4, 5, 4).transpose(0, 2, 1, 3).reshape(25, 4, 4)
        )
        start = (h @ start) % PRIME
    return seq


def constraint_rows(seq, n):
    lhs = np.zeros((25, 4, 4, 4, 4), dtype="int64")
    rhs = np.zeros_like(lhs)
    for i in range(n + 1):
        x, y = seq[i], seq[n - i]
        coefficient = (-1) ** i * comb(n, i) % PRIME
        lhs += coefficient * np.einsum("vai,vbj->vijab", x, y)
        rhs += coefficient * np.einsum("via,vjb->vijab", x, y)
    common = (lhs - rhs).reshape(-1, 16) % PRIME
    ports = np.zeros((25, 16, 10, 16), dtype="int64")
    for j in range(5):
        for i in range(5):
            pair = 5 * j + i
            ports[pair, :, 5 + j, :] = lhs[pair].reshape(16, 16)
            ports[pair, :, i, :] = -rhs[pair].reshape(16, 16)
    return common, ports.reshape(-1, 160) % PRIME


def row_basis(rows):
    """An echelon basis with row entries always reduced before multiplication."""
    a = rows.copy()
    rank = 0
    for col in range(a.shape[1]):
        eligible = np.flatnonzero(a[rank:, col])
        if not len(eligible):
            continue
        pivot = rank + int(eligible[0])
        a[[rank, pivot]] = a[[pivot, rank]]
        a[rank] = a[rank] * pow(int(a[rank, col]), -1, PRIME) % PRIME
        a[rank + 1 :] = (a[rank + 1 :] - a[rank + 1 :, col, None] * a[rank]) % PRIME
        rank += 1
    return a[:rank]


def certificate():
    board = camera.front.optics.pairing.formation.federation(1)
    route = camera.front.optics.pairing.formation.region(board)
    rings = camera.screens(board, route)
    records = {}
    for odd in (False, True):
        both = {}
        for q in (-1, 1):
            h = camera.front.native.generator(board, route, q, odd)
            seq = powers(h, rings)
            both[q] = seq
            common = np.zeros((0, 16), dtype="int64")
            separate = np.zeros((0, 160), dtype="int64")
            changes = []
            for n in range(6, ORDER + 1):
                rows, ports = constraint_rows(seq, n)
                common = row_basis(np.vstack([common, rows]))
                separate = row_basis(np.vstack([separate, ports]))
                changes.append([n, len(common), len(separate)])
                if len(separate) == 160:
                    break
            records[str((odd, q))] = changes
        shared = np.zeros((0, 160), dtype="int64")
        tagged = np.zeros((0, 320), dtype="int64")
        pooled = []
        for n in range(6, ORDER + 1):
            pos = constraint_rows(both[1], n)[1]
            neg = constraint_rows(both[-1], n)[1]
            shared = row_basis(np.vstack([shared, (pos + neg) % PRIME]))
            tagged = row_basis(np.vstack([tagged, np.column_stack([pos, neg])]))
            pooled.append([n, len(shared), len(tagged)])
            if len(tagged) == 320:
                break
        records[f"pooled_{odd}"] = pooled
    return dict(prime=PRIME, derivative_limit=ORDER, ranks=records)


if __name__ == "__main__":
    print(json.dumps(certificate(), indent=2))
