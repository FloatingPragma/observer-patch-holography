"""Exact finite-sector reduction and variational floor for connected neutral words."""

from itertools import combinations

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh

from model.candidate_laws.interval_echo import candidate
from model.executable_model import cell, federation


def construct():
    # Independent basis: first choose all four sites, then choose their two positive signs.
    basis = []
    connected = []
    links = cell().seams
    for support in combinations(range(12), 4):
        seen = {support[0]}
        while True:
            neighbors = {p for p, q in links if p in support and q in seen} | {
                q for p, q in links if q in support and p in seen
            }
            if neighbors <= seen:
                break
            seen.update(neighbors)
        for up in combinations(support, 2):
            basis.append(
                tuple(1 if p in up else -1 if p in support else 0 for p in range(12))
            )
            connected.append(len(seen) == 4)
    index = {w: i for i, w in enumerate(basis)}
    r = []
    c = []
    vals = []
    references = np.flatnonzero(federation(0).degree[:12] == 5)
    k = candidate.KAPPA
    for i, w in enumerate(basis):
        diagonal = k * sum((1 + w[p]) ** 2 - 1 for p in references)
        for p, q in links:
            if w[p] != w[q]:
                changed = list(w)
                changed[p], changed[q] = changed[q], changed[p]
                r.append(i)
                c.append(index[tuple(changed)])
                vals.append(-0.5)
                diagonal += 0.5 + k * (w[p] - w[q]) ** 2
        r.append(i)
        c.append(i)
        vals.append(diagonal)
    matrix = coo_matrix((vals, (r, c)), shape=(len(basis), len(basis))).tocsr()
    return matrix, np.array(connected), np.array(basis)


def measure():
    matrix, connected, words = construct()
    compact = matrix[connected][:, connected]
    # Work directly in doubled energy units. Off-diagonal entries are nonpositive.
    if not np.array_equal(2 * matrix.data, (2 * matrix.data).astype(int)):
        raise ArithmeticError("Generator entries must be integral in doubled units")
    compact_integer = (2 * compact).astype("int64")
    off = compact_integer.copy()
    off.setdiag(0)
    if max(off.data) > 0:
        raise ArithmeticError("Positive off-diagonal invalidates the row-sum bound")
    lower = int(compact_integer.sum(axis=1).min()) / 2
    witness = int(np.argmin(matrix.diagonal()))
    upper = matrix.diagonal()[witness]

    def ground(op):
        e, v = eigsh(op, k=1, which="SA", tol=1e-13, v0=np.ones(op.shape[0]))
        residual = np.linalg.norm(op @ v[:, 0] - e[0] * v[:, 0])
        return float(e[0]), float(residual)

    return {
        "all": ground(matrix),
        "connected_only": ground(matrix[connected][:, connected]),
        "other_only": ground(matrix[~connected][:, ~connected]),
        "connected_words": int(sum(connected)),
        "connected_exact_lower": lower,
        "full_variational_upper": float(upper),
        "witness_word": words[witness].tolist(),
    }


if __name__ == "__main__":
    print(measure())
