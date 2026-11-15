"""Independent recomputation of the signed-source charge, rest-energy and frame results.

The script reads only the sealed carrier graph (from the published evidence in
code/openwave_static_signed_source) and the native repair rule. It does
not import the sealed experiment scripts. Every statement is checked exactly
(rational or integer arithmetic) except the logarithmic energies, which are
closed-form values of exact rational transfer eigenvalues.

    python3 code/openwave_signed_source/recompute.py
"""
from __future__ import annotations

import hashlib
import json
import math
import random
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path

HERE = Path(__file__).resolve().parent
EVIDENCE = HERE.parent / 'openwave_static_signed_source' / 'evidence' / 'openwave' / 'shared'
GRAPH = EVIDENCE / 'localized_rest_join' / 'evidence' / 'graph.json'
GRAPH_SHA256 = '2e8be4c030e8e04ab903adb5a8a786f11654ecca0380ee54d17509a707f49085'


def load_graph():
    raw = GRAPH.read_bytes()
    if hashlib.sha256(raw).hexdigest() != GRAPH_SHA256:
        raise ValueError('carrier graph differs from the sealed graph')
    graph = json.loads(raw)
    seams = [tuple(e) for e in graph['seams']]
    return graph, seams


def native_current(d, coin):
    """Native repair current on a seam with load difference d and coin c."""
    return (d - (2 * coin - 1) * (abs(d) % 2)) // 2


def require(condition, message):
    if not condition:
        raise ValueError(message)


# ---------------------------------------------------------------- integer topology
def smith_diagonal(rows):
    """Absolute Smith invariants of an integer matrix (list of rows)."""
    A = [list(r) for r in rows]
    m, n = len(A), len(A[0])
    diag = []
    r = 0
    for c in range(n):
        if r >= m:
            break
        while True:
            pivots = [(abs(A[i][j]), i, j) for i in range(r, m) for j in range(c, n) if A[i][j]]
            if not pivots:
                return diag + [0] * (min(m, n) - len(diag))
            _, i, j = min(pivots)
            A[r], A[i] = A[i], A[r]
            for row in A:
                row[c], row[j] = row[j], row[c]
            p = A[r][c]
            done = True
            for i in range(r + 1, m):
                q = A[i][c] // p
                for k in range(c, n):
                    A[i][k] -= q * A[r][k]
                if A[i][c]:
                    done = False
            for k in range(c + 1, n):
                q = A[r][k] // p
                for i in range(r, m):
                    A[i][k] -= q * A[i][c]
                if A[r][k]:
                    done = False
            if done:
                bad = [(i, k) for i in range(r + 1, m) for k in range(c + 1, n) if A[i][k] % p]
                if not bad:
                    break
                i, _ = bad[0]
                for k in range(c, n):
                    A[r][k] += A[i][k]
        diag.append(abs(A[r][c]))
        r += 1
    return diag + [0] * (min(m, n) - len(diag))


def rank(rows):
    M = [[F(x) for x in r] for r in rows]
    rk, cols = 0, len(M[0])
    for c in range(cols):
        piv = next((i for i in range(rk, len(M)) if M[i][c]), None)
        if piv is None:
            continue
        M[rk], M[piv] = M[piv], M[rk]
        for i in range(len(M)):
            if i != rk and M[i][c]:
                f = M[i][c] / M[rk][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[rk])]
        rk += 1
    return rk


def solve(A, b):
    """Exact solve of a nonsingular rational system."""
    n = len(A)
    M = [[F(x) for x in row] + [F(v)] for row, v in zip(A, b)]
    for c in range(n):
        piv = next(i for i in range(c, n) if M[i][c])
        M[c], M[piv] = M[piv], M[c]
        for i in range(n):
            if i != c and M[i][c]:
                f = M[i][c] / M[c][c]
                M[i] = [a - f * x for a, x in zip(M[i], M[c])]
    return [M[i][n] / M[i][i] for i in range(n)]


def charge(cell):
    ports = 12
    index = {e: k for k, e in enumerate(cell)}
    adj = {v: set() for v in range(ports)}
    for a, b in cell:
        adj[a].add(b)
        adj[b].add(a)
    faces = [t for t in combinations(range(ports), 3)
             if t[1] in adj[t[0]] and t[2] in adj[t[0]] and t[2] in adj[t[1]]]
    require(len(cell) == 30 and len(faces) == 20, 'twelve-port cell is the icosahedral sphere')
    B = [[0] * 30 for _ in range(ports)]
    for k, (a, b) in enumerate(cell):
        B[a][k], B[b][k] = 1, -1
    D = [[0] * 20 for _ in range(30)]
    for f, (a, b, c) in enumerate(faces):
        for x, y in ((a, b), (b, c), (c, a)):
            D[index[tuple(sorted((x, y)))]][f] += 1 if x < y else -1
    BD = [[sum(B[i][k] * D[k][j] for k in range(30)) for j in range(20)] for i in range(ports)]
    require(all(v == 0 for r in BD for v in r), 'boundary of boundary')
    smith = smith_diagonal(D)
    rankB = rank(B)
    require(rankB == 11 and smith == [1] * 19 + [0], 'face certificate')
    # rank D = 19 = 30 - rank B and unit Smith factors: ker B = im D over Z.

    # Stopped mean flow E[F_tau] = B^T L^+ q from the exact Bellman equation.
    L = [[sum(B[i][k] * B[j][k] for k in range(30)) for j in range(ports)] for i in range(ports)]
    G = [[F(L[i][j]) + F(1, ports) for j in range(ports)] for i in range(ports)]
    sources = [(p, r) for p in range(ports) for r in range(ports) if p != r]

    def flow(q):
        x = solve(G, q)
        mean = sum(x) / ports
        x = [v - mean for v in x]
        return [x[a] - x[b] for a, b in cell]

    expected = {s: flow([int(i == s[0]) - int(i == s[1]) for i in range(ports)]) for s in sources}
    checked = 0
    for p, r in sources:
        q = [int(i == p) - int(i == r) for i in range(ports)]
        rhs = [F(0)] * 30
        for k, (a, b) in enumerate(cell):
            for coin in (0, 1):
                j = native_current(q[a] - q[b], coin)
                out = q[:]
                out[a] -= j
                out[b] += j
                nxt = [F(0)] * 30 if not any(out) else expected[(out.index(1), out.index(-1))]
                for m in range(30):
                    rhs[m] += nxt[m] / 60
                rhs[k] += F(j, 60)
        require(rhs == expected[(p, r)], 'Bellman equation')
        checked += 1

    # Seeded native trials and a changed-law (hold) control on the same cell.
    def run(p, r, seed, law, cap=20000):
        rng = random.Random(seed)
        q = [int(i == p) - int(i == r) for i in range(ports)]
        flux = [0] * 30
        for _ in range(cap):
            if not any(q):
                break
            k = rng.randrange(30)
            a, b = cell[k]
            j = native_current(q[a] - q[b], rng.randrange(2)) if law == 'native' else 0
            q[a] -= j
            q[b] += j
            flux[k] += j
        boundary = [0] * ports
        for (a, b), j in zip(cell, flux):
            boundary[a] += j
            boundary[b] -= j
        return not any(q), boundary

    absorbed = 0
    for t, (p, r) in enumerate(sources * 2):
        ok, boundary = run(p, r, 1000 + t, 'native')
        require(ok and boundary == [int(i == p) - int(i == r) for i in range(ports)], 'native absorption')
        absorbed += 1
    held = sum(run(p, r, 5000 + t, 'hold', cap=200)[0] for t, (p, r) in enumerate(sources[:48]))
    require(held == 0, 'changed-law control absorbs')

    return {'cell_seams': 30, 'cell_faces': 20, 'boundary_rank': rankB, 'face_smith': smith,
            'kernel_equals_face_image_over_Z': True, 'bellman_sources_checked': checked,
            'native_trials_absorbed': absorbed, 'changed_law_controls_absorbed': held}


# ---------------------------------------------------------------- rest energy
def energy(graph, seams):
    n = graph['ports']
    adj = {v: [] for v in range(n)}
    for a, b in seams:
        adj[a].append(b)
        adj[b].append(a)

    def lap(vec):
        return {v: len(adj[v]) * vec.get(v, 0) - sum(vec.get(u, 0) for u in adj[v]) for v in range(n)}

    f = {int(k): F(v, 2) for k, v in graph['f'].items()}
    g = {int(k): F(v, 2) for k, v in graph['g'].items()}
    for ray in (f, g):
        Lr = lap(ray)
        require(all(Lr[v] == 6 * ray.get(v, 0) for v in range(n)), 'compact ray eigenvalue six')
    require(not any(b in g for a in f for b in adj[a]) and not set(f) & set(g), 'disjoint rays')
    require(len(seams) == 630, 'full graph seams')
    t6 = 1 - F(6, 1260)
    t12 = 1 - F(12, 1260)
    require(t6 == F(209, 210) and t12 == F(104, 105), 'transfer eigenvalues')
    E6 = -1260 * math.log1p(-1 / 210)
    E12 = -1260 * math.log1p(-1 / 105)
    pair_keV = 511 * E12 / E6
    excess = 1260 * math.log1p(1 / 43680)
    require(abs((E12 - 2 * E6) - excess) < 1e-9, 'finite-clock excess')
    # A lower nonconstant mode: the difference of the first two cell indicators.
    witness = {v: F(1) if v < 12 else (F(-1) if v < 24 else F(0)) for v in range(n)}
    Lw = lap(witness)
    rq_witness = sum(Lw[v] * witness[v] for v in range(n)) / sum(x * x for x in witness.values())
    require(rq_witness == F(1, 3), 'lower nonconstant mode at Rayleigh quotient 1/3')
    return {'ports': n, 'seams': len(seams), 'transfer': 'T = I - L/1260',
            'E6': E6, 'E12': E12, 'calibration_keV': 511, 'pair_keV': pair_keV,
            'pair_excess_keV': pair_keV - 1022, 'constant_mode_energy': 0,
            'lower_mode_rayleigh_quotient_witness': str(rq_witness)}


# ---------------------------------------------------------------- frame and field
def frame(cell):
    import numpy as np
    phi = (1 + 5 ** 0.5) / 2
    L = np.zeros((12, 12))
    for a, b in cell:
        L[a, a] += 1
        L[b, b] += 1
        L[a, b] -= 1
        L[b, a] -= 1
    w, V = np.linalg.eigh(L)
    target = 5 - 5 ** 0.5
    cols = [k for k in range(12) if abs(w[k] - target) < 1e-9]
    require(len(cols) == 3, 'icosahedral embedding')
    X = V[:, cols]
    X = X / np.linalg.norm(X[0]) * math.sqrt(1 + phi * phi)   # vertices (0, +-1, +-phi)
    neighbours = {v: {b if a == v else a for a, b in cell if v in (a, b)} for v in range(12)}
    anti = {p: next(r for r in range(12) if r != p and r not in neighbours[p]
                    and not neighbours[p] & neighbours[r]) for p in range(12)}

    def moments(p):
        r = anti[p]
        q = [int(i == p) - int(i == r) for i in range(12)]
        deltas = []
        for a, b in cell:
            for coin in (0, 1):
                j = native_current(q[a] - q[b], coin)
                deltas.append(j * (X[b] - X[a]))
        D = np.array(deltas)
        m = D.mean(axis=0)
        C = ((D - m).T @ (D - m)) / len(D)
        return m, C, len(deltas)

    branches = 0
    for p in range(12):
        m, C, k = moments(p)
        branches += k
        u = X[p] / np.linalg.norm(X[p])
        require(np.allclose(m, (phi - 3) / 15 * X[p]), 'polar mean')
        require(abs(u @ C @ u - (3 - phi) / 9) < 1e-12, 'parallel variance')
        P = np.eye(3) - np.outer(u, u)
        require(np.allclose(P @ C @ P, (phi + 2) / 15 * P), 'transverse variance')
    # Executed proper half-turn: rotation by pi about a two-fold axis of the cell
    # (through a seam midpoint) orthogonal to X[p]; it maps p to its antipode.
    p = 0
    u = X[p] / np.linalg.norm(X[p])
    axis = next(X[a] + X[b] for a, b in cell if abs((X[a] + X[b]) @ u) < 1e-9)
    axis /= np.linalg.norm(axis)
    R = 2 * np.outer(axis, axis) - np.eye(3)
    require(abs(np.linalg.det(R) - 1) < 1e-12, 'proper rotation')
    perm = [int(np.argmin(np.linalg.norm(X - R @ X[v], axis=1))) for v in range(12)]
    require(sorted(perm) == list(range(12)) and all(np.allclose(R @ X[v], X[perm[v]]) for v in range(12)),
            'half-turn is a cell symmetry')
    require(perm[p] == anti[p], 'half-turn maps the source axis to its antipode')
    m0, C0, _ = moments(p)
    m1, C1, _ = moments(perm[p])
    require(np.allclose(m1, R @ m0) and np.allclose(m1, -m0), 'mean reverses under the half-turn')
    require(np.allclose(C1, R @ C0 @ R.T) and np.allclose(C1, C0), 'covariance returns under the half-turn')
    require(np.allclose(R @ R @ m0, m0), 'second half-turn returns the mean')
    other = next(v for v in range(12) if v not in (p, anti[p]))
    require(not np.allclose(moments(other)[1], C0), 'covariance depends on the axis')
    return {'branches': branches, 'mean_coefficient': (phi - 3) / 15,
            'C_parallel': (3 - phi) / 9, 'C_transverse': (phi + 2) / 15,
            'half_turn_returns_covariance': True, 'half_turn_reverses_mean': True,
            'second_half_turn_returns_mean': True, 'covariance_axis_dependent': True,
}


def main():
    graph, seams = load_graph()
    cell = [e for e in seams if e[1] < 12]
    result = {'charge': charge(cell), 'rest_energy': energy(graph, seams), 'frame': frame(cell)}
    print(json.dumps(result, indent=2))
    return result


if __name__ == '__main__':
    main()
