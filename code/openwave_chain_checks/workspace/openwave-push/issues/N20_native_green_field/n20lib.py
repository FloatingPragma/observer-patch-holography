#!/usr/bin/env python3
"""n20lib: graphs, exact rational algebra, the absorbing chain of the literal rule, samplers and
Bloch sums for worker N20 (the native Green field).

Conventions (those of the task and of I05):
  seam e = (a, b) is oriented a -> b;  B[a, e] = +1, B[b, e] = -1;  L = B B^T.
  event (e, coin):  d = n_a - n_b,  j = (d - (2 coin - 1)(|d| mod 2)) // 2,  n -> n - j B_e,  F_e -> F_e + j.
  Hence n + B F is constant.  The record field of a history is F; (B F)_x is the flux out of port x.

Nothing here writes outside this directory.  The reference kernel `kernel/ophk.py` is imported, not copied.
"""
from __future__ import annotations

import itertools
import math
import sys
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = Path(__file__).resolve().parent
PUSH = HERE.parents[1]
if str(PUSH / "kernel") not in sys.path:
    sys.path.insert(0, str(PUSH / "kernel"))
import ophk  # noqa: E402  (read-only reference kernel)


# ----------------------------------------------------------------------------------------------
# graphs
# ----------------------------------------------------------------------------------------------
class Graph:
    def __init__(self, name, n, seams, meta=None):
        self.name = name
        self.N = int(n)
        self.seams = [tuple(int(v) for v in s) for s in seams]
        self.M = len(self.seams)
        self.meta = meta or {}

    def incidence_fr(self):
        b = [[Fr(0)] * self.M for _ in range(self.N)]
        for e, (a, c) in enumerate(self.seams):
            b[a][e] = Fr(1)
            b[c][e] = Fr(-1)
        return b

    def laplacian_int(self):
        lap = [[0] * self.N for _ in range(self.N)]
        for a, b in self.seams:
            lap[a][a] += 1
            lap[b][b] += 1
            lap[a][b] -= 1
            lap[b][a] -= 1
        return lap

    def incidence_sparse(self):
        s = np.asarray(self.seams)
        m = self.M
        rows = np.concatenate([s[:, 0], s[:, 1]])
        cols = np.concatenate([np.arange(m), np.arange(m)])
        vals = np.concatenate([np.ones(m), -np.ones(m)])
        return sp.csr_matrix((vals, (rows, cols)), shape=(self.N, m))

    def degrees(self):
        deg = [0] * self.N
        for a, b in self.seams:
            deg[a] += 1
            deg[b] += 1
        return deg

    def connected(self):
        adj = [[] for _ in range(self.N)]
        for a, b in self.seams:
            adj[a].append(b)
            adj[b].append(a)
        seen = {0}
        stack = [0]
        while stack:
            v = stack.pop()
            for w in adj[v]:
                if w not in seen:
                    seen.add(w)
                    stack.append(w)
        return len(seen) == self.N


def carrier_graph():
    c = ophk.carrier()
    return Graph("carrier", 12, c.seams)


def k4_graph():
    return Graph("K4", 4, list(itertools.combinations(range(4), 2)))


def path_graph(n=4):
    return Graph(f"path{n}", n, [(i, i + 1) for i in range(n - 1)])


def paw_graph():
    return Graph("paw", 4, [(0, 1), (1, 2), (0, 2), (2, 3)])


def house_graph():
    return Graph("house", 5, [(0, 1), (1, 2), (2, 3), (3, 0), (0, 4), (1, 4)])


def chang_graph():
    """Seidel switching of the triangular graph T(8) with respect to the eight edges of an 8-cycle.
    Strongly regular with parameters (28, 12, 6, 4)."""
    verts = list(itertools.combinations(range(8), 2))
    idx = {v: i for i, v in enumerate(verts)}
    switch = {idx[tuple(sorted((i, (i + 1) % 8)))] for i in range(8)}
    seams = []
    for i, u in enumerate(verts):
        for j in range(i + 1, 28):
            w = verts[j]
            adj = len(set(u) & set(w)) == 1
            if (i in switch) != (j in switch):
                adj = not adj
            if adj:
                seams.append((i, j))
    return Graph("chang", 28, seams)


def triangular_graph8():
    verts = list(itertools.combinations(range(8), 2))
    seams = [(i, j) for i in range(28) for j in range(i + 1, 28) if len(set(verts[i]) & set(verts[j])) == 1]
    return Graph("T8", 28, seams)


# ---- the cubic routing of RER model/native_spatial_refinement, rebuilt from its definition -------------
# Numbers a + b*phi are pairs (a, b); phi^2 = phi + 1.
def _qmul(x, y):
    return (x[0] * y[0] + x[1] * y[1], x[0] * y[1] + x[1] * y[0] + x[1] * y[1])


def _qsub(x, y):
    return (x[0] - y[0], x[1] - y[1])


def _qsign(x):
    return (x[0] + x[1] * (1 + 5 ** 0.5) / 2 > 0) - (x[0] + x[1] * (1 + 5 ** 0.5) / 2 < 0)


_Z, _ONE, _MIN, _PHI, _NPHI = (0, 0), (1, 0), (-1, 0), (0, 1), (0, -1)
# port order of model.common_source.slow_geometry.raw_frame
RAW_FRAME = ((_MIN, _PHI, _Z), (_ONE, _PHI, _Z), (_MIN, _NPHI, _Z), (_ONE, _NPHI, _Z),
             (_Z, _MIN, _PHI), (_Z, _ONE, _PHI), (_Z, _MIN, _NPHI), (_Z, _ONE, _NPHI),
             (_PHI, _Z, _MIN), (_PHI, _Z, _ONE), (_NPHI, _Z, _MIN), (_NPHI, _Z, _ONE))


def frame_seams():
    """Thirty internal seams: pairs of frame vectors at squared distance four (exact in Q(sqrt 5))."""
    out = []
    for i in range(12):
        for j in range(i + 1, 12):
            d = (0, 0)
            for a, b in zip(RAW_FRAME[i], RAW_FRAME[j]):
                t = _qsub(a, b)
                sq = _qmul(t, t)
                d = (d[0] + sq[0], d[1] + sq[1])
            if d == (4, 0):
                out.append((i, j))
    return tuple(out)


def cubic_links(control=False):
    """Ordered port pairs (p, m): port p of a cell is joined to port m of the next cell along the axis.

    Declared routing (PROOF.md): in each coordinate-zero group the two ports with positive product of the other
    two coordinates are the positive ports; each is joined to the port obtained by flipping the cyclic last
    nonzero coordinate.  control=True is the alternative axis-zero partition of the repository (sign of coordinate
    one, flip of coordinate one, on axis 0 only)."""
    out = []
    for axis in range(3):
        row = []
        for i, p in enumerate(RAW_FRAME):
            if p[axis] != _Z:
                continue
            positive = _qsign(_qmul(p[(axis + 1) % 3], p[(axis + 2) % 3])) > 0
            if control and axis == 0:
                positive = _qsign(p[1]) > 0
            if positive:
                q = list(p)
                flip = (axis + 1) % 3 if (control and axis == 0) else (axis + 2) % 3
                q[flip] = (-q[flip][0], -q[flip][1])
                row.append((i, RAW_FRAME.index(tuple(q))))
        assert len(row) == 2
        out.append(tuple(row))
    assert sorted(x for row in out for e in row for x in e) == list(range(12))
    return tuple(out)


def cubic_graph(size, control=False):
    """Periodic cubic routing: 12 L^3 ports, 36 L^3 seams.  Port id = 12 * cell + p, cell = i L^2 + j L + k.
    Seam order: the 30 internal seams of every cell first, then the inter-cell seams (cell, axis, link);
    an inter-cell seam is oriented along the positive axis."""
    assert size >= 2
    coords = list(itertools.product(range(size), repeat=3))
    index = {c: i for i, c in enumerate(coords)}
    internal = frame_seams()
    links = cubic_links(control)
    seams = []
    for cell in range(len(coords)):
        seams += [(12 * cell + a, 12 * cell + b) for a, b in internal]
    inter_start = len(seams)
    face = []          # (cell, axis) for every inter-cell seam
    for cell, c in enumerate(coords):
        for axis in range(3):
            nxt = list(c)
            nxt[axis] = (nxt[axis] + 1) % size
            other = index[tuple(nxt)]
            for a, b in links[axis]:
                seams.append((12 * cell + a, 12 * other + b))
                face.append((cell, axis))
    g = Graph(f"cubic{size}{'_control' if control else ''}", 12 * len(coords), seams,
              meta={"size": size, "coords": coords, "index": index, "inter_start": inter_start,
                    "face": face, "links": links, "internal": internal})
    assert len(set(tuple(sorted(s)) for s in seams)) == len(seams)
    return g


def cubic_symbol_zero_int(control=False):
    h = [[0] * 12 for _ in range(12)]
    for i in range(12):
        h[i][i] = 6
    for a, b in frame_seams():
        h[a][b] -= 1
        h[b][a] -= 1
    for row in cubic_links(control):
        for a, b in row:
            h[a][b] -= 1
            h[b][a] -= 1
    return h


def cubic_quadratic_coefficient(control=False):
    """Exact second-order coefficient kappa_ab of lambda_0(k) = sum kappa_ab k_a k_b + O(k^4):
    kappa = (links per axis / 12) delta_ab - j_a^T G j_b / 12 with G the pseudo-inverse of H(0)."""
    h0 = [[Fr(v) for v in row] for row in cubic_symbol_zero_int(control)]
    g = pinv_laplacian_fr(h0)
    js = []
    for row in cubic_links(control):
        j = [Fr(0)] * 12
        for a, b in row:
            j[a] -= 1
            j[b] += 1
        js.append(j)
    gram = [[sum(js[a][p] * g[p][q] * js[b][q] for p in range(12) for q in range(12)) for b in range(3)]
            for a in range(3)]
    kappa = [[(Fr(2, 12) if a == b else Fr(0)) - gram[a][b] / 12 for b in range(3)] for a in range(3)]
    eig8 = all(sum(h0[p][q] * js[a][q] for q in range(12)) == 8 * js[a][p] for a in range(3) for p in range(12))
    return kappa, gram, eig8


def ring_graph(lx, ly):
    f = ophk.flat(lx, ly, "ring")
    g = Graph(f"ring{lx}x{ly}", f.ports, [tuple(map(int, s)) for s in f.seams],
              meta={"lx": lx, "ly": ly, "intra": f.intra, "points": f.cell_points, "carriers": f.carriers})
    return g


# ----------------------------------------------------------------------------------------------
# exact rational linear algebra
# ----------------------------------------------------------------------------------------------
def lu_fr(a):
    """P A = L U with Fractions.  Returns (lu, perm)."""
    n = len(a)
    lu = [row[:] for row in a]
    perm = list(range(n))
    for k in range(n):
        piv = next((i for i in range(k, n) if lu[i][k] != 0), None)
        if piv is None:
            raise ZeroDivisionError("singular matrix")
        if piv != k:
            lu[k], lu[piv] = lu[piv], lu[k]
            perm[k], perm[piv] = perm[piv], perm[k]
        rk = lu[k]
        inv = 1 / rk[k]
        nzcols = [c for c in range(k + 1, n) if rk[c] != 0]
        for i in range(k + 1, n):
            ri = lu[i]
            if ri[k] != 0:
                f = ri[k] * inv
                ri[k] = f
                for c in nzcols:
                    ri[c] -= f * rk[c]
    return lu, perm


def lu_solve_fr(lu, perm, b):
    """Solve A x = b (b a vector of Fractions)."""
    n = len(lu)
    y = [b[perm[i]] for i in range(n)]
    for i in range(n):
        ri = lu[i]
        s = y[i]
        for c in range(i):
            if ri[c] != 0 and y[c] != 0:
                s -= ri[c] * y[c]
        y[i] = s
    for i in range(n - 1, -1, -1):
        ri = lu[i]
        s = y[i]
        for c in range(i + 1, n):
            if ri[c] != 0 and y[c] != 0:
                s -= ri[c] * y[c]
        y[i] = s / ri[i]
    return y


def lu_solve_transpose_fr(lu, perm, b):
    """Solve A^T x = b."""
    n = len(lu)
    y = list(b)
    for i in range(n):                     # U^T y = b
        s = y[i]
        for c in range(i):
            if lu[c][i] != 0 and y[c] != 0:
                s -= lu[c][i] * y[c]
        y[i] = s / lu[i][i]
    for i in range(n - 1, -1, -1):         # L^T z = y
        s = y[i]
        for c in range(i + 1, n):
            if lu[c][i] != 0 and y[c] != 0:
                s -= lu[c][i] * y[c]
        y[i] = s
    x = [Fr(0)] * n
    for i in range(n):
        x[perm[i]] = y[i]
    return x


def inverse_fr(a):
    n = len(a)
    lu, perm = lu_fr(a)
    cols = []
    for j in range(n):
        e = [Fr(0)] * n
        e[j] = Fr(1)
        cols.append(lu_solve_fr(lu, perm, e))
    return [[cols[j][i] for j in range(n)] for i in range(n)]


def pinv_laplacian_fr(lap):
    """Pseudo-inverse of a connected-graph Laplacian: (L + J/N)^(-1) - J/N, exact."""
    n = len(lap)
    a = [[Fr(lap[i][j]) + Fr(1, n) for j in range(n)] for i in range(n)]
    inv = inverse_fr(a)
    return [[inv[i][j] - Fr(1, n) for j in range(n)] for i in range(n)]


def weighted_laplacian_fr(g, w):
    lap = [[Fr(0)] * g.N for _ in range(g.N)]
    for (a, b), we in zip(g.seams, w):
        lap[a][a] += we
        lap[b][b] += we
        lap[a][b] -= we
        lap[b][a] -= we
    return lap


# ----------------------------------------------------------------------------------------------
# the absorbing Markov chain of the literal rule
# ----------------------------------------------------------------------------------------------
def build_chain(g, starts, weights=None):
    """Breadth-first closure of the literal rule `ophk.native` from the given load configurations.

    Returns (states, trans): states is a list of load tuples; trans[i] is a list of
    (probability, next index, seam, coin, current) over ALL (seam, coin) events of state i.
    weights: seam draw probabilities (Fractions summing to one); default uniform.  Coins are fair."""
    w = weights if weights is not None else [Fr(1, g.M)] * g.M
    index, states, queue = {}, [], []
    for s in starts:
        s = tuple(s)
        if s not in index:
            index[s] = len(states)
            states.append(s)
            queue.append(s)
    trans = []
    head = 0
    while head < len(queue):
        s = queue[head]
        head += 1
        row = []
        for e, seam in enumerate(g.seams):
            for coin in (0, 1):
                out, (_, cur) = ophk.native(s, seam, coin)
                if out not in index:
                    index[out] = len(states)
                    states.append(out)
                    queue.append(out)
                row.append((w[e] / 2, index[out], e, coin, cur))
        trans.append(row)
    return states, trans


def pair_start(g, x, y, q):
    s = [q] * g.N
    s[x] += 1
    s[y] -= 1
    return tuple(s)


class PairChain:
    """Exact statistics of an excess and a deficit over the filling q up to their annihilation."""

    def __init__(self, g, q=0, weights=None, record=None):
        self.g, self.q = g, q
        starts = [pair_start(g, x, y, q) for x in range(g.N) for y in range(g.N) if x != y]
        self.states, self.trans = build_chain(g, starts, weights)
        uniform = tuple([q] * g.N)
        self.absorbing = self.states.index(uniform)
        self.tr = [i for i in range(len(self.states)) if i != self.absorbing]
        self.pos = {s: k for k, s in enumerate(self.tr)}
        n = len(self.tr)
        self.n = n
        # record(e, coin, current) -> recorded current (mutation hook); default: the literal C record
        self.record = record if record is not None else (lambda e, coin, cur: cur)
        a = [[Fr(0)] * n for _ in range(n)]
        for k, s in enumerate(self.tr):
            a[k][k] += 1
            for p, nxt, e, coin, cur in self.trans[s]:
                if nxt != self.absorbing:
                    a[k][self.pos[nxt]] -= p
        self.lu, self.perm = lu_fr(a)
        self.xy = []
        for s in self.tr:
            st = self.states[s]
            self.xy.append((st.index(q + 1), st.index(q - 1)))
        self.state_of = {xy: k for k, xy in enumerate(self.xy)}

    def mean_current_matrix(self):
        """c[k][e] = one-step mean recorded current on seam e from transient state k."""
        c = [[Fr(0)] * self.g.M for _ in range(self.n)]
        for k, s in enumerate(self.tr):
            for p, nxt, e, coin, cur in self.trans[s]:
                r = self.record(e, coin, cur)
                if r:
                    c[k][e] += p * r
        return c

    def visits(self, x, y):
        """Expected number of events spent in each transient state before annihilation."""
        b = [Fr(0)] * self.n
        b[self.state_of[(x, y)]] = Fr(1)
        return lu_solve_transpose_fr(self.lu, self.perm, b)

    def mean_field(self, x, y, c=None):
        c = c if c is not None else self.mean_current_matrix()
        v = self.visits(x, y)
        return [sum(v[k] * c[k][e] for k in range(self.n) if v[k] != 0 and c[k][e] != 0) for e in range(self.g.M)], v

    def lifetimes(self):
        """h[k] = E[tau] from transient state k."""
        return lu_solve_fr(self.lu, self.perm, [Fr(1)] * self.n)

    def all_mean_fields(self, c=None):
        """m[k][e] = E[F_tau | start k], solved column by column."""
        c = c if c is not None else self.mean_current_matrix()
        cols = [lu_solve_fr(self.lu, self.perm, [c[k][e] for k in range(self.n)]) for e in range(self.g.M)]
        return [[cols[e][k] for e in range(self.g.M)] for k in range(self.n)]

    def second_moments(self, m):
        """S[k] = E[|F_tau|^2 | start k] from the one-step recursion (m = all_mean_fields)."""
        rhs = [Fr(0)] * self.n
        for k, s in enumerate(self.tr):
            acc = Fr(0)
            for p, nxt, e, coin, cur in self.trans[s]:
                r = self.record(e, coin, cur)
                if r:
                    acc += p * r * r
                    if nxt != self.absorbing:
                        acc += 2 * p * r * m[self.pos[nxt]][e]
            rhs[k] = acc
        return lu_solve_fr(self.lu, self.perm, rhs)

    def crossing_counts(self, v):
        """Expected number of unit currents on each seam before and at annihilation."""
        out = [Fr(0)] * self.g.M
        for k, s in enumerate(self.tr):
            if v[k] == 0:
                continue
            for p, nxt, e, coin, cur in self.trans[s]:
                if cur:
                    out[e] += v[k] * p * abs(cur)
        return out


def one_load_objects(g, weights=None):
    """Exact one-load objects from the literal rule: T (column = start), mean current per start port."""
    w = weights if weights is not None else [Fr(1, g.M)] * g.M
    n, m = g.N, g.M
    t = [[Fr(0)] * n for _ in range(n)]
    cur = [[Fr(0)] * m for _ in range(n)]       # cur[z][e] = E[j_e | load at z]
    for z in range(n):
        loads = [0] * n
        loads[z] = 1
        for e, seam in enumerate(g.seams):
            for coin in (0, 1):
                out, (_, j) = ophk.native(tuple(loads), seam, coin)
                t[out.index(1)][z] += w[e] / 2
                if j:
                    cur[z][e] += w[e] / 2 * j
    return t, cur


# ----------------------------------------------------------------------------------------------
# samplers
# ----------------------------------------------------------------------------------------------
def mc_literal(g, x, y, q, nhist, rng, weights=None, unrecorded_coin=None, max_events=10 ** 9):
    """The literal rule, event by event, on integer loads, vectorised over independent histories.

    Every event of every history draws a seam (uniformly, or with `weights`) and a fair coin and applies
    d = n_a - n_b, j = (d - (2c - 1)(|d| mod 2)) // 2.  A history is followed until its loads are uniform.
    unrecorded_coin: mutation control; currents of events with this coin are not written to the record.
    Returns tau (events), F (nhist x M integer records), final loads."""
    seams = np.asarray(g.seams, dtype=np.int64)
    m_, n_ = g.M, g.N
    loads = np.full((nhist, n_), q, dtype=np.int32)
    loads[:, x] += 1
    loads[:, y] -= 1
    rec = np.zeros((nhist, m_), dtype=np.int32)
    tau = np.zeros(nhist, dtype=np.int64)
    exc = np.full(nhist, 2, dtype=np.int32)          # sum over ports of |n - q|
    alive = np.arange(nhist)
    cum = None if weights is None else np.cumsum(np.asarray(weights, float))
    t = 0
    while alive.size and t < max_events:
        t += 1
        k = alive.size
        if cum is None:
            e = rng.integers(0, m_, size=k)
        else:
            e = np.minimum(np.searchsorted(cum, rng.random(k), side="right"), m_ - 1)
        c = rng.integers(0, 2, size=k)
        a = seams[e, 0]
        b = seams[e, 1]
        na = loads[alive, a]
        nb = loads[alive, b]
        d = na - nb
        j = (d - (2 * c - 1) * (np.abs(d) % 2)) // 2
        nz = np.nonzero(j)[0]
        if nz.size:
            idx = alive[nz]
            jj = j[nz]
            an, bn = a[nz], b[nz]
            loads[idx, an] = na[nz] - jj
            loads[idx, bn] = nb[nz] + jj
            if unrecorded_coin is None:
                rec[idx, e[nz]] += jj
            else:
                rec[idx, e[nz]] += jj * (c[nz] != unrecorded_coin)
            exc[idx] += (np.abs(na[nz] - jj - q) + np.abs(nb[nz] + jj - q)
                         - np.abs(na[nz] - q) - np.abs(nb[nz] - q))
            done = idx[exc[idx] == 0]
            if done.size:
                tau[done] = t
                alive = alive[exc[alive] != 0]
    return tau, rec, loads


def neighbour_tables(g):
    """Padded tables: NB[port, slot] neighbour (-1 = no seam), SE seam index, SG = +1 if port is the first
    endpoint of the seam."""
    lists = [[] for _ in range(g.N)]
    for e, (a, b) in enumerate(g.seams):
        lists[a].append((b, e, 1))
        lists[b].append((a, e, -1))
    dmax = max(len(v) for v in lists)
    nb = -np.ones((g.N, dmax), dtype=np.int64)
    se = np.zeros((g.N, dmax), dtype=np.int64)
    sg = np.zeros((g.N, dmax), dtype=np.int64)
    for p, v in enumerate(lists):
        for s, (b, e, sign) in enumerate(v):
            nb[p, s], se[p, s], sg[p, s] = b, e, sign
    return nb, se, sg, dmax


def mc_jump(g, x, y, nhist, rng, tables=None, max_trials=10 ** 8):
    """Exact sampler of the literal rule for one excess (at x) and one deficit (at y) on a uniform filling.

    Events on seams with zero contrast change neither loads nor records (j = 0), so they are not drawn one by
    one.  A trial selects one of the 2 dmax (walker, slot) alternatives uniformly; every literal event is a
    trial with probability dmax / M, independently of the state, and a trial moves a given unit across a given
    seam with probability 1 / (2 dmax), which is 1 / (2M) per literal event.  The event count of a history is
    drawn from the exact negative-binomial law of its trial count.
    Returns dict(tau, trials, currents, F)."""
    nb, se, sg, dmax = tables if tables is not None else neighbour_tables(g)
    xs = np.full(nhist, x, dtype=np.int64)
    ys = np.full(nhist, y, dtype=np.int64)
    rec = np.zeros((nhist, g.M), dtype=np.int16)
    trials = np.zeros(nhist, dtype=np.int64)
    currents = np.zeros(nhist, dtype=np.int64)
    alive = np.arange(nhist)
    k = 0
    while alive.size and k < max_trials:
        k += 1
        r = rng.integers(0, 2 * dmax, size=alive.size)
        w = r >= dmax                      # False: the excess moves; True: the deficit moves
        s = r - dmax * w
        xa = xs[alive]
        ya = ys[alive]
        mover = np.where(w, ya, xa)
        other = np.where(w, xa, ya)
        tgt = nb[mover, s]
        ok = np.nonzero(tgt >= 0)[0]
        if ok.size == 0:
            continue
        idx = alive[ok]
        mv, sv, tg, wv = mover[ok], s[ok], tgt[ok], w[ok]
        # excess moving mover -> target carries one unit along that direction; a deficit moving
        # mover -> target is one unit moving target -> mover
        rec[idx, se[mv, sv]] += (sg[mv, sv] * np.where(wv, -1, 1)).astype(np.int16)
        currents[idx] += 1
        hit = tg == other[ok]
        xs[idx] = np.where(wv, xa[ok], tg)
        ys[idx] = np.where(wv, tg, ya[ok])
        if hit.any():
            done = idx[hit]
            trials[done] = k
            keep = np.ones(alive.size, dtype=bool)
            keep[ok[hit]] = False
            alive = alive[keep]
    rho = dmax / g.M
    tau = trials + rng.negative_binomial(np.maximum(trials, 1), rho)
    return {"tau": tau, "trials": trials, "currents": currents, "F": rec, "rho": rho}


# ----------------------------------------------------------------------------------------------
# float helpers on explicit graphs
# ----------------------------------------------------------------------------------------------
def green_solve(g, rho):
    """phi = L^+ rho for a neutral rho on a connected graph (sparse LU on the grounded Laplacian)."""
    b = g.incidence_sparse()
    lap = (b @ b.T).tocsc()
    n = g.N
    rho = np.asarray(rho, float)
    assert abs(rho.sum()) < 1e-9
    keep = np.arange(1, n)
    lu = spla.splu(lap[keep][:, keep].tocsc())
    phi = np.zeros(n)
    phi[1:] = lu.solve(rho[1:])
    phi -= phi.mean()
    return phi, b


def coulomb_field(g, x, y):
    rho = np.zeros(g.N)
    rho[x] += 1
    rho[y] -= 1
    phi, b = green_solve(g, rho)
    return b.T @ phi, phi, float(phi[x] - phi[y])


# ----------------------------------------------------------------------------------------------
# Bloch sums on the cubic routing
# ----------------------------------------------------------------------------------------------
def cubic_symbol_batch(kx, ky, kz, control=False):
    """H(k) for arrays of momenta (shape (n,)); returns (n, 12, 12) complex."""
    n = len(kx)
    base = np.array(cubic_symbol_zero_int(control), dtype=complex)
    links = cubic_links(control)
    for row in links:                       # remove the k = 0 link entries, they are re-added with phases
        for a, b in row:
            base[a, b] += 1
            base[b, a] += 1
    h = np.broadcast_to(base, (n, 12, 12)).copy()
    for axis, kk in enumerate((kx, ky, kz)):
        ph = np.exp(1j * kk)
        for a, b in links[axis]:
            h[:, a, b] -= ph
            h[:, b, a] -= ph.conj()
    return h


def cubic_green_cell_source(size, control=False, want_diag=False):
    """Phi[p, R] = sum_q L^+[(R, p), (0, q)] / 12: potential at port p of cell R of a unit load spread uniformly
    over the twelve ports of cell 0, on the periodic cubic routing of the given size (Bloch sum, FFT).
    With want_diag also returns diag[p] = L^+[(0, p), (0, p)]."""
    grid = 2 * np.pi * np.arange(size) / size
    ky, kz = np.meshgrid(grid, grid, indexing="ij")
    ky, kz = ky.ravel(), kz.ravel()
    ghat = np.zeros((12, size, size, size), dtype=complex)
    diag = np.zeros(12)
    ones = np.ones((len(ky), 12, 1), dtype=complex)
    h0 = np.array(cubic_symbol_zero_int(control), dtype=float)
    p0 = np.linalg.pinv(h0)
    for i, kx in enumerate(grid):
        h = cubic_symbol_batch(np.full(len(ky), kx), ky, kz, control)
        if i == 0:
            h[0] = np.eye(12)                # placeholder at k = 0, replaced below
        if want_diag:
            inv = np.linalg.inv(h)
            x = inv.sum(axis=2)
            d = np.einsum("nii->ni", inv).real
            if i == 0:
                d[0] = np.diag(p0)
            diag += d.sum(axis=0)
        else:
            x = np.linalg.solve(h, ones)[:, :, 0]
        if i == 0:
            x[0] = 0.0                       # H(0)^+ applied to the uniform vector
        ghat[:, i] = (x.T / 12).reshape(12, size, size)
    phi = np.fft.ifftn(ghat, axes=(1, 2, 3)).real
    if want_diag:
        return phi, diag / size ** 3
    return phi


def cubic_face_flux(phi, axis, control=False):
    """Mean record flux through the face between cell R and cell R + e_axis (sum over its two seams)."""
    out = 0.0
    for a, b in cubic_links(control)[axis]:
        out = out + phi[a] - np.roll(phi[b], -1, axis=axis)
    return out


def ring_green_carrier_source(size):
    """Carrier-averaged potential on the flat ring screen with size x size unit cells (two carriers each):
    phibar[R] = (1/144) sum over ports of A(R) and of A(0) of L^+ (Bloch sum, FFT);
    also returns Phi[p, R] for the 24 ports of the unit cell (unit load spread over the 12 ports of A(0))."""
    b1 = 2 * np.pi * np.array([1 / math.sqrt(3), -1 / 3])
    b2 = 2 * np.pi * np.array([0.0, 2 / 3])
    src = np.zeros(24)
    src[:12] = 1 / 12
    ghat = np.zeros((24, size, size), dtype=complex)
    for i in range(size):
        for j in range(size):
            h = ophk.bloch_laplacian(i / size * b1 + j / size * b2, "ring")
            if i == 0 and j == 0:
                ghat[:, i, j] = np.linalg.pinv(h) @ src
            else:
                ghat[:, i, j] = np.linalg.solve(h, src)
    phi = np.fft.ifft2(ghat, axes=(1, 2)).real
    return phi[:12].mean(axis=0), phi


# ----------------------------------------------------------------------------------------------
# small utilities
# ----------------------------------------------------------------------------------------------
def frs(x):
    return str(x)


def frl(v):
    return [str(x) for x in v]


def matvec_fr(a, v):
    return [sum(a[i][j] * v[j] for j in range(len(v)) if a[i][j] != 0 and v[j] != 0) for i in range(len(a))]


def reff_fr(lp, x, y):
    return lp[x][x] + lp[y][y] - 2 * lp[x][y]
