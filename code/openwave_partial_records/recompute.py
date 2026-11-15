"""Independent recomputation of the OpenWave partial results on native records.

The script reads the sealed carrier graph published with the signed-source
package, the archived evidence in ``evidence/`` and the native repair rule.
It imports none of the sealed experiment scripts. The only RER module it
uses is the published CAR field law (``model/defect_field_bridge/full_domain``),
whose branch columns are compared with an independent Jordan-Wigner count.
Every statement about the source is checked in integer or rational
arithmetic; the harmonic, the well spectra, the radial Coulomb levels and
the residence fit are floating-point checks with stated tolerances.

    python3 code/openwave_partial_records/recompute.py
"""
from __future__ import annotations

import csv
import hashlib
import itertools
import json
import math
import random
import sys
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
EVIDENCE = HERE / 'evidence'
GRAPH = (ROOT / 'code' / 'openwave_static_signed_source' / 'evidence' / 'openwave' / 'shared'
         / 'localized_rest_join' / 'evidence' / 'graph.json')
GRAPH_SHA256 = '2e8be4c030e8e04ab903adb5a8a786f11654ecca0380ee54d17509a707f49085'
RECORDS = EVIDENCE / 'row20_residence' / 'records.csv'

V, P, K, E = 240, 256, 1260, 2048      # ports, padded ports, native events, padded events
DIM = P * E


def require(condition, message):
    if not condition:
        raise ValueError(message)


def native_current(d, coin):
    """Native repair current on a seam with load difference d and coin c."""
    return (d - (2 * coin - 1) * (abs(d) % 2)) // 2


def load_graph():
    raw = GRAPH.read_bytes()
    if hashlib.sha256(raw).hexdigest() != GRAPH_SHA256:
        raise ValueError('carrier graph differs from the sealed graph')
    graph = json.loads(raw)
    seams = [tuple(e) for e in graph['seams']]
    require(graph['ports'] == V and len(seams) == 630, 'full carrier graph')
    return graph, seams


def compact_ray(graph):
    """Integer 2f for f = (|4> - |5> - |6> + |7>)/2."""
    fi = np.zeros(P, dtype=np.int64)
    for k, v in graph['f'].items():
        fi[int(k)] = v
    return fi


# ------------------------------------------------------------ rows 3 and 4
def source_arcs(seams):
    """Native event arcs of one positive unit: arc destinations and transfer counts.

    Arc (p, 2*seam + coin) holds the unit at port p before that native event.
    A port off the seam has load difference 0 and current 0, so its arc is held.
    """
    require(native_current(0, 0) == native_current(0, 1) == 0, 'zero current off the seam')
    dest = np.arange(DIM, dtype=np.int64).reshape(P, E)
    counts = np.zeros((V, V), dtype=np.int64)
    counts[np.arange(V), np.arange(V)] = K
    for ei, (a, b) in enumerate(seams):
        for c in (0, 1):
            for p in (a, b):
                j = native_current(int(p == a) - int(p == b), c)
                if j:
                    q = b if j == 1 else a
                    dest[p, 2 * ei + c] = q * E + 2 * ei + (1 - c)
                    counts[p, p] -= 1
                    counts[q, p] += 1
    return dest.ravel(), counts


def clock(graph, seams):
    dest, counts = source_arcs(seams)
    idx = np.arange(DIM)
    require(np.array_equal(dest[dest], idx) and np.array_equal(np.sort(dest), idx),
            'native reversal S is an involution of the complete padded domain')
    fi = compact_ray(graph)
    require(np.array_equal(counts @ fi[:V], 1254 * fi[:V]), 'A*SA f = (209/210) f')
    moved = dest.reshape(P, E)[:V, :K] // E
    contrast = int(((fi[moved] - fi[:V, None]) ** 2).sum())
    lam = 1 - F(contrast, 2 * K * int(fi @ fi))
    require(contrast == 48 and lam == F(209, 210), 'D/Var = 1/210 from the native currents')

    def S(x):
        return x[dest]

    def RK(x):
        """K (2AA* - I) on valid positions and K I on invalid ones, in integers."""
        y = x.reshape(P, E)
        out = -K * y
        out[:V, :K] += 2 * y[:V, :K].sum(axis=1)[:, None]
        out[V:] = K * y[V:]
        return out.ravel()

    # Exact integer invariant span: U = 2 sqrt(K) u, SU, and the triangle cycle v.
    U = np.zeros((P, E), dtype=np.int64)
    U[:, :K] = fi[:, None]
    U = U.ravel()
    SU = S(U)
    require(np.array_equal(S(RK(U)), K * SU), 'W u = S u')
    require(np.array_equal(S(RK(SU)), 2 * 1254 * SU - K * U), 'W S u = 2 lambda S u - u')
    require(int(U @ U) == 4 * K and int(U @ SU) == 4 * 1254 and int(SU @ SU) == 4 * K,
            '<u,u> = 1, <u,Su> = lambda')

    edge_id = {frozenset(e): i for i, e in enumerate(seams)}
    triangle = next(t for t in itertools.combinations(range(12), 3)
                    if any(fi[x] for x in t)
                    and all(frozenset(e) in edge_id for e in itertools.combinations(t, 2)))
    vi = np.zeros(DIM, dtype=np.int64)
    for a, b in zip(triangle, triangle[1:] + triangle[:1]):
        ei = edge_id[frozenset((a, b))]
        c = 0 if a == seams[ei][0] else 1
        z = a * E + 2 * ei + c
        require(dest[z] // E == b, 'actual native arc of the triangle')
        vi[z], vi[dest[z]] = 1, -1
    require(np.array_equal(S(vi), -vi), 'S v = -v')
    require(np.array_equal(RK(vi), -K * vi), 'A* v = 0')
    require(int(U @ vi) == 0 and int(SU @ vi) == 0, 'v is orthogonal to u and Su')

    # Floating harmonic on all position/event outcomes.
    lamf = 209 / 210
    s = math.sqrt(1 - lamf * lamf)
    theta = math.acos(lamf)
    u = U / (2 * math.sqrt(K))
    w = (SU - lamf * U) / (2 * math.sqrt(K) * s)
    v = vi / math.sqrt(6)
    kept = np.flatnonzero((U != 0) | (SU != 0) | (vi != 0))
    witness = next(int(z) for z in np.flatnonzero(vi) if fi[z // E])
    a, b, d = u[kept], w[kept], v[kept]
    coef = np.array([d * d / 2 + (a * a + b * b) / 4, d * a, d * b, (a * a - b * b) / 4, a * b / 2])
    k = int(np.searchsorted(kept, witness))
    c1_closed = -1 / (2 * math.sqrt(7560))
    require(abs(coef[1, k] - c1_closed) < 1e-15, 'witness C1 = v_z u_z = -1/(2 sqrt 7560)')

    def W(x):
        y = x.reshape(P, E)
        out = -y.copy()
        out[:V, :K] += 2 * y[:V, :K].sum(axis=1)[:, None] / K
        out[V:] = y[V:]
        return out.ravel()[dest]

    chi = (u + v) / math.sqrt(2)
    harmonic = 0.0
    total = 0.0
    outside = 0.0
    for n in range(129):
        trig = np.array([1, math.cos(n * theta), math.sin(n * theta),
                         math.cos(2 * n * theta), math.sin(2 * n * theta)])
        harmonic = max(harmonic, float(np.max(np.abs(chi[kept] ** 2 - trig @ coef))))
        total = max(total, abs(float(chi @ chi) - 1))
        outside = max(outside, float(np.sum(chi.reshape(P, E)[12:] ** 2)))
        chi = W(chi)
    require(harmonic < 1e-12 and total < 1e-12, 'theta harmonic on all outcomes')

    support = set((kept // E).tolist())
    envelope = set(np.flatnonzero(fi).tolist()) | set(triangle)
    for p, q in seams:
        if fi[p]:
            envelope.add(q)
        if fi[q]:
            envelope.add(p)
    require(support == envelope and max(envelope) < 12, 'invariant span inside one twelve-port cell')
    return {'lambda': str(lam), 'D_over_Var': str(1 - lam), 'source_current_sum_squares': contrast,
            'full_arc_dimension': DIM, 'body_arcs': V * K, 'charged_columns_checked': V * K,
            'involution_columns': DIM, 'theta': 'acos(209/210)', 'triangle': list(triangle),
            'witness_arc': witness, 'witness_C1': '-1/(2 sqrt(7560))',
            'all_position_event_labels_retained': DIM, 'explicit_nonzero_support_labels': int(len(kept)),
            'walk_relations': ['W u = lambda u + s w', 'W w = -s u + lambda w', 'W v = v'],
            'position_envelope': sorted(envelope), 'within_one_twelve_port_cell': True,
            'harmonic_words_checked': 129, '_theta_float': theta, '_witness': coef[:, k].tolist(),
            '_outside_probability': outside}


# ------------------------------------------------------------ row 6
def jw_hop(mask, i, j):
    """c_j^dagger c_i on an occupation mask, Jordan-Wigner order by mode index."""
    require(mask >> i & 1 and (i == j or not mask >> j & 1), 'hop domain')
    sign = (-1) ** bin(mask & ((1 << i) - 1)).count('1')
    mask ^= 1 << i
    sign *= (-1) ** bin(mask & ((1 << j) - 1)).count('1')
    return mask | 1 << j, sign


def spin(graph, seams):
    if str(ROOT) not in sys.path:
        sys.path.insert(0, str(ROOT))
    from model.defect_field_bridge.full_domain import branch_action
    raw_code = {0: 0, 1: 1, -1: 2}

    def column(p, s, seam, coin, sites):
        a, b = seam
        j = native_current(int(p == a) - int(p == b), coin)
        dest = p if not j else (b if j == 1 else a)
        mask = 1 << (2 * p + s)
        out, sign = jw_hop(mask, 2 * p + s, 2 * dest + s) if j else (mask, 1)
        got = branch_action(mask, 'C', coin, sites=sites, flavors=2, seam=seam)
        require(got == (raw_code[j], out, sign, j), 'C-law CAR column equals J(K_h x I2)')
        return dest, j

    # Full graph: every seam, coin, endpoint and spin (all moving and held endpoint columns).
    endpoint_columns = 0
    for seam in seams:
        for coin in (0, 1):
            for p in seam:
                for s in (0, 1):
                    column(p, s, seam, coin, V)
                    endpoint_columns += 1

    # Small case: the twelve-port cell as a complete CAR space, every column and the operators.
    cell = [e for e in seams if e[1] < 12]
    n = 24
    sigma = [np.array([[0, 1], [1, 0]], dtype=complex), np.array([[0, -1j], [1j, 0]]),
             np.array([[1, 0], [0, -1]], dtype=complex)]
    # Spin generators S^a = sum_p c^dagger_(p,t) (sigma^a/2)_(t,s) c_(p,s) on the one-particle sector;
    # one-particle mode 2p+s is basis index 2p+s.
    spin_ops = []
    for sg in sigma:
        op = np.zeros((n, n), dtype=complex)
        for p in range(12):
            for s in (0, 1):
                for t in (0, 1):
                    if sg[t, s]:
                        out, sign = jw_hop(1 << (2 * p + s), 2 * p + s, 2 * p + t)
                        op[out.bit_length() - 1, 2 * p + s] += sign * sg[t, s] / 2
        spin_ops.append(op)
    casimir = sum(op @ op for op in spin_ops)
    require(np.allclose(casimir, 0.75 * np.eye(n)), 'one-particle CAR spin Casimir 3/4')
    cell_columns = 0
    for seam in cell:
        completeness = np.zeros((n, n))
        for coin in (0, 1):
            branches = defaultdict(lambda: np.zeros((n, n)))
            Ks = defaultdict(lambda: np.zeros((12, 12)))
            for p in range(12):
                for s in (0, 1):
                    dest, j = column(p, s, seam, coin, 12)
                    branches[j][2 * dest + s, 2 * p + s] = 1
                    Ks[j][dest, p] = 1
                    cell_columns += 1
            for j, Vh in branches.items():
                require(np.array_equal(Vh, np.kron(Ks[j], np.eye(2))), 'V_h J = J (K_h x I2)')
                for op in spin_ops:
                    require(np.allclose(Vh @ op, op @ Vh), 'internal field spin is conserved')
                completeness += Vh.T @ Vh / 2
        require(np.array_equal(completeness, np.eye(n)), 'branch effects sum to the identity')

    # Same state: f x chi_s keeps Lf = 6f and the transfer eigenvalue 209/210.
    fi = compact_ray(graph)[:V]
    L = np.zeros((V, V), dtype=np.int64)
    for a, b in seams:
        L[a, a] += 1
        L[b, b] += 1
        L[a, b] -= 1
        L[b, a] -= 1
    require(np.array_equal(L @ fi, 6 * fi), 'L f = 6 f')
    E6 = -1260 * math.log1p(-1 / 210)
    require(abs(E6 - 6.014331228348673) < 1e-12, 'E6 = -1260 log(209/210)')
    return {'full_graph_endpoint_columns': endpoint_columns, 'cell_columns': cell_columns,
            'crossing_sign': 1, 'intertwiner': 'V_h J = J(K_h tensor I_2)', 'spin_casimir': '3/4',
            'spin_commutes_with_every_branch': True, 'branch_effects_sum_to_identity': True,
            'compact_ray_L_eigenvalue': 6, 'transfer_eigenvalue': '209/210',
            'E6': '-1260 log(209/210)'}


# ------------------------------------------------------------ row 8
def contact(seams, sealed_energy):
    removed = 0
    for a, b in seams:
        for j in (1, -1):
            for coin in (0, 1):
                current = native_current(2 * j, coin)
                require(current == j and (j - current, -j + current) == (0, 0), 'contact removes the pair')
                removed += 1
    nb = [[] for _ in range(V)]
    for a, b in seams:
        nb[a].append(b)
        nb[b].append(a)
    seed = seams[0]
    vec = {seed: 1}
    moments, sizes = [1], [1]
    for _ in range(16):
        out = defaultdict(int)
        for (p, r), value in vec.items():
            out[p, r] += (len(nb[p]) + len(nb[r])) * value
            for q in nb[p]:
                if q != r:
                    out[q, r] -= value
            for q in nb[r]:
                if q != p:
                    out[p, q] -= value
        vec = {key: x for key, x in out.items() if x}
        moments.append(vec.get(seed, 0))
        sizes.append(len(vec))
    lower = sum((F(moments[n], n * 1260 ** (n - 1)) for n in range(1, 17)), F(0))
    upper = lower + F(1260) * F(24, 1260) ** 17 / (17 * (1 - F(24, 1260)))
    law = []
    for k in range(17):
        value = F(math.comb(16, k), 24 ** 16) * sum(
            (F(math.comb(16 - k, r) * (-1) ** r * 24 ** (16 - k - r) * moments[k + r])
             for r in range(17 - k)), F(0))
        require(value >= 0, 'pointer effect is positive')
        law.append(value)
    require(sum(law) == 1, 'pointer outcomes sum to one')
    require(moments == sealed_energy['LD_moments'] and sizes == sealed_energy['krylov_support_sizes'],
            'Dirichlet moments')
    require([str(lower), str(upper)] == sealed_energy['mean_energy_interval_exact'], 'energy enclosure')
    require(list(map(str, law)) == sealed_energy['bernstein_energy_pointer_probabilities'], 'energy pointer')
    require(float(lower) == sealed_energy['mean_energy_decimal_diagnostic'] == 11.051928582384567,
            'mean contact energy')
    require(upper - lower < F(1, 10 ** 27), 'enclosure width')

    # Reciprocal-work ring: two swap layers, characteristic field, continuity, contact emission.
    rng = random.Random(20261115)
    wave_checks = 0
    for L in (8, 10, 12):
        u = [rng.choice((-1, 0, 1)) for _ in range(L)]
        v = [rng.choice((-1, 0, 1)) for _ in range(L)]
        frames = []
        for _ in range(3 * L):
            frames.append((u[:], v[:]))
            s0, s1 = v[:], u[:]                     # swap the two slots at every site
            nxt0, nxt1 = s0[:], s1[:]
            for i in range(L):                     # swap lane1 at i with lane0 at i+1
                nxt1[i], nxt0[(i + 1) % L] = s0[(i + 1) % L], s1[i]
            require(nxt0 == [u[i - 1] for i in range(L)] and nxt1 == [v[(i + 1) % L] for i in range(L)],
                    'u moves right and v moves left by one site per cycle')
            for i in range(L):
                flux_in = u[i - 1] ** 2 - v[i] ** 2
                flux_out = u[i] ** 2 - v[(i + 1) % L] ** 2
                require(nxt0[i] ** 2 + nxt1[i] ** 2 - (u[i] ** 2 + v[i] ** 2) == flux_in - flux_out,
                        'local flux continuity')
            u, v = nxt0, nxt1
        for t in range(1, len(frames) - 1):
            for i in range(L):
                for sign in (1, -1):
                    field = [x + sign * y for x, y in zip(*frames[t])]
                    before = [x + sign * y for x, y in zip(*frames[t - 1])]
                    after = [x + sign * y for x, y in zip(*frames[t + 1])]
                    require(after[i] + before[i] == field[i - 1] + field[(i + 1) % L], 'wave equation')
                    wave_checks += 1
        for j in (1, -1):
            u, v = [0] * L, [0] * L
            u[0] = v[0] = native_current(2 * j, 0)
            require((u[0] + v[0], u[0] - v[0], u[0] ** 2 + v[0] ** 2) == (2 * j, 0, 2),
                    'contact writes E = 2j, B = 0, field norm 2 = released load square')
            for _ in range(L // 2):
                u = [u[i - 1] for i in range(L)]
                v = [v[(i + 1) % L] for i in range(L)]
            require(u[L // 2] == v[L // 2] == j, 'records coincide at the antipodal site')
    return {'contact_branches_removed': removed, 'contact_current': 'j for both coins',
            'pair_preparation': list(seed), 'dirichlet_moments': 16, 'pointer_outcomes': len(law),
            'mean_energy_interval_exact': [str(lower), str(upper)],
            'ring_wave_equalities_checked': wave_checks,
            'ring_field': 'u_i(t+1) = u_(i-1)(t), v_i(t+1) = v_(i+1)(t), E = u+v, B = u-v'}


# ------------------------------------------------------------ row 20
def residence(sealed_fit):
    with RECORDS.open() as handle:
        rows = list(csv.DictReader(handle))
    native = [r for r in rows if r['mode'] == '0']
    require(len(rows) == 6048 and len(native) == 5040, 'trajectory count')
    require(all(r['censored'] == '0' for r in rows), 'no censored trajectory')
    shells = [Counter() for _ in range(18)]
    for x, y, z in itertools.product(range(-17, 18), repeat=3):
        r2 = x * x + y * y + z * z
        if r2 < 324:
            shells[math.isqrt(r2)][r2] += 1
    cells = np.array([sum(s.values()) for s in shells[3:15]], dtype=float)
    counts = np.array([[int(r['n' + str(k)]) for k in range(3, 15)] for r in native], dtype=float)
    ports = np.array([int(r['initial_port']) for r in native])
    Y = counts / (12 * cells)
    groups = [Y[ports == p] for p in range(12)]
    require(all(len(g) == 420 for g in groups), 'balanced initial ports')
    mean = sum(g.mean(axis=0) for g in groups) / 12
    cov = sum(np.cov(g, rowvar=False, ddof=1) / len(g) for g in groups) / 144
    Wt = np.linalg.inv(cov)
    squared = np.array(sorted({s for shell in shells[3:15] for s in shell}), dtype=float)
    weights = np.array([[shell.get(int(s), 0) / sum(shell.values()) for s in squared]
                        for shell in shells[3:15]])

    def score(p):
        X = np.column_stack([np.ones(12), weights @ squared ** (-p / 2)])
        ab = np.linalg.solve(X.T @ Wt @ X, X.T @ Wt @ mean)
        if ab[1] <= 0 or np.any(X @ ab <= 0):
            return math.inf, ab
        r = mean - X @ ab
        return float(r @ Wt @ r), ab

    grid = np.linspace(0.3, 2.0, 1701)
    p0 = float(grid[int(np.argmin([score(p)[0] for p in grid]))])
    lo, hi = p0 - 0.002, p0 + 0.002
    g = (math.sqrt(5) - 1) / 2
    while hi - lo > 1e-12:
        x1, x2 = hi - g * (hi - lo), lo + g * (hi - lo)
        if score(x1)[0] < score(x2)[0]:
            hi = x2
        else:
            lo = x1
    p = (lo + hi) / 2
    best, ab = score(p)

    def edge(a, b):
        fa = score(a)[0] - best - 3.841458820694124
        for _ in range(200):
            m = (a + b) / 2
            fm = score(m)[0] - best - 3.841458820694124
            if (fm > 0) == (fa > 0):
                a, fa = m, fm
            else:
                b = m
        return (a + b) / 2

    interval = [edge(0.3, p), edge(2.0, p)]
    require(abs(p - sealed_fit['exponent']) < 1e-7, 'residence exponent')
    require(all(abs(x - y) < 1e-7 for x, y in zip(interval, sealed_fit['profile95_approx'])),
            'profile interval')
    require(abs(ab[1] - sealed_fit['coefficients_a_b'][1]) < 1e-6 and ab[1] > 0, 'positive amplitude')
    return {'native_trajectories': len(native), 'native_currents': int(sum(int(r['jumps']) for r in native)),
            'fit_shells': [3, 14], 'model': 'a + b mean(r^-p)', 'exponent': round(p, 6),
            'profile95': [round(interval[0], 6), round(interval[1], 6)]}


# ------------------------------------------------------------ row 31
def orbital(seams, sealed_probe):
    # Internal bands of the twelve-port cell, exactly: L(L-6)(L^2-10L+20) = 0.
    Lc = np.zeros((12, 12), dtype=np.int64)
    for a, b in (e for e in seams if e[1] < 12):
        Lc[a, a] += 1
        Lc[b, b] += 1
        Lc[a, b] -= 1
        Lc[b, a] -= 1
    I12 = np.eye(12, dtype=np.int64)
    quad = Lc @ Lc - 10 * Lc + 20 * I12
    require(not np.any(Lc @ (Lc - 6 * I12) @ quad), 'cell spectrum 0, 6, 5 +- sqrt 5')
    P6 = -Lc @ quad
    require(np.array_equal(P6 @ P6, 24 * P6) and np.trace(P6) == 120, 'rank-five L = 6 projector')

    # Native first hit on a 5x5x5 box of twelve-port cells: the six external moves are equal.
    cells = list(itertools.product(range(-2, 3), repeat=3))
    ci = {x: i for i, x in enumerate(cells)}
    cell_seams = [e for e in seams if e[1] < 12]
    edges = []
    for x in cells:
        edges.extend((12 * ci[x] + a, 12 * ci[x] + b) for a, b in cell_seams)
        for axis in range(3):
            if x[axis] < 2:
                y = list(x)
                y[axis] += 1
                edges.extend((12 * ci[x] + p, 12 * ci[tuple(y)] + p) for p in range(12))
    nbr = defaultdict(list)
    for a, b in edges:
        for c in (0, 1):                               # the one moving coin of each endpoint
            if native_current(1, c):
                nbr[a].append(b)
            if native_current(-1, c):
                nbr[b].append(a)
    interior = [12 * ci[x] + p for x in cells if 0 < max(map(abs, x)) < 2 for p in range(12)]
    index = {q: i for i, q in enumerate(interior)}
    A = np.eye(len(interior))
    rhs = np.zeros(len(interior))
    for q in interior:
        for r in nbr[q]:
            x = cells[r // 12]
            if x == (0, 0, 0):
                rhs[index[q]] += 1 / len(nbr[q])
            elif max(map(abs, x)) < 2:
                A[index[q], index[r]] -= 1 / len(nbr[q])
    h = np.linalg.solve(A, rhs)
    per_cell = defaultdict(set)
    for q in interior:
        per_cell[cells[q // 12]].add(round(float(h[index[q]]), 12))
    require(all(len(vals) == 1 for vals in per_cell.values()), 'hit probability is constant on each cell')
    exact = {1: F(5, 22), 2: F(1, 11), 3: F(1, 22)}
    for x, vals in per_cell.items():
        require(abs(vals.pop() - float(exact[sum(c != 0 for c in x)])) < 1e-12, 'simple cubic first hit')
    require([str(exact[k]) for k in (1, 2, 3)] == sealed_probe['source_hit_probabilities'],
            'sealed first-hit probabilities')

    # Finite feedback well S = L_periodic(3^3) + lambda (1 - h) below each band threshold.
    outer = list(itertools.product(range(-1, 2), repeat=3))
    oi = {x: i for i, x in enumerate(outer)}
    Ls = np.zeros((27, 27))
    for x in outer:
        for axis in range(3):
            y = list(x)
            y[axis] = (y[axis] + 2) % 3 - 1
            i, j = oi[x], oi[tuple(y)]
            Ls[i, i] += 1
            Ls[j, j] += 1
            Ls[i, j] -= 1
            Ls[j, i] -= 1
    accept = np.array([0.0 if x == (0, 0, 0) else 1 - float(exact[sum(c != 0 for c in x)]) for x in outer])
    wells = []
    for lam, sealed in zip((5 - math.sqrt(5), 6.0, 5 + math.sqrt(5)), sealed_probe['finite_wells']):
        low = float(np.linalg.eigvalsh(Ls + lam * np.diag(accept))[0])
        off = float(np.linalg.eigvalsh(Ls + lam * np.eye(27))[0])
        require(abs(low - sealed['lowest_squared_energy']) < 1e-9 and low < lam and abs(off - lam) < 1e-9,
                'attractive well below its band, removed without feedback')
        wells.append(round(low, 9))

    # Coulomb shells: radial levels of -Laplacian - beta/r at -beta^2/(4 n^2), n^2 states per shell.
    from scipy.linalg import eigh_tridiagonal
    beta, R, N = 1.0, 200.0, 8000
    r = np.arange(1, N + 1) * (R / (N + 1))
    step = R / (N + 1)
    levels = {}
    for l in range(4):
        diag = 2 / step ** 2 + l * (l + 1) / r ** 2 - beta / r
        off = -np.ones(N - 1) / step ** 2
        vals = eigh_tridiagonal(diag, off, select='v', select_range=(-1.0, -0.0125),
                                eigvals_only=True)
        levels[l] = sorted(vals.tolist())
    shells = {}
    for n in range(1, 5):
        target = -beta ** 2 / (4 * n * n)
        members = [l for l in range(4) for x in levels[l] if abs(x / target - 1) < 2e-3]
        require(members == list(range(n)), 'shell n carries l = 0 .. n-1')
        shells[n] = sum(2 * l + 1 for l in members)
        require(shells[n] == n * n, 'n^2 orbital states per shell')
    W = math.sqrt(6) / (32 * math.pi ** 3) * math.gamma(1 / 24) * math.gamma(5 / 24) \
        * math.gamma(7 / 24) * math.gamma(11 / 24)
    c_star = 1 / (4 * math.pi * W / 6)
    require(abs(W - 1.516386059151978) < 1e-12, 'Watson integral')
    return {'cell_bands': ['0', '5-sqrt(5)', '6', '5+sqrt(5)'], 'L6_projector_rank': 5,
            'first_hit': {'(1,0,0)': '5/22', '(1,1,0)': '1/11', '(1,1,1)': '1/22'},
            'hit_constant_on_cells': True, 'well_band_thresholds': ['5-sqrt(5)', '6', '5+sqrt(5)'],
            'shell_multiplicity': {str(n): shells[n] for n in shells},
            'coulomb_coefficient': 'c_* = 1/(4 pi G(0)), G(0) = W/6',
            '_c_star': c_star, '_wells': wells}


def main():
    graph, seams = load_graph()
    energy = json.loads((EVIDENCE / 'neutral_source_completion' / 'RECEIPT.json').read_text())
    fit = json.loads((EVIDENCE / 'row20_residence' / 'RESULT.json').read_text())['fit']
    probe = json.loads((EVIDENCE / 'row31_probe_binding' / 'RECEIPT.json').read_text())
    result = {'clock': clock(graph, seams), 'spin': spin(graph, seams),
              'contact': contact(seams, energy['contact_energy']), 'residence': residence(fit),
              'orbital': orbital(seams, probe)}
    print(json.dumps(result, indent=2))
    return result


if __name__ == '__main__':
    main()
