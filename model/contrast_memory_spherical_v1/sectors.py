"""The law on sectors of fixed load number, the changed laws, the emitter and the reader.

The law on the loads of a carrier with seams `E` and port degrees `d_v`:

    H_kappa = (1/2) sum_e (I - SWAP_e) + kappa C,   C = sum_e (n_a - n_b)^2 + sum_v (6 - d_v) n_v^2
                                                      = 6 sum_v n_v^2 - 2 sum_e n_a n_b.

It conserves the number of loads and the multiset of loads, so it splits into sectors. On loads 0 and 1 with `N`
loads, `C = 6 N - 2 (agreed seams)`: the sector Hamiltonians below are written as `H_kappa - 6 N kappa`, so that
`N` free loads at rest sit at 0 and a pair on a seam lies `2 kappa` lower in the diagonal.

Changed laws used as controls, with the same alphabet, routing and symmetry: the exchange alone (`kappa = 0`) and
the agreement energy reversed (`-kappa`).

Emitter: one fuel state at a seam `e = (a, b)` with energy `E_f` equal to the diagonal of the law at the pair on
`e`, coupled to that pair with strength `g`, a pulse of duration `pi/(2g)` with the law acting throughout.

Reader: every seam has two port observers that report the agreed record "both ends hold one load". On the sector
of two loads the outcomes are one record per seam and "no pair"; fuel left in the emitter counts as "no pair".
On three loads the face reading has one record per face and "no agreed face".
"""
from __future__ import annotations

import itertools
import math

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

G_EMIT = 64.0


def pair_index(n, lo, hi):
    """Index of the binary pair {lo < hi} among the n(n-1)/2 pairs of n ports."""
    return lo * n - lo * (lo + 1) // 2 + (hi - lo - 1)


class PairSector:
    """Two loads on distinct ports of a carrier (binary two-load sector)."""

    def __init__(self, car):
        self.car = car
        n = self.n = car.n
        self.aa, self.bb = np.triu_indices(n, 1)
        self.P = len(self.aa)
        keys = np.sort(np.array([a * n + b for a, b in car.seams], dtype=np.int64))
        self.adjacent = np.isin(self.aa.astype(np.int64) * n + self.bb, keys)
        self.seam_states = pair_index(n, np.array(car.seams)[:, 0], np.array(car.seams)[:, 1])
        self._exchange = None

    def index(self, a, b):
        return int(pair_index(self.n, min(a, b), max(a, b)))

    def exchange(self):
        """`(1/2) sum_e (I - SWAP_e)` on the sector."""
        if self._exchange is None:
            t, deg, n = self.car.table, self.car.degree, self.n
            lines, cols = [], []
            here = np.arange(self.P, dtype=np.int64)
            for mover, other in ((self.aa, self.bb), (self.bb, self.aa)):
                for j in range(t.shape[1]):
                    dest = t[mover, j]
                    keep = (dest >= 0) & (dest != other)
                    lines.append(pair_index(n, np.minimum(dest[keep], other[keep]), np.maximum(dest[keep], other[keep])))
                    cols.append(here[keep])
            lines, cols = np.concatenate(lines), np.concatenate(cols)
            diag = (deg[self.aa] + deg[self.bb] - 2 * self.adjacent) / 2
            off = sparse.csr_matrix((np.full(len(lines), -0.5), (lines, cols)), shape=(self.P, self.P))
            del lines, cols
            self._exchange = (off + sparse.diags(diag)).tocsr()
        return self._exchange

    def hamiltonian(self, kappa):
        """`H_kappa - 12 kappa` on the sector: exchange minus `2 kappa` on every pair on a seam."""
        return (self.exchange() - sparse.diags(2.0 * kappa * self.adjacent)).tocsr()

    def cell_weights(self, prob):
        w = np.zeros(self.car.cells)
        np.add.at(w, self.aa // 12, prob / 2)
        np.add.at(w, self.bb // 12, prob / 2)
        return w


def doubled_hamiltonian(car, kappa):
    """`H_kappa - 12 kappa` on the doubled loads (one port with load 2): `L/2 + 12 kappa`. The exchange moves the
    doubled load as one unit, and the doubled loads form a reducing sector."""
    return (car.laplacian() / 2 + 12.0 * kappa * sparse.identity(car.n)).tocsr()


class BinarySector:
    """`N` loads on distinct ports of a small carrier, with bit-mask states."""

    def __init__(self, car, N):
        self.car, self.N = car, N
        self.states = list(itertools.combinations(range(car.n), N))
        self.index = {s: i for i, s in enumerate(self.states)}
        dim = len(self.states)
        K = np.zeros((dim, dim))
        agreed = np.zeros(dim)
        for i, s in enumerate(self.states):
            occ = set(s)
            for a, b in car.seams:
                if a in occ and b in occ:
                    agreed[i] += 1
                elif (a in occ) != (b in occ):
                    K[i, i] += 0.5
                    K[self.index[tuple(sorted(occ ^ {a, b}))], i] -= 0.5
        self.K, self.agreed = K, agreed

    def hamiltonian(self, kappa, shifted=True):
        """`H_kappa` (or `H_kappa - 6 N kappa` when shifted) as a dense matrix."""
        cost = (0.0 if shifted else 6.0 * self.N) - 2.0 * self.agreed
        return self.K + kappa * np.diag(cost)


class MultisetSector:
    """`N` loads with every port holding 0, 1 or 2 (doubled loads included), dense; for small carriers."""

    def __init__(self, car, N):
        self.car, self.N = car, N
        states = []
        for combo in itertools.combinations_with_replacement(range(car.n), N):
            occ = [0] * car.n
            for v in combo:
                occ[v] += 1
            if max(occ) <= 2:
                states.append(tuple(occ))
        self.states = states
        self.index = {s: i for i, s in enumerate(states)}
        dim = len(states)
        K = np.zeros((dim, dim))
        C = np.zeros(dim)
        deg = car.degree
        for i, s in enumerate(states):
            C[i] = sum((s[a] - s[b]) ** 2 for a, b in car.seams) + sum((6 - deg[v]) * s[v] ** 2 for v in range(car.n))
            for a, b in car.seams:
                if s[a] != s[b]:
                    t = list(s)
                    t[a], t[b] = t[b], t[a]
                    K[i, i] += 0.5
                    K[self.index[tuple(t)], i] -= 0.5
        self.K, self.C = K, C
        self.histogram = [tuple(sorted((x for x in s if x), reverse=True)) for s in states]

    def hamiltonian(self, kappa):
        return self.K + kappa * np.diag(self.C)


class TripleSector:
    """Three loads on distinct ports, vectorized; `H_kappa - 18 kappa` with three free loads at 0."""

    def __init__(self, car):
        n = self.n = car.n
        self.car = car
        tri = np.array(list(itertools.combinations(range(n), 3)), dtype=np.int64)
        keys = np.sort(np.array([a * n + b for a, b in car.seams], dtype=np.int64))
        a, b, c = tri.T
        self.agreed = (np.isin(a * n + b, keys).astype(np.int64) + np.isin(a * n + c, keys) + np.isin(b * n + c, keys))
        self.tri = tri
        self.dim = len(tri)
        t, deg = car.table, car.degree
        lines, cols = [], []
        here = np.arange(self.dim)
        for k in range(3):
            mover, o1, o2 = tri[:, k], tri[:, (k + 1) % 3], tri[:, (k + 2) % 3]
            for j in range(t.shape[1]):
                dest = t[mover, j]
                keep = (dest >= 0) & (dest != o1) & (dest != o2)
                new = np.sort(np.stack([dest[keep], o1[keep], o2[keep]], axis=1), axis=1)
                lines.append(self.rank(new[:, 0], new[:, 1], new[:, 2]))
                cols.append(here[keep])
        lines, cols = np.concatenate(lines), np.concatenate(cols)
        diag = (deg[a] + deg[b] + deg[c] - 2 * self.agreed) / 2
        self.K = (sparse.csr_matrix((np.full(len(lines), -0.5), (lines, cols)), shape=(self.dim,) * 2)
                  + sparse.diags(diag)).tocsr()
        self.faces = self.agreed == 3

    def rank(self, x, y, z):
        n = self.n
        c3 = lambda m: m * (m - 1) * (m - 2) // 6
        c2 = lambda m: m * (m - 1) // 2
        return c3(n) - c3(n - x) + c2(n - x - 1) - c2(n - y) + (z - y - 1)

    def hamiltonian(self, kappa):
        return (self.K - sparse.diags(2.0 * kappa * self.agreed)).tocsr()


# ----------------------------------------------------------------------------------------------------- emitter

def emit(H, p, g=G_EMIT):
    """Fuel at the pair state `p`: pulse `pi/(2g)` of `[[E_f, g], [g, H]] - E_f` with the law acting throughout.

    Returns the state on (fuel, sector) and `E_f` (the diagonal of `H` at `p`)."""
    Ef = float(H[p, p])
    dim = H.shape[0]
    coupling = sparse.csr_matrix(([g], ([0], [p])), shape=(1, dim))
    prep = sparse.bmat([[None, coupling], [coupling.T, H - Ef * sparse.identity(dim, format='csr')]], format='csr')
    seed = np.zeros(dim + 1, complex)
    seed[0] = 1.0
    return expm_multiply(-1j * (math.pi / (2 * g)) * prep, seed), Ef


def emitter_error_bound(da, db, g=G_EMIT):
    """Norm error of the emitter on a seam with end degrees `da`, `db`: `sqrt(da + db - 2)/(2g)`."""
    return math.sqrt(da + db - 2) / (2 * g)


def record_bound(kappa, amax, beta, da, db, g=G_EMIT):
    """`Pr(no pair, doubled load or fuel) <= [2 beta/(2 kappa - amax) + emitter error]^2` (capped at 1)."""
    if 2 * kappa <= amax:
        return 1.0
    return min(1.0, (2 * beta / (2 * kappa - amax) + emitter_error_bound(da, db, g)) ** 2)


def theorem_bound(kappa, g=None):
    """The carrier-free bound `40/(2 kappa - 6)^2`, with the emitter error for degrees up to six when `g` is set."""
    if kappa <= 3:
        return 1.0
    eps = 0.0 if g is None else math.sqrt(10) / (2 * g)
    return min(1.0, (2 * math.sqrt(10) / (2 * kappa - 6) + eps) ** 2)


# ------------------------------------------------------------------------------------------------------ reader

def pair_outcomes(sector, psi):
    """Outcome law of the seam reader: `[P(no pair), P(record on seam 0), ...]` in the order of `car.seams`.
    `psi` is the sector part; any weight outside it (fuel) counts as "no pair"."""
    prob = np.abs(psi) ** 2
    seams = prob[sector.seam_states]
    return np.concatenate([[1.0 - seams.sum()], seams])


def free_run(sector, H, psi, times, first=None, shift=0.0):
    """Unconditional record probabilities at the listed values of the evolution parameter (audit)."""
    Hs = (H - shift * sparse.identity(H.shape[0], format='csr')).tocsr() if shift else H
    lines, prev = [], 0.0
    psi = psi.copy()
    for t in times:
        if t > prev:
            psi = expm_multiply(-1j * (t - prev) * Hs, psi)
            prev = t
        out = pair_outcomes(sector, psi)
        prob = np.abs(psi) ** 2
        on_first = float(prob[first]) if first is not None else 0.0
        cells = sector.cell_weights(prob / max(prob.sum(), 1e-300))
        lines.append({'t': float(t), 'no_pair': float(out[0]), 'pair_on_first_seam': on_first,
                     'pair_on_another_seam': float(out[1:].sum() - on_first),
                     'cells_covered': float(1.0 / np.sum(cells ** 2))})
    return lines, psi


def start_vector(dim, seed=20261118):
    """A fixed generic start vector for Lanczos (a symmetric start would stay in one symmetry sector)."""
    return np.random.default_rng(seed).standard_normal(dim)


# ------------------------------------------------------------------------------------------------- seam events

class SeamEvents:
    """One seam acting once on the sector of two loads: `exp(-i tau [(I - SWAP_e)/2 - 2 kappa n_a n_b])`.

    The event touches only the pairs with exactly one load on the seam (they are exchanged with amplitude
    `(1 -+ e^{-i tau})/2`) and the pair on the seam (phase `e^{2 i kappa tau}`); all other pairs are unchanged."""

    def __init__(self, sector, kappa):
        self.sector, self.kappa = sector, kappa
        n = sector.n
        self.seams = sector.car.seams
        self.tables = []
        for a, b in self.seams:
            xs = np.array([x for x in range(n) if x not in (a, b)])
            ia = pair_index(n, np.minimum(a, xs), np.maximum(a, xs))
            ib = pair_index(n, np.minimum(b, xs), np.maximum(b, xs))
            self.tables.append((ia, ib, int(pair_index(n, a, b))))

    def apply(self, psi, k, tau):
        """Apply the event of seam number `k` with angle `tau` in place."""
        ia, ib, both = self.tables[k]
        ph = np.exp(-1j * tau)
        c, s = (1 + ph) / 2, (1 - ph) / 2
        va, vb = psi[ia], psi[ib]
        psi[ia] = c * va + s * vb
        psi[ib] = c * vb + s * va
        psi[both] *= np.exp(2j * self.kappa * tau)

    def matrix(self, k, tau):
        """The event of seam `k` as a dense matrix (small sectors)."""
        P = self.sector.P
        U = np.zeros((P, P), complex)
        for i in range(P):
            col = np.zeros(P, complex)
            col[i] = 1
            self.apply(col, k, tau)
            U[:, i] = col
        return U

    def start(self, seam=None):
        psi = np.zeros(self.sector.P, complex)
        a, b = self.seams[0] if seam is None else seam
        psi[pair_index(self.sector.n, min(a, b), max(a, b))] = 1
        return psi

    def no_pair(self, psi):
        return float(1 - (np.abs(psi[self.sector.seam_states]) ** 2).sum())


def layers(car):
    """Layers of disjoint seams: a five-colouring of the 30 seams of the cell (by search), used in every cell,
    and the external seams as a sixth layer."""
    cell = sorted({(a % 12, b % 12) for a, b in car.seams if a // 12 == b // 12})
    colour = {}

    def place(k):
        if k == len(cell):
            return True
        a, b = cell[k]
        used = {colour[e] for e in colour if a in e or b in e}
        for c in range(5):
            if c not in used:
                colour[cell[k]] = c
                if place(k + 1):
                    return True
                del colour[cell[k]]
        return False

    if not place(0):
        raise ValueError('no five-colouring of the seams of the cell')
    out = [[] for _ in range(6)]
    for i, (a, b) in enumerate(car.seams):
        out[colour[(a % 12, b % 12)] if a // 12 == b // 12 else 5].append(i)
    out = [layer for layer in out if layer]
    for layer in out:
        ports = [p for i in layer for p in car.seams[i]]
        if len(ports) != len(set(ports)):
            raise ValueError('a layer is no matching')
    return out


def pulse_period(events, layer_list, tau):
    """One symmetric period: layers 1..m-1 for tau/2, layer m for tau, then back down for tau/2."""
    order = [(layer, tau / 2) for layer in layer_list[:-1]] + [(layer_list[-1], tau)] + \
            [(layer, tau / 2) for layer in reversed(layer_list[:-1])]

    def apply(psi):
        for layer, t in order:
            for k in layer:
                events.apply(psi, k, t)
        return psi

    return apply
