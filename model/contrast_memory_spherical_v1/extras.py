"""Further values the contrast-memory law gives on its carriers (audit on the operators and the spectra).

(a) one load: the levels of `L/2` on one cell and on the federations of 20, 80 and 320 cells;
(b) three and four loads on one cell, doubled loads included: the lowest level, its binding against every split
    into smaller objects, the fragmentation threshold, and the number of levels below it;
(c) the pair on the federation of 20 cells at kappa 2, 8, 128: levels below two free loads, the band of the
    lowest 630 levels, and the effective hop of a pair between an internal seam and the external seam beside it,
    against `1/(8 kappa)`;
(d) one cell: the average over the evolution parameter of the probability that a pair started on a seam is read
    on a seam, for kappa 0 to 8;
(e) a ring of 24 ports: the lowest two-load level of every total momentum `K` against the branch
    `E(K) = Delta - cos^2(K/2)/Delta`, `Delta = 1 + 2 kappa`.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.sparse.linalg import eigsh

from . import carrier as cr
from .common import clean
from .sectors import MultisetSector, PairSector, start_vector


def grouped(w, tol=1e-8):
    """Distinct levels with multiplicities, from sorted levels."""
    out = []
    for x in w:
        if out and abs(x - out[-1][0]) < tol:
            out[-1][1] += 1
        else:
            out.append([float(x), 1])
    return out


# (a) ---------------------------------------------------------------------------------------------------------

def one_load(levels=(0, 1, 2), shown=12):
    out = {}
    for key in ('cell',) + tuple(levels):
        car = cr.carrier(key)
        w = np.linalg.eigvalsh(car.laplacian().toarray() / 2)
        g = grouped(w)
        nonzero = [x for x in g if x[0] > 1e-9]
        out['one cell' if key == 'cell' else cr.NAMES[key]] = {
            'ports': car.n, 'lowest_distinct_levels': g[:shown], 'gap': nonzero[0][0],
            'levels_below_1': int((w < 1).sum()), 'largest_level': float(w[-1])}
    s5 = math.sqrt(5)
    out['one_cell_closed_form'] = {'levels': [0.0, (5 - s5) / 2, 3.0, (5 + s5) / 2], 'multiplicities': [1, 3, 5, 3]}
    return out


# (b) ---------------------------------------------------------------------------------------------------------

def partitions(n, largest=None):
    largest = n if largest is None else largest
    if n == 0:
        yield ()
        return
    for k in range(min(n, largest), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest


def few_loads(kappas=(0.5, 1.0, 2.0, 4.0, 8.0, 128.0), Nmax=4):
    car = cr.cell()
    secs = {N: MultisetSector(car, N) for N in range(1, Nmax + 1)}
    out = {'states': {str(N): len(s.states) for N, s in secs.items()},
           'binary_states': {str(N): sum(1 for h in s.histogram if max(h) == 1) for N, s in secs.items()}}
    for kappa in kappas:
        spectra = {N: np.linalg.eigh(s.hamiltonian(kappa)) for N, s in secs.items()}
        E = {N: float(sp[0][0]) for N, sp in spectra.items()}
        entry = {}
        for N in range(3, Nmax + 1):
            w, U = spectra[N]
            sec = secs[N]
            by_hist = {}
            for h in sorted(set(sec.histogram)):
                idx = [i for i, x in enumerate(sec.histogram) if x == h]
                by_hist['+'.join(map(str, h))] = float(np.linalg.eigvalsh(sec.hamiltonian(kappa)[np.ix_(idx, idx)])[0])
            splits = {'+'.join(map(str, p)): sum(E[k] for k in p) - E[N] for p in partitions(N) if len(p) > 1}
            cheapest = min(splits, key=splits.get)
            threshold = E[N] + splits[cheapest]
            entry[str(N)] = {'lowest_level': E[N], 'lowest_by_load_multiset': by_hist,
                           'binding_against_splits': splits, 'fragmentation_channel': cheapest,
                           'fragmentation_threshold': splits[cheapest],
                           'levels_below_threshold': int((w < threshold - 1e-9).sum()),
                           'doubled_load_weight_in_lowest': float(
                               (U[[i for i, x in enumerate(sec.histogram) if max(x) == 2], 0] ** 2).sum())}
        entry['E1'], entry['E2'] = E[1], E[2]
        out[str(kappa)] = entry
    return out


# (c) ---------------------------------------------------------------------------------------------------------

def pair_band(kappas=(2.0, 8.0, 128.0), extra=10):
    ps = PairSector(cr.federation(0))
    car = ps.car
    m = len(car.seams)
    a, b = car.external[0]
    c = next(v for v in car.neighbours[a] if v // 12 == a // 12)
    i_ext, i_int = car.seams.index((a, b)), car.seams.index((min(a, c), max(a, c)))
    out = {'external_seam': [a, b], 'internal_seam': [min(a, c), max(a, c)], 'band_size': m}
    for kappa in kappas:
        H = ps.hamiltonian(kappa)
        w, V = eigsh(H, k=m + extra, which='SA', v0=start_vector(H.shape[0]))
        order = np.argsort(w)
        w, V = w[order], V[:, order]
        below = int((w < 0).sum())
        if not below < len(w):
            raise ValueError('more levels below two free loads than computed')
        entry = {'levels_below_two_free_loads': below, 'lowest': float(w[0]), 'level_630': float(w[m - 1]),
               'level_631': float(w[m]), 'gap_above_lowest_630': float(w[m] - w[m - 1])}
        isolated = entry['gap_above_lowest_630'] > 0.1
        entry['lowest_630_isolated'] = isolated
        if isolated:
            # des Cloizeaux effective Hamiltonian of the isolated band on the seam states
            B = V[ps.seam_states, :m]
            ev, evec = np.linalg.eigh(B @ B.T)
            S_inv_half = (evec / np.sqrt(ev)) @ evec.T
            Heff = S_inv_half @ (B * w[:m]) @ B.T @ S_inv_half
            hop = float(Heff[i_ext, i_int])
            entry.update(band_width=float(w[m - 1] - w[0]), smallest_seam_overlap=float(np.sqrt(ev[0])),
                       effective_hop_external_internal=hop, one_over_8kappa=1 / (8 * kappa),
                       hop_over_minus_one_over_8kappa=hop * (-8 * kappa))
        out[str(kappa)] = entry
    return out


# (d) ---------------------------------------------------------------------------------------------------------

def average_adjacency(kappas=(0.0, 0.5, 1.0, 2.0, 4.0, 8.0)):
    ps = PairSector(cr.cell())
    start = ps.index(0, 1)
    out = {}
    for kappa in kappas:
        w, U = np.linalg.eigh(ps.hamiltonian(kappa).toarray())
        total, i = 0.0, 0
        while i < len(w):
            j = i
            while j < len(w) and abs(w[j] - w[i]) < 1e-9:
                j += 1
            phi = U[:, i:j] @ U[start, i:j]            # projection of the start state on the level
            total += float((phi[ps.adjacent] ** 2).sum())
            i = j
        out[str(kappa)] = total
    return out


# (e) ---------------------------------------------------------------------------------------------------------

def ring_branch(L=24, kappas=(0.5, 1.0, 2.0, 8.0)):
    ring = cr.Carrier('ring', L, [(x, (x + 1) % L) for x in range(L)])
    ps = PairSector(ring)
    blocks = []
    for ell in range(L):
        K = 2 * np.pi * ell / L
        vecs = []
        for r in range(1, L // 2 + 1):
            period = L // 2 if 2 * r == L else L
            if (ell * period) % L:
                continue
            v = np.zeros(ps.P, complex)
            for x in range(period):
                v[ps.index(x, (x + r) % L)] = np.exp(1j * K * x) / np.sqrt(period)
            vecs.append(v)
        blocks.append((K, np.column_stack(vecs)))
    Q = np.column_stack([q for _, q in blocks])
    iso = float(np.abs(Q.conj().T @ Q - np.eye(ps.P)).max())
    out = {'ports': L, 'binary_pair_states': ps.P, 'momentum_basis_isometry_error': iso}
    for kappa in kappas:
        H = ps.hamiltonian(kappa).toarray() + 4 * kappa * np.eye(ps.P)   # chain law without degree term
        delta = 1 + 2 * kappa
        worst, resid, lines = 0.0, 0.0, []
        for K, q in blocks:
            h = q.conj().T @ H @ q
            resid = max(resid, float(np.abs(H @ q - q @ h).max()))
            e = float(np.linalg.eigvalsh(h)[0])
            formula = delta - math.cos(K / 2) ** 2 / delta
            worst = max(worst, abs(e - formula))
            lines.append([K, e, formula])
        out[str(kappa)] = {'largest_difference': worst, 'block_residual': resid,
                           'rest_level_K0': lines[0][1], 'formula_K0': lines[0][2]}
    return out


def run(full=True):
    out = {'one_load': one_load((0, 1, 2) if full else (0,)), 'few_loads': few_loads(),
           'average_adjacency': average_adjacency(), 'ring_branch': ring_branch()}
    if full:
        out['pair_band'] = pair_band()
    return clean(out)


if __name__ == '__main__':
    import json
    print(json.dumps(run(), indent=1))
