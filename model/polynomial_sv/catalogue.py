"""Bound objects of the contrast-memory law and their bindings.

A load costs 6 in `C`, and every agreed seam (both ends loaded) takes 2 off: an object of `N` loads costs
`6 N - 2 (agreed seams)`. Part 1 lists the cheapest objects of `N` loads on one cell and the least gain of each
against every split into two parts. Part 2 gives the binding `B_N = N E_1 - E_N` of the lowest level of `N` loads
on one cell against `N` free loads at finite `kappa` (`E_1 = 6 kappa`). Part 3 places the pair on the federation of
20 cells: the pair level cell by cell, the twenty lowest levels, and the hop between cells.
"""
from __future__ import annotations

import itertools

import numpy as np
from scipy.sparse.linalg import eigsh

from . import carrier as cr
from .common import clean
from .sectors import BinarySector, PairSector, start_vector

SHAPES = {1: 'one port', 2: 'one seam', 3: 'one face', 4: 'two faces on a seam', 5: 'three faces round a port',
          6: 'a port and its five neighbours', 7: 'that cap and one more port', 8: 'the cell less a face and a port',
          9: 'the cell less a face', 10: 'the cell less a seam', 11: 'the cell less a port', 12: 'the whole cell'}
KAPPAS = (0.5, 1.0, 2.0, 4.0, 8.0, 128.0)


def strong_binding():
    car = cr.cell()
    adj = car.adjacency().toarray().astype(bool)
    best = {}
    for N in range(1, car.n + 1):
        top, count = -1, 0
        for S in itertools.combinations(range(car.n), N):
            e = int(adj[np.ix_(S, S)].sum() // 2)
            if e > top:
                top, count = e, 1
            elif e == top:
                count += 1
        best[N] = (top, count, 6 * N - 2 * top)
    lines = {}
    for N, (e, count, C) in best.items():
        gain = min((best[k][2] + best[N - k][2] - C for k in range(1, N // 2 + 1)), default=None)
        lines[str(N)] = {'shape': SHAPES[N], 'how_many': count, 'agreed_seams': e, 'cost': C,
                        'least_gain_against_a_split': gain}
    return lines


def finite_kappa(kappas=KAPPAS, Nmax=6):
    car = cr.cell()
    secs = {N: BinarySector(car, N) for N in range(1, Nmax + 1)}
    out = {}
    for kappa in kappas:
        E1 = 6 * kappa
        entry = {}
        for N, sec in secs.items():
            w, U = np.linalg.eigh(sec.hamiltonian(kappa, shifted=False))
            entry[str(N)] = {'E': float(w[0]), 'multiplicity': int((np.abs(w - w[0]) < 1e-8).sum()),
                           'agreed_seams_in_it': float((U[:, 0] ** 2) @ sec.agreed), 'binding': float(N * E1 - w[0])}
        out[str(kappa)] = entry
    return out


def federation_pair(kappas=(1.0, 2.0, 4.0, 8.0, 16.0)):
    ps = PairSector(cr.federation(0))
    car = ps.car
    a, b = car.external[0]
    c = next(v for v in car.neighbours[a] if v // 12 == a // 12)
    K = ps.exchange()
    two_step = float((K @ K)[ps.index(a, b), ps.index(a, c)])
    cellof = ps.aa // 12
    same_cell = ps.aa // 12 == ps.bb // 12
    lines = {}
    for kappa in kappas:
        H = ps.hamiltonian(kappa)
        local = []
        for cnum in range(car.cells):
            sel = np.flatnonzero(same_cell & (cellof == cnum))
            local.append(float(np.linalg.eigvalsh(H[sel][:, sel].toarray())[0]))
        w = np.sort(eigsh(H, k=21, which='SA', v0=start_vector(H.shape[0]), return_eigenvectors=False))
        lines[str(kappa)] = {'lowest': float(w[0]), 'in_cell_spread': max(local) - min(local),
                            'kinds_of_cell': len({round(x, 6) for x in local}),
                            'width_of_twenty_lowest': float(w[19] - w[0]), 'gap_above_twenty': float(w[20] - w[19]),
                            'two_step_hop_1_over_8kappa': two_step / (2 * kappa)}
    return {'external_seam': [a, b], 'internal_seam': [min(a, c), max(a, c)], 'two_step_element': two_step,
            'by_kappa': lines}


def run(full=False):
    out = {'strong_binding': strong_binding(), 'finite_kappa': finite_kappa()}
    pair = {k: v['2']['agreed_seams_in_it'] for k, v in out['finite_kappa'].items()}
    out['adjacency_in_lowest_pair_level'] = pair
    out['every_object_gains_against_every_split'] = all(
        r['least_gain_against_a_split'] is None or r['least_gain_against_a_split'] >= 2
        for r in out['strong_binding'].values())
    if full:
        out['federation_pair'] = federation_pair()
    return clean(out)


if __name__ == '__main__':
    import json
    print(json.dumps(run(True), indent=1))
