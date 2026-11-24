"""Grow routes from local cell geometry and test native one-load propagation."""
import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.sparse import diags
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo.openwave.row_30_autonomous.channel import (
    K, column_checks, matrices,
)
from model.candidate_laws.interval_echo.openwave.row_30_route.experiment import (
    portable, shortest_internal,
)
from model.executable_model import federation


HERE = Path(__file__).resolve().parent
SHOTS = 262144


def native_routing(car):
    gateways = defaultdict(list)
    across = {}
    for a, b in car.external:
        across[a] = b
        across[b] = a
        gateways[a // 12].append(a)
        gateways[b // 12].append(b)
    options = {}
    for ports in gateways.values():
        assert len(ports) == 3
        for entry in ports:
            options[entry] = [shortest_internal(car, entry, exit_port)
                              for exit_port in sorted(ports) if exit_port != entry]
    return options, across


def grow(options, across, entry, rule, high):
    history = []
    seen = set()
    while entry // 12 not in seen:
        seen.add(entry // 12)
        choices = options[entry]
        if rule == 'antipodal':
            choices = [path for path in choices if len(path) == 4]
            if not choices:
                return history, 'no antipodal exit'
        length = (min if rule == 'nearest' else max)(map(len, choices))
        tied = [p for p in choices if len(p) == length]
        selected = (max if high else min)(tied, key=lambda p: p[-1])
        history.append(selected)
        entry = across[selected[-1]]
    return history, 'revisited cell'


def bloch(k):
    h = np.diag([-float(K), 0., -float(K)]).astype(complex)
    h[0,1] = h[1,0] = h[1,2] = h[2,1] = -K / 2
    h[0,2] = -K / 2 * np.exp(-1j*k)
    h[2,0] = h[0,2].conjugate()
    return np.linalg.eigvalsh(h)


def gaps():
    edges = np.array([bloch(0.), bloch(np.pi)])
    return [(float(edges[:,i].max()), float(edges[:,i+1].min())) for i in (0,1)]


def close_to_trimer(car, route):
    return all(tuple(car.degree[p]) == (6,5,6) for p in route)


def timed_clicks(diag, off, e, vec, sites, cells, seed):
    source = np.zeros(len(sites))
    source[0] = 1
    op = diags((off, diag, off), (-1,0,1), format='csr')
    rng = np.random.default_rng(seed)
    reads = []
    for time in (0.,1.,4.):
        psi = expm_multiply(-1j*time*op, source)
        error = np.linalg.norm(psi - vec @ (np.exp(-1j*time*e)*vec[0]))
        assert error < 1e-10
        probabilities = abs(psi)**2
        assert abs(probabilities.sum()-1.) < 1e-10
        clicks = rng.multinomial(SHOTS, probabilities / probabilities.sum())
        reads.append(dict(time=time, solver_error=error, total_shots=SHOTS,
                          cell_detectors=[dict(offset=o, cell=cells[o],
                            probability=probabilities[sites//12==cells[o]].sum(),
                            clicks=clicks[sites//12==cells[o]].sum())
                            for o in (1,3) if o < len(cells)]))
    return reads


def spectrum_gap(e, vec, intervals):
    rows = []
    for a,b in intervals:
        interior = (e > a+1e-9) & (e < b-1e-9)
        rows.append(dict(inside=int(interior.sum()), total_source_weight=float(
            (vec[0, interior]**2).sum()), middle_half_mean_projector=float(
            np.sum(vec[len(e)//4:len(e)-len(e)//4, interior]**2,axis=1).mean())))
    return rows


def test_route(car, route, seed):
    if not route:
        return None
    sites = np.concatenate(route)
    assert len(sites) == len(set(sites))
    assert all(b in car.neighbours[a] for a,b in zip(sites,sites[1:]))
    cells = [p[0]//12 for p in route]
    out = dict(cells=cells, sites=sites, segment_words=[car.degree[p].tolist() for p in route])
    for control, degree in [('native', True), ('degree_reference_removed', False)]:
        d, off = matrices(car,sites,degree)
        e, vec = eigh_tridiagonal(d.astype(float),off.astype(float))
        r = dict(interior_gaps=spectrum_gap(e,vec,gaps()), trace=timed_clicks(
            d,off,e,vec,sites,cells,seed+int(degree)))
        if degree and close_to_trimer(car,route) and len(route)>1:
            # Every interior k=n*pi/M has three standing modes of the symmetric,
            # finite trimer path. This check sets k using only the observed cell count.
            k = np.pi*np.arange(1,len(route))/len(route)
            bands = np.stack([bloch(q) for q in k])
            mismatch = np.min(abs(bands[:,:,None]-e), axis=2)
            r['standing_spectrum'] = dict(interior_wavenumbers=k,
                predicted_lines=bands, max_nearest_line_error=mismatch.max(),
                distinct_matches=len(set(np.argmin(abs(bands[:,:,None]-e),axis=2).flat)))
        out[control]=r
    return out


def run():
    all_boards = {}
    for level in (0,1,2):
        car = federation(level)
        options, across = native_routing(car)
        rows = {}
        for rule,high in [('antipodal',False), ('farthest',False), ('farthest',True),
                          ('nearest',False), ('nearest',True)]:
            lengths, words, reasons, best = Counter(), Counter(), Counter(), []
            trimer = Counter()
            for start in sorted(options):
                route, stop = grow(options,across,start,rule,high)
                lengths[len(route)] += 1
                reasons[stop] += 1
                for p in route:
                    words[','.join(map(str, car.degree[p]))] += 1
                if route and close_to_trimer(car,route):
                    trimer[len(route)] += 1
                if len(route)>len(best):
                    best = route
            rows[rule+('_high' if high else '_low')] = dict(start_count=len(options),
                length_histogram=lengths, stop_histogram=reasons, segment_visits=words,
                all_trimer_histogram=trimer, longest=test_route(car,best,2712+level))
        all_boards[car.cells] = rows
    return portable(dict(native_column_errors=column_checks(), shots=SHOTS,
                         trimer_bulk_gaps=gaps(), boards=all_boards))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    result=run()
    if args.output:
        with args.output.open('x') as out:
            json.dump(result,out,separators=(',',':'),allow_nan=False)
            out.write('\n')
    elif result != json.loads((HERE/'results.json').read_text()):
        raise SystemExit('Local-route experiment differs from receipt')
    for size,rows in result['boards'].items():
        for name,sample in rows.items():
            r=sample['longest']
            print(size,name, max(map(int,sample['length_histogram'])), sum(sample['all_trimer_histogram'].values()),
                  [i['inside'] for i in r['native']['interior_gaps']] if r else None)


if __name__ == '__main__':
    main()
