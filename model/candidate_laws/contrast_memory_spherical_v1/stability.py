"""Particle stability: a pair of loads written on one seam persists in free evolution, read on agreed records.

Observers and records. Every seam has two port observers; they report one agreed record, "both ends hold one
load". The complete outcome set on the sector of two loads is one record per seam and "no pair"; fuel left in the
emitter counts as "no pair". The claim is on these records: a pair written on one seam by the local emitter is
read as a pair on some seam with probability at least one minus the bound at the first reading after free
evolution. It is read on other seams than the first, and the law with the agreement energy removed loses it.

Each sample uses a fresh preparation and an uninterrupted free flight. The sampled parameter values, the levels
and the state vectors below are audit.

    python -m model.candidate_laws.contrast_memory_spherical_v1.stability          # one cell and 20 cells
    python -m model.candidate_laws.contrast_memory_spherical_v1.stability --full   # also 80 and 320 cells, the face object on 20 cells
"""
from __future__ import annotations

import argparse
import json
import math
from operator import index

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh, expm_multiply

from model.executable_model import carrier as cr
from . import theorem
from .common import clean
from .sectors import (G_EMIT, PairSector, SeamEvents, TripleSector, emit, layers, pulse_period,
                      record_bound, start_vector, theorem_bound)

TIMES = (0, 4, 8, 12, 16, 20, 24, 28, 32)
LONG_TIMES = (0, 64, 256, 1024)
SEED = 20261118
FACE = (0, 1, 5)
APART = (0, 130)          # two ports of different cells of the federation of 20 cells


def law(sector, kappa, kind='law'):
    """The law and the changed laws with the same alphabet, routing and symmetry."""
    if kind == 'law':
        return sector.hamiltonian(kappa)
    if kind == 'exchange alone':
        return sector.hamiltonian(0.0)
    if kind == 'agreement energy reversed':
        return sector.hamiltonian(-kappa)
    raise ValueError(kind)


def lowest_level(H, adjacent, dense_limit=5000, tol=1e-8):
    """Lowest level of the sector (two free loads at rest at 0), its multiplicity, and the probability of a pair on
    a seam in it, averaged over its eigenspace (independent of the basis). Audit on the spectrum.

    Small sectors use dense diagonalization. Large ones use Lanczos for the smallest algebraic level, which needs
    no shift. The off-diagonal elements of the law on two loads are -1/2 and its moves connect every pair of
    ports, so by the Perron-Frobenius theorem the lowest level is simple with a one-signed eigenvector; the sparse
    path checks that sign (to rounding) and fails otherwise."""
    if H.shape[0] <= dense_limit:
        w, U = np.linalg.eigh(H.toarray())
        deg = int((w - w[0] < tol).sum())
    else:
        w, U = eigsh(H, k=1, which='SA', v0=start_vector(H.shape[0]))
        v = U[:, 0] * np.sign(U[:, 0].sum())
        if not v.min() > -1e-8 * v.max():
            raise ValueError('the lowest Lanczos vector is not one-signed')
        deg = 1
    weight = float((U[adjacent, :deg] ** 2).sum() / deg)
    return float(w[0]), weight, deg


def summary(lines):
    return {'max_no_pair': max(r['no_pair'] for r in lines),
            'max_pair_on_another_seam': max(r['pair_on_another_seam'] for r in lines),
            'cells_covered_at_last_reading': lines[-1]['cells_covered']}


_SECTORS = {}
_RUNS = {}


def sector(key):
    if key not in _SECTORS:
        _SECTORS.clear()
        _SECTORS[key] = PairSector(cr.carrier(key))
    return _SECTORS[key]


def law_run(key, kappa, kind='law'):
    """Emitter on the carrier's emitter seam, then the readings at TIMES; cached per law."""
    tag = (key, float(kappa), kind)
    if tag not in _RUNS:
        ps = sector(key)
        _RUNS[tag] = emitted_run(ps, law(ps, kappa, kind), cr.EMITTER_SEAM[key])
    return _RUNS[tag]


def emitted_run(ps, H, seam, times=TIMES):
    """Read fresh copies at each specified exposure; keep the fuel as its own refusal."""
    a, b = sorted(index(v) for v in seam)
    if not (0 <= a < b < ps.n and b in ps.car.neighbours[a]):
        raise ValueError("Emitter needs a seam between two distinct ports of the carrier")
    times = tuple(float(t) for t in times)
    if any(not math.isfinite(t) or t < 0 for t in times):
        raise ValueError("Free exposures must be finite and nonnegative")
    p = ps.index(a, b)
    state, Ef = emit(H, p)
    fuel = float(abs(state[0])**2)
    psi = last_state = state[1:]
    flight = H - Ef*sparse.eye(ps.P, format="csr")
    other_seams = ps.seam_states[ps.seam_states != p]
    previous = 0.
    records = {}
    for t in sorted(set(times)):
        if t > previous:
            psi = expm_multiply((-1j*(t-previous))*flight, psi)
        previous = t
        prob = np.abs(psi)**2
        material = prob.sum()
        if not math.isclose(material + fuel, 1., rel_tol=0, abs_tol=1e-8):
            raise ArithmeticError("Pair flight lost preparation probability")
        cells = ps.cell_weights(prob/material) if material else None
        records[t] = {
            "t": t, "no_pair": float(fuel + prob[~ps.adjacent].sum()),
            "pair_on_first_seam": float(prob[p]),
            "pair_on_another_seam": float(prob[other_seams].sum()),
            "cells_covered": float(1/(cells @ cells)) if cells is not None else None,
        }
        if t == times[-1]:
            last_state = psi
    return [records[t].copy() for t in times], fuel, last_state, Ef


def law_table(key, kappas, consts=None):
    """The headline table on one carrier: bound, largest P(no pair), largest P(pair on another seam)."""
    ps = sector(key)
    car = ps.car
    c = consts or theorem.constants(car)
    amax, beta = c['adjacent_block_largest_level'], c['coupling_largest_singular_value']
    seam = cr.EMITTER_SEAM[key]
    da, db = (int(car.degree[v]) for v in seam)
    out = {'pair_states': ps.P, 'emitter_seam': list(seam), 'emitter_seam_degrees': [da, db],
           'adjacent_block_largest_level': amax, 'coupling_largest_singular_value': beta, 'by_kappa': {}}
    for kappa in kappas:
        lines, fuel, _, Ef = law_run(key, kappa)
        s = summary(lines)
        s.update(bound_with_carrier_constants=record_bound(kappa, amax, beta, da, db),
                 bound_of_the_theorem=theorem_bound(kappa, G_EMIT), fuel_left=fuel,
                 readings=[{k: r[k] for k in ('t', 'no_pair', 'pair_on_another_seam')} for r in lines])
        s['within_bound'] = s['max_no_pair'] <= s['bound_with_carrier_constants']
        out['by_kappa'][str(kappa)] = s
    return out


def long_parameter(kappa=8.0):
    ps = sector(0)
    lines, _, _, _ = emitted_run(ps, law(ps, kappa), cr.EMITTER_SEAM[0], LONG_TIMES)
    return [{k: r[k] for k in ('t', 'no_pair', 'pair_on_another_seam', 'cells_covered')} for r in lines]


def without_the_property(kappas=(8.0, 128.0)):
    """Two loads on ports of different cells: the largest probability of an agreed pair record."""
    ps = sector(0)
    out = {}
    for kappa in kappas:
        H = law(ps, kappa)
        psi = np.zeros(ps.P, complex)
        psi[ps.index(*APART)] = 1
        worst, prev = 0.0, 0
        for t in TIMES:
            if t > prev:
                psi = expm_multiply(-1j * (t - prev) * H, psi)
                prev = t
            worst = max(worst, float((np.abs(psi[ps.seam_states]) ** 2).sum()))
        out[str(kappa)] = worst
    return {'ports': list(APART), 'largest_P_agreed_pair': out}


def changed_laws(kappa=8.0, carriers=('cell', 0)):
    """The law, the exchange alone and the reversed agreement energy: free persistence and the lowest level."""
    out = {}
    for key in carriers:
        ps = sector(key)
        rows_k = {}
        for kind in ('law', 'exchange alone', 'agreement energy reversed'):
            H = law(ps, kappa, kind)
            lines, _, _, _ = law_run(key, kappa, kind)
            s = summary(lines)
            s['lowest_level'], s['pair_in_lowest_level'], s['lowest_level_multiplicity'] = lowest_level(H, ps.adjacent)
            s['level_below_two_free_loads'] = s['lowest_level'] < -1.0 and s['pair_in_lowest_level'] > 0.5
            rows_k[kind] = s
        out['one cell' if key == 'cell' else cr.NAMES[key]] = rows_k
    return out


def face_object(keys=('cell',), kappa=8.0):
    """Three loads on the three ports of one face: P(no agreed face) and P(agreed record on another face)."""
    out = {}
    for key in keys:
        car = cr.carrier(key)
        ts = TripleSector(car)
        start = int(ts.rank(*FACE))
        if not ts.faces[start]:
            raise ValueError('the start is no face')
        res = {'states': ts.dim}
        for kind, k in (('law', kappa), ('exchange alone', 0.0)):
            H = ts.hamiltonian(k)
            psi = np.zeros(ts.dim, complex)
            psi[start] = 1
            prev, no_face, other = 0, [], []
            for t in TIMES:
                if t > prev:
                    psi = expm_multiply(-1j * (t - prev) * H, psi)
                    prev = t
                p = np.abs(psi) ** 2
                no_face.append(float(1 - p[ts.faces].sum()))
                other.append(float(p[ts.faces].sum() - p[start]))
            res[kind] = {'max_no_face': max(no_face), 'max_other_face': max(other)}
        out['one cell' if key == 'cell' else cr.NAMES[key]] = res
        del ts
    return out


def single_readings(kappa=8.0, runs=100000, t=32.0):
    """Sampled single readings on 20 cells at one value of the parameter, each from a fresh run."""
    ps = sector(0)
    seam = cr.EMITTER_SEAM[0]
    if t != TIMES[-1]:
        raise ValueError('single-reading time differs from final evolution time')
    lines, _, psi, _ = law_run(0, kappa)
    prob = np.abs(psi) ** 2
    outcomes = np.clip(np.concatenate([[1 - prob[ps.seam_states].sum()], prob[ps.seam_states]]), 0, None)
    rng = np.random.default_rng(SEED)
    draws = rng.choice(len(outcomes), size=runs, p=outcomes / outcomes.sum())
    counts = np.bincount(draws, minlength=len(outcomes))
    first = 1 + ps.car.seams.index(tuple(seam))
    return {'runs': runs, 'seed': SEED, 't': t, 'no_pair_runs': int(counts[0]), 'exact_no_pair': float(outcomes[0]),
            'standard_error': math.sqrt(outcomes[0] * (1 - outcomes[0]) / runs),
            'distinct_seams_read': int((counts[1:] > 0).sum()), 'runs_on_first_seam': int(counts[first]),
            'exact_first_seam': float(outcomes[first])}


def layered_pulses(kappa=8.0, settings=(('cell', 32.0, (1 / 64, 1 / 256)), (0, 8.0, (1 / 64,)))):
    """The law run as pulses on layers of disjoint seams, both orders of the layers, from a pair on a seam."""
    out = {}
    for key, T, taus in settings:
        car = cr.carrier(key)
        ps = PairSector(car)
        ev = SeamEvents(ps, kappa)
        lay = layers(car)
        psi0 = ev.start()
        cont = float(1 - (np.abs(expm_multiply(-1j * T * law(ps, kappa), psi0)[ps.seam_states]) ** 2).sum())
        lines = {'T': T, 'continuous': cont}
        for tau in taus:
            for name, order in (('order 1', lay), ('order 2', lay[::-1])):
                step = pulse_period(ev, order, tau)
                psi = psi0.copy()
                for _ in range(int(round(T / tau))):
                    psi = step(psi)
                lines[f'1/{round(1 / tau)} {name}'] = ev.no_pair(psi)
        out['one cell' if key == 'cell' else cr.NAMES[key]] = lines
    return out


def run(full=False):
    _RUNS.clear()
    consts = {'cell': None}
    out = {'law': {}}
    plan = [('cell', (8.0, 128.0)), (0, (8.0, 128.0) if full else (8.0,))]
    if full:
        plan += [(1, (8.0, 128.0)), (2, (8.0,))]
    for key, kappas in plan:
        out['law']['one cell' if key == 'cell' else cr.NAMES[key]] = law_table(key, kappas, consts.get(key))
        if key == 0:
            out['single_readings'] = single_readings()
    out['changed_laws'] = changed_laws()
    out['layered_pulses'] = layered_pulses(settings=(('cell', 32.0, (1 / 64, 1 / 256)), (0, 8.0, (1 / 64,)))
                                           if full else (('cell', 32.0, (1 / 64, 1 / 256)),))
    out['face_object'] = face_object(('cell', 0) if full else ('cell',))
    if full:
        out['long_parameter'] = long_parameter()
        out['without_the_property'] = without_the_property()
    return clean(out)


def verdict(result) -> dict:
    """The four parts of the claim, from a result of `run`."""
    law_rows = [r for c in result['law'].values() for r in c['by_kappa'].values()]
    ch = result['changed_laws']
    return {
        'pair_kept_within_the_bound': all(r['within_bound'] and r['max_no_pair'] < 0.05 for r in law_rows),
        'pair_read_on_other_seams': all(r['max_pair_on_another_seam'] > 0.05 for r in law_rows),
        'level_below_two_free_loads_under_the_law': all(c['law']['level_below_two_free_loads'] for c in ch.values()),
        'exchange_alone_loses_the_pair': all(c['exchange alone']['max_no_pair'] > 0.5 for c in ch.values()),
        'reversed_sign_has_no_level_below_two_free_loads': all(
            not c['agreement energy reversed']['level_below_two_free_loads'] for c in ch.values()),
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--full', action='store_true')
    args = ap.parse_args(argv)
    res = run(args.full)
    res['verdict'] = verdict(res)
    print(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
