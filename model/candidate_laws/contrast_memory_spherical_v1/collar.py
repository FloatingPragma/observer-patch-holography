"""The unrecovered collar: a seam whose private memory is not blank, under the complete local law.

The collar of a seam is its private memory of eleven states; recovery is its return to blank. One cell, the seam
`e0 = (0, 1)` with its complete memory, every other seam recovered (blank memory, so it acts through its seam term);
0, 1 and 2 loads. The local law of `e0` on loads x memory is `h_kappa = U G_kappa U` from `law.py`.

Probe (audit on the operators): the law has no element between "memory blank" and "memory not blank"; the levels
of a collar with no load anywhere; one load moves on the carrier without the seam; beside a record of size one a
load on a port of `e0` lies below a load elsewhere, by about `2 kappa` where the record fits the load and by about
`kappa` where it does not.

Reading: a local fuel at level 0 is converted into a contrast record of the collar by a pulse of area `pi/4` with
the law acting; after `q` cycles a second pulse with one of four phase settings reads the fuel. The record is
"fuel present". With `Z(q) = [p_0 - p_pi - i (p_pi/2 - p_3pi/2)]/4` the observer reads the phase rate per cycle of
the symmetric and of the antisymmetric record; their difference is the cycle in exchange units and their ratio
returns the coupling, with no use of the duration of a cycle.
"""
from __future__ import annotations

import itertools
import math

import numpy as np

from model.executable_model import carrier as cr
from . import law as lw
from .common import clean

G, CYCLE, QUERIES, TRIALS, SEED = 64.0, 1 / 16, 33, 10 ** 6, 20261118
SETTINGS = (0.0, math.pi / 2, math.pi, 3 * math.pi / 2)
M = 1 + lw.N_LAB
SEAMS = cr.CELL_SEAMS
E0, OTHERS = SEAMS[0], SEAMS[1:]
H_EXCHANGE, H_CONTRAST = lw.local_term_parts()


def configs(loads):
    out = []
    for combo in itertools.combinations_with_replacement(range(12), loads):
        occ = [0] * 12
        for v in combo:
            occ[v] += 1
        out.append(tuple(occ))
    return out


def build(loads, kappa, without_e0=False):
    """The complete law with the memory of `e0` kept; `without_e0`: the carrier with `e0` taken out."""
    cf = configs(loads)
    ci = {c: i for i, c in enumerate(cf)}
    m_states = 1 if without_e0 else M
    H = np.zeros((len(cf) * m_states,) * 2)
    idx = lambda c, m: ci[c] * m_states + m
    a0, b0 = E0
    for c in cf:
        diag = kappa * sum(x * x for x in c)             # degree term of one cell: every port has five seams
        for a, b in OTHERS:
            diag += kappa * (c[a] - c[b]) ** 2
            if c[a] != c[b]:
                diag += 0.5
                s = list(c)
                s[a], s[b] = s[b], s[a]
                for m in range(m_states):
                    H[idx(tuple(s), m), idx(c, m)] -= 0.5
        for m in range(m_states):
            H[idx(c, m), idx(c, m)] += diag
        if without_e0:
            continue
        src = lw.SRC.index((c[a0], c[b0]))
        for m in range(M):
            col = src if m == 0 else lw.N_SRC + (m - 1) * lw.N_SRC + src
            for entry in np.flatnonzero(np.abs(H_EXCHANGE[:, col]) + np.abs(H_CONTRAST[:, col]) > 1e-14):
                m2, src2 = (0, lw.SRC[entry]) if entry < lw.N_SRC else (1 + (entry - lw.N_SRC) // lw.N_SRC,
                                                                     lw.SRC[(entry - lw.N_SRC) % lw.N_SRC])
                s = list(c)
                s[a0], s[b0] = src2
                H[idx(tuple(s), m2), idx(c, m)] += H_EXCHANGE[entry, col] + kappa * H_CONTRAST[entry, col]
    if not np.allclose(H, H.T):
        raise ValueError('collar Hamiltonian is not symmetric')
    return cf, ci, H


def beside_a_record(cf, H, kappa, free):
    """One load beside a seam that holds a record of size one: the levels below a load elsewhere.

    A load elsewhere starts at `free + kappa`: the lowest level of a free load plus the level of the record with
    empty ports. The levels below it fall into two groups, split at their largest gap. In the lower group the
    record fits the load on its port and the seam charges nothing; in the upper group it does not fit. The weight
    is the mean probability, over the levels of a group, that the load sits on a port of the seam.
    """
    one = np.array([m > 0 and abs(lw.LABELS[m - 1][1]) == 1 for _ in cf for m in range(M)])
    on_port = np.array([float(c[E0[0]] + c[E0[1]]) for c in cf for _ in range(M)])[one]
    w, v = np.linalg.eigh(H[np.ix_(one, one)])
    start = free + kappa
    below = np.flatnonzero(w < start - 1e-9)
    out = {'a_load_elsewhere_starts_at': float(start), 'levels_below': int(len(below))}
    if len(below) < 2:
        return out
    cut = int(np.argmax(np.diff(w[below]))) + 1
    weight = lambda group: float(np.mean((v[:, group] ** 2).T @ on_port))
    out.update({'levels_where_the_record_fits': [float(x) for x in w[below[:cut]]],
                'levels_where_the_record_does_not_fit': [float(x) for x in w[below[cut:]]],
                'weight_on_a_port_of_the_seam': [weight(below[:cut]), weight(below[cut:])],
                'lowest_level_of_a_load_elsewhere': float(w[len(below)]),
                'weight_on_a_port_of_the_seam_there': weight(np.array([len(below)]))})
    return out


def probe(kappas=(8.0, 0.0)):
    out = {}
    for kappa in kappas:
        lines = {}
        for loads in (0, 1, 2):
            cf, ci, H = build(loads, kappa)
            blank = np.array([m == 0 for _ in cf for m in range(M)])
            entry = {'states': int(H.shape[0]),
                   'largest_element_between_blank_and_not_blank': float(np.abs(H[np.ix_(blank, ~blank)]).max()),
                   'lowest_level_blank': float(np.linalg.eigvalsh(H[np.ix_(blank, blank)])[0]),
                   'lowest_level_not_blank': float(np.linalg.eigvalsh(H[np.ix_(~blank, ~blank)])[0])}
            if loads == 0:
                w = np.linalg.eigvalsh(H[np.ix_(~blank, ~blank)])
                entry['levels_with_no_load'] = sorted({round(float(x), 9) + 0.0 for x in w})
            if loads == 1:
                _, _, Href = build(1, kappa, without_e0=True)
                ref = np.linalg.eigvalsh(Href)
                levels = np.linalg.eigvalsh(H)
                entry['load_moves_on_the_carrier_without_the_seam'] = bool(
                    all(np.abs(levels - (x + kappa)).min() < 1e-9 for x in ref))
                entry['beside_a_record_of_size_one'] = beside_a_record(cf, H, kappa, entry['lowest_level_blank'])
            lines[str(loads)] = entry
        out[str(kappa)] = lines
    return out


def reader(loads, kappa, record_sign, start_config=None):
    """Fuel probabilities `p[setting, q]` for the record `(|c=0,+1> + sign |c=0,-1>)/sqrt2` of `e0`."""
    cf, ci, H = build(loads, kappa)
    n = H.shape[0]
    T = np.zeros((n, n))
    for c in cf:
        for lab, amp in (((0, 1), 1.0), ((0, -1), float(record_sign))):
            T[ci[c] * M + 1 + lw.LABELS.index(lab), ci[c] * M] = amp / math.sqrt(2)
    Z0 = np.zeros((n, n))
    wp, Vp = np.linalg.eigh(np.block([[H, G * T], [G * T.T, H]]))       # blocks: fuel absent, fuel present
    wf, Vf = np.linalg.eigh(np.block([[H, Z0], [Z0, H]]))
    pulse = (Vp * np.exp(-1j * math.pi / (4 * G) * wp)) @ Vp.conj().T
    step = (Vf * np.exp(-1j * CYCLE * wf)) @ Vf.conj().T
    c0 = start_config if start_config is not None else cf[0]
    psi = np.zeros(2 * n, complex)
    psi[n + ci[c0] * M] = 1
    psi = pulse @ psi
    p = np.zeros((4, QUERIES))
    for q in range(QUERIES):
        for s, theta in enumerate(SETTINGS):
            x = psi.copy()
            x[n:] *= np.exp(-1j * theta)
            p[s, q] = float((np.abs((pulse @ x)[n:]) ** 2).sum())
        psi = step @ psi
    return p


def interpret(p):
    z = (p[0] - p[2] - 1j * (p[1] - p[3])) / 4
    phase = np.unwrap(np.angle(z / z[0]))
    return z, -float(np.polyfit(np.arange(len(phase)), phase, 1)[0])


def reading(loads, kappa, rng=None, start_config=None):
    rates, contrast = {}, {}
    for name, sign in (('symmetric', +1), ('antisymmetric', -1)):
        p = reader(loads, kappa, sign, start_config)
        if rng is not None:
            p = rng.binomial(TRIALS, np.clip(p, 0, 1)) / TRIALS
        z, rates[name] = interpret(p)
        contrast[name] = float(np.abs(z).max() - np.abs(z).min())
    cycle = rates['antisymmetric'] - rates['symmetric']
    return {'rate_symmetric': rates['symmetric'], 'rate_antisymmetric': rates['antisymmetric'],
            'cycle_in_exchange_units': cycle,
            'coupling_read': rates['symmetric'] / cycle if abs(cycle) > 1e-12 else None,
            'largest_change_of_contrast': max(contrast.values())}


def run(full=False):
    out = {'probe': probe(),
           'reading': {'probability_law': reading(0, 8.0),
                       'sampled': reading(0, 8.0, rng=np.random.default_rng(SEED)),
                       'without_agreement_energy': reading(0, 0.0)},
           'cycle': CYCLE, 'trials_per_setting': TRIALS, 'queries': QUERIES}
    if full:
        cf = configs(1)
        far = next(c for c in cf if c[E0[0]] == 0 and c[E0[1]] == 0)
        near = next(c for c in cf if c[E0[0]] == 1)
        out['reading']['one_load_off_the_seam'] = reading(1, 8.0, start_config=far)
        out['reading']['one_load_on_a_port_of_the_seam'] = reading(1, 8.0, start_config=near)
    return clean(out)


if __name__ == '__main__':
    import json
    print(json.dumps(run(True), indent=1))
