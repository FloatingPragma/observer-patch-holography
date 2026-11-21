"""The law as local events with no common clock.

An event is one seam acting once, `exp(-i tau [(I - SWAP_e)/2 - 2 kappa n_a n_b])`. A program is a sequence of
events. A port observer takes part in the events of its own seams and counts them; that count is its clock. Two
events that share no port commute, so a program fixes the records only through the order of events that each
port observer sees. The position in the sequence is audit; P(no pair) is the largest value over all readings
along the program, one reading after every block of as many events as there are seams.

Programs (fixed seeds):
  rounds            every seam once per round, in a fresh random order each round
  rounds, two rates every round: each seam of the fast half acts twice and every other seam once, fresh order
  independent       every event draws its seam at random, all seams alike
  two rates         every event draws its seam at random, seams of the fast half twice as often

The member with committed contrast: after every event the contrast of the acting seam (which end holds a load)
is committed to a public record; exact density matrices on one cell.
"""
from __future__ import annotations

import numpy as np

from model.executable_model import carrier as cr
from .common import clean
from .sectors import PairSector, SeamEvents

KAPPA = 8.0
SEED = 20261118
KINDS = ('rounds', 'independent', 'two rates')


def fast_half(ev):
    n = ev.sector.n
    return np.array([(a * 2 < n) and (b * 2 < n) for a, b in ev.seams])


def program(kind, ev, rounds, rng):
    m = len(ev.seams)
    if kind == 'rounds':
        return np.concatenate([rng.permutation(m) for _ in range(rounds)])
    if kind == 'independent':
        return rng.integers(0, m, size=rounds * m)
    if kind == 'two rates':
        w = np.where(fast_half(ev), 2.0, 1.0)
        return rng.choice(m, size=rounds * m, p=w / w.sum())
    if kind == 'rounds, two rates':
        base = np.concatenate([np.arange(m), np.flatnonzero(fast_half(ev))])
        return np.concatenate([rng.permutation(base) for _ in range(rounds)])
    raise ValueError(kind)


def execute(ev, seq, tau, every):
    """Apply the events; return the final state, the largest P(no pair) over the readings, the observers' counts."""
    psi = ev.start()
    worst = 0.0
    clock = np.zeros(ev.sector.n, np.int64)
    for i, k in enumerate(seq):
        ev.apply(psi, int(k), tau)
        if (i + 1) % every == 0:
            worst = max(worst, ev.no_pair(psi))
    ends = np.array(ev.seams)[np.asarray(seq, dtype=np.int64)]
    np.add.at(clock, ends[:, 0], 1)
    np.add.at(clock, ends[:, 1], 1)
    return psi, worst, clock


def programs(settings=(('cell', ((1 / 64, 32), (1 / 256, 32))), (0, ((1 / 64, 8),)))):
    lines = {}
    for key, pairs in settings:
        ev = SeamEvents(PairSector(cr.carrier(key)), KAPPA)
        m = len(ev.seams)
        name = 'one cell' if key == 'cell' else cr.NAMES[key]
        for tau, T in pairs:
            rounds = int(round(T / tau))
            for kind in KINDS:
                rng = np.random.default_rng([SEED, m, int(round(1 / tau)), KINDS.index(kind)])
                seq = program(kind, ev, rounds, rng)
                psi, worst, clock = execute(ev, seq, tau, every=m)
                lines[f'{name} | 1/{round(1 / tau)} | {kind}'] = {
                    'events': int(len(seq)), 'readings': int(len(seq) // m), 'largest_no_pair': worst,
                    'norm': float(np.vdot(psi, psi).real), 'smallest_observer_count': int(clock.min()),
                    'largest_observer_count': int(clock.max())}
            if tau == 1 / 64:
                rng = np.random.default_rng([SEED, m, 99])
                seq = program('rounds, two rates', ev, rounds, rng)
                psi, worst, clock = execute(ev, seq, tau, every=m)
                lines[f'{name} | 1/64 | rounds, two rates'] = {
                    'events': int(len(seq)), 'largest_no_pair': worst, 'norm': float(np.vdot(psi, psi).real),
                    'smallest_observer_count': int(clock.min()), 'largest_observer_count': int(clock.max())}
    return lines


def angle_scan(inverse_angles=(256, 1024), total=32):
    """Independent seams on one cell at smaller event angles, the same total angle."""
    ev = SeamEvents(PairSector(cr.cell()), KAPPA)
    out = {}
    for inv in inverse_angles:
        rng = np.random.default_rng([SEED, 30, inv, 5])
        seq = program('independent', ev, total * inv, rng)
        _, worst, _ = execute(ev, seq, 1 / inv, every=30)
        out[f'1/{inv}'] = worst
    return out


def same_order_for_every_observer(key, seed, tau=1 / 64, swaps=20000):
    """Two sequences in which every port observer sees its events in the same order give the same records.

    The first sequence is one round of independent events (one event per seam, seams drawn at random). The
    second sequence exchanges neighbouring events that share no port. The comparison is of the seam records:
    the probabilities of every outcome of the reader."""
    ev = SeamEvents(PairSector(cr.carrier(key)), KAPPA)
    rng = np.random.default_rng(seed)
    seq = list(program('independent', ev, 1, rng))
    other = seq.copy()
    done = 0
    for _ in range(swaps):
        i = int(rng.integers(0, len(other) - 1))
        if not set(ev.seams[other[i]]) & set(ev.seams[other[i + 1]]):
            other[i], other[i + 1] = other[i + 1], other[i]
            done += 1
    orders_equal = all(
        [k for k in seq if p in ev.seams[k]] == [k for k in other if p in ev.seams[k]] for p in range(ev.sector.n))
    p1 = np.abs(execute(ev, seq, tau, 10 ** 9)[0]) ** 2
    p2 = np.abs(execute(ev, other, tau, 10 ** 9)[0]) ** 2
    r1, r2 = p1[ev.sector.seam_states], p2[ev.sector.seam_states]
    return {'exchanges_of_events': done, 'sequences_differ': seq != other, 'every_observer_order_equal': orders_equal,
            'largest_seam_record_difference': float(np.abs(r1 - r2).max()),
            'largest_state_probability_difference': float(np.abs(p1 - p2).max())}


def committed_member(rounds=2048, tau=1 / 64, marks=(1, 8, 64, 512, 2048)):
    """Unread memory against the contrast of the acting seam committed at every event (one cell, rounds)."""
    ps = PairSector(cr.cell())
    ev = SeamEvents(ps, KAPPA)
    P = ps.P
    gates = [ev.matrix(k, tau) for k in range(len(ev.seams))]
    labels = []
    for a, b in ev.seams:
        lab = np.zeros(P, int)
        for i in range(P):
            x, y = int(ps.aa[i]), int(ps.bb[i])
            lab[i] = (1 if (a in (x, y) and b not in (x, y)) else 0) - (1 if (b in (x, y) and a not in (x, y)) else 0)
        labels.append(lab[:, None] == lab[None, :])
    rng = np.random.default_rng([SEED, 30, 77])
    seq = np.concatenate([rng.permutation(len(ev.seams)) for _ in range(rounds)])
    out = {}
    for member in ('unread', 'committed'):
        rho = np.zeros((P, P), complex)
        s = ps.index(*ev.seams[0])
        rho[s, s] = 1
        worst, trace = 0.0, {}
        for i, k in enumerate(seq):
            U = gates[int(k)]
            rho = U @ rho @ U.conj().T
            if member == 'committed':
                rho = rho * labels[int(k)]
            if (i + 1) % len(ev.seams) == 0:
                lost = float(1 - rho.diagonal().real[ps.seam_states].sum())
                worst = max(worst, lost)
                r = (i + 1) // len(ev.seams)
                if r in marks:
                    trace[str(r)] = lost
        out[member] = {'largest_no_pair': worst, 'by_round': trace, 'trace': float(np.trace(rho).real)}
    out['no_pair_for_uniform_pairs_on_one_cell'] = 36 / 66
    return out


def run(full=True):
    out = {'same_order_for_every_observer': {
        'one cell': same_order_for_every_observer('cell', [SEED, 7]),
        '20 cells': same_order_for_every_observer(0, [SEED, 8])}}
    if full:
        out['programs'] = programs()
        out['independent_seams_by_angle'] = angle_scan()
        out['committed_member'] = committed_member()
    return clean(out)


if __name__ == '__main__':
    import json
    print(json.dumps(run(), indent=1))
