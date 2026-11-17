"""Proposed common positive instrument, not a compiled observer apparatus.

All fourteen loads are in [-8,8]. Port 12/13 is the reference seam. An
experiment supplies its ticket list and commit/erase observation setting.
"""
from math import sqrt

PORTS, BOUND = 14, 8
ALLOWED = frozenset(((0,1),(0,5),(0,7),(0,10),(0,11),(1,5),(1,7),
    (1,8),(1,9),(2,3),(2,4),(2,6),(2,10),(2,11),(3,4),(3,6),
    (3,8),(3,9),(4,5),(4,9),(4,11),(5,9),(5,11),(6,7),(6,8),
    (6,10),(7,8),(7,10),(8,9),(10,11),(12,13)))


def event(loads, edge, coin):
    if len(loads) != PORTS or any(type(x) is not int or abs(x) > BOUND for x in loads):
        raise ValueError('Outside declared common source domain')
    a, b = edge
    if not (0 <= a < PORTS and 0 <= b < PORTS and a != b) or type(coin) is not int or coin not in (0, 1):
        raise ValueError('Invalid source operation')
    if tuple(sorted(edge)) not in ALLOWED:
        raise ValueError('Nonadjacent common-source ports')
    difference = loads[a] - loads[b]
    current = (difference - (2*coin-1)*(abs(difference) % 2)) // 2
    out = list(loads)
    out[a] -= current
    out[b] += current
    return tuple(out), current


def instrument(amplitudes, tickets):
    """V|x> = sum_r K_r|x>|r>, with every ticket/coin/current flag.

    `amplitudes` is a source ket, and no path is sampled or discarded here.
    Fine-record partial injectivity makes this a complete isometry. Complex
    input amplitudes do not alter the positive native branch coefficients.
    """
    if not tickets:
        raise ValueError('Empty scheduler')
    gain = 1 / sqrt(2*len(tickets))
    result = {}
    for q, amplitude in amplitudes.items():
        for ticket, edge in enumerate(tickets):
            for coin in (0, 1):
                out, current = event(q, edge, coin)
                branch = result.setdefault((ticket, coin, current), {})
                branch[out] = branch.get(out, 0) + gain*amplitude
    return result


def norm2(vector):
    return sum(abs(a)**2 for a in vector.values())


def observe(amplitudes, tickets, mode):
    fine = instrument(amplitudes, tickets)
    if mode == 'commit':
        return {'fine': fine, 'outcome_probabilities': {r: norm2(v) for r, v in fine.items()}}
    if mode != 'optical_erase':
        raise ValueError('Unregistered observation mode')
    # Extend the original row-7 Householder analyzer by identity on higher
    # current flags. The anchor is the same (ticket0,coin0,current-1) flag.
    support = [(t, c, j) for t in range(len(tickets))
               for c in (0, 1) for j in (-1, 0, 1)]
    anchor = support[0]
    gain = 1 / sqrt(len(support))
    norm = sqrt(2*(1-gain))
    u = {r: (gain-int(r == anchor))/norm for r in support}
    contraction = {}
    for r, coefficient in u.items():
        for q, amplitude in fine.get(r, {}).items():
            contraction[q] = contraction.get(q, 0) + coefficient*amplitude
    outcomes = {}
    for r in set(fine) | set(support):
        vector = dict(fine.get(r, {}))
        if r in u:
            for q, amplitude in contraction.items():
                vector[q] = vector.get(q, 0)-2*u[r]*amplitude
        outcomes[r] = vector
    success = outcomes[anchor]
    failures = {r: v for r, v in outcomes.items() if r != anchor}
    return {'fine': fine, 'success': success,
            'success_probability': norm2(success),
            'failure_states': failures,
            'all_failure_probability': sum(norm2(v) for v in failures.values())}
