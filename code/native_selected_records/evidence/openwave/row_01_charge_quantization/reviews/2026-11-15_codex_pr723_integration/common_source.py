"""Proposed common positive instrument, not a compiled observer apparatus.

All fourteen loads are in [-8,8]. Port 12/13 is the reference seam. An
experiment supplies its ticket list and commit/erase observation setting.
"""
from math import sqrt

PORTS, BOUND = 14, 8


def event(loads, edge, coin):
    if len(loads) != PORTS or any(type(x) is not int or abs(x) > BOUND for x in loads):
        raise ValueError('Outside declared common source domain')
    a, b = edge
    if not (0 <= a < PORTS and 0 <= b < PORTS and a != b) or type(coin) is not int or coin not in (0, 1):
        raise ValueError('Invalid source operation')
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
    # The success vector has support only on j=-1,0,1. Higher currents in
    # the shared load domain are retained as failures, not renormalized away.
    gain = 1 / sqrt(6*len(tickets))
    success = {}
    for (_, _, j), vector in fine.items():
        if j in (-1, 0, 1):
            for q, a in vector.items():
                success[q] = success.get(q, 0) + gain*a
    return {'fine': fine, 'success': success,
            'success_probability': norm2(success),
            'all_failure_probability': sum(norm2(v) for v in fine.values())-norm2(success)}
