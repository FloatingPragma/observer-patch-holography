"""Postprocess the complete local E instrument without changing its branches."""
from fractions import Fraction as F

from model.native_fair_local_readout import readout as r


def profile(p, weights):
    r.authenticate(p)
    if (type(weights) is not tuple or len(weights) != len(p['ports'])
            or any(type(w) is not F for w in weights)):
        raise ValueError('one exact rational weight per retained port required')
    return dict(zip(p['ports'], weights, strict=True))


def current_increments(p, events, attempts, weights, *, selector='uniform'):
    """Delta of the weighted occupation, including every boundary current.

    Completeness of incident events is a caller premise at this algebraic
    interface. Use query() to bind it to the original continuation kernel.
    A zero raw word has no normalized current or particle trajectory.
    """
    phi = profile(p, weights)
    r.integer(attempts)
    if type(events) is not tuple or selector not in ('uniform', 'fixed'):
        raise ValueError('literal chronological events and selector required')
    previous = 0
    increments = []
    for event in events:
        if type(event) is not tuple or len(event) != 2:
            raise ValueError('tick/record tuple required')
        tick, rec = event
        r.integer(tick, previous + 1, attempts)
        previous = tick
        if not r.valid_record(rec, p, selector):
            return None
        a, b = rec.seam
        increments.append((tick, (phi.get(b, F()) - phi.get(a, F())) * rec.current))
    return tuple(increments)


def observe(p, events, attempts, *, selector='uniform', population=None, total=None):
    """One complete history label; erasures and zero words remain explicit."""
    decoded = r.decode(p, events, attempts, selector=selector,
                       population=population, total=total)
    result = {'decoded': decoded, 'snapshots': (), 'attempts': attempts,
              'current_events': events, 'absolute_history_reconstructed': False}
    if decoded['status'] != 'DECODED':
        return result
    initial = decoded['initial_mask']
    state = {port: (initial >> i) & 1 for i, port in enumerate(p['ports'])}
    snapshots = [(0, initial)]
    for tick, rec in events:
        a, b = rec.seam
        if a in state:
            state[a] -= rec.current
        if b in state:
            state[b] += rec.current
        if any(v not in (0, 1) for v in state.values()):
            raise ArithmeticError('authenticated decoded history left binary occupation')
        snapshots.append((tick, sum(state[v] << i for i, v in enumerate(p['ports']))))
    final = snapshots[-1][1]
    if final != decoded['terminal_mask']:
        raise ArithmeticError('terminal history differs from the original decoder')
    if snapshots[-1][0] != attempts:
        snapshots.append((attempts, final))
    result.update(snapshots=tuple(snapshots), absolute_history_reconstructed=True)
    return result


def ward(p, observed, weights):
    """Verify the discrete charge balance on every reconstructed prefix.

    The increments are observable even without an absolute occupation census.
    No cut current is called stress, energy, or a relativistic gauge current.
    """
    phi = profile(p, weights)
    events = observed['current_events']
    selector = events[0][1].selector if events else 'uniform'
    increments = current_increments(p, events, observed['attempts'], weights,
                                    selector=selector)
    if observed['decoded']['status'] == 'ZERO_RAW_WORD':
        return {'status': 'ZERO_RAW_WORD', 'increments': None, 'net_change': None}
    if increments is None:
        raise ValueError('a nonzero readout cannot contain an invalid E slot')
    net = sum((delta for _, delta in increments), F())
    values = ()
    if observed['absolute_history_reconstructed']:
        values = tuple((tick, sum((phi[v] * ((mask >> i) & 1)
                                  for i, v in enumerate(p['ports'])), F()))
                       for tick, mask in observed['snapshots'])
        cursor = 0
        accumulated = F()
        for tick, value in values:
            while cursor < len(increments) and increments[cursor][0] <= tick:
                accumulated += increments[cursor][1]
                cursor += 1
            if value - values[0][1] != accumulated:
                raise ArithmeticError('prefix charge/current Ward identity failed')
    return {'status': observed['decoded']['status'], 'increments': increments,
            'net_change': net, 'weighted_occupations': values}


def query(value, p, slots, density, *, use_global_total=False):
    """Retain the original CP output, FIFO, archive, and every zero outcome."""
    actual = r.raw_query(value, p, slots, density, use_global_total=use_global_total)
    spec = ({'population': 240 * 4**value['carrier']['level'],
             'total': value['carrier']['total']} if use_global_total else {})
    selector = value['commands'][0]['selector']['kind']
    observed = observe(p, actual['observed_events'], len(slots), selector=selector, **spec)
    if observed['decoded'] != actual['decoded']:
        raise ArithmeticError('history postprocessing changed the complete readout')
    return {**actual, 'charge_history': observed}
