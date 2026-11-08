"""Exact bounds for local raw records, disturbance and declared binned clocks.

The Poisson comparison is a mathematical target. The admitted scheduler
continues to make one independently uniform seam selection per attempt.
"""
from collections import defaultdict
from fractions import Fraction as F
from math import comb, factorial

from model.native_local_response import bounds, cone
from model.native_capacity_restriction import restriction


def integer(value, lower=0):
    if type(value) is not int or value < lower:
        raise ValueError('an exact integer in the stated domain required')
    return value


def rational(value, lower=F(0), upper=None):
    if type(value) is not F or value < lower or (upper is not None and value > upper):
        raise ValueError('an exact rational in the stated domain required')
    return value


def aligned_clock(seams, partition):
    """One declared observation partition, with boundaries on the native grid.

    Ordering within each bin is retained; exact global tick addresses are
    traced by this observation map. No physical clock identification follows.
    """
    integer(seams, 1)
    if (type(partition) is not tuple or len(partition) < 2 or
            any(type(t) is not F for t in partition) or partition[0] != 0 or
            any(a >= b for a, b in zip(partition, partition[1:]))):
        raise ValueError('a strictly increasing rational partition starting at zero required')
    ticks = tuple(seams*t for t in partition)
    if any(t.denominator != 1 for t in ticks):
        raise ValueError('clock boundaries must align; unproved rounding is rejected')
    return tuple(int(t) for t in ticks)


def tick_bin(tick, seams, partition):
    ticks = aligned_clock(seams, partition)
    integer(tick, 1)
    if tick > ticks[-1]:
        raise ValueError('the native tick lies beyond the declared horizon')
    return next(i for i, right in enumerate(ticks[1:]) if tick <= right)


def poissonization_bound(seams, internal_seams, partition):
    """Half-diamond allowance after complete local observation.

    The exact one-tick marked selector TV is lambda*(1-exp(-lambda)).
    Its rational upper bound lambda**2 is propagated through the actual
    complete CP instruments, without conditioning on their outcomes.
    """
    integer(internal_seams)
    ticks = aligned_clock(seams, partition)
    if internal_seams > seams:
        raise ValueError('the comparison window cannot contain more than all seams')
    return min(F(1), F(ticks[-1]*internal_seams**2, seams**2))


def radius_errors(root_size, degree, radius, seams, partition):
    integer(root_size, 1); integer(degree, 1); integer(radius)
    attempts = aligned_clock(seams, partition)[-1]
    length = radius+1
    discrete = (root_size*comb(attempts, length)*F(degree, seams)**length
                if attempts >= length else F())
    continuous = root_size*(degree*partition[-1])**length/factorial(length)
    return {'native_cone': min(F(1), discrete),
            'Poisson_cone': min(F(1), continuous)}


def canonical_state_bound(population, occupied, retained, density):
    integer(population, 1); integer(occupied); integer(retained)
    rational(density, upper=F(1))
    if occupied > population or retained > population:
        raise ValueError('occupation and comparison window must fit the source')
    return min(F(1), restriction.bernoulli_bound(population, retained)
               + retained*abs(F(occupied, population)-density))


def local_budget(population, occupied, seams, retained, internal_seams,
                 root_size, degree, radius, partition, density):
    """Native full local instrument to the common finite Poisson window.

    Geometry, endpoint orientation, observed seam addresses and initial
    joint window/reference reductions must match. The canonical specialization
    here requires an independent, common reference, as stated in PROOF.md.
    """
    state = canonical_state_bound(population, occupied, retained, density)
    clocks = poissonization_bound(seams, internal_seams, partition)
    cone_error = radius_errors(root_size, degree, radius, seams, partition)
    return {'initial_half_trace_upper': state,
            'complete_window_instrument_halfdiamond_upper': clocks,
            **cone_error,
            'native_to_common_finite_window_upper': min(F(1), state+clocks+cone_error['native_cone']),
            'native_to_matching_infinite_local_process_upper':
                min(F(1), state+clocks+sum(cone_error.values()))}


def poisson_tail(mean, cutoff):
    rational(mean); integer(cutoff)
    return min(F(1), mean**(cutoff+1)/factorial(cutoff+1))


def exp_minus_interval(value, order=20):
    """Rigorous alternating-series enclosure on [0,1], no float evidence."""
    rational(value, upper=F(1)); integer(order, 1)
    even = 2*order
    upper = sum(((-value)**k/F(factorial(k)) for k in range(even+1)), F())
    lower = upper-value**(even+1)/factorial(even+1)
    return lower, upper


def selector_law_control(seams, internal_seams, partition):
    integer(internal_seams)
    attempts = aligned_clock(seams, partition)[-1]
    if not 0 < internal_seams <= seams:
        raise ValueError('a nonempty window selector required')
    rate = F(internal_seams, seams)
    lo, hi = exp_minus_interval(rate)
    # Per tick: all common singleton mass is lambda*exp(-lambda), and
    # all common empty mass is 1-lambda. All two-or-more marks are retained.
    tv = (rate*(1-hi), rate*(1-lo))
    if not 0 <= tv[0] <= tv[1] <= rate**2:
        raise ArithmeticError('marked selector coupling exceeds its rational allowance')
    return {'one_tick_exact_TV': 'lambda*(1-exp(-lambda))',
            'lambda': rate, 'one_tick_TV_interval': tv,
            'attempts': attempts,
            'complete_word_TV_upper': min(F(1), attempts*rate**2),
            'multiple_marks_retained': True}


def observe_schedules(family, density, schedules, final, observed_edges,
                      seams, partition, *, law='mixture'):
    """Exact finite control of the full binned local CQ observation.

    Input schedules include every global tick. Their rational weights are
    declared selector laws. Every raw slot, zero block and residual quantum
    state is retained, with arbitrary external-reference matrix units.
    """
    attempts = aligned_clock(seams, partition)[-1]
    if (type(schedules) is not dict or not schedules or
            any(type(w) is not F or w < 0 for w in schedules.values()) or
            sum(schedules.values()) != 1):
        raise ValueError('a normalized exact selector distribution required')
    if (type(observed_edges) is not tuple or
            any(type(e) is not tuple or len(e) != 2 or
                any(type(v) is not int for v in e) or not e[0] < e[1]
                for e in observed_edges) or
            observed_edges != tuple(sorted(set(observed_edges))) or
            any(e not in family.seams for e in observed_edges)):
        raise ValueError('canonical native observed seam addresses required')
    out = {}
    for schedule, weight in schedules.items():
        if type(schedule) is not tuple or len(schedule) != attempts:
            raise ValueError('every schedule must retain the common global attempt count')
        observed = tuple(i for i, e in enumerate(schedule) if e in observed_edges)
        blocks = cone.instrument(family, density, schedule, final, observed, law=law)
        for history, block in blocks.items():
            key = tuple((tick_bin(i+1, seams, partition), slot, record)
                        for i, slot, record in history)
            target = out.setdefault(key, defaultdict(F))
            for index, value in block.items():
                target[index] += weight*value
    return {key: {i: v for i, v in block.items() if v} for key, block in out.items()}


def halo_signature(level, cell=0):
    """Explicit orientation-preserving source halo isomorphism.

    Root ports retain their native local number. Each external endpoint is
    named by its incident root port, so every observed raw address has a
    reversible transport between the displayed source windows.
    """
    family = cone.native.federation_family(level, (1,)*12)
    data = bounds.halo_geometry(level, cell)
    root = set(data['root']); window = set(data['window'])
    names = {v: v-12*cell for v in root}
    attachment = []
    for a, b in family.seams:
        if (a in root) != (b in root):
            inside, outside = (a, b) if a in root else (b, a)
            names[outside] = 12+names[inside]
            attachment.append(names[inside])
    internal = tuple(sorted(tuple(sorted((names[a], names[b])))
                            for a, b in family.seams if {a, b} <= window))
    # Reordering an endpoint would also require relabeling the coin/current.
    # Cell zero has all root addresses below its external endpoints.
    orientation_preserved = all(names[a] < names[b]
                                for a, b in family.seams if {a, b} <= window)
    return {'root': data['root'], 'window': data['window'],
            'address_map': tuple(sorted(names.items())),
            'attachment_ports': tuple(sorted(attachment)),
            'internal_edges': internal,
            'endpoint_orientation_preserved': orientation_preserved}
