"""Recover the continuing product clock from complete endpoint-coin counts.

The decoder never receives the apparatus's port weights, source loads or
reference.  Its fixed training tree determines the estimate; the other nineteen
seams only test compatibility.  Conditional on any past quantum/classical
history, the next placement coin has probability theta_a/(theta_a+theta_b).
Consequently a fixed count budget has a simultaneous Hoeffding guarantee even
on one continuing, disturbed carrier with an arbitrary untouched reference.

This is a statistical guarantee for the declared coin law, not authentication
of a seeded sampler.  A stopped protocol must retain the identical training
tape and terminal carrier in both comparisons.  The complete error bound is
not conditioned on admission.  Subsequent source compilation and preparation
have their own errors; this module does not implement either service.
"""
from __future__ import annotations

from collections import deque
from fractions import Fraction as F
import hashlib
from math import gcd, lcm
from pathlib import Path

from ..native_modular_reference import native_census
from ..observer import record_recovery

ROOT = Path(__file__).resolve().parents[2]
TRIALS_PER_SEAM = 1 << 20
TREE = native_census.TREE
SEAMS = record_recovery.seams()


def _fraction(value):
    if type(value) is int:
        return F(value)
    if type(value) is not F:
        raise ValueError('an original exact integer or canonical Fraction required')
    n, d = value.numerator, value.denominator
    if (type(n) is not int or type(d) is not int or d <= 0
            or gcd(n, d) != 1):
        raise ValueError('canonical original rational slots required')
    return F(n, d)


def _parameter(value):
    result = F(1) if value is None else _fraction(value)
    if abs(result) > 1:
        raise ValueError('the fixed inference contract requires |parameter| <= 1')
    return result


def _design():
    seams = tuple(record_recovery.seams())
    tree = tuple(native_census.TREE)
    if (len(seams) != 30 or len(set(seams)) != 30
            or any(type(e) is not tuple or len(e) != 2
                   or any(type(i) is not int for i in e) or not 0 <= e[0] < e[1] < 12
                   for e in seams)
            or len(tree) != 11 or len({tuple(sorted(e)) for e in tree}) != 11
            or any(tuple(sorted(e)) not in seams for e in tree)):
        raise ArithmeticError('the native thirty seams and census spanning tree changed')
    adjacency = [[] for _ in range(12)]
    for k, (a, b) in enumerate(tree):
        # Edge coordinate is log(theta_a/theta_b), even if the tree is oriented.
        adjacency[a].append((b, k, -1))
        adjacency[b].append((a, k, 1))
    paths = {0: (0,)*11}
    queue = deque([0])
    while queue:
        a = queue.popleft()
        for b, k, sign in adjacency[a]:
            if b not in paths:
                row = list(paths[a]); row[k] += sign
                paths[b] = tuple(row); queue.append(b)
    if len(paths) != 12:
        raise ArithmeticError('the source census tree is disconnected')
    rows = tuple(paths[i] for i in range(12))
    diameter = max(sum(abs(a-b) for a, b in zip(x, y, strict=True))
                   for x in rows for y in rows)
    if diameter != 5:
        raise ArithmeticError('the fixed census tree diameter changed')
    return seams, tree, rows, diameter


def _counts(counts):
    seams, _, _, _ = _design()
    if (type(counts) is not dict or len(counts) != 30
            or any(type(e) is not tuple or len(e) != 2
                   or any(type(i) is not int for i in e) for e in counts)
            or set(counts) != set(seams)):
        raise ValueError('all thirty canonical ascending native seam keys required')
    result = {}
    for edge in seams:
        value = counts[edge]
        if (type(value) is not tuple or len(value) != 2
                or any(type(x) is not int for x in value)
                or value[0] != TRIALS_PER_SEAM or not 0 <= value[1] <= value[0]):
            raise ValueError('each seam needs the fixed integer trials and complete coin1 count')
        result[edge] = (value[0], value[1])
    return result


def _outward(value, *, lower=False, scale=10**12):
    numerator = value.numerator*scale
    integer = numerator//value.denominator if lower else -((-numerator)//value.denominator)
    return F(integer, scale)


def fixed_budget(parameter=None):
    """Fresh exact, transcript-independent inference and failure allowances.

    e^x >= sum_(j=0)^96 x^j/j! supplies an exact rational upper bound for
    60 exp(-2*N*epsilon^2).  No rounded exponential enters the certificate.
    """
    s = _parameter(parameter)
    _, _, _, diameter = _design()
    epsilon = F(3, 1000)
    exponent = 2*TRIALS_PER_SEAM*epsilon**2
    term = total = F(1)
    for j in range(1, 97):
        term *= exponent/j
        total += term
    failure = _outward(F(60)/total)
    phase = abs(s)*F(25, 4)*diameter*epsilon
    if not failure < F(4, 10**7) or phase > F(3, 32):
        raise ArithmeticError('the fixed inference budget changed')
    return {'trials_per_seam': TRIALS_PER_SEAM, 'seam_count': 30,
            'native_attempts': 30*TRIALS_PER_SEAM,
            'probability_radius': epsilon, 'tree_diameter': diameter,
            'hoeffding_exponent': exponent, 'exponential_lower_sum_terms': 97,
            'statistical_failure_upper': failure,
            'parameter': s, 'clock_halfdiamond_upper_on_confidence_event': phase,
            'complete_inference_halfdiamond_upper': min(F(1), phase+failure),
            'bound_conditioned_on_admission': False}


def _odds(p):
    return p/(1-p)


def _path_interval(coefficients, intervals):
    lo = hi = F(1)
    for coefficient, (a, b) in zip(coefficients, intervals, strict=True):
        if coefficient >= 0:
            lo *= a**coefficient; hi *= b**coefficient
        else:
            lo *= b**coefficient; hi *= a**coefficient
    return lo, hi


def _feasible(intervals):
    """Exact multiplicative difference constraints, with an attained witness.

    For L <= theta_a/theta_b <= U the edges are b->a with multiplier U
    and a->b with multiplier 1/L.  A product-negative cycle is equivalent
    to infeasibility.  Bellman--Ford starts from a super-source of weight1.
    The returned positive rational potentials satisfy ALL input constraints;
    they are only a feasibility witness, never the clock estimate.
    """
    edges = [(b, a, hi) for (a, b), (lo, hi) in intervals.items()]
    edges += [(a, b, 1/lo) for (a, b), (lo, hi) in intervals.items()]
    witness = [F(1)]*12
    for iteration in range(12):
        changed = False
        for a, b, multiplier in edges:
            candidate = witness[a]*multiplier
            if candidate < witness[b]:
                witness[b] = candidate; changed = True
        if not changed:
            if not all(lo <= witness[a]/witness[b] <= hi
                       for (a, b), (lo, hi) in intervals.items()):
                raise ArithmeticError('shared interval feasibility witness failed')
            return tuple(w/witness[0] for w in witness)
        if iteration == 11:
            return None
    raise AssertionError('unreachable')


def exact_identification(edge_probabilities):
    """Identify primitive bounded integer weights from exact tree probabilities.

    This oracle is NOT applied to empirical frequencies.  Exact tree odds fix
    every weight up to a common factor; integer denominators and the primitive
    gcd then fix that factor.  All thirty supplied probabilities are checked.
    """
    seams, tree, paths, _ = _design()
    if (type(edge_probabilities) is not dict
            or any(type(e) is not tuple or len(e) != 2
                   or any(type(i) is not int for i in e) for e in edge_probabilities)
            or set(edge_probabilities) != set(seams)):
        raise ValueError('all thirty exact canonical seam probabilities required')
    probabilities = {e: _fraction(edge_probabilities[e]) for e in seams}
    if any(not F(1, 5) <= p <= F(4, 5) for p in probabilities.values()):
        raise ValueError('endpoint probabilities outside the declared family')
    odds = tuple(_odds(probabilities[tuple(sorted(e))])**(1 if e[0] < e[1] else -1)
                 for e in tree)
    relative = tuple(_path_interval(row, tuple((q, q) for q in odds))[0] for row in paths)
    denominator = lcm(*(q.denominator for q in relative))
    integers = tuple(int(q*denominator) for q in relative)
    divisor = gcd(*integers)
    weights = tuple(n//divisor for n in integers)
    if (max(weights) > 32 or min(weights) < 1 or max(weights) > 4*min(weights)
            or any(F(weights[a], weights[a]+weights[b]) != probabilities[a, b]
                   for a, b in seams)):
        raise ValueError('exact odds have no common primitive bounded-family realization')
    return weights


def infer_counts(counts, *, parameter=None):
    """Fit only the fixed tree; validate every held-out edge without refitting.

    Malformed data raises ValueError. Statistically possible data incompatible
    with the declared model returns an explicit stop, with no usable clock.
    Neither stop nor acceptance is postselection in the complete protocol.
    """
    counts = _counts(counts)
    budget = fixed_budget(parameter)
    seams, tree, paths, _ = _design()
    epsilon = budget['probability_radius']
    rows, probabilities, centers, intervals = [], {}, {}, {}
    for edge in seams:
        n, k = counts[edge]
        empirical = F(k, n)
        lo, hi = max(F(1, 5), empirical-epsilon), min(F(4, 5), empirical+epsilon)
        center = min(F(4, 5), max(F(1, 5), empirical))
        probabilities[edge] = (lo, hi)
        rows.append({'seam': edge, 'trials': n, 'coin1_count': k,
                     'empirical_probability': empirical, 'probability_interval': (lo, hi)})
        centers[edge] = center
        if lo <= hi:
            intervals[edge] = (_odds(lo), _odds(hi))
    result = {'status': 'stop', 'stop_reason': 'probability_domain_incompatible',
              'design': {'tree': tree, 'held_out_seams': tuple(e for e in seams if e not in
                         {tuple(sorted(x)) for x in tree}),
                         'held_out_used_for_fitting': False},
              'observations': tuple(rows), 'budget': budget,
              'probability_intervals': probabilities,
              'relative_weights': None, 'weight_intervals': None,
              'held_out_checks': (), 'shared_interval_feasible': False,
              'shared_feasibility_witness': None,
              'exact_integer_weights_identified': False}
    if len(intervals) != 30:
        return result
    tree_intervals = tuple(intervals[tuple(sorted(e))] if e[0] < e[1]
                           else (1/intervals[tuple(sorted(e))][1],
                                 1/intervals[tuple(sorted(e))][0]) for e in tree)
    tree_centers = tuple(_odds(centers[tuple(sorted(e))])**(1 if e[0] < e[1] else -1)
                         for e in tree)
    relative = tuple(_path_interval(row, tuple((p, p) for p in tree_centers))[0]
                     for row in paths)
    weight_intervals = tuple(_path_interval(row, tree_intervals) for row in paths)
    checks = []
    for edge in result['design']['held_out_seams']:
        a, b = edge
        coefficients = tuple(x-y for x, y in zip(paths[a], paths[b], strict=True))
        lo, hi = _path_interval(coefficients, tree_intervals)
        observed_lo, observed_hi = intervals[edge]
        checks.append({'seam': edge, 'tree_odds_interval': (lo, hi),
                       'observed_odds_interval': (observed_lo, observed_hi),
                       'compatible': max(lo, observed_lo) <= min(hi, observed_hi)})
    result['held_out_checks'] = tuple(checks)
    if not all(row['compatible'] for row in checks):
        result['stop_reason'] = 'held_out_cycle_incompatible'
        return result
    witness = _feasible(intervals)
    if witness is None:
        result['stop_reason'] = 'no_shared_interval_realization'
        return result
    # These rational coordinate enclosures concern the SAME true theta/theta_0,
    # not independently selected references. The fixed Lipschitz theorem
    # controls their modular action without evaluating a logarithm here.
    result.update(status='accept', stop_reason=None, relative_weights=relative,
                  weight_intervals=weight_intervals,
                  shared_interval_feasible=True, shared_feasibility_witness=witness)
    return result


def exact_certificate():
    seams, tree, _, diameter = _design()
    menu = sorted({F(a, a+b) for a in range(1, 33) for b in range(1, 33)
                   if max(a, b) <= 4*min(a, b)})
    gap = min(b-a for a, b in zip(menu, menu[1:]))
    if len(menu) != 489 or gap != F(1, 3843):
        raise ArithmeticError('the exact bounded-family endpoint menu changed')
    return {'schema': 'oph.local_repair_family.inference.v1',
            'native_seams': [list(e) for e in seams],
            'training_tree': [list(e) for e in tree], 'tree_diameter': diameter,
            'held_out_edges': 19, 'trials_per_seam': TRIALS_PER_SEAM,
            'probability_radius': '3/1000', 'endpoint_probability_domain': ['1/5', '4/5'],
            'exact_endpoint_probabilities': len(menu), 'minimum_distinct_probability_gap': str(gap),
            'exact_identification': 'Connected exact edge odds fix theta up to scale; primitive integer normalization fixes scale.',
            'empirical_counts_treated_as_exact_probabilities': False,
            'fixed_parameter_one_phase_bound': '3/32',
            'arbitrary_untouched_reference': True,
            'source_state_independent_coin_probability_required': True,
            'all_training_records_and_stopped_carriers_retained_by_consumer': True,
            'bound_conditioned_on_admission': False,
            'calibration_schedule': 'Fixed requested seams; distinct from subsequent positive full-seam burn-in.',
            'physical_time_or_port_weight_selection_derived': False}


def source_bindings():
    names = sorted(set((*native_census.SOURCE_PATHS,
                        'model/observer/record_recovery.py',
                        'model/local_repair_family/inference.py',
                        'model/local_repair_family/test_inference.py')))
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names}
