"""Exact occupation-inclusion operators for the native exclusion transfer.

These operators act on the population transfer Hilbert space. Their
intertwiners do not intertwine every retained quantum outcome branch.
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations
from math import comb

from ..weighted_source_family import family as source


def integer(value, lower=0):
    if type(value) is not int or value < lower:
        raise ValueError('an exact integer in the declared domain required')
    return value


def degree_data(ports, degree):
    integer(ports, 1); integer(degree)
    if degree > ports:
        raise ValueError('occupation exceeds the port count')
    return ports, degree


def coefficients(values, ports, degree=None):
    integer(ports, 1)
    if type(values) is not dict:
        raise ValueError('an exact sparse coefficient dictionary required')
    if degree is not None:
        degree_data(ports, degree)
    result = {}
    for subset, value in values.items():
        if (type(subset) is not tuple or tuple(sorted(set(subset))) != subset
                or any(type(i) is not int or not 0 <= i < ports for i in subset)
                or type(value) not in (int, F)):
            raise ValueError('canonical subsets and exact real coefficients required')
        if degree is None:
            degree = len(subset)
        if len(subset) != degree:
            raise ValueError('one homogeneous subset degree required')
        if value:
            result[subset] = value
    return result, degree


def lower(values, ports, degree=None):
    values, _ = coefficients(values, ports, degree)
    result = defaultdict(int)
    for subset, value in values.items():
        for i in subset:
            result[tuple(j for j in subset if j != i)] += value
    return {s: v for s, v in result.items() if v}


def raise_once(values, ports, degree=None):
    values, _ = coefficients(values, ports, degree)
    result = defaultdict(int)
    for subset, value in values.items():
        for i in range(ports):
            if i not in subset:
                result[tuple(sorted((*subset, i)))] += value
    return {s: v for s, v in result.items() if v}


def inclusion(values, ports, target, degree):
    """Evaluate squarefree fields on a bounded complete target sector."""
    degree_data(ports, target); degree_data(ports, degree)
    values, _ = coefficients(values, ports, degree)
    if target < degree or comb(ports, target) > 10000:
        raise ValueError('ordered degrees and at most10000 evaluated target states required')
    return {state: sum(value for subset, value in values.items() if set(subset) <= set(state))
            for state in combinations(range(ports), target)}


def norm_factor(ports, degree, occupied):
    """Uniform-sector norm of an included primitive unit Euclidean vector."""
    degree_data(ports, degree); degree_data(ports, occupied)
    if not degree <= min(occupied, ports-occupied):
        raise ValueError('the primitive chain must reach the target occupation')
    return F(comb(ports-2*degree, occupied-degree), comb(ports, occupied))


def dimensions(ports, occupied):
    degree_data(ports, occupied)
    return tuple(comb(ports, r)-(comb(ports, r-1) if r else 0)
                 for r in range(min(occupied, ports-occupied)+1))


def contrast(first, second, ports):
    """(n_a-n_b)(n_c-n_d), as four degree-two coefficients."""
    for edge in (first, second):
        if (type(edge) is not tuple or len(edge) != 2
                or any(type(i) is not int for i in edge)
                or not 0 <= edge[0] < edge[1] < ports):
            raise ValueError('two ordered original endpoint pairs required')
    if set(first) & set(second):
        raise ValueError('the two probe seams must be disjoint')
    a, b = first; c, d = second
    return {tuple(sorted((i, j))): sign for i, j, sign in
            ((a,c,1),(a,d,-1),(b,c,-1),(b,d,1))}


def fair_graph(family):
    if type(family) is not source.Family or set(family.weights) != {1} or set(family.seam_weights) != {1}:
        raise ValueError('the native fair vertex and seam weights required')
    return family


def counts_step(family, values, degree):
    """Apply the exact2M-count transfer, without building its sector matrix."""
    family = fair_graph(family)
    values, _ = coefficients(values, len(family.weights), degree)
    neighbours = [[] for _ in family.weights]
    for a, b in family.seams:
        neighbours[a].append(b); neighbours[b].append(a)
    result = defaultdict(int); scale = 2*len(family.seams)
    for subset, value in values.items():
        selected = set(subset); crossing = 0
        for a in subset:
            for b in neighbours[a]:
                if b not in selected:
                    image = tuple(sorted((selected-{a}) | {b}))
                    result[image] += value; crossing += 1
        result[subset] += (scale-crossing)*value
    return {s: v for s, v in result.items() if v}


def dot(first, second):
    return sum(value*second.get(subset, 0) for subset, value in first.items())


def probe_data(level):
    graph = source.federation_family(level, (1,)*12)
    inter = next(e for e in graph.seams if e[0]//12 != e[1]//12)
    cell = inter[0]//12
    # The second seam is local to the first endpoint's cell and disjoint.
    internal = next(e for e in graph.seams if e[0]//12 == e[1]//12 == cell
                    and not set(e) & set(inter))
    return graph, inter, internal, contrast(inter, internal, len(graph.weights))


def autocorrelation(level, steps=6):
    integer(steps)
    if steps > 12:
        raise ValueError('the exact sparse evaluator permits at most12 attempts')
    family, first, second, h = probe_data(level)
    if lower(h, len(family.weights), 2):
        raise ArithmeticError('probe contains a one-body component')
    values = h; denominator = 1; scale = 2*len(family.seams); rows=[]
    for t in range(steps+1):
        rows.append({'attempts':t,'normalized_covariance':str(F(dot(h, values),4*denominator)),
                     'supported_pair_coefficients':len(values)})
        if t != steps:
            values=counts_step(family,values,2);denominator*=scale
            if lower(values,len(family.weights),2):
                raise ArithmeticError('native transfer left the primitive degree-two sector')
    return {'level':level,'ports':len(family.weights),'seams':len(family.seams),
            'probe_seams':[list(first),list(second)],'primitive_degree_two_dimension':comb(len(family.weights),2)-len(family.weights),
            'rows':rows}
