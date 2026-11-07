"""Exact weighted binary repair on a declared finite graph.

Endpoint translation and raw C/E codes are imported from the native source.
The connected weight assignment is an extension choice, not an admitted
weighted federation program or a physical capacity law.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations
from math import prod

from ..federated_observer import source as geometry
from ..local_repair_family import local_bias as local
from ..native_program import interpreter
from ..observer import native


def integer(value, lower=0):
    if type(value) is not int or value < lower:
        raise ValueError('exact integer in the stated domain required')
    return value


@dataclass(frozen=True)
class Family:
    weights: tuple[int, ...]
    seams: tuple[tuple[int, int], ...]
    seam_weights: tuple[int, ...]

    def __post_init__(self):
        if (type(self.weights) is not tuple or len(self.weights) < 2
                or any(type(w) is not int or w <= 0 for w in self.weights)
                or type(self.seams) is not tuple or not self.seams
                or any(type(e) is not tuple or len(e) != 2
                       or any(type(i) is not int for i in e)
                       or not 0 <= e[0] < e[1] < len(self.weights) for e in self.seams)
                or tuple(sorted(set(self.seams))) != self.seams
                or type(self.seam_weights) is not tuple
                or len(self.seam_weights) != len(self.seams)
                or any(type(w) is not int or w <= 0 for w in self.seam_weights)):
            raise ValueError('positive vertex/seam weights and a canonical finite graph required')
        reached = {0}
        while True:
            following = reached | {v for e in self.seams if reached.intersection(e) for v in e}
            if following == reached:
                break
            reached = following
        if len(reached) != len(self.weights):
            raise ValueError('a connected declared graph required')

    def selection(self, edge):
        return F(self.seam_weights[self.seams.index(edge)], sum(self.seam_weights))

    def coin(self, edge, coin):
        if edge not in self.seams or type(coin) is not int or coin not in (0, 1):
            raise ValueError('an existing ordered seam and binary coin required')
        a, b = edge
        return F(self.weights[a if coin else b], self.weights[a]+self.weights[b])


def local_family(parameters=None):
    if parameters is None:
        parameters = local.Parameters(port_weights=(25, 25)+(8,)*10)
    if type(parameters) is not local.Parameters:
        raise ValueError('existing local Parameters required')
    return Family(parameters.port_weights, tuple(local.SEAMS), parameters.seam_weights)


def federation_family(level, port_word=(25, 25)+(8,)*10):
    """Choose a repeated weight word on the literal retained federation graph."""
    parameters = local.Parameters(port_weights=port_word)
    graph = geometry.federation(level)
    return Family(parameters.port_weights*graph.cells, graph.seams, (1,)*len(graph.seams))


def sector_masks(ports, occupied):
    integer(ports); integer(occupied)
    if occupied > ports:
        raise ValueError('occupation exceeds the carrier')
    return tuple(sum(1 << i for i in chosen) for chosen in combinations(range(ports), occupied))


def image(mask, edge, coin, ports):
    integer(mask); integer(ports, 2)
    if (type(edge) is not tuple or len(edge) != 2
            or any(type(i) is not int for i in edge)
            or not 0 <= edge[0] < edge[1] < ports
            or type(coin) is not int or coin not in (0, 1)):
        raise ValueError('ordered endpoints and a binary native coin required')
    if mask.bit_length() > ports:
        raise ValueError('binary source mask exceeds the carrier')
    a, b = edge
    endpoints = ((mask >> a) & 1, (mask >> b) & 1)
    output, (_, current) = native.native(endpoints, (0, 1), coin)
    target = (mask & ~((1 << a) | (1 << b))) | (output[0] << a) | (output[1] << b)
    return target, current, endpoints[0]-endpoints[1]


def apply_branch(family, density, edge, law, coin, label, weight, *, hold=False):
    """One sparse source/reference operator block; complex-linearity is implicit."""
    if not weight or label is None:
        return {}
    mapping = {}
    for i in {key[k] for key in density for k in (0, 2)}:
        target, current, contrast = image(i, edge, coin, len(family.weights))
        if hold:
            if label == 0:
                mapping[i] = i
        elif label == (current if law == 'C' else contrast):
            mapping[i] = target
    return {(mapping[i], r, mapping[j], s): weight*v
            for (i, r, j, s), v in density.items() if i in mapping and j in mapping and v}


def boundary(family, window):
    if (type(window) is not tuple or tuple(sorted(set(window))) != window
            or any(type(i) is not int or not 0 <= i < len(family.weights) for i in window)):
        raise ValueError('a canonical retained vertex subset required')
    inside = set(window)
    return tuple(e for e in family.seams if len(inside.intersection(e)) == 1)


def attempt(family, density, *, law='mixture', c_weight=F(1, 2), repair=F(1), cut=None):
    """Complete fixed-selector instrument, including hold and all raw abort slots.

    The optional comparison control replaces crossing seams by the identity
    at raw code zero while retaining the seam, law, coin and all zero slots.
    It is not an alternative admitted native program.
    """
    if type(family) is not Family or type(density) is not dict:
        raise ValueError('a declared family and finite sparse operator required')
    if law not in ('C', 'E', 'mixture') or type(c_weight) is not F or not 0 <= c_weight <= 1:
        raise ValueError('an exact C/E mixture required')
    if type(repair) is not F or not 0 <= repair <= 1:
        raise ValueError('an exact attempt probability required')
    for key, value in density.items():
        if (type(key) is not tuple or len(key) != 4
                or any(type(i) is not int or i < 0 for i in key)
                or max(key[0], key[2]).bit_length() > len(family.weights)
                or type(value) not in (int, F)):
            raise ValueError('bounded masks, reference indices and exact real coefficients required')
    crossing = () if cut is None else boundary(family, cut)
    output = [(interpreter.Record('hold'), {k: (1-repair)*v for k, v in density.items() if (1-repair)*v})]
    mix = c_weight if law == 'mixture' else F(law == 'C')
    for edge in family.seams:
        for kind, fraction in (('C', mix), ('E', 1-mix)):
            for coin in (0, 1):
                decoding = interpreter.circuit.record_values(kind, coin)
                weight = repair*family.selection(edge)*fraction*family.coin(edge, coin)
                for raw in range(1 << interpreter.circuit.RECORD_BITS[kind]):
                    label = decoding.get(raw)
                    current = (None if label is None else label if kind == 'C'
                               else (label-(2*coin-1)*(abs(label) % 2))//2)
                    record = interpreter.Record('repair', 'uniform', edge, kind, coin, raw, label, current)
                    output.append((record, apply_branch(family, density, edge, kind, coin, label,
                                                       weight, hold=edge in crossing)))
    return tuple(output)


def population_column(family, mask):
    result = defaultdict(F)
    for edge in family.seams:
        for coin in (0, 1):
            target, _, _ = image(mask, edge, coin, len(family.weights))
            result[target] += family.selection(edge)*family.coin(edge, coin)
    return dict(result)


def modular_ratio(family, edge, current):
    """Exact reference probability ratio across one recorded native current."""
    if edge not in family.seams or type(current) is not int or current not in (-1, 0, 1):
        raise ValueError('an existing seam and native binary current required')
    a, b = edge
    return F(family.weights[b], family.weights[a])**current


def elementary(weights, degree):
    """Coefficient of z**degree in product_i (1+w_i*z), without state enumeration."""
    if type(weights) is not tuple or any(type(w) is not int or w <= 0 for w in weights):
        raise ValueError('finite positive integer weights required')
    if type(degree) is not int:
        raise ValueError('integer coefficient degree required')
    if not 0 <= degree <= len(weights):
        return 0
    values = [1]+[0]*degree
    for i, w in enumerate(weights):
        for k in range(min(i+1, degree), 0, -1):
            values[k] += w*values[k-1]
    return values[degree]


def canonical(weights, occupied):
    z = elementary(weights, occupied)
    if not z:
        raise ValueError('nonempty fixed-number carrier required')
    return {mask: F(prod(weights[i] for i in range(len(weights)) if mask >> i & 1), z)
            for mask in sector_masks(len(weights), occupied)}


def window_law(weights, occupied, window):
    """Full weighted number mixture under the literal occupation partial trace."""
    # The graph is irrelevant to a static diagonal state; validate the subset directly.
    if (type(window) is not tuple or tuple(sorted(set(window))) != window
            or any(type(i) is not int or not 0 <= i < len(weights) for i in window)):
        raise ValueError('a canonical retained vertex subset required')
    z = elementary(weights, occupied)
    if not z:
        raise ValueError('valid total occupation required')
    inside = tuple(weights[i] for i in window)
    outside = tuple(w for i, w in enumerate(weights) if i not in set(window))
    return tuple(F(elementary(inside, k)*elementary(outside, occupied-k), z)
                 for k in range(len(window)+1))


def restricted_state(weights, occupied, window):
    law = window_law(weights, occupied, window)
    retained = tuple(weights[i] for i in window)
    return {mask: weight*p for k, weight in enumerate(law) if weight
            for mask, p in canonical(retained, k).items()}


def grand_canonical(weights, fugacity):
    if type(fugacity) is not F or fugacity <= 0:
        raise ValueError('declared positive rational fugacity required')
    elementary(weights, 0)
    z = prod(1+fugacity*w for w in weights)
    return tuple(fugacity**k*elementary(weights, k)/z for k in range(len(weights)+1))


def cut_bounds(family, window, *, attempts, internal_repairs):
    integer(attempts); integer(internal_repairs)
    crossing = boundary(family, window)
    inside = set(window)
    internal = tuple(e for e in family.seams if set(e) <= inside)
    b = sum((family.selection(e) for e in crossing), F())
    a = sum((family.selection(e) for e in internal), F())
    if not a:
        raise ValueError('at least one internal seam required for the local repair clock')
    return {'internal_probability': a, 'boundary_probability': b,
            'global_attempt_error_upper': 1-(1-b)**attempts,
            'no_boundary_through_kth_internal': (a/(a+b))**internal_repairs,
            'internal_clock_error_upper': 1-(a/(a+b))**internal_repairs}


def local_generator(first, second):
    """Detailed-balance conjugate I-S on |00>,|10>,|01>,|11>."""
    import sympy as s
    integer(first, 1); integer(second, 1)
    block = s.zeros(4)
    block[1, 1], block[2, 2] = s.Rational(second, first+second), s.Rational(first, first+second)
    block[1, 2] = block[2, 1] = -s.sqrt(first*second)/(first+second)
    return s.ImmutableMatrix(block)
