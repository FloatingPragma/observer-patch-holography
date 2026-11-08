"""Exact exterior-path Schur response, without a global pair matrix."""
from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations

import sympy as sp

from model.native_contact_response import contact
from model.native_spatial_scaling import geometry


@dataclass(frozen=True)
class Problem:
    neighbours: tuple
    core: tuple


def validate(problem):
    if type(problem) is not Problem or type(problem.neighbours) is not tuple:
        raise ValueError('an exact finite adjacency problem required')
    v = len(problem.neighbours)
    if (type(problem.core) is not tuple or not 2 <= len(problem.core) <= 12
            or tuple(sorted(set(problem.core))) != problem.core
            or any(type(i) is not int or not 0 <= i < v for i in problem.core)):
        raise ValueError('two through twelve distinct canonical core ports required')
    for i, row in enumerate(problem.neighbours):
        if (type(row) is not tuple or len(row) > 6 or tuple(sorted(set(row))) != row
                or any(type(j) is not int or not 0 <= j < v or j == i for j in row)):
            raise ValueError('canonical simple graph with degree at most six required')
    if any(i not in problem.neighbours[j] for i, row in enumerate(problem.neighbours) for j in row):
        raise ValueError('symmetric source adjacency required')
    return problem


@lru_cache(maxsize=12)
def problem(level, cell=0):
    contact.integer(level, 0, 2)
    g = geometry.graph(level); geometry.halo(g, cell)
    return validate(Problem(g.neighbours, tuple(range(12*cell, 12*(cell+1)))))


def sparse(values, ports):
    return contact.source.coefficients(values, ports, 2)[0]


def apply(problem, values, exterior=False):
    """G or QGQ on exact sparse pair coefficients; no scheduler replacement."""
    if type(exterior) is not bool:
        raise ValueError('a literal exterior flag required')
    values = sparse(values, len(problem.neighbours)); core = set(problem.core)
    if exterior and any(set(pair) <= core for pair in values):
        raise ValueError('exterior input cannot contain a core pair')
    result = defaultdict(int)
    for pair, value in values.items():
        occupied = set(pair)
        for a in pair:
            for b in problem.neighbours[a]:
                if b not in occupied:
                    image = tuple(sorted((occupied-{a}) | {b}))
                    result[pair] += value
                    if not exterior or not set(image) <= core:
                        result[image] -= value
    return {pair: value for pair, value in result.items() if value}


def blocks(problem):
    validate(problem)
    pairs = tuple(combinations(problem.core, 2)); where = {p: i for i, p in enumerate(pairs)}
    a = sp.zeros(len(pairs)); columns = []
    for j, pair in enumerate(pairs):
        column = apply(problem, {pair: 1})
        for image in set(column) & set(where):
            a[where[image], j] = column[image]
        columns.append({p: v for p, v in column.items() if p not in where})
    gram = sp.Matrix([[sum(v*columns[j].get(p, 0) for p, v in column.items())
                       for j in range(len(pairs))] for column in columns])
    core = set(problem.core)
    if any(sum(b not in core for b in problem.neighbours[a]) > 1 for a in problem.core):
        raise ValueError('at most one outward seam per core port required')
    boundary_count = [sum(b not in core for a in pair for b in problem.neighbours[a]) for pair in pairs]
    expected = sp.diag(*boundary_count)
    if gram != expected or max(boundary_count, default=0) > 2:
        raise ValueError('distinct outward images and at most one outward seam per core port required')
    if a != a.T:
        raise ArithmeticError('core compression is not symmetric')
    return pairs, sp.ImmutableMatrix(a), tuple(columns), sp.ImmutableMatrix(gram)


def moments(problem, degree=2):
    """B*(24I-D)^k B; only reachable exterior pair configurations are visited."""
    contact.integer(degree, 0, 6)
    pairs, a, initial, gram = blocks(problem)
    columns = initial; matrices = []; support = []
    for k in range(degree+1):
        matrix = sp.Matrix([[sum(v*columns[j].get(p, 0) for p, v in column.items())
                             for j in range(len(pairs))] for column in initial])
        if matrix != matrix.T or (k == 0 and matrix != gram):
            raise ArithmeticError('exterior moments lost exact symmetry or coupling')
        matrices.append(sp.ImmutableMatrix(matrix)); support.append(len(set().union(*(set(c) for c in columns))))
        if k < degree:
            following = []
            for column in columns:
                image = apply(problem, column, exterior=True)
                keys = set(column) | set(image)
                following.append({p: 24*column.get(p, 0)-image.get(p, 0)
                                  for p in keys if 24*column.get(p, 0)-image.get(p, 0)})
            columns = tuple(following)
    return {'pairs': pairs, 'core_action': a, 'coupling_gram': gram,
            'moments': tuple(matrices), 'reachable_exterior_counts': tuple(support)}


def errors(s, degree):
    s = contact.rational(s); contact.integer(degree, 0)
    if s <= 0:
        raise ValueError('strictly positive off-spectrum shift required')
    ratio = F(24)/(24+s); tail = ratio**(degree+1)/s
    return {'shift': s, 'degree': degree, 'geometric_ratio': ratio,
            'self_energy_error_upper': 2*tail,
            'core_resolvent_error_upper': 2*tail/(s*s),
            'normalized_Abel_error_upper': 2*tail/s}


def approximation(problem, s=F(24), degree=2):
    bounds = errors(s, degree); data = moments(problem, degree); n = len(data['pairs'])
    s = bounds['shift']; sigma = sp.zeros(n)
    for k, matrix in enumerate(data['moments']):
        sigma += matrix/sp.Rational((24+s)**(k+1))
    f = data['core_action']+sp.Rational(s)*sp.eye(n)-sigma
    return {**data, **bounds, 'self_energy': sp.ImmutableMatrix(sigma),
            'schur': sp.ImmutableMatrix(f)}


def probe(problem):
    if len(problem.core) != 12:
        raise ValueError('the native twelve-port cell required for this record probe')
    offset = problem.core[0]
    edges = tuple(tuple(offset+i for i in edge) for edge in contact.PROBE)
    if any(b not in problem.neighbours[a] for a, b in edges):
        raise ValueError('the original native contrast seams are required')
    values = contact.source.contrast(*edges, len(problem.neighbours))
    return sp.Matrix([values.get(pair, 0) for pair in combinations(problem.core, 2)])


def response(problem, s=F(24), degree=2):
    data = approximation(problem, s, degree); h = probe(problem)
    inverse = data['schur'].inv(method='DM')
    lower = contact.fraction((h.T*inverse*h)[0]/4)*data['shift']
    upper = min(F(1), lower+data['normalized_Abel_error_upper'])
    if not 0 < lower <= upper <= 1:
        raise ArithmeticError('the connected native Abel enclosure failed')
    return {'shift': data['shift'], 'degree': degree,
            'normalized_Abel_interval': (lower, upper),
            'reachable_exterior_counts': data['reachable_exterior_counts'],
            'self_energy_error_upper': data['self_energy_error_upper'],
            'core_resolvent_error_upper': data['core_resolvent_error_upper'],
            'normalized_Abel_error_upper': data['normalized_Abel_error_upper']}


def certificate():
    rows = []
    isolated = -24*contact.response(-24)
    for level, cell in ((0, 0), (0, 1), (1, 0), (2, 0)):
        p = problem(level, cell); data = moments(p, 2)
        boundary = tuple(a for a in p.core if any(b not in p.core for b in p.neighbours[a]))
        local = contact.data()['two']
        if data['core_action'] != local+data['coupling_gram'] or len(boundary) != 3:
            raise ArithmeticError('the actual routed compression differs from Gcell+boundary occupation')
        result = response(p)
        if cell == 0 and isolated-result['normalized_Abel_interval'][1] <= F(1, 120):
            raise ArithmeticError('the exact connected/isolated response separation failed')
        rows.append({'level': level, 'cell': cell, 'ports': len(p.neighbours),
                     'seams': sum(map(len, p.neighbours))//2, 'boundary_ports': boundary,
                     'core_pairs': 66, 'coupling_rank': data['coupling_gram'].rank(),
                     **result})
    return {'finite_source_controls': rows, 'core_projection_inserted_inside_self_energy': False,
            'isolated_normalized_Abel_value_at_shift24': isolated,
            'cell0_connected_discrepancy_lower': F(1, 120),
            'operator_bound': 24, 'coupling_squared_norm_upper': 2,
            'negative_spectral_parameter_only': True,
            'positive_energy_pole_continuation_certified': False}
