"""Exact contact-mode leakage when an isolated native cell is connected."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from math import lcm

import sympy as sp

from model.native_spatial_scaling import geometry
from . import contact as c


@lru_cache(maxsize=1)
def cubic_basis():
    g = c.data()['two']
    annihilator = c.polynomial(g, c.CUBIC.all_coeffs())
    columns = annihilator.nullspace()
    if len(columns) != 15:
        raise ArithmeticError('the full cubic eigenspace dimension changed')
    integers = []
    for column in columns:
        scale = lcm(*(int(sp.denom(x)) for x in column))
        integers.append(tuple(int(scale*x) for x in column))
    result = sp.ImmutableMatrix.hstack(*(sp.ImmutableMatrix(x) for x in integers))
    if annihilator*result != sp.zeros(66, 15) or c.data()['primitive']*result != result:
        raise ArithmeticError('the cubic basis is not native primitive degree two')
    return result


def leading_minors(matrix):
    """Fraction-free Sylvester test, preserving the exact integer matrix."""
    if (type(matrix) is not tuple or not matrix or any(type(row) is not tuple for row in matrix)
            or any(len(row) != len(matrix) for row in matrix)
            or any(type(x) is not int for row in matrix for x in row)):
        raise ValueError('a nonempty exact square integer matrix is required')
    n = len(matrix)
    if any(matrix[i][j] != matrix[j][i] for i in range(n) for j in range(n)):
        raise ValueError('the Sylvester matrix must be symmetric')
    work = [list(row) for row in matrix]; previous = 1; minors = []
    for k in range(n):
        pivot = work[k][k]; minors.append(pivot)
        if pivot <= 0:
            return tuple(minors)
        for i in range(k+1, n):
            for j in range(i, n):
                numerator = pivot*work[i][j]-work[i][k]*work[k][j]
                value, remainder = divmod(numerator, previous)
                if remainder:
                    raise ArithmeticError('fraction-free principal-minor division failed')
                work[i][j] = work[j][i] = value
        for i in range(k+1, n):
            work[i][k] = work[k][i] = 0
        previous = pivot
    return tuple(minors)


@lru_cache(maxsize=1)
def port_grams():
    n = cubic_basis(); size = n.cols
    grams = []
    for port in range(12):
        rows = [i for i, pair in enumerate(c.data()['pairs']) if port in pair]
        grams.append(tuple(tuple(sum(int(n[k, i])*int(n[k, j]) for k in rows)
                                 for j in range(size)) for i in range(size)))
    norm_matrix = n.T*n
    norm = tuple(tuple(int(norm_matrix[i, j]) for j in range(size)) for i in range(size))
    return tuple(grams), norm


def boundary_matrix(ports, lower=F(1, 100)):
    if (type(ports) is not tuple or len(ports) != 3 or tuple(sorted(set(ports))) != ports
            or any(type(i) is not int or not 0 <= i < 12 for i in ports)):
        raise ValueError('three distinct ordered native boundary ports are required')
    lower = c.rational(lower)
    if lower < 0:
        raise ValueError('a nonnegative candidate leakage bound is required')
    grams, norm = port_grams()
    return tuple(tuple(lower.denominator*sum(grams[p][i][j] for p in ports)
                       -lower.numerator*norm[i][j] for j in range(15)) for i in range(15))


def certificate():
    checks = []
    for ports in combinations(range(12), 3):
        minors = leading_minors(boundary_matrix(ports))
        if len(minors) != 15 or min(minors) <= 0:
            raise ArithmeticError('a native cubic mode violates the uniform leakage bound')
        checks.append(ports)
    rows = []
    for level in (0, 1, 2):
        graph = geometry.graph(level)
        patterns = set()
        for cell in range(graph.cells):
            halo = geometry.halo(graph, cell)
            # The owning geometry validates the three distinct endpoint ports,
            # distinct outside neighbors, and every original internal seam.
            ports = tuple(sorted(a % 12 if a//12 == cell else b % 12
                                 for a, b in graph.seams if (a//12 == cell) != (b//12 == cell)))
            if ports not in checks:
                raise ArithmeticError('actual federation has a different cell boundary')
            patterns.add(ports)
        rows.append({'level': level, 'cells': graph.cells, 'boundary_patterns': sorted(patterns)})
    return {'cubic_subspace_dimension': 15, 'distinct_boundary_triplets': len(checks),
            'positive_principal_minors_per_triplet': 15,
            'outside_action_norm_squared_lower': F(1, 100),
            'outside_action_norm_lower': F(1, 10),
            'evaluated_federations': rows,
            'all_level_attachment_requires_returned_routing_postconditions': True}


def endpoint_lower(time, attempts=None):
    """Lower bound on a fixed isolated-mode comparison after connection.

    Heat evolution is a comparison operator; native attempts are separate.
    """
    time = c.rational(time)
    if time < 0:
        raise ValueError('a nonnegative exchange parameter is required')
    value = max(F(), time/10-288*time*time)
    if attempts is not None:
        c.integer(attempts, 1)
        value = max(F(), value-F(3, 5*attempts))
    return value
