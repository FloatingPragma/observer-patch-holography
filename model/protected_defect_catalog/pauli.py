"""All Pauli-product charges of the actual selected exchange source."""
from fractions import Fraction as F
from itertools import product

import numpy as np

from model.native_ordered_ground import source


def graph_data(sites, edges, fields=()):
    if type(sites) is not int or not 1 <= sites <= 16:
        raise ValueError('one through sixteen sites required')
    normalized = []
    for edge in edges:
        if len(edge) != 2 or any(type(x) is not int or not 0 <= x < sites for x in edge):
            raise ValueError('valid graph endpoints required')
        a, b = sorted(edge)
        if a == b or (a, b) in normalized:
            raise ValueError('simple graph required')
        normalized.append((a, b))
    fields = tuple(fields)
    if len(set(fields)) != len(fields) or any(type(x) is not int or not 0 <= x < sites for x in fields):
        raise ValueError('distinct field sites required')
    return tuple(sorted(normalized)), tuple(sorted(fields))


def components(sites, edges):
    todo = set(range(sites)); out = []
    while todo:
        found = {min(todo)}
        while True:
            enlarged = found | {v for edge in edges if set(edge) & found for v in edge}
            if enlarged == found:
                break
            found = enlarged
        todo -= found; out.append(tuple(sorted(found)))
    return tuple(out)


def charges(sites, edges, fields=()):
    """Complete classification modulo scalar phases, not all conserved operators."""
    edges, fields = graph_data(sites, edges, fields)
    blocks = components(sites, edges)
    generators = []
    for block in blocks:
        mask = sum(1 << x for x in block)
        generators.append((mask, 0))
        if not set(block) & set(fields):
            generators.append((0, mask))
    if len(generators) > 16:
        raise ValueError('explicit charge list exceeds the declared enumeration budget')
    result = []
    for selectors in product((0, 1), repeat=len(generators)):
        x = z = 0
        for selected, (gx, gz) in zip(selectors, generators, strict=True):
            if selected:
                x ^= gx; z ^= gz
        result.append((x, z))
    return tuple(sorted(result))


def commutes_with_terms(x, z, edges, fields):
    # Distinct nonzero Pauli terms are linearly independent. Conjugation by
    # one Pauli word only changes their signs, so cancellation cannot hide
    # a noncommuting term in [P,H].
    for a, b in edges:
        edge = (1 << a) | (1 << b)
        if ((z & edge).bit_count() | (x & edge).bit_count()) & 1:
            return False
    return all(not ((z >> a) & 1) for a in fields)


def exhaustive(sites, edges, fields=()):
    if type(sites) is not int or not 1 <= sites <= 8:
        raise ValueError('exhaustive census supports at most eight sites')
    edges, fields = graph_data(sites, edges, fields)
    return tuple((x, z) for x in range(1 << sites) for z in range(1 << sites)
                 if commutes_with_terms(x, z, edges, fields))


def word_matrix(sites, x, z):
    size = 1 << sites
    if any(type(a) is not int or not 0 <= a < size for a in (x, z)):
        raise ValueError('Pauli words must fit their sites')
    matrix = np.zeros((size, size), dtype=complex)
    phase = (1j)**((x & z).bit_count())
    for column in range(size):
        matrix[column ^ x, column] = phase * (-1)**((z & column).bit_count())
    return matrix


def source_polynomial():
    """Exact integer numerator in the source's Marshall frame, no fitted H."""
    edges, _ = source.cube()
    actual, scale = source.matrix()
    coefficients = {(0, 0): -F(len(edges), 4)}
    for a, b in edges:
        edge = (1 << a) | (1 << b)
        for word, coefficient in (((edge, 0), -F(1, 4)),
                                  ((edge, edge), -F(1, 4)),
                                  ((0, edge), F(1, 4))):
            coefficients[word] = coefficients.get(word, F(0)) + coefficient
    for a in range(8):
        coefficients[1 << a, 0] = -F(1, 8)
    reproduced = np.zeros((256, 256), dtype=np.int64)
    for (x, z), coefficient in coefficients.items():
        scaled = coefficient * scale
        if scaled.denominator != 1:
            raise ArithmeticError('source scale does not exactly represent Pauli coefficients')
        word = word_matrix(8, x, z)
        if np.any(word.imag) or not np.array_equal(word.real, word.real.astype(np.int64)):
            raise ArithmeticError('source polynomial ceased to be real integral')
        reproduced += scaled.numerator * word.real.astype(np.int64)
    if not np.array_equal(actual, reproduced):
        raise ArithmeticError('actual source Hamiltonian differs from exact Pauli polynomial')
    return edges, actual, scale, coefficients


def census():
    edges, actual, scale, polynomial = source_polynomial()
    selected = charges(8, edges, range(8))
    zero_field = charges(8, edges)
    if exhaustive(8, edges, range(8)) != selected or exhaustive(8, edges) != zero_field:
        raise ArithmeticError('graph theorem differs from exhaustive Pauli census')
    for x, z in selected:
        p = word_matrix(8, x, z)
        if not np.array_equal(p @ actual, actual @ p):
            raise ArithmeticError('classified charge does not commute with source')
    if selected != ((0, 0), (255, 0)):
        raise ArithmeticError('selected source charge census changed')
    return {'sites': 8, 'edges': len(edges), 'Pauli_words_exhausted_per_field': 65536,
            'exact_source_matrix_scale': scale,
            'nonzero_Pauli_terms_including_identity': len(polynomial),
            'selected_field_charge_words_modulo_phase': [list(p) for p in selected],
            'zero_field_charge_words_modulo_phase': [list(p) for p in zero_field],
            'selected_field_charge_group': 'Z2 generated by global X',
            'zero_field_charge_group_modulo_phase': 'Z2 x Z2 generated by global X and Z',
            'nonidentity_strictly_local_Pauli_charge_count': 0,
            'all_conserved_operator_algebra_classified': False,
            'all_topological_defects_classified': False,
            'physical_particle_identification_count': 0}
