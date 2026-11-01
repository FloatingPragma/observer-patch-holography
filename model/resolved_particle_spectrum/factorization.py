"""Exact eta-preserving grouping of the original chronological Strang word.

Only the ideal word is compressed. Complete native source errors are bound
separately; no equality of a rounded native word with this ideal word is used.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations_with_replacement

import numpy as np
from scipy import sparse

from ..cube_eta_protection import cube, spectrum
from ..field_symmetry import symmetry
from ..field_stress.energy import GaussianInteger as GI
from ..certified_reference import interval as iv


STEP = F(39, 2**19)
PAIRS = tuple(combinations_with_replacement(range(8), 2))


def ordered_edges(terms=None):
    terms = tuple(term for group in cube.GEOMETRY.groups for term in group) if terms is None else tuple(terms)
    if len(terms) != 24: raise ValueError('all24 original chronological spin hoppings required')
    edges = []
    for j in range(0, len(terms), 2):
        first, second = terms[j:j+2]
        x, y = first[0]//2, first[1]//2
        if ((second[0]//2, second[1]//2) != (x, y) or
                len({first[0], first[1], second[0], second[1]}) != 4):
            raise ValueError('consecutive spin-complete disjoint edge pair required')
        delta = [b-a for a, b in zip(cube.GEOMETRY.sites[x], cube.GEOMETRY.sites[y], strict=True)]
        if sorted(delta) != [0, 0, 1]: raise ValueError('nearest positive coordinate edge required')
        axis = delta.index(1)
        block = np.zeros((2, 2), complex)
        for p, q, z in (first, second): block[p % 2, q % 2] = 16*z
        if not np.array_equal(block, 1j*np.asarray(symmetry.SIGMA[axis])):
            raise ValueError('literal edge does not equal the original Weyl link')
        edges.append((x, y, axis, symmetry.link_sign(cube.GEOMETRY.sites[x], axis)))
    return tuple(edges)


def symmetric_generator(x, y, sign):
    """Integer generator in the inherited UNnormalized singlet basis."""
    indices = {pair: j for j, pair in enumerate(PAIRS)}
    result = np.zeros((36, 36), dtype=np.int64)
    for col, pair in enumerate(PAIRS):
        coefficient = 1 if pair[0] == pair[1] else 2
        for slot, site in enumerate(pair):
            if site not in (x, y): continue
            changed = list(pair); changed[slot] = y if site == x else x
            target = tuple(sorted(changed))
            result[indices[target], col] += coefficient*sign
    for row, (a, b) in enumerate(PAIRS):
        denominator = 1 if a == b else 2
        if np.any(result[row] % denominator): raise ArithmeticError('nonintegral singlet generator')
        result[row] //= denominator
    return result


@lru_cache(maxsize=1)
def certificate():
    ops = cube.operators(); t, _, v_eta = cube.parent_lift()
    frame = GI(cube.restrict(ops.frame.real, 2), cube.restrict(ops.frame.imag, 2))
    _, singlets, gram, _, _, q = spectrum.singlet_data()
    bonds = {}; index = {site: j for j, site in enumerate(cube.GEOMETRY.sites)}; n = 0
    for x, site in enumerate(cube.GEOMETRY.sites):
        for axis in range(3):
            neighbor = list(site); neighbor[axis] += 1
            y = index.get(tuple(neighbor))
            if y is not None: bonds[x, y] = ops.bonds[n]; n += 1
    reference = np.zeros(len(cube.weights(6)), dtype=np.int64)
    # eta^3 vacuum is an integer vector; its normalization is irrelevant here.
    reference[:] = t@v_eta
    for x, y, _, sign in ordered_edges():
        bond = bonds[x, y]
        for part in (bond.real, bond.imag):
            h2, h6 = cube.restrict(part, 2), cube.restrict(part, 6)
            if (h6@t-t@h2).nnz or np.any(h6@reference):
                raise ArithmeticError('individual complete edge lost the exact eta lift')
        h2 = GI(cube.restrict(bond.real, 2), cube.restrict(bond.imag, 2))
        transformed = frame@h2@frame.adjoint()
        if transformed.imag.nnz: raise ArithmeticError('Hubbard edge is not real')
        target = symmetric_generator(x, y, sign)
        if np.any(transformed.real@singlets-singlets@target):
            raise ArithmeticError('symmetric-square singlet generator differs from actual source edge')
    occupations = lambda number: np.array([sum((int(mask) >> (2*j)) & 3 == 3 for j in range(8))
                                           for mask in cube.weights(number)], dtype=np.int64)
    d2, d6 = (sparse.diags(occupations(n), format='csr', dtype=np.int64) for n in (2, 6))
    if (d6@t-t@d2-2*t).nnz or np.any(d6@reference-3*reference):
        raise ArithmeticError('complete contact layer lost its eta phase')
    if list(gram) != [1 if x == y else 2 for x, y in PAIRS] or int(q@(gram*q)) != 56:
        raise ArithmeticError('inherited singlet chart or local probe normalization changed')
    return {'original_spin_hopping_terms': 24, 'consecutive_complete_edges': 12,
            'all_three_coordinate_axes_retained': True,
            'each_complete_edge_commutes_with_eta': True,
            'contact_lift_identity': 'D_N6 T=T D_N2+2T',
            'one_step_intertwiner': 'S_N6(h) T=exp(-2ih) T S_N2(h)',
            'reference_one_step_phase': 'exp(-3ih)',
            'orthonormal_singlet_chart_dimension': 36,
            'eta3_singlet_subdimension': 35,
            'local_probe_integer_squared_norm': 56,
            'compression_applies_to': 'ideal chronological Strang word; native complete outputs have separate norm bounds',
            'projection_of_the_rounded_native_word_claimed': False}


def add(a, b): return a[0]+b[0], a[1]+b[1]
def multiply(a, b): return a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]
def scale(a, factor): return a[0]*factor, a[1]*factor
def zero(): return iv.Interval.point(0), iv.Interval.point(0)
def one(): return iv.Interval.point(1), iv.Interval.point(0)


@lru_cache(maxsize=1)
def step_intervals():
    """Outward160-bit entry enclosures of the exact36x36 unitary step."""
    certificate()
    sine, cosine = iv.sine_cosine(STEP/32)
    forward = [[one() if i == j else zero() for j in range(8)] for i in range(8)]
    for x, y, _, sign in ordered_edges():
        left, right = forward[x], forward[y]
        phase = iv.Interval.point(0), -sign*sine
        forward[x] = [add(scale(a, cosine), multiply(phase, b)) for a, b in zip(left, right, strict=True)]
        forward[y] = [add(multiply(phase, a), scale(b, cosine)) for a, b in zip(left, right, strict=True)]
    root_two = iv.sqrt_fraction(F(2))
    w = []
    for a, b in PAIRS:
        row = []
        for x, y in PAIRS:
            if a == b and x == y: value = multiply(forward[a][x], forward[a][x])
            elif a == b: value = scale(multiply(forward[a][x], forward[a][y]), root_two)
            elif x == y: value = scale(multiply(forward[a][x], forward[b][x]), root_two)
            else: value = add(multiply(forward[a][x], forward[b][y]), multiply(forward[a][y], forward[b][x]))
            row.append(value)
        w.append(row)
    sine_d, cosine_d = iv.sine_cosine(STEP)
    phase_d = cosine_d, -sine_d
    dw = [[multiply(phase_d, entry) if a == b else entry for entry in row]
          for (a, b), row in zip(PAIRS, w, strict=True)]
    # Reverse chronological hopping is W^T because each individual
    # Hubbard-frame hopping unitary is complex symmetric.
    step = [[sum_complex(multiply(w[k][i], dw[k][j]) for k in range(36))
             for j in range(36)] for i in range(36)]
    return ([[entry[0] for entry in row] for row in step],
            [[entry[1] for entry in row] for row in step])


def sum_complex(entries):
    result = zero()
    for entry in entries: result = add(result, entry)
    return result


def probe_vector():
    q = spectrum.singlet_data()[-1]
    if any(q[j] for j, (x, y) in enumerate(PAIRS) if x != y):
        raise ArithmeticError('integer normalized-chart probe no longer diagonal')
    return [int(x) for x in q]
