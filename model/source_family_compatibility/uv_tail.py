"""Volume-independent original CAR energy tails from local source commutators."""
from fractions import Fraction as F
from functools import lru_cache
from math import factorial
import numpy as np

from . import composite


def _detached(value):
    if isinstance(value, F):
        return F(value.numerator, value.denominator)
    if isinstance(value, dict):
        return {key: _detached(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_detached(item) for item in value]
    if isinstance(value, tuple):
        return tuple(_detached(item) for item in value)
    return value


def parameters(c, u, degree):
    if type(c) not in (int, F) or type(u) not in (int, F) or c <= 0 or u < 0:
        raise ValueError('exact positive hopping and nonnegative contact required')
    if type(degree) is not int or not 0 <= degree <= 6:
        raise ValueError('cubic source maximum degree in 0..6 required')
    c, u = F(c), F(u)
    return (u+c*degree)/2, u/2+2*c*degree, (c*degree/2)**2+u*u/4


def ball(radius):
    if type(radius) is not int or radius < 0:
        raise ValueError('nonnegative integer lattice radius required')
    return (4*radius**3+6*radius**2+8*radius+3)//3


def local_derivative_bound(c, u, degree, order):
    if type(order) is not int or order < 0:
        raise ValueError('nonnegative integer derivative order required')
    a, b, _ = parameters(c, u, degree)
    return F(1) if order == 0 else a*b**(order-1)*factorial(order-1)


def even_moment_bound(c, u, degree, order):
    """Signed retarded moment 2*order of any unit lattice smearing."""
    if type(order) is not int or order < 0:
        raise ValueError('nonnegative integer half-moment order required')
    _, _, second = parameters(c, u, degree)
    if order == 0:
        return F(1)
    if order == 1:
        return second
    return 4*ball(2*order)*local_derivative_bound(c, u, degree, order)**2


def tail_budget(c, u, degree, cutoff, max_order=64):
    if type(cutoff) not in (int, F) or cutoff <= 0:
        raise ValueError('exact positive energy cutoff required')
    if type(max_order) is not int or not 1 <= max_order <= 256:
        raise ValueError('finite moment search in 1..256 required')
    e = F(cutoff)
    rows = [(even_moment_bound(c, u, degree, n)/e**(2*n), n)
            for n in range(1, max_order+1)]
    value, order = min(rows)
    return {'energy_cutoff': e, 'moment_order': 2*order,
            'original_total_charge_tail_upper': min(F(1), value),
            'finite_search_orders': max_order,
            'conditional_probability_used': False,
            'uniform_in_volume_and_unit_smearing': True}


def exponential_tail(c, u, degree, cutoff):
    if type(cutoff) not in (int, F) or cutoff <= 0:
        raise ValueError('exact positive energy cutoff required')
    a, b, _ = parameters(c, u, degree); e = F(cutoff)
    if not b:
        return F(0)
    m = int(e//(2*b))
    if m < 1:
        return F(1)
    return min(F(1), 4*ball(2*m+2)*(a/e)**2/F(4)**m)


def exponential_moment(c, u, degree, slope):
    """Bound integral exp(slope*E) over both original positive charge measures."""
    if type(slope) not in (int, F) or slope < 0:
        raise ValueError('exact nonnegative exponential slope required')
    a, b, _ = parameters(c, u, degree); x = F(slope)*b
    if x >= 1:
        raise ValueError('slope must lie inside the proved source analytic strip')
    return F(1) if not b else 1+23*a*F(slope)/(1-x)**2


def _positive(matrix):
    """Exact small rational positive-semidefinite elimination."""
    a = [[F(x) for x in row] for row in matrix]; size = len(a)
    if any(len(row) != size for row in a) or any(a[i][j] != a[j][i]
            for i in range(size) for j in range(size)):
        raise ValueError('symmetric square rational Gram required')
    for k in range(size):
        pivot = a[k][k]
        if pivot < 0:
            raise ArithmeticError('source moment exceeds its exact bound')
        if not pivot:
            if any(a[k][j] for j in range(k, size)):
                raise ArithmeticError('zero Gram pivot has a nonzero row')
            continue
        for i in range(k+1, size):
            q = a[i][k]/pivot
            for j in range(i, size):
                a[i][j] -= q*a[k][j]; a[j][i] = a[i][j]


def _sum_entries(matrix):
    return sum((int(x) for x in matrix.flat), 0)


@lru_cache(maxsize=2)
def _finite_control(sites):
    if type(sites) is not int or sites not in (2, 3):
        raise ValueError('complete two- or three-site source control required')
    c, u, e0 = ((F(3, 2), F(4), F(-9, 2)) if sites == 2
                else (F(2), F(3, 7), F(-45, 14)))
    data = composite.source((sites, 1, 1), c, u)
    numerator, denominator, rank = composite.ground(data, e0)
    rho = numerator.astype(object); h = data['H'].astype(object)
    scale = data['scale']; modes = 2*sites; degree = 1 if sites == 2 else 2
    derivatives = [a.astype(object) for a in data['a']]
    rows = []
    for n in range(1, 7):
        derivatives = [d@h-h@d for d in derivatives]
        local = local_derivative_bound(c, u, degree, n)
        for d in derivatives:
            row_norm = max(sum(abs(int(x)) for x in row) for row in d)
            col_norm = max(sum(abs(int(x)) for x in row) for row in d.T)
            if F(max(row_norm, col_norm), scale**n) > local:
                raise ArithmeticError('actual local derivative exceeds its source bound')
        weighted = [d@rho+rho@d for d in derivatives]
        gram = [[F(_sum_entries(derivatives[i]*weighted[j]), denominator*scale**(2*n))
                 for j in range(modes)] for i in range(modes)]
        _positive(gram)
        bound = even_moment_bound(c, u, degree, n)
        _positive([[bound*(i == j)-gram[i][j] for j in range(modes)]
                   for i in range(modes)])
        if n == 1:
            one = data['one'].astype(object)
            expected = one@one
            if any(gram[i][j] != F(int(expected[i, j]), scale**2)+(u*u/4)*(i == j)
                   for i in range(modes) for j in range(modes)):
                raise ArithmeticError('actual source second-moment identity differs')
        rows.append({'moment_order': 2*n, 'uniform_operator_bound': bound,
                     'actual_original_moment_trace': sum(gram[i][i] for i in range(modes)),
                     'local_operator_derivative_bound_checked': True,
                     'full_original_Gram_bound_checked': True})
    return {'sites': sites, 'full_Fock_dimension': 1 << modes, 'c': c, 'U': u,
            'normalized_full_ground_rank': len(h)-rank, 'moments': rows,
            'arbitrary_precision_integer_products': True,
            'spectral_pole_or_mass_inferred': False}


def finite_control(sites):
    return _detached(_finite_control(sites))


def certificate():
    c, u, degree = F(1, 8), F(1), 6
    return {'source': 'same centered full invariant cubic ground and original CAR fields',
            'local_derivative_constants': parameters(c, u, degree)[:2],
            'default_original_charge_tails': [tail_budget(c, u, degree, e)
                                              for e in (8, 16, 20, 24, 32)],
            'explicit_exponential_tail_controls': [
                {'energy_cutoff': e, 'original_tail_upper': exponential_tail(c, u, degree, e)}
                for e in (24, 48, 96)],
            'default_exponential_moment_half_strip': exponential_moment(c, u, degree, F(1, 4)),
            'full_source_controls': [finite_control(sites) for sites in (2, 3)],
            'uniform_in_volume_and_normalized_wavepacket': True,
            'separate_charge_measures_determined_by_all_their_moments': True,
            'finite_moment_list_identifies_a_pole': False,
            'positive_shrinking_window_weight_proved': False,
            'physical_species': [], 'physical_masses': []}
