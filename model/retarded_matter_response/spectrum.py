"""Complete local ground-down spectral measure, with exact admission.

Modular reconstruction and a floating Hermitian eigensolver only propose
polynomial coefficients and brackets. Integer annihilation, a nonzero
modular Hankel determinant, and outward integer interval arithmetic admit
the complete78-line measure of c_origin,down^dagger G. No raw source
preparation, sampling, infinite-volume pole or relativistic mass is claimed.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from math import floor, isqrt

import numpy as np
from sympy import nextprime
from sympy.ntheory.residue_ntheory import sqrt_mod

from .response import dilation


DEGREE = 78
BITS = 1024
SCALE = 1 << BITS
ROOT_BITS = 384
NORM2 = 36**2


def _mul(a, b):
    return a[0]*b[0]+3*a[1]*b[1], a[0]*b[1]+a[1]*b[0]


def _add(a, b):
    return a[0]+b[0], a[1]+b[1]


def _bm(sequence, prime):
    """Modular Berlekamp--Massey proposer; no exact rank is inferred here."""
    c, old, degree, offset, previous = [1], [1], 0, 1, 1
    for n in range(len(sequence)):
        discrepancy = sum(c[i]*int(sequence[n-i]) for i in range(degree+1)) % prime
        if not discrepancy:
            offset += 1
            continue
        saved = c.copy()
        factor = discrepancy*pow(previous, -1, prime) % prime
        c += [0]*max(0, len(old)+offset-len(c))
        for j, value in enumerate(old):
            c[j+offset] = (c[j+offset]-factor*value) % prime
        if 2*degree <= n:
            degree, old, previous, offset = n+1-degree, saved, discrepancy, 1
        else:
            offset += 1
    return c, degree


def _modular_moments(h, a, b, prime, root, count):
    initial = (a+root*b) % prime
    value = initial.copy()
    # All reductions use int64, with explicit scalar preflights.
    if len(a)*(prime-1)**2 >= np.iinfo(np.int64).max:
        raise ArithmeticError('modular scalar product would overflow')
    row_norm = int(np.max(np.asarray(abs(h).sum(axis=1)), initial=0))
    if row_norm*(prime-1) >= np.iinfo(np.int64).max:
        raise ArithmeticError('modular source product would overflow')
    result = []
    for _ in range(count):
        result.append(int(initial@value) % prime)
        value = h@value % prime
    return result


def _det_mod(matrix, prime):
    a = np.array(matrix, dtype=np.int64, copy=True) % prime
    determinant = 1
    for column in range(len(a)):
        found = next((j for j in range(column, len(a)) if a[j, column]), None)
        if found is None:
            return 0
        if found != column:
            a[[column, found]] = a[[found, column]]
            determinant = -determinant
        pivot = int(a[column, column])
        determinant = determinant*pivot % prime
        inverse = pow(pivot, -1, prime)
        for j in range(column+1, len(a)):
            a[j] = (a[j]-(int(a[j, column])*inverse % prime)*a[column]) % prime
    return determinant % prime


def _source():
    data = dilation('ground', 'down')
    a, b = (x[:, 0].copy() for x in data.vectors)
    if data.h16.shape != (560, 560) or data.norm_square != 1:
        raise ArithmeticError('the original ground-down local source changed')
    rows = tuple(tuple(zip(data.h16.indices[data.h16.indptr[i]:data.h16.indptr[i+1]].tolist(),
                           data.h16.data[data.h16.indptr[i]:data.h16.indptr[i+1]].tolist(), strict=True))
                 for i in range(560))
    return data.h16, a, b, rows


def _apply(rows, value):
    """Unbounded Python integers; no overflow or floating admission."""
    return [sum(coefficient*value[j] for j, coefficient in row) for row in rows]


def _annihilates(coefficients, rows, a, b):
    value, other = [0]*len(a), [0]*len(a)
    for c, d in coefficients:
        value, other = _apply(rows, value), _apply(rows, other)
        value = [x+c*y+3*d*z for x, y, z in zip(value, a, b, strict=True)]
        other = [x+c*z+d*y for x, y, z in zip(other, a, b, strict=True)]
    return not any(value) and not any(other)


@lru_cache(maxsize=1)
def exact_data():
    h, a, b, rows = _source()
    # This finite reconstruction is a proposal. Its height estimate is not
    # an acceptance premise: the exact vector equation below is mandatory.
    modulus, rational, radical, prime, primes = 1, [0]*79, [0]*79, 10**6, []
    while modulus <= 2*26**DEGREE:
        prime = int(nextprime(prime))
        roots = sqrt_mod(3, prime, all_roots=True)
        if len(roots) != 2:
            continue
        polynomials = []
        for root in roots:
            proposed, degree = _bm(_modular_moments(h, a, b, prime, int(root), 160), prime)
            if degree != DEGREE or len(proposed) != DEGREE+1:
                raise ArithmeticError('modular proposal did not retain degree78')
            polynomials.append(proposed)
        inv2, inv2root = pow(2, -1, prime), pow(2*int(roots[0]), -1, prime)
        ca = [(x+y)*inv2 % prime for x, y in zip(*polynomials, strict=True)]
        cb = [(x-y)*inv2root % prime for x, y in zip(*polynomials, strict=True)]
        inverse = pow(modulus, -1, prime)
        rational = [x+modulus*((y-x)*inverse % prime) for x, y in zip(rational, ca, strict=True)]
        radical = [x+modulus*((y-x)*inverse % prime) for x, y in zip(radical, cb, strict=True)]
        modulus *= prime
        primes.append(prime)
    rational = [x if x <= modulus//2 else x-modulus for x in rational]
    radical = [x if x <= modulus//2 else x-modulus for x in radical]
    polynomial = tuple(zip(rational, radical, strict=True))
    aa, bb = a.tolist(), b.tolist()
    if polynomial[0] != (1, 0) or not _annihilates(polynomial, rows, aa, bb):
        raise ArithmeticError('proposed spectral polynomial fails exact source annihilation')
    # 97 is prime and10^2=3 mod97. The nonzero Hankel determinant is
    # a lower bound on the characteristic-zero cyclic dimension.
    moments97 = _modular_moments(h, a, b, 97, 10, 2*DEGREE-1)
    hankel = [[moments97[i+j] for j in range(DEGREE)] for i in range(DEGREE)]
    determinant = _det_mod(hankel, 97)
    if not determinant:
        raise ArithmeticError('the exact cyclic-rank lower witness vanished')
    moments, value, other = [], aa, bb
    for _ in range(DEGREE):
        moments.append((sum(x*y+3*z*w for x, y, z, w in zip(aa, value, bb, other, strict=True)),
                        sum(x*w+z*y for x, y, z, w in zip(aa, value, bb, other, strict=True))))
        value, other = _apply(rows, value), _apply(rows, other)
    if moments[0] != (NORM2, 0):
        raise ArithmeticError('local charged spectral measure is not normalized')
    # <v,(z-H16)^-1 v>=numerator(z)/(1296 polynomial(z)).
    numerator = []
    for j in range(DEGREE):
        value = (0, 0)
        for k in range(j+1):
            value = _add(value, _mul(polynomial[k], moments[j-k]))
        numerator.append(value)
    derivative = tuple(((DEGREE-i)*x, (DEGREE-i)*y)
                       for i, (x, y) in enumerate(polynomial[:-1]))
    # The previously symmetry-resolved lines are exactly1/10 each; the
    # complete minimal polynomial now excludes any coincident residual atom.
    for point in ((0, -3), (16, -3)):
        def evaluate(coefficients):
            value = (0, 0)
            for coefficient in coefficients:
                value = _add(_mul(value, point), coefficient)
            return value
        num, den = evaluate(numerator), evaluate(derivative)
        if evaluate(polynomial) != (0, 0) or den == (0, 0) or tuple(10*x for x in num) != tuple(NORM2*x for x in den):
            raise ArithmeticError('known exact charged line or its weight changed')
    return {'polynomial': polynomial, 'numerator': tuple(numerator), 'derivative': derivative,
            'moments': tuple(moments), 'rank_prime': 97, 'rank_root': 10,
            'hankel_determinant': determinant, 'reconstruction_primes': tuple(primes)}


@dataclass(frozen=True)
class _Interval:
    """Fixed1024-bit outward dyadic intervals, with integer-only operations."""
    lo: int
    hi: int

    def __post_init__(self):
        if type(self.lo) is not int or type(self.hi) is not int or self.lo > self.hi:
            raise ValueError('ordered integer interval endpoints required')

    def __add__(self, other):
        return _Interval(self.lo+other.lo, self.hi+other.hi)

    def __mul__(self, other):
        products = (self.lo*other.lo, self.lo*other.hi, self.hi*other.lo, self.hi*other.hi)
        return _Interval(min(products)//SCALE, -((-max(products))//SCALE))

    def fractions(self):
        return F(self.lo, SCALE), F(self.hi, SCALE)


def _coefficient(pair):
    a, b = pair
    root = isqrt(3*SCALE*SCALE)
    if not root*root < 3*SCALE*SCALE < (root+1)**2:
        raise ArithmeticError('strict sqrt3 integer enclosure failed')
    return _Interval(a*SCALE+min(b*root, b*(root+1)), a*SCALE+max(b*root, b*(root+1)))


def _evaluate(coefficients, point):
    value = _Interval(0, 0)
    for coefficient in coefficients:
        value = value*point+coefficient
    return value


def _point_sign(coefficients, numerator, bits):
    value = numerator << (BITS-bits)
    interval = _evaluate(coefficients, _Interval(value, value))
    return -1 if interval.hi < 0 else 1 if interval.lo > 0 else 0


def _brackets(coefficients):
    # An eigensolver proposes intervals only. EXACT78 disjoint sign changes
    # certify that every root of the degree78 polynomial is accounted for.
    eigenvalues = np.linalg.eigvalsh(_source()[0].toarray())
    if not np.all(np.isfinite(eigenvalues)):
        raise ArithmeticError('nonfinite root proposal')
    bits, result = 24, []
    for value in eigenvalues:
        left = floor(float(value)*(1 << bits))-1
        right = left+3
        if result and left <= result[-1][1]:
            continue
        sl, sr = _point_sign(coefficients, left, bits), _point_sign(coefficients, right, bits)
        if sl*sr == -1:
            result.append((left, right))
    if len(result) != DEGREE:
        raise ArithmeticError('proposed root brackets do not exhaust78 exact sign changes')
    return tuple(result), bits


def _refine(coefficients, bracket, bits, target=ROOT_BITS):
    left, right = bracket
    sign = _point_sign(coefficients, left, bits)
    if sign*_point_sign(coefficients, right, bits) != -1:
        raise ValueError('strict sign-changing bracket required')
    while bits < target:
        middle = left+right
        sm = _point_sign(coefficients, middle, bits+1)
        if not sm:
            raise ArithmeticError('fixed exact interval precision cannot decide root bisection')
        if sm == sign:
            left, right = middle, 2*right
        else:
            left, right = 2*left, middle
        bits += 1
    return _Interval(left << (BITS-bits), right << (BITS-bits))


def _ratio(numerator, denominator):
    a, b = numerator.fractions(), denominator.fractions()
    if b[0] <= 0 <= b[1]:
        raise ArithmeticError('residue denominator enclosure contains zero')
    values = [x/y/NORM2 for x in a for y in b]
    return min(values), max(values)


def _outward(interval, places=24):
    scale = 10**places
    low, high = (x*scale for x in interval)
    return F(low.numerator//low.denominator, scale), F(-((-high.numerator)//high.denominator), scale)


@lru_cache(maxsize=1)
def certificate():
    data = exact_data()
    pol, num, derivative = (tuple(_coefficient(x) for x in data[key])
                            for key in ('polynomial', 'numerator', 'derivative'))
    brackets, bits = _brackets(pol)
    rows, total_low, total_high = [], F(), F()
    for index, bracket in enumerate(brackets):
        root = _refine(pol, bracket, bits)
        weight = _ratio(_evaluate(num, root), _evaluate(derivative, root))
        if not 0 < weight[0] <= weight[1] < 1 or weight[1]-weight[0] >= F(1, 10**20):
            raise ArithmeticError('positive narrow spectral residue enclosure failed')
        shifted = root+_coefficient((0, 4))
        energy = tuple(x/16 for x in shifted.fractions())
        if energy[0] <= 0:
            raise ArithmeticError('resolved charged frequency is not positive above G')
        energy, weight = _outward(energy), _outward(weight)
        if weight[0] <= 0:
            raise ArithmeticError('published outward weight lost strict positivity')
        total_low += weight[0]
        total_high += weight[1]
        rows.append({'line': index+1, 'gap_above_ground_interval': list(map(str, energy)),
                     'normalized_weight_interval': list(map(str, weight)),
                     'diagnostic_gap': float(sum(energy)/2), 'diagnostic_weight': float(sum(weight)/2),
                     'H16_root_bracket': [str(F(bracket[0], 1 << bits)), str(F(bracket[1], 1 << bits))]})
    if not total_low <= 1 <= total_high or total_high-total_low >= F(1, 10**20):
        raise ArithmeticError('complete spectral normalization enclosure failed')
    return {'source':'original cube H, c=1/8,U=1, groundG at E_G=-sqrt3/4',
            'probe':'c_origin,down^dagger G, norm1, physical frame at origin unchanged',
            'full_source_sector_dimension':560, 'exact_local_cyclic_dimension':DEGREE,
            'resolved_distinct_positive_lines':DEGREE,
            'spectral_variable':'x is an eigenvalue of16H_N5; measured gap=(x+4sqrt3)/16',
            'polynomial_descending_Qsqrt3':[list(x) for x in data['polynomial']],
            'resolvent_numerator_descending_Qsqrt3':[list(x) for x in data['numerator']],
            'resolvent_formula':'<v,(z-16H_N5)^(-1)v>=numerator(z)/(1296 polynomial(z))',
            'exact_vector_annihilation':True, 'rank_witness_prime':data['rank_prime'],
            'rank_witness_sqrt3_image':data['rank_root'],
            'rank_witness_78_Hankel_determinant':data['hankel_determinant'],
            'reconstruction_is_proposer_only':True,
            'reconstruction_primes':list(data['reconstruction_primes']),
            'floating_eigenvalues_are_bracket_proposals_only':True,
            'integer_interval_bits':BITS, 'root_bisection_bits':ROOT_BITS,
            'rows':rows, 'weight_sum_interval':[str(total_low),str(total_high)],
            'known_exact_lines':[{'gap':'sqrt3/16','weight':'1/10'},
                                 {'gap':'1+sqrt3/16','weight':'1/10'}],
            'remaining_symmetry_component_lines':76, 'remaining_symmetry_component_total_weight':'4/5',
            'unresolved_local_probe_weight':'0',
            'entire560_state_spectrum_or_degeneracies_resolved':False,
            'native_preparation_detector_or_finite_shot_resolution_executed':False,
            'infinite_volume_poles_or_relativistic_masses_identified':False}
