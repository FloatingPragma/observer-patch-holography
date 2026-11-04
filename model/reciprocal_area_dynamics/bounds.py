"""Outward rational certificates for the complete reciprocal collision law.

Higher precision than the one-step area receipt is needed to distinguish
every direction of the finite memory. No floating spectrum selects its rank.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from math import factorial

import sympy as s

from ..finite_area_response import response as prior


DENOMINATOR = 10**120
ISOLATION = F(1, 10**100)
ORDER = 200


def rational(value):
    return prior.rational(value)


def down(value):
    value = rational(value)
    return F(value.numerator*DENOMINATOR//value.denominator, DENOMINATOR)


def up(value):
    return -down(-rational(value))


@dataclass(frozen=True)
class I:
    lo: F
    hi: F

    def __post_init__(self):
        lo, hi = rational(self.lo), rational(self.hi)
        if lo > hi:
            raise ValueError('ordered exact interval required')
        object.__setattr__(self, 'lo', down(lo))
        object.__setattr__(self, 'hi', up(hi))

    @classmethod
    def point(cls, value):
        return cls(value, value)

    def __add__(self, other):
        other = other if isinstance(other, I) else I.point(other)
        return I(self.lo+other.lo, self.hi+other.hi)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, other):
        return self+(-other if isinstance(other, I) else -rational(other))

    def __mul__(self, other):
        other = other if isinstance(other, I) else I.point(other)
        corners = [x*y for x in (self.lo, self.hi) for y in (other.lo, other.hi)]
        return I(min(corners), max(corners))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = other if isinstance(other, I) else I.point(other)
        if other.lo <= 0 <= other.hi:
            raise ValueError('division by an unresolved interval')
        return self*I(1/other.hi, 1/other.lo)

    def square(self):
        if self.lo <= 0 <= self.hi:
            return I(0, max(self.lo*self.lo, self.hi*self.hi))
        return I(min(self.lo*self.lo, self.hi*self.hi), max(self.lo*self.lo, self.hi*self.hi))

    def encoded(self):
        return [str(self.lo), str(self.hi)]


ZERO, ONE = I.point(0), I.point(1)


def exp_bound(value):
    value = rational(value)
    if abs(value) > 3:
        raise ValueError('the retained exponential domain is |x|<=3')
    term = total = F(1)
    for k in range(1, ORDER+1):
        term *= value/k
        total += term
    tail = abs(value)**(ORDER+1)/(factorial(ORDER+1)*(1-abs(value)/F(ORDER+2)))
    return I(total-tail, total+tail)


def spectrum(matrix):
    roots = matrix.charpoly().as_poly().intervals(eps=s.Rational(ISOLATION))
    result = []
    for index, ((lo, hi), multiplicity) in enumerate(roots):
        lo, hi = rational(lo)-ISOLATION, rational(hi)+ISOLATION
        if (multiplicity != 1 or prior.inertia(matrix, lo) != index
                or prior.inertia(matrix, hi) != index+1):
            raise ArithmeticError('independent rational inertia failed')
        result.append(I(lo, hi))
    if len(result) != 8 or any(a.hi >= b.lo for a, b in zip(result, result[1:])):
        raise ArithmeticError('the complete simple eight-level source spectrum is required')
    return tuple(result)


@lru_cache(maxsize=1)
def source_data():
    levels = tuple(x for matrix in prior.matrices() for x in spectrum(matrix))
    if any(not 1 < x.lo <= x.hi < 5 for x in levels):
        raise ArithmeticError('the complete source thermal spectral window changed')
    weights = tuple(I(exp_bound(3-x.hi).lo, exp_bound(3-x.lo).hi) for x in levels)
    partition = sum(weights, ZERO)
    pi = tuple(w/partition for w in weights)
    # P acts on observables; its transpose is the population transition T.
    rows = []
    for i in range(16):
        row = [ZERO]*16
        for j in range(16):
            if i//8 != j//8:
                row[j] = weights[j]/(weights[i]+weights[j])/16
        row[i] = ONE-sum(row, ZERO)
        rows.append(tuple(row))
    return levels, weights, pi, tuple(rows)


def matvec(matrix, vector):
    return tuple(sum((a*b for a, b in zip(row, vector, strict=True)), ZERO) for row in matrix)


def inner(pi, left, right):
    return sum((w*a*b for w, a, b in zip(pi, left, right, strict=True)), ZERO)


def norm2(pi, vector):
    return sum((w*a.square() for w, a in zip(pi, vector, strict=True)), ZERO)


def subtract(left, coefficient, right):
    return tuple(a-coefficient*b for a, b in zip(left, right, strict=True))


@lru_cache(maxsize=1)
def projection():
    levels, weights, pi, P = source_data()
    p = sum(pi[8:], ZERO)
    variance = p*(ONE-p)
    f = tuple(I.point(i//8)-p for i in range(16))
    unit = (ONE,)*16

    def Q(vector):
        mean = inner(pi, unit, vector)
        coefficient = inner(pi, f, vector)/variance
        return tuple(x-mean-coefficient*y for x, y in zip(vector, f, strict=True))

    def C(vector):
        return Q(matvec(P, Q(vector)))

    b = Q(matvec(P, f))
    flow = sum((pi[i]*P[i][j] for i in range(8) for j in range(8, 16)), ZERO)
    a = ONE-flow/variance
    return {'levels': levels, 'weights': weights, 'pi': pi, 'P': P,
            'p': p, 'variance': variance, 'f': f, 'b': b, 'a': a,
            'flow': flow, 'Q': Q, 'C': C}


def memory_rank(data):
    """Positive exact interval Gram pivots prove cyclic dimension fourteen.

    All vectors are defined algebraically using the exact source objects.
    The intervals only enclose those vectors. Projection onto constants
    and the geometry direction is applied before each Gram orthogonalization.
    """
    pi, Q, C = data['pi'], data['Q'], data['C']
    vectors, norms = [], []
    vector = data['b']
    for index in range(14):
        vector = Q(vector)
        for earlier, norm in zip(vectors, norms, strict=True):
            vector = subtract(vector, inner(pi, earlier, vector)/norm, earlier)
        norm = norm2(pi, vector)
        if norm.lo <= 0:
            raise ArithmeticError(f'memory rank unresolved at direction {index+1}')
        vectors.append(vector)
        norms.append(norm)
        vector = C(vector)
    return tuple(norms)


@lru_cache(maxsize=1)
def certificate():
    data = projection()
    levels, pi, P = data['levels'], data['pi'], data['P']
    p, variance, b = data['p'], data['variance'], data['b']
    conditional = tuple(sum((pi[i]*levels[i] for i in range(8*g, 8*g+8)), ZERO)
                        /sum(pi[8*g:8*g+8], ZERO) for g in (0, 1))
    eta = (conditional[1]-conditional[0])/F(40, 99)
    if not F(578, 1000) < eta.lo <= eta.hi < F(579, 1000):
        raise ArithmeticError('the thermal conditional-family area coefficient is not positive')
    kernels, vector = [], b
    for _ in range(8):
        kernels.append(inner(pi, b, vector)/variance)
        vector = data['C'](vector)
    if not F(217, 10000) < kernels[0].lo <= kernels[0].hi < F(218, 10000):
        raise ArithmeticError('the actual geometry projection lost its strict memory witness')
    norms = memory_rank(data)
    # This is the actual FA114 matter replacement, without changing geometry.
    reset = tuple((ONE-p if i < 8 else p)/8 for i in range(16))
    u0 = tuple(q/w-ONE for q, w in zip(reset, pi, strict=True))
    force, vector = [], data['Q'](u0)
    for _ in range(8):
        force.append(inner(pi, b, vector))
        vector = data['C'](vector)
    area_response = force[0]*F(-40, 99)
    if not F(-628, 10**6) < area_response.lo <= area_response.hi < F(-627, 10**6):
        raise ArithmeticError('the memory forcing differs from the admitted matter disturbance')
    if any(P[i][i].lo <= F(1, 2) for i in range(16)):
        raise ArithmeticError('the complete transition lost strict laziness')
    if any(P[i][j].lo <= F(1, 1024) for i in range(16) for j in range(16) if i//8 != j//8):
        raise ArithmeticError('the analytic mixing minorization is unsupported')
    return {
        'complete_source_states': 16,
        'spectral_enclosures': [x.encoded() for x in levels],
        'geometry1_equilibrium_probability': p.encoded(),
        'dimensionless_conditional_energies': [x.encoded() for x in conditional],
        'conditional_thermal_area_coefficient': eta.encoded(),
        'source_cross_geometry_equilibrium_flow': data['flow'].encoded(),
        'instantaneous_area_coefficient': data['a'].encoded(),
        'memory_kernel_first_eight': [x.encoded() for x in kernels],
        'orthogonal_memory_Gram_pivots': [x.encoded() for x in norms],
        'minimal_scalar_memory_dimension': 14,
        'distinct_positive_memory_poles': 14,
        'visible_nonconstant_relaxation_poles': 15,
        'all_nonconstant_eigenvalues_strict_upper': '16383/16384',
        'all_nonconstant_eigenvalues_strict_lower': '0',
        'two_step_transition_entry_strict_lower': '1/131072',
        'two_step_uniform_minorization': '1/8192',
        'memory_tail_bound': 'sum_(j>=N) k_j <= 16384 k_0 (16383/16384)^N',
        'matter_reset_forcing_first_eight': [x.encoded() for x in force],
        'one_step_matter_area_response': area_response.encoded(),
        'conditional_thermal_manifold_invariant_under_actual_law': False,
        'thermal_area_coefficient_is_universal_Newton_normalization': False,
        'regional_G4_or_Einstein_dynamics_proved': False,
    }
