"""Rational matrix and thermal bounds for the unchanged reciprocal source.

The primary response changes the matter state, then applies one unchanged
thermal update. A separate equilibrium susceptibility changes the matter
Hamiltonian and requires the corresponding bath and energy-readout family.
Neither response is an Einstein or continuum-field construction.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from math import factorial, isqrt

import sympy as s

from ..native_equilibration import reciprocal_vacuum as source


DENOMINATOR = 10**30
ORDER = 60
GEOMETRIES = (F(9, 10), F(11, 10))


def rational(value):
    if type(value) in (int, F):
        return F(value)
    if isinstance(value, s.Rational):
        return F(int(value.p), int(value.q))
    raise ValueError('exact real rational input required')


def floor(value):
    value = rational(value)
    return F(value.numerator*DENOMINATOR//value.denominator, DENOMINATOR)


def ceil(value):
    return -floor(-rational(value))


@dataclass(frozen=True)
class Interval:
    lo: F
    hi: F

    def __post_init__(self):
        lo, hi = rational(self.lo), rational(self.hi)
        if lo > hi:
            raise ValueError('ordered rational enclosure required')
        object.__setattr__(self, 'lo', floor(lo))
        object.__setattr__(self, 'hi', ceil(hi))

    @classmethod
    def point(cls, value):
        value = rational(value)
        return cls(value, value)

    def __add__(self, other):
        other = other if isinstance(other, Interval) else Interval.point(other)
        return Interval(self.lo+other.lo, self.hi+other.hi)

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self+-other if isinstance(other, Interval) else self+(-rational(other))

    def __mul__(self, other):
        other = other if isinstance(other, Interval) else Interval.point(other)
        corners = [a*b for a in (self.lo, self.hi) for b in (other.lo, other.hi)]
        return Interval(min(corners), max(corners))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = other if isinstance(other, Interval) else Interval.point(other)
        if other.lo <= 0 <= other.hi:
            raise ValueError('division by an unresolved or zero interval')
        return self*Interval(1/other.hi, 1/other.lo)

    def encoded(self):
        return [str(self.lo), str(self.hi)]


def exponential(point):
    point = rational(point)
    if abs(point) > 3:
        raise ValueError('the exponential certificate requires |x| <= 3')
    term, total = F(1), F(1)
    for n in range(1, ORDER+1):
        term *= point/n
        total += term
    bound = abs(point)
    tail = bound**(ORDER+1)/(factorial(ORDER+1)*(1-bound/F(ORDER+2)))
    return Interval(total-tail, total+tail)


def beta_interval():
    root = isqrt(14*DENOMINATOR**2)
    return Interval(F(DENOMINATOR, root+1), F(DENOMINATOR, root))


@lru_cache(maxsize=1)
def matrices():
    """Independent rational beta*H, checked against the source action."""
    x = s.Matrix([[0, 1], [1, 0]])
    eye = s.eye(2)
    coupling = s.kronecker_product(x, x, eye)+s.kronecker_product(eye, x, x)
    diagonal = s.diag(*(2*mask.bit_count()+3 for mask in range(8)))
    result = []
    for r, actual in zip(GEOMETRIES, source.default_hamiltonians(), strict=True):
        a = r/4+1/(7*r)+3*r*r/28
        c = -1/(7*r)
        exact = s.ImmutableMatrix(s.Rational(a)*diagonal+s.Rational(c)*coupling)
        if exact != (source.BETA*actual).applyfunc(s.simplify):
            raise ArithmeticError('independent rational action differs from the reciprocal source')
        result.append(exact)
    if tuple(source.GEOMETRIES) != GEOMETRIES or source.BETA != 1/s.sqrt(14):
        raise ArithmeticError('the retained geometry/temperature specification changed')
    return tuple(result)


def inertia(matrix, point):
    """Rational LDL congruence independently certifies eigenvalue counts.

    Nonzero scalar pivots are required. The root-isolation endpoints are
    expanded rationally, and every admitted endpoint is checked here.
    """
    point = rational(point)
    if matrix.rows != matrix.cols or matrix != matrix.T:
        raise ValueError('exact symmetric square matrix required')
    a = [[rational(matrix[i, j])-(point if i == j else 0)
          for j in range(matrix.cols)] for i in range(matrix.rows)]
    negative = 0
    for k in range(len(a)):
        pivot = a[k][k]
        if pivot == 0:
            raise ValueError('zero LDL pivot is not a certified spectral endpoint')
        negative += pivot < 0
        for i in range(k+1, len(a)):
            for j in range(i, len(a)):
                updated = a[i][j]-a[i][k]*a[k][j]/pivot
                a[i][j] = a[j][i] = updated
    return negative


def spectral_enclosures(matrix):
    epsilon = F(1, 10**22)
    isolated = matrix.charpoly().as_poly().intervals(eps=s.Rational(epsilon))
    rows = []
    for index, ((lo, hi), multiplicity) in enumerate(isolated):
        lo, hi = rational(lo)-epsilon, rational(hi)+epsilon
        if multiplicity != 1 or inertia(matrix, lo) != index or inertia(matrix, hi) != index+1:
            raise ArithmeticError('independent matrix inertia did not certify a simple source level')
        rows.append(Interval(lo, hi))
    if len(rows) != 8 or any(a.hi >= b.lo for a, b in zip(rows, rows[1:])):
        raise ArithmeticError('the full eight-state spectrum is not isolated')
    return tuple(rows)


def spectral_weights(matrix):
    spectrum = spectral_enclosures(matrix)
    # The common exp(3) factor cancels from every probability and rate.
    weights = tuple(Interval(exponential(3-e.hi).lo, exponential(3-e.lo).hi)
                    for e in spectrum)
    return spectrum, weights


def matrix_exponential_diagonal(matrix):
    """Spectral-free rational Taylor enclosure of exp(3I-beta H)."""
    size = matrix.rows
    centered = [[rational(matrix[i, j])-(3 if i == j else 0)
                 for j in range(size)] for i in range(size)]
    bound = max(sum(abs(x) for x in row) for row in centered)
    if bound > 3:
        raise ValueError('the retained rational matrix-exponential domain was exceeded')
    rows = [[(k, value) for k, value in enumerate(row) if value] for row in centered]
    term = [[F(int(i == j)) for j in range(size)] for i in range(size)]
    total = [row.copy() for row in term]
    for n in range(1, ORDER+1):
        term = [[-sum((value*term[k][j] for k, value in rows[i]), F())/n
                 for j in range(size)] for i in range(size)]
        total = [[total[i][j]+term[i][j] for j in range(size)] for i in range(size)]
    tail = bound**(ORDER+1)/(factorial(ORDER+1)*(1-bound/F(ORDER+2)))
    return tuple(Interval(total[i][i]-tail, total[i][i]+tail) for i in range(size)), bound, tail


def fixed_law_response(weights, area=(F(20, 9), F(20, 11))):
    if len(weights) != 2 or any(len(block) != 8 for block in weights):
        raise ValueError('both complete eight-state field spectra are required')
    if any(w.lo <= 0 for block in weights for w in block):
        raise ValueError('faithful positive thermal weights required')
    partition = [sum(block, Interval.point(0)) for block in weights]
    p1 = partition[1]/sum(partition, Interval.point(0))
    cross = Interval.point(0)
    for w0 in weights[0]:
        for w1 in weights[1]:
            cross += Interval(w1.lo/(w0.hi+w1.lo), w1.hi/(w0.lo+w1.hi))
    probability_response = (cross-64*p1)/128
    return {'p1': p1, 'cross': cross, 'dp1': probability_response,
            'area': probability_response*(rational(area[1])-rational(area[0]))}


def equilibrium_response():
    """Same local occupation perturbation in both Hamiltonian blocks.

    The Gibbs law and energy-dependent reservoir bank follow the perturbed
    Hamiltonians. This is separate from the fixed-law state intervention.
    """
    diagonals, proofs = [], []
    occupation = s.diag(0, 0, 0, 0, 1, 1, 1, 1)
    for matrix in matrices():
        if matrix*occupation == occupation*matrix:
            raise ArithmeticError('the noncommuting matter perturbation control disappeared')
        diagonal, bound, tail = matrix_exponential_diagonal(matrix)
        diagonals.append(diagonal)
        proofs.append({'centered_norm_bound': str(bound), 'Taylor_tail_upper': str(ceil(tail))})
    z = [sum(block, Interval.point(0)) for block in diagonals]
    population = [value/(z[0]+z[1]) for value in z]
    mean = [sum(block[4:], Interval.point(0))/weight for block, weight in zip(diagonals, z, strict=True)]
    covariance = population[0]*population[1]*F(40, 99)*(mean[0]-mean[1])
    response = -beta_interval()*covariance
    return {'population': population, 'occupation': mean, 'response': response, 'proof': proofs}


def metric_certificate():
    r = s.Symbol('r', positive=True)
    density = s.diag(-1/r, 1/r, r*r, r*r)
    metric = (s.sqrt(-density.det())*density.inv()).applyfunc(s.simplify)
    area = s.simplify(2*s.sqrt(metric[2:4, 2:4].det()))
    volume = s.simplify(s.sqrt(metric[1:4, 1:4].det()))
    if area != 2/r or volume != 1:
        raise ArithmeticError('independent fixed-volume metric determinant changed')
    return {'metric': [[str(x) for x in row] for row in metric.tolist()],
            'two_face_area': str(area), 'normalized_proper_volume': str(volume),
            'branch_lapse': 'r',
            'clock_dictionary': 'd tau=r dt; proper generator H/r; proper inverse temperature r beta',
            'one_update_duration_or_common_proper_temperature_calibrated': False,
            'branch_areas': ['20/9', '20/11'],
            'branch_area_difference': '40/99',
            'internal_geometry_block_dimensions': [1, 1],
            'internal_edge_log_dimensions': ['0', '0'],
            'positive_affine_area_coefficient_from_these_internal_dimensions': False,
            'geometry_Shannon_entropy_is_separate': True}


def randomization_certificate():
    total = s.zeros(8)
    for i in range(8):
        for j in range(8):
            matrix_unit = s.zeros(8)
            matrix_unit[i, j] = 1
            total += matrix_unit.T*matrix_unit/8
    if total != s.eye(8):
        raise ArithmeticError('the complete randomization instrument lost an outcome')
    return {'replacement_Kraus_records': 64, 'identity_mixture_record': 1,
            'all_65_records_retained_at_zero_and_unit_amplitudes': True,
            'complete_trace_preservation_exact': True,
            'geometry_marginal_preserved_before_bath_update': True,
            'possible_dilation': 'SWAP three field bits with three supplied maximally mixed bits; retain displaced input and mixture flag',
            'initialization_memory_and_intervention_energy_are_supplied': True}


def certificate():
    spectra, weights = zip(*(spectral_weights(matrix) for matrix in matrices()), strict=True)
    primary = fixed_law_response(weights)
    secondary = equilibrium_response()
    if not (F(-628, 10**6) < primary['area'].lo <= primary['area'].hi < F(-627, 10**6)):
        raise ArithmeticError('the complete fixed-law area response sign was not certified')
    if not (F(-675, 10**6) < secondary['response'].lo <= secondary['response'].hi < F(-673, 10**6)):
        raise ArithmeticError('the noncommuting equilibrium susceptibility sign was not certified')
    identical = fixed_law_response((weights[0], weights[0]))
    if not identical['area'].lo <= 0 <= identical['area'].hi:
        raise ArithmeticError('equal-Hamiltonian deletion control lost its zero response')
    tracial = fixed_law_response(((Interval.point(1),)*8,)*2)
    if tracial['area'] != Interval.point(0):
        raise ArithmeticError('the infinite-temperature control is not exactly zero')
    flat_area = fixed_law_response(weights, (F(1), F(1)))
    if flat_area['area'] != Interval.point(0):
        raise ArithmeticError('a constant geometric observable responded')
    partitions = [sum(block, Interval.point(0)) for block in weights]
    total_weight = sum(partitions, Interval.point(0))
    equilibrium_energy = sum((level*weight for levels, block in zip(spectra, weights, strict=True)
                              for level, weight in zip(levels, block, strict=True)), Interval.point(0))/total_weight
    reset_energy = sum((partition/total_weight*rational(s.trace(matrix)/8)
                        for matrix, partition in zip(matrices(), partitions, strict=True)), Interval.point(0))
    work = (reset_energy-equilibrium_energy)/beta_interval()
    if not F(2697, 1000) < work.lo <= work.hi < F(1349, 500):
        raise ArithmeticError('the matter randomization energy injection was not certified')
    return {'geometry': metric_certificate(),
            'dimensionless_source_matrices': [[[str(x) for x in row] for row in matrix.tolist()]
                                            for matrix in matrices()],
            'complete_spectra': [[row.encoded() for row in block] for block in spectra],
            'LDL_eigenvalue_counts': [[i, i+1] for i in range(8)],
            'exact_detailed_balance': source.exact_balance_identity()['population_detailed_balance_exact'],
            'complete_randomization_instrument': randomization_certificate(),
            'fixed_law_state_response': {
                'intervention': '(1-epsilon) id + epsilon replacement_by_I8/8 on matter; geometry unchanged',
                'unchanged_Hamiltonian_bath_gaps_readout_and_transition': True,
                'all_64_cross_geometry_energy_edges_and_rejections_retained': True,
                'source_native_randomization_and_full_controller_executed': False,
                'geometry1_equilibrium_probability': primary['p1'].encoded(),
                'cross_rate_sum': primary['cross'].encoded(),
                'geometry1_probability_slope': primary['dp1'].encoded(),
                'area_slope': primary['area'].encoded(),
                'mean_intervention_energy_injection_per_unit_epsilon': work.encoded(),
                'apparatus_work_heat_decomposition_derived': False,
                'simplified_area_slope_enclosure': ['-157/250000', '-627/1000000'],
                'valid_mixture_amplitudes': ['0', '1'],
                'response_exactly_linear_in_epsilon': True},
            'equilibrium_Hamiltonian_family_response': {
                'perturbation': 'H_g(lambda)=H_g+lambda n_left in both branches',
                'requires_matching_bath_gaps_and_energy_readout_family': True,
                'commuting_perturbation_assumed': False,
                'conditional_occupation': [x.encoded() for x in secondary['occupation']],
                'area_susceptibility': secondary['response'].encoded(),
                'matrix_exponential_proofs': secondary['proof']},
            'controls': {'identical_conditional_Hamiltonians_zero_response': True,
                         'constant_area_zero_response': True,
                         'beta_zero_zero_response': True,
                         'equilibrium_input_unchanged_by_full_thermal_law': True},
            'complete_interacting_QFT_or_Einstein_geometry_established': False,
            'same_reference_as_Weyl_particle_source': False}
