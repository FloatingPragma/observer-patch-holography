"""Rational reconstruction from the certified complete eight-spin ground."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product

import numpy as np
import sympy as sp

from model.native_ordered_ground import certificate, experiment, source


def upper(value, bits=64):
    value = F(value)
    scale = 1 << bits
    return F(-((-value.numerator*scale)//value.denominator), scale)


def nearest(numerator, denominator):
    if denominator <= 0:
        raise ValueError('positive integer denominator required')
    return (2*numerator+denominator)//(2*denominator)


def partial_outer(left, right):
    if len(left) != 256 or len(right) != 256:
        raise ValueError('complete eight-spin vectors required')
    return [[sum(int(left[4*k+i])*int(right[4*k+j]) for k in range(64))
             for j in range(4)] for i in range(4)]


def integer_energy(site):
    """Eight times the actual half-bond plus selection-field density."""
    if type(site) is not int or not 0 <= site < 8:
        raise ValueError('a declared cube site is required')
    edges, odd = source.cube()
    out = np.zeros((256, 256), dtype=np.int64)
    for mask in range(256):
        for x, y in edges:
            if site not in (x, y):
                continue
            if ((mask >> x) ^ (mask >> y)) & 1:
                out[mask, mask] -= 2
                out[mask ^ (1 << x) ^ (1 << y), mask] += 2
        out[mask ^ (1 << site), mask] += 1 if site in odd else -1
    return out


@lru_cache(maxsize=1)
def state():
    witness = experiment.witness()
    checked = certificate.verify(witness)
    h, scale = source.matrix()
    if scale != 8:
        raise ArithmeticError('the declared eight-site source scale changed')
    odd = source.cube()[1]
    sign = np.array([(-1)**sum((mask >> j) & 1 for j in odd)
                     for mask in range(256)], dtype=np.int64)
    h = sign[:, None]*h*sign[None, :]
    v = tuple(int(x)*int(s) for x, s in zip(witness['vectors'][:, 0], sign, strict=True))
    norm = sum(x*x for x in v)
    e = [integer_energy(x) for x in range(8)]
    if not np.array_equal(sum(e), h):
        raise ArithmeticError('source energy densities do not sum to its action')
    raw = partial_outer(v, v)
    rho = tuple(tuple(F(x, norm) for x in row) for row in raw)
    if sum(rho[j][j] for j in range(4)) != 1:
        raise ArithmeticError('regional restriction lost normalization')
    # Sylvester's criterion is exact; no numerical eigenvalue supplies a floor.
    floor = F(1, 40)
    shifted = sp.Matrix([[sp.Rational(x.numerator, x.denominator)
                          for x in row] for row in rho])-sp.eye(4)/40
    minors = tuple(shifted[:k, :k].det() for k in range(1, 5))
    if any(x <= 0 for x in minors):
        raise ArithmeticError('the source regional spectral floor failed')
    error = F(checked['normalized_ground_vector_error'])
    if not error < floor-F(1, 41):
        raise ArithmeticError('actual selected marginal lacks the retained floor')
    # Products use a proved int64 range; the actual matrix remains unreduced.
    if 256*int(np.max(abs(h)))*int(np.max(abs(e[0]))) >= 1 << 60:
        raise ArithmeticError('exact energy-current integer product could overflow')
    current = h@e[0]-e[0]@h
    if not np.array_equal(current.T, -current):
        raise ArithmeticError('the source current is not anti-Hermitian before i')
    cv = tuple(sum(int(current[i, j])*v[j] for j in range(256)) for i in range(256))
    left = partial_outer(cv, v)
    right = partial_outer(v, cv)
    tangent = tuple(tuple(F(left[i][j]+right[i][j], 64*norm)
                          for j in range(4)) for i in range(4))
    if sum(tangent[j][j] for j in range(4)) != 0:
        raise ArithmeticError('source current tangent lost trace zero')
    norm_current = F(max(sum(abs(int(x)) for x in row) for row in current), 64)
    return {'rho': rho, 'tangent': tangent, 'vector': v, 'vector_norm_squared': norm,
            'ground_vector_error': error, 'reference_floor': floor,
            'actual_floor': F(1, 41), 'principal_minors': minors,
            'energy_integer': h, 'density_integer': e[0],
            'current_imaginary_integer': current, 'current_norm_upper': norm_current,
            'full_source_certificate': checked}


@lru_cache(maxsize=8, typed=True)
def logarithm(order=1024, bits=96):
    """Fixed-point polynomial -log(rho), with every rounding cost charged."""
    if type(order) is not int or order < 1 or type(bits) is not int or bits < 16:
        raise ValueError('positive series order and at least 16 fractional bits required')
    data = state()
    scale = 1 << bits
    r = np.array([[nearest(x.numerator*scale, x.denominator) for x in row]
                  for row in data['rho']], dtype=object)
    if not np.array_equal(r, r.T):
        raise ArithmeticError('regional quantization is not symmetric')
    q = scale*np.eye(4, dtype=object)-r
    rounding = F(2, scale)  # four entries per row, each within half a grid unit
    contraction = 1-data['reference_floor']+rounding
    if not 0 < contraction < 1:
        raise ValueError('quantized logarithm lacks a convergent series')
    power = scale*np.eye(4, dtype=object)
    result = np.zeros((4, 4), dtype=object)
    for n in range(1, order+1):
        # Symmetrization preserves Hermiticity and contracts the accumulated error.
        multiplication = power@q+q@power
        power = np.array([[nearest(int(x), 2*scale) for x in row]
                          for row in multiplication], dtype=object)
        term = np.array([[nearest(int(x), n) for x in row]
                         for row in power], dtype=object)
        result += term
    tail = upper(contraction**(order+1)/((order+1)*(1-contraction)))
    polynomial_rounding = 2*order*rounding
    reference_error = data['ground_vector_error']+rounding
    log_reference_error = reference_error/(data['reference_floor']-reference_error)
    error = upper(tail+polynomial_rounding+log_reference_error)
    matrix = tuple(tuple(F(int(x), scale) for x in row) for row in result)
    return {'matrix': matrix, 'operator_error_upper': error, 'order': order, 'bits': bits,
            'tail_upper': tail, 'fixed_point_error_upper': polynomial_rounding,
            'reference_log_error_upper': upper(log_reference_error)}


PAULI = {'I': sp.eye(2), 'X': sp.Matrix([[0, 1], [1, 0]]),
         'Y': sp.Matrix([[0, -sp.I], [sp.I, 0]]), 'Z': sp.diag(1, -1)}


def pauli(label):
    if len(label) != 2 or any(x not in PAULI for x in label):
        raise ValueError('two declared Pauli factors required')
    return sp.kronecker_product(*(PAULI[x] for x in reversed(label)))


def pairing(matrix, operator):
    value = sum(matrix[i][j]*operator[j, i] for i in range(4) for j in range(4))
    value = sp.simplify(value)
    if not value.is_Rational:
        raise ArithmeticError('Hermitian regional pairing is not real rational')
    return F(int(value.p), int(value.q))


def coefficients(matrix):
    result = {''.join(label): pairing(matrix, pauli(''.join(label)))/4
              for label in product('IXYZ', repeat=2)}
    rebuilt = sum((sp.Rational(v.numerator, v.denominator)*pauli(k)
                   for k, v in result.items()), sp.zeros(4))
    target = sp.Matrix([[sp.Rational(v.numerator, v.denominator) for v in row]
                        for row in matrix])
    if rebuilt != target:
        raise ArithmeticError('complete Pauli reconstruction lost a coefficient')
    return result


def source_response(log=None):
    data = state()
    log = logarithm() if log is None else log
    coef = coefficients(log['matrix'])
    eta, epsilon, v = log['operator_error_upper'], data['ground_vector_error'], data['current_norm_upper']
    # ||K_actual|| <= log(41) <4, derived from the actual positive floor.
    slope = sum(data['tangent'][i][j]*log['matrix'][j][i]
                for i in range(4) for j in range(4))
    slope_error = upper(2*v*eta+16*v*epsilon)
    difference = coef['XX']-coef['YY']
    if not difference+eta < F(-13, 100):
        raise ArithmeticError('physical-energy boost exclusion lost its strict witness')
    if not slope+slope_error < F(-53, 1000):
        raise ArithmeticError('source entropy/area countercontrol lost its strict witness')
    fisher = F(0)
    actual_fisher_lower, actual_fisher_upper = F(0), F(0)
    outcomes = []
    for label in coef:
        if label == 'II':
            outcomes.append({'pauli': label, 'probabilities': (F(1, 16), F(0)),
                             'derivatives': (F(0), F(0))})
            continue
        mean = pairing(data['rho'], pauli(label))
        rate = pairing(data['tangent'], pauli(label))
        if abs(mean) >= 1:
            raise ArithmeticError('complete Pauli readout lacks a positive outcome')
        fisher += rate*rate/(16*(1-mean*mean))
        outcomes.append({'pauli': label, 'probabilities': ((1+mean)/32, (1-mean)/32),
                         'derivatives': (rate/32, -rate/32)})
        def square_interval(lo, hi):
            return (F(0) if lo <= 0 <= hi else min(lo*lo, hi*hi), max(lo*lo, hi*hi))
        rlo, rhi = square_interval(mean-2*epsilon, mean+2*epsilon)
        alo, ahi = square_interval(rate-4*v*epsilon, rate+4*v*epsilon)
        if rhi >= 1:
            raise ArithmeticError('actual complete readout lacks a certified support floor')
        actual_fisher_lower += alo/(16*(1-rlo))
        actual_fisher_upper += ahi/(16*(1-rhi))
    hs = sum(x*x for row in data['tangent'] for x in row)
    mu = data['reference_floor']
    lower, higher = hs/4, hs/(32*mu*(1-2*mu))
    if not 0 < lower <= fisher <= higher:
        raise ArithmeticError('complete detector response metric failed')
    return {'pauli_coefficients': coef, 'XX_minus_YY': difference,
            'coefficient_error_upper': eta,
            'regional_entropy_slope': slope, 'entropy_slope_error_upper': slope_error,
            'exact_ground_total_energy_derivative': F(0), 'static_cut_area_derivative': F(0),
            'complete_pauli_response_records': outcomes, 'rational_reference_fisher_metric': fisher,
            'rational_reference_fisher_lower': lower, 'rational_reference_fisher_upper': higher,
            'actual_source_fisher_lower': actual_fisher_lower,
            'actual_source_fisher_upper': actual_fisher_upper,
            'tangent_hilbert_schmidt_squared': hs}


def finite_pulse(parameter=F(1, 10000)):
    """Strict finite entropy change under the same source-current control."""
    if type(parameter) is not F or parameter <= 0:
        raise ValueError('positive exact current-pulse parameter required')
    data, response = state(), source_response()
    v, mu = data['current_norm_upper'], data['actual_floor']
    distance = parameter*v
    if distance >= mu:
        raise ValueError('pulse leaves the certified faithful neighborhood')
    remainder = 8*parameter**2*v**2+distance**2/(mu-distance)
    slope_upper = response['regional_entropy_slope']+response['entropy_slope_error_upper']
    entropy_upper = parameter*slope_upper+remainder
    return {'parameter': parameter, 'regional_halftrace_distance_upper': distance,
            'entropy_linearization_error_upper': remainder,
            'actual_regional_entropy_change_upper': entropy_upper,
            'fixed_cut_generalized_entropy_change_upper': entropy_upper,
            'independent_positive_Newton_constant_required_for_zero_area_term': True}


def current_filter(time=F(16)):
    """A time-averaged derivative current cannot yield nonzero limiting momentum."""
    if type(time) is not F or time <= 0:
        raise ValueError('positive exact source time required')
    data = state()
    h, e, c = data['energy_integer'], data['density_integer'], data['current_imaginary_integer']
    force = h@c-c@h
    if not np.any(force):
        raise ArithmeticError('local source current unexpectedly became conserved')
    density_bound = F(max(sum(abs(int(x)) for x in row) for row in e), 8)
    nonzero = tuple(int(x) for x in np.argwhere(force != 0)[0])
    return {'time': time, 'density_norm_upper': density_bound,
            'current_norm_upper': data['current_norm_upper'],
            'averaged_current_norm_upper': 2*density_bound/time,
            'averaged_current_derivative_norm_upper': 2*data['current_norm_upper']/time,
            'force_nonzero_matrix_position': nonzero,
            'force_matrix_entry': F(-int(force[nonzero]), 512),
            'averaging_this_derivative_current_has_zero_norm_limit': True,
            'nonzero_spatial_translation_generator_constructed': False}
