"""Deterministic rational reconstruction on the certified eight-line support.

No floating inverse is trusted. Rational Lagrange polynomials propose a
96-bit inverse; outward dyadic arithmetic certifies its residual against
the exact source frequency intervals. Counts, rather than supplied spectral
weights, determine the reconstructed amplitudes.
"""
from fractions import Fraction as F
from functools import lru_cache

from ..certified_reference import interval as iv
from ..particle_identification import vacuum

BITS = 96
SCALE = 1 << BITS
TIMES = tuple(39*n for n in range(8))
ZERO = (F(), F())
ONE = (F(1), F())


def add(a, b): return a[0]+b[0], a[1]+b[1]
def neg(a): return -a[0], -a[1]
def mul(a, b): return a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]
def div(a, b):
    d = b[0]**2+b[1]**2
    if d == 0: raise ValueError('distinct rational sampling nodes required')
    c = mul(a, (b[0], -b[1]))
    return c[0]/d, c[1]/d


def dyadic(value):
    scaled = F(value)*SCALE+F(1, 2)
    return F(scaled.numerator//scaled.denominator, SCALE)


def trig_interval(low, high):
    """Complex exp(-i*x), enclosing every x in a rational interval."""
    if low > high: raise ValueError('ordered exact angle interval required')
    interval = iv.Interval(iv.Interval.point(low).lo, iv.Interval.point(high).hi)
    pi = iv.pi()
    midpoint = F(interval.lo+interval.hi, 2*iv.SCALE)
    pi_midpoint = F(pi.lo+pi.hi, 2*iv.SCALE)
    quotient = midpoint/(2*pi_midpoint)+F(1, 2)
    turns = quotient.numerator//quotient.denominator
    reduced = interval-2*turns*pi
    center = F(reduced.lo+reduced.hi, 2*iv.SCALE)
    radius = F(reduced.hi-reduced.lo, 2*iv.SCALE)
    if abs(center)+radius > 4: raise ValueError('reduced angle outside certified Taylor domain')
    sine, cosine = iv.sine_cosine(center)
    error = iv.Interval.point(radius).hi
    sine = iv.Interval(sine.lo-error, sine.hi+error)
    cosine = iv.Interval(cosine.lo-error, cosine.hi+error)
    return cosine, -sine


def interval_product(point, rectangle):
    a, b = point; x, y = rectangle
    return a*x-b*y, a*y+b*x


def norm_upper(real, imag):
    bound = iv.sqrt_fraction(real.absolute_upper()**2+imag.absolute_upper()**2)
    return F(bound.hi, iv.SCALE)


@lru_cache(maxsize=1)
def source_rows():
    return tuple(sorted(vacuum.certificate()['spectral_rows'],
                        key=lambda row: F(row['omega_interval'][0])))


@lru_cache(maxsize=1)
def certified_inverse():
    frequencies = [tuple(map(F, row['omega_interval'])) for row in source_rows()]
    nodes = []
    for low, high in frequencies:
        x, y = trig_interval(39*low, 39*high)
        nodes.append((dyadic(F(x.lo+x.hi, 2*iv.SCALE)),
                      dyadic(F(y.lo+y.hi, 2*iv.SCALE))))
    inverse = []
    for j, node in enumerate(nodes):
        coefficients = [ONE]; denominator = ONE
        for k, other in enumerate(nodes):
            if k == j: continue
            polynomial = [ZERO]*(len(coefficients)+1)
            for n, coefficient in enumerate(coefficients):
                polynomial[n] = add(polynomial[n], mul(neg(other), coefficient))
                polynomial[n+1] = add(polynomial[n+1], coefficient)
            coefficients = polynomial
            denominator = mul(denominator, add(node, neg(other)))
        inverse.append(tuple(tuple(dyadic(x) for x in div(coefficient, denominator))
                             for coefficient in coefficients))
    inverse = tuple(inverse)
    matrix = [[trig_interval(t*low, t*high) for low, high in frequencies] for t in TIMES]
    row_residuals = []
    for j, row in enumerate(inverse):
        bounds = []
        for k in range(8):
            real = iv.Interval.point(-int(j == k)); imag = iv.Interval.point(0)
            for n, coefficient in enumerate(row):
                x, y = interval_product(coefficient, matrix[n][k]); real += x; imag += y
            bounds.append(norm_upper(real, imag))
        row_residuals.append(sum(bounds, F()))
    residual = max(row_residuals)
    inverse_norm = max(sum((F(iv.sqrt_fraction(x*x+y*y).hi, iv.SCALE) for x, y in row), F())
                       for row in inverse)
    if residual >= F(1, 10**18) or inverse_norm >= 2:
        raise ArithmeticError('rational inverse lost its certified residual or norm margin')
    return inverse, inverse_norm, residual


def matrix_certificate():
    inverse, norm, residual = certified_inverse()
    return {'times': list(TIMES), 'source_support_count': 8,
            'construction': 'exact rational Lagrange inverse at dyadic nodes, rounded to96bits and certified against exact source frequencies',
            'inverse_entries_real_imag': [[[str(x),str(y)] for x,y in row] for row in inverse],
            'inverse_infinity_norm_upper': str(norm),
            'inverse_infinity_norm_simple_upper': '2',
            'inverse_times_true_V_minus_identity_infinity_norm_upper': str(residual),
            'residual_simple_upper': '1/1000000000000000000',
            'frequency_intervals': [row['omega_interval'] for row in source_rows()],
            'source_weights_used_to_construct_inverse': False,
            'floating_linear_solver_used': False,
            'unknown_frequency_support_inferred_from_counts': False}


def reconstruct(counts_zero, counts_quadrature, trials, systematic_error, sampling_error):
    if type(trials) is not int or trials <= 0: raise ValueError('positive integer trial count required')
    for counts in (counts_zero, counts_quadrature):
        if len(counts) != 8 or any(type(n) is not int or not 0 <= n <= trials for n in counts):
            raise ValueError('eight valid integer counts per pointer setting required')
    for error in (systematic_error, sampling_error):
        if type(error) not in (int, F) or error < 0: raise ValueError('nonnegative exact probability errors required')
    inverse, norm, residual = certified_inverse()
    observations = [((32*F(a,trials)-17)/15, (16-32*F(b,trials))/15)
                    for a,b in zip(counts_zero, counts_quadrature, strict=True)]
    # Complex quadrature noise <=2sqrt(2)e<3e; elastic removal adds16/15.
    error = norm*F(16,5)*(systematic_error+sampling_error)+residual
    estimates = []
    for index, row in enumerate(inverse):
        total = ZERO
        for coefficient, observation in zip(row, observations, strict=True):
            total = add(total, mul(coefficient, observation))
        low, high = total[0]-error, total[0]+error
        estimates.append({'line': index+1, 'real_estimate': str(total[0]),
                          'imaginary_residual': str(total[1]),
                          'weight_interval': [str(low),str(high)],
                          'strictly_positive': low > 0})
    return {'trials_per_setting': trials,
            'counts_phase_zero': list(counts_zero), 'counts_phase_pi_over_two': list(counts_quadrature),
            'total_probability_error_budget': str(systematic_error+sampling_error),
            'weight_error_upper': str(error), 'weights': estimates,
            'all_eight_weights_strictly_positive': all(row['strictly_positive'] for row in estimates),
            'source_weights_used_in_reconstruction': False,
            'coverage_condition': 'the stated simultaneous sampling event and systematic probability bounds hold',
            'identification_domain': 'the eight source-specified frequencies, not arbitrary unknown spectral support'}
