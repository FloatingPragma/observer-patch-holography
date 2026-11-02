"""Exact constants and finite bounds; the infinite-lattice proof is in README."""
from fractions import Fraction as F
from math import isqrt


def rational(value, name, *, zero=False):
    if isinstance(value, bool) or not isinstance(value, (int, F)):
        raise ValueError(name + ' must be an exact rational')
    value = F(value)
    if value < 0 or (value == 0 and not zero):
        raise ValueError(name + ' has an invalid sign')
    return value


def sqrt_upper(value, bits=80):
    value = rational(value, 'radicand', zero=True)
    if type(bits) is not int or bits < 1:
        raise ValueError('positive integer precision required')
    scale = 1 << bits
    n = isqrt(value.numerator * scale * scale // value.denominator)
    lower = F(n, scale)
    return lower if lower * lower == value else F(n + 1, scale)


def constants():
    # Per-row weighted Parseval, including all three second derivatives.
    moment = 1 + 6 + 9
    # Infinity shells: inner cube <=27 L^3, exterior <=52 L^3.
    shell = 27 + 52
    schur = 64
    assert 2 * moment * shell < schur**2
    # sqrt(11/8) <6/5, using the centered-ground occupation certificate.
    assert F(11, 8) < F(6, 5)**2
    microscopic = 3 + F(6, 5) * 312 * schur
    assert microscopic < 2**15
    # E(1+|Z|)^2=2+2 sqrt(2/pi) <4.
    return {'weighted_row_square': moment, 'lattice_weight_sum': shell,
            'schur_coefficient': schur, 'schur_time_power': '3/2',
            'microscopic_error_coefficient': 2**15,
            'microscopic_time_power': '3/2',
            'gaussian_moment_coefficient': 4,
            'smoothed_error_coefficient': 2**17,
            'free_phase_space_coefficient': 1536}


def spectral_bound(c, coupling, radius, width_root):
    """Upper bound for each unit-smearing composite low-energy weight.

    epsilon=width_root**2. The trivial upper bound is used if the simple
    rational exponential majorant 1/(1-x) is unavailable. No spectral atom
    or interacting phase is inferred from a failure to obtain a small bound.
    """
    c = rational(c, 'kinetic scale')
    coupling = rational(coupling, 'contact coupling', zero=True)
    radius = rational(radius, 'energy radius', zero=True)
    q = rational(width_root, 'square root of Gaussian width')
    if q > 1:
        raise ValueError('Gaussian width must be at most one')
    k = constants()
    epsilon = q*q
    exponent = radius**2 / (2*c*c*epsilon**2)
    terms = (k['free_phase_space_coefficient'] * q**12
             + k['smoothed_error_coefficient'] * (coupling/c) / q**3)
    bound = F(1) if exponent >= 1 else min(F(1), terms / (1-exponent))
    return {'epsilon': epsilon, 'exponential_argument': exponent,
            'unmultiplied_bound': terms, 'weight_upper': bound,
            'coupling_squared_weight_upper': coupling**2 * bound}


def field_error(c, coupling, radius, width_root, *, time=F(1), momentum=F(1)):
    radius = rational(radius, 'positive splitting radius')
    time = rational(time, 'absolute time', zero=True)
    momentum = rational(momentum, 'kinetic window bound', zero=True)
    row = spectral_bound(c, coupling, radius, width_root)
    # Take the root after multiplication by U^2: an absolute root error
    # must not be amplified by a possibly enormous growing coupling.
    error = (time*sqrt_upper(row['coupling_squared_weight_upper'])
             + F(coupling)*(2+time*momentum)/radius)/2
    return min(F(2), error)


def power_exponents(alpha):
    if isinstance(alpha, bool) or not isinstance(alpha, (int, F)):
        raise ValueError('exact rational power required')
    alpha = F(alpha)
    exponents = {'weighted_composite': (14*alpha-4)/5,
                 'width_floor': 2*alpha-F(3, 2),
                 'high_energy_tail': alpha-F(1, 2)}
    return exponents, all(x < 0 for x in exponents.values())


def example(base):
    if type(base) is not int or base < 2:
        raise ValueError('integer base at least two required')
    # c=b^80, U=b^20, R=b^40, epsilon=b^-8. Thus U=c^(1/4).
    c, coupling, radius, q = F(base**80), F(base**20), F(base**40), F(1, base**4)
    row = spectral_bound(c, coupling, radius, q)
    # This independent closed form checks every normalization and power.
    direct = F(1536+2**17, base**8) / (1-F(1, 2*base**64))
    assert row['coupling_squared_weight_upper'] == min(coupling**2, direct)
    return {'base': base, 'c': str(c), 'U': str(coupling), 'R': str(radius),
            **{key: str(value) for key, value in row.items()},
            'unit_time_unit_window_field_error_upper':
                str(field_error(c, coupling, radius, q))}
