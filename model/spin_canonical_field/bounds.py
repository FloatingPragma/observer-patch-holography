"""Exact finite constants used in the analytic continuum construction."""
from fractions import Fraction as F
from math import isqrt


def rational(value, name, *, positive=False):
    if type(value) not in (int, F):
        raise ValueError(name + ' must be an exact rational')
    value = F(value)
    if value < 0 or (positive and value == 0):
        raise ValueError(name + ' has an invalid sign')
    return value


def sqrt_upper(value, bits=80):
    value = rational(value, 'radicand')
    if type(bits) is not int or bits < 1:
        raise ValueError('positive integer precision required')
    scale = 1 << bits
    n = isqrt(value.numerator * scale * scale // value.denominator)
    result = F(n, scale)
    return result if result*result == value else F(n+1, scale)


def site_weyl_error(a, b):
    if type(a) not in (int, F) or type(b) not in (int, F):
        raise ValueError('real rational angles required')
    a, b = abs(F(a)), abs(F(b))
    return (a*a*b+a*b*b)/4


def energy_kick(exchange, edge_square_sum):
    """Squared increment of the energy square root, C_A/2."""
    return rational(exchange, 'exchange', positive=True)*rational(edge_square_sum, 'edge square sum')/4


def energy_after(energy, kick_squared):
    energy = rational(energy, 'energy')
    kick_squared = rational(kick_squared, 'squared kick')
    return energy+kick_squared+2*sqrt_upper(energy*kick_squared)


def transport_error_squared(time, energy, kick_squared):
    time = rational(time, 'absolute time')
    energy = rational(energy, 'energy')
    kick_squared = rational(kick_squared, 'squared kick')
    # min(4, ...) is the distance-squared bound for two unit vectors.
    return min(F(4), 4*time*sqrt_upper(energy*kick_squared)+10*time*kick_squared)


def scaling_example(n):
    if type(n) is not int or n < 1:
        raise ValueError('positive integer scale required')
    p = F(1, n*n)
    # J=1, E=p and the macroscopic-order kick has delta^2=p^3.
    raw = 4*p+10*p*p
    bound = transport_error_squared(1/p, p, p**3)
    if bound != min(F(4), raw):
        raise ArithmeticError('scaled transport normalization changed')
    return {'scale': n, 'p': str(p), 'time': str(1/p),
            'history_energy_upper': str(p), 'order_kick_squared_upper': str(p**3),
            'transport_distance_squared_upper': str(bound),
            'untruncated_transport_bound': str(raw),
            'order_self_averaging_rate_claimed': False}


def constants():
    # Uniform transverse density <=|k|/sqrt(20)<|k|/4.
    if not F(1, 20) < F(1, 16):
        raise ArithmeticError('uniform-field bound')
    # Staggered density <=sqrt(6/(2|k|^2/3))=3/|k|.
    return {'magnetization_strict_lower': '1/1000',
            'uniform_second_moment_coefficient': '1/4',
            'staggered_second_moment_coefficient': '3',
            'uniform_energy_coefficient': '1/4',
            'staggered_energy_coefficient': '3',
            'source_L1_cone_speed_upper': '72',
            'one_site_weyl_remainder_coefficient': '1/4',
            'transport_mixed_coefficient': 4,
            'transport_quadratic_coefficient': 10}
