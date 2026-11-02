"""Exact compact energy splines and certified Fourier resource bounds."""
from fractions import Fraction as F


def exact(x, *, positive=False):
    if type(x) not in (int, F):
        raise ValueError('an exact rational, not a Boolean or float, is required')
    x = F(x)
    if positive and x <= 0:
        raise ValueError('a positive exact rational is required')
    return x


def parameters(a, b, d):
    a, b, d = (exact(x, positive=True) for x in (a, b, d))
    if not d < a <= b:
        raise ValueError('require 0<d<a<=b')
    return a, b, d


def triangular_cdf(x, d):
    x, d = exact(x), exact(d, positive=True)
    if x <= -d/2:
        return F()
    if x <= 0:
        return 2*(x+d/2)**2/d**2
    if x < d/2:
        return 1-2*(d/2-x)**2/d**2
    return F(1)


def energy_filter(energy, a, b, d):
    a, b, d = parameters(a, b, d)
    e = exact(energy)
    return triangular_cdf(e-a+d/2, d)-triangular_cdf(e-b-d/2, d)


def fourier_bounds(a, b, d):
    a, b, d = parameters(a, b, d)
    length = b-a+d
    ratio = 4*length/d
    power = 0
    while 2**power < ratio:
        power += 1
    return {'support_lower': a-d, 'support_upper': b+d,
            'plateau_lower': a, 'plateau_upper': b,
            'box_length': length, 'dyadic_log_ceiling': power,
            'kernel_L1_upper': F(2+2*power, 3),
            'kernel_first_absolute_moment_upper': 16/(3*d),
            'kernel_derivative_upper': (b+d)*length/6,
            'tail_numerator': 16/(3*d*d),
            'scalar_Lipschitz_upper': 2/d}


def tail_upper(a, b, d, time_cutoff):
    t = exact(time_cutoff, positive=True)
    return min(fourier_bounds(a, b, d)['kernel_L1_upper'],
               fourier_bounds(a, b, d)['tail_numerator']/t**2)


def goldstone_filter(m=F(1, 1000)):
    m = exact(m, positive=True)
    if m > F(1, 2):
        raise ValueError('spin-half magnetization cannot exceed one half')
    return fourier_bounds(m/8, 32/m**2, m/16)


def encode(value):
    if type(value) is F:
        return str(value)
    if type(value) is dict:
        return {k: encode(v) for k, v in value.items()}
    if type(value) in (tuple, list):
        return [encode(v) for v in value]
    return value
