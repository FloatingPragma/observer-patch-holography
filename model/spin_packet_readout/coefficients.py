"""Rational enclosures for the actual finite Fourier-LCU coefficients."""
from fractions import Fraction as F
from math import factorial
from .filters import exact, parameters


def mul(a, b):
    values = [x*y for x in a for y in b]
    return min(values), max(values)


def div(a, b):
    if b[0] <= 0 <= b[1]:
        raise ValueError('interval denominator contains zero')
    return mul(a, (1/b[1], 1/b[0]))


def pi_interval(bits):
    if type(bits) is not int or bits < 8:
        raise ValueError('at least eight integer precision bits required')
    def arctan(q):
        total = F()
        n = 0
        while True:
            total += F((-1)**n, (2*n+1)*q**(2*n+1))
            n += 1
            term = F((-1)**n, (2*n+1)*q**(2*n+1))
            if abs(term) <= F(1, 1 << (bits+8)):
                return min(total, total+term), max(total, total+term)
    a, b = arctan(5), arctan(239)
    return 16*a[0]-4*b[1], 16*a[1]-4*b[0]


def trig_interval(x, bits=64):
    x = exact(x)
    if type(bits) is not int or bits < 8:
        raise ValueError('at least eight integer precision bits required')
    extra = max(0, abs(x.numerator).bit_length()-x.denominator.bit_length())+8
    lo, hi = pi_interval(bits+extra)
    center = (lo+hi)/2
    turns = x/(2*center)+F(1, 2)
    turns = turns.numerator//turns.denominator
    reduced = x-2*turns*center
    uncertainty = abs(turns)*(hi-lo)
    if abs(reduced) > 4:
        raise ArithmeticError('rational range reduction failed')
    n = 4
    tolerance = F(1, 1 << (bits+4))
    while max(abs(reduced)**(2*n+2)/factorial(2*n+2),
              abs(reduced)**(2*n+3)/factorial(2*n+3)) > tolerance:
        n += 1
    cosine = sum(((-1)**j*reduced**(2*j)/factorial(2*j) for j in range(n+1)), F())
    sine = sum(((-1)**j*reduced**(2*j+1)/factorial(2*j+1) for j in range(n+1)), F())
    ec = abs(reduced)**(2*n+2)/factorial(2*n+2)+uncertainty
    es = abs(reduced)**(2*n+3)/factorial(2*n+3)+uncertainty
    return (cosine-ec, cosine+ec), (sine-es, sine+es)


def kernel_rectangle(time, a, b, d, bits=64):
    a, b, d = parameters(a, b, d)
    t = exact(time)
    pi = pi_interval(bits+16)
    length = b-a+d
    if t == 0:
        return div((length, length), (2*pi[0], 2*pi[1])), (F(), F())
    _, first = trig_interval(length*t/2, bits)
    _, second = trig_interval(d*t/4, bits)
    amp = mul(div(first, mul(pi, (t, t))),
              mul(div(second, (d*t/4, d*t/4)), div(second, (d*t/4, d*t/4))))
    co, si = trig_interval((a+b)*t/2, bits)
    return mul(amp, co), mul(amp, (-si[1], -si[0]))


def midpoint_coefficient(time, step, a, b, d, error):
    """Return rational complex coordinates within the requested modulus error."""
    step, error = exact(step, positive=True), exact(error, positive=True)
    bits = 16
    while True:
        real, imag = kernel_rectangle(time, a, b, d, bits)
        radius = step*((real[1]-real[0])+(imag[1]-imag[0]))/2
        if radius <= error:
            return {'real': step*(real[0]+real[1])/2,
                    'imag': step*(imag[0]+imag[1])/2,
                    'complex_error_upper': radius, 'precision_bits': bits}
        bits += 16


def sqrt_interval(value, bits=64):
    from math import isqrt
    value = exact(value)
    if value < 0 or type(bits) is not int or bits < 8:
        raise ValueError('nonnegative rational and at least eight precision bits required')
    scale = 1 << bits
    square = value*scale*scale
    lower = isqrt(square.numerator//square.denominator)
    upper = lower+(lower*lower*square.denominator < square.numerator)
    return F(lower, scale), F(upper, scale)


def contractive_coefficients(rectangles, error, bits=64):
    """Round a known l1-contractive complex vector without losing contractivity.

    Input rectangles must enclose that vector; this norm premise is explicit,
    not inferred from possibly overlapping rectangles.
    """
    error = exact(error, positive=True)
    if not rectangles or type(bits) is not int or bits < 8:
        raise ValueError('nonempty coefficient enclosures and precision required')
    centers = []
    uncertainty = F()
    for rectangle in rectangles:
        if len(rectangle) != 2 or any(len(interval) != 2 for interval in rectangle):
            raise ValueError('real and imaginary intervals required')
        real, imag = tuple(tuple(exact(x) for x in interval) for interval in rectangle)
        if real[0] > real[1] or imag[0] > imag[1]:
            raise ValueError('ordered interval endpoints required')
        centers.append(((real[0]+real[1])/2, (imag[0]+imag[1])/2))
        uncertainty += (real[1]-real[0]+imag[1]-imag[0])/2
    if uncertainty > error/4:
        raise ValueError('refine input intervals before coefficient normalization')
    while True:
        radii = [sqrt_interval(r*r+i*i, bits) for r, i in centers]
        radial_error = sum((hi-lo for lo, hi in radii), F())
        if radial_error <= error/4:
            break
        bits += 16
    upper = sum((hi for _, hi in radii), F())
    denominator = max(F(1), upper)
    result = [(real/denominator, imag/denominator) for real, imag in centers]
    return {'coefficients': result,
            'l1_norm_upper': upper/denominator,
            'l1_error_upper': 2*uncertainty+radial_error,
            'input_vector_l1_at_most_one_required': True,
            'precision_bits': bits}
