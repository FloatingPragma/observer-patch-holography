"""Rational admissibility proof at the declared decimal source coordinates.

The input decimals are exact conditional model coordinates. These enclosures
do not certify that a physical endpoint selects those coordinates.
"""
from fractions import Fraction as F
from math import factorial

from model.native_flavor_orbit.intervals import I, pi_interval, rounded_out, sqrt_rational
from .construction import source


def rounded(value):
    return rounded_out(value, 28)


def cosine(value):
    if max(abs(value.lo),abs(value.hi)) > 5:
        raise ValueError('cosine enclosure domain is [-5,5]')
    total = sum(((-1)**k*value**(2*k)/factorial(2*k) for k in range(40)),I(0))
    error = F(5)**80/factorial(80)
    return rounded(total+I(-error,error))


def exponential(value):
    if max(abs(value.lo),abs(value.hi)) > 1:
        raise ValueError('exponential enclosure domain is [-1,1]')
    total = sum((value**k/factorial(k) for k in range(30)),I(0))
    error = F(3, factorial(30))  # e^|x| < 3 on this domain.
    return rounded(total+I(-error,error))


def source_shape_boxes(branch):
    pixel,alpha = (I(x) for x in source(branch))
    pi = pi_interval()
    root2,root5 = sqrt_rational(2),sqrt_rational(5)
    # alpha_P^2 avoids an unnecessary interval square root of pi.
    alpha_p_squared = rounded((pixel-(1+root5)/2)**2/pi)
    chi = rounded((alpha**2/512)/(1-alpha/77))
    zeta = rounded((alpha_p_squared/21)/(1-alpha/27-alpha**2/135))
    rho = rounded(root2*exponential(-chi))
    squares = [rounded((1+rho*cosine(rounded(F(2,9)+zeta+2*pi*j/3)))**2)
               for j in range(3)]
    squares.sort(key=lambda x:x.lo)
    if any(a.hi >= b.lo for a,b in zip(squares,squares[1:])):
        raise ArithmeticError('source shape ordering is not certified')
    return squares


def normalized_mass_boxes(branch):
    squares = source_shape_boxes(branch)
    total = sum(squares,I(0))
    return [rounded(x/total) for x in squares]


def mass_norm_bound(branch):
    """Prove sum(m_i/v)<1 without numerical cube-root evaluation."""
    squares = source_shape_boxes(branch)
    alpha = I(source(branch)[1])
    kappa = (alpha/50)/(1+alpha/31+alpha**2/310)
    numerator_cubed = exponential(3*kappa)*sum(squares,I(0))**3
    denominator_cubed = 2*6**14*squares[0]*squares[1]*squares[2]
    cube = numerator_cubed/denominator_cubed
    if not 0 < cube.lo <= cube.hi < 1:
        raise ArithmeticError('MCPR total mass/v bound is not certified')
    return {'sum_mass_over_v_cubed':rounded_out(cube,18).strings(),
            'sum_mass_over_v_strict_upper_bound':'1',
            'operator_norm_mass_matrix_over_v_strict_upper_bound':'1'}


def _sign(value):
    if value.lo > 0:
        return 1
    if value.hi < 0:
        return -1
    raise ArithmeticError('polynomial endpoint has no certified sign')


def _root_box(polynomial,lo,hi):
    lo,hi = F(lo),F(hi)
    left,right = _sign(polynomial(I(lo))),_sign(polynomial(I(hi)))
    if left == right:
        raise ArithmeticError('root interval has no sign change')
    return I(lo,hi)


def admissibility(branch):
    """Three disjoint cubic roots and two disjoint admissible inverse roots.

    Endpoint signs plus the polynomial degrees exhaust their real roots.
    There is no floating point step in this certificate.
    """
    lam = normalized_mass_boxes(branch)
    e2 = sum((lam[i]*lam[j] for i in range(3) for j in range(i)),I(0))
    e3 = lam[0]*lam[1]*lam[2]
    determinant = I(F(2,5))-2*sqrt_rational(5)/25
    cubic = lambda z:z**3-z**2+F(5,4)*e2*z-e3/determinant
    w = _root_box(cubic,'.00099154','.00099157')
    u = _root_box(cubic,'.07058200','.07058206')
    v = _root_box(cubic,'.92842638','.92842646')
    A = v*(v+w)
    B = v*(1-v+w)+(v+w)*w
    quadratic = lambda p:-A*p*p+B*p-w
    inverse = []
    for low,high in (('.01914','.01917'),('.05997','.06001')):
        p = _root_box(quadratic,low,high)
        q = 1-v/(1-p*(v+w))
        if not 0 < p.lo < p.hi < F(1,12) or not 0 < q.lo < q.hi < F(1,6):
            raise ArithmeticError('inverse branch is outside native gate bounds')
        inverse.append({'p':p.strings(),'q':rounded_out(q,12).strings()})
    return {'source_interpretation':'exact supplied decimal P and alpha_U; conditional endpoint',
            'normalized_mass_intervals':[x.strings() for x in lam],
            'projector_weight_intervals':[x.strings() for x in (u,v,w)],
            'inverse_branches':inverse,'absolute_mass_bound':mass_norm_bound(branch),
            'cubic_root_count':3,'inverse_root_count':2,
            'method':'rational Taylor remainders and polynomial endpoint signs with degree exhaustion'}
