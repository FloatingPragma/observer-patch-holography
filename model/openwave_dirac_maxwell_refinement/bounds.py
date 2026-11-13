"""Rational weighted-flux and functional-calculus error bounds."""
from fractions import Fraction as F
from .operators import ELECTRIC, MAGNETIC, SPEED_SQUARED, root_tail


def dyadic_upper(value, bits=100):
    value = F(value)
    denominator = 1 << bits
    return F(-((-value.numerator*denominator)//value.denominator),denominator)


def exp_upper(x, terms=80):
    """An exact rational upper bound; no floating exponential is an input."""
    x = F(x)
    if x < 0:
        # exp(x) <= 1 / partial_sum(exp(-x)), all omitted terms positive.
        term = total = F(1)
        for k in range(1, terms+1):
            term *= -x/k
            total += term
        return dyadic_upper(1/total)
    if x >= terms+2:
        raise ValueError('increase Taylor order for this exponent')
    term = total = F(1)
    for k in range(1, terms+1):
        term *= x/k
        total += term
    next_term = term*x/(terms+1)
    return dyadic_upper(total+next_term/(1-x/F(terms+2)))


def infinite_flux_bound(cap, charge, mass, time, degree=6, regulator=F(1,16), weight=F(1,2)):
    """Source state error before finite eigensystem/rounding errors.

    The electric diagonal is unbounded and treated exactly. Its interaction
    picture preserves the weighted operator norms used here.
    """
    m, t, r, s = F(mass), abs(F(time)), F(regulator), F(weight)
    if cap < 1 or m < 1 or r < 0 or s <= 0:
        raise ValueError('positive domain and charged family m>=1 required')
    rho = SPEED_SQUARED/(m*m)
    rho_s = rho*exp_upper(2*s*abs(charge))
    if rho_s >= 1:
        raise ValueError('weighted square-root series not certified in this domain')
    b_s = dyadic_upper(m*rho_s/(1-rho_s)+MAGNETIC*exp_upper(s)+r*exp_upper(s*abs(charge)))
    b_0 = m*rho/(1-rho)+MAGNETIC+r
    tail = exp_upper(-s*(cap+1)+t*b_s)
    galerkin = (1+t*b_0)*tail
    polynomial = t*root_tail(m,degree)
    return {'weight_exact': str(s), 'weighted_root_ratio_upper_exact': str(rho_s),
            'weighted_bounded_generator_norm_upper_exact': str(b_s),
            'bounded_generator_norm_upper_exact': str(b_0),
            'full_flux_tail_norm_upper_exact': str(tail),
            'full_to_Galerkin_state_norm_upper_exact': str(galerkin),
            'root_polynomial_state_norm_upper_exact': str(polynomial),
            'state_norm_error_before_numerical_diagonalization_exact': str(galerkin+polynomial),
            'finite_numerical_errors_included': False}
