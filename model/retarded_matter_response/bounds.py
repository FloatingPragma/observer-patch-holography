"""Rational finite-frequency error bounds for the exact ground response.

The charged-sector gap comes from the retained original full-Fock ground
certificate. It is not a fitted pole, thermodynamic gap or continuum mass.
"""
from fractions import Fraction as F

from .response import GAP


def _rational(value, name, positive=False):
    if type(value) not in (int, F) or value < 0 or (positive and value == 0):
        raise ValueError('exact '+('positive ' if positive else 'nonnegative ')+name+' required')
    return F(value)


def _order(order):
    if type(order) is not int or not 0 <= order <= 32:
        raise ValueError('integer expansion order0..32 required')
    return order


def memory_taylor_remainder(radius, order):
    """For |z|<=radius<delta, after powers0..order of the exact memory.

    Sigma(z)=-sum z^n W^*B^(-n-1)W, B=QKQ|RanQ>=delta.
    ||W||^2=1/4, including every eliminated state of this source sector.
    """
    radius = _rational(radius, 'frequency radius'); order = _order(order)
    if radius >= GAP: raise ValueError('frequency radius must be below7/1024')
    return radius**(order+1)/(4*GAP**(order+1)*(GAP-radius))


def response_taylor_remainder(radius, order):
    """Same disk, expansion -sum z^n V^*K^(-n-1)V of G itself."""
    return 4*memory_taylor_remainder(radius, order)


def large_frequency_remainder(modulus_lower, order):
    """G=sum_(n=0)^order M_n/z^(n+1)+R, ||R|| bounded below.

    On this560-state ground-down sector0<K<2; M0..M4 are computed exactly.
    The analytic bound is valid at every complex |z|>2, not only Imz>0.
    """
    modulus = _rational(modulus_lower, 'frequency modulus', True); order = _order(order)
    if modulus <= 2: raise ValueError('large-frequency modulus must exceed2')
    return F(2)**(order+1)/(modulus**(order+1)*(modulus-2))


def certificate():
    # lambda=sqrt3/16, ||K||<=1+9lambda<2 since243<256.
    # N5 charged gap from original global ground is strictly7/1024.
    if not 9**2*3 < 16**2: raise ArithmeticError('strict kinetic upper bound failed')
    return {'channel':'ground down only', 'strict_generator_interval':['7/1024','2'],
            'upper_bound_proof':'D<=1 on Nup4,Ndown1; kinetic norm<=5lambda; -E_G=4lambda; 1+9lambda<2',
            'memory_operator_squared_norm':'1/4',
            'compressed_generator':'B=QKQ restricted to RanQ; same lower and upper bounds',
            'low_frequency_coefficients':'C_n=W^*B^(-n-1)W; Sigma(z)=-sum_(n>=0)z^n C_n',
            'coefficient_order':'0<=C_(n+1)<=C_n/delta as quadratic forms',
            'low_frequency_remainder':'r^(p+1)/[4 delta^(p+1)(delta-r)] for|z|<=r<delta',
            'low_frequency_examples':[{'radius':str(GAP/d),'order':p,
                                        'memory_remainder':str(memory_taylor_remainder(GAP/d,p)),
                                        'response_remainder':str(response_taylor_remainder(GAP/d,p))}
                                       for d,p in ((2,4),(4,4),(4,8))],
            'large_frequency_remainder':'2^(p+1)/[R^(p+1)(R-2)] for|z|>=R>2',
            'large_frequency_examples':[{'modulus':str(r),'order':4,
                                          'remainder':str(large_frequency_remainder(r,4))}
                                         for r in (4,8,16)],
            'low_frequency_coefficients_defined_by_exact_source_matrices':True,
            'low_frequency_inverse_matrices_numerically_evaluated':False,
            'spatial_long_wavelength_or_continuum_limit_proved':False,
            'gap_extended_to_excited_reference_or_infinite_volume':False}
