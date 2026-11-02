"""Exact rational current fractions and fractional-energy tail controls."""
from fractions import Fraction as F


def rational(x, name):
    if isinstance(x,bool) or not isinstance(x,(int,F)):
        raise ValueError(name+' must be an exact integer or Fraction')
    return F(x)


def soft_fraction(m0,kappa,radius=1):
    m,k,r=[rational(x,n) for x,n in zip((m0,kappa,radius),('m0','kappa','radius'))]
    if not 0<m<=F(1,2) or not m*m/2<=k<=F(1,2) or r<=0:
        raise ValueError('invalid order, kappa or momentum radius')
    a,b=m*m*r/8,24*r/(m*m)
    weight=m*m*r/48
    return {'field_band_lower':a,'field_band_upper':b,
            'open_source_band_lower':a/2,'open_source_band_upper':2*b,
            'field_weight_lower':weight,
            'soft_first_moment_lower':a*weight,
            'total_first_moment_upper':k*r*r/2,
            'soft_current_fraction_lower':m**4/(192*k),
            'universal_fraction_lower':m**4/96}


def fractional_tail(first,degree,scale):
    """For alpha=1-1/degree, cutoff scale**degree: tail <= first/scale."""
    d=rational(first,'first')
    if type(degree) is not int or degree<2 or type(scale) is not int or scale<1 or d<0:
        raise ValueError('invalid fractional tail parameters')
    return {'alpha':F(degree-1,degree),'cutoff':F(scale**degree),
            'fractional_tail_upper':d/scale}


def escaping_current(n):
    """Positive measures obeying the P moment bounds; not a Heisenberg state.

    Low part uniform on [1,2], mass 1/24; high atom at n, mass 1/(16n).
    First moment is exactly 1/8; its escaping half is always 1/16.
    """
    if type(n) is not int or n<3:raise ValueError('n must be an integer at least three')
    low=F(1,24);high=F(1,16*n)
    inverse_upper=low+high/n
    assert inverse_upper<=F(1,12) and low+high<=F(1,4)
    return {'low_continuous_mass':low,'high_energy':F(n),'high_mass':high,
            'first_moment':F(1,8),'escaping_first_moment':F(1,16),
            'inverse_moment_upper':inverse_upper,'total_mass':low+high,
            'weak_limit_has_a_positive_energy_atom':False,
            'is_a_Heisenberg_ground_state_construction':False}


def serialize(value):
    if isinstance(value,F):return str(value)
    if isinstance(value,dict):return {k:serialize(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [serialize(v) for v in value]
    return value
