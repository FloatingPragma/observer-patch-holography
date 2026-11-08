"""Exact common-ground bounds and a faithful-reference vacuum diagnostic."""
from fractions import Fraction as F
from math import comb, prod


def integer(value, low, name):
    if type(value) is not int or value < low:
        raise ValueError(f'{name} must be an integer >= {low}')
    return value


def fraction(value, name):
    if type(value) is not F:
        raise ValueError(f'{name} must be an exact Fraction')
    return value


def dicke_weights(vertices, occupation, region):
    integer(vertices, 1, 'vertices'); integer(occupation, 0, 'occupation')
    integer(region, 0, 'region')
    if occupation > vertices or region > vertices or region > 64:
        raise ValueError('occupation and region must lie inside the source; at most 64 regional ports')
    denominator=prod(vertices-j for j in range(region))
    return tuple(F(comb(region,k)*prod(occupation-j for j in range(k))
                   *prod(vertices-occupation-j for j in range(region-k)),denominator)
                 if 0 <= occupation-k <= vertices-region else F()
                 for k in range(region+1))


def local_ground_witness(vertices, occupation):
    integer(vertices,2,'vertices');integer(occupation,0,'occupation')
    if occupation>vertices:raise ValueError('occupation outside the source')
    return F(occupation*(vertices-occupation),2*vertices*(vertices-1))


def cooling_length(vertices,gap,seams,error):
    """A sufficient length without constructing binomial(V,q) or a global state."""
    integer(vertices,1,'vertices');integer(seams,1,'seams')
    fraction(gap,'exchange gap');fraction(error,'error')
    if not 0<gap<=2*seams or not 0<error<1:
        raise ValueError('positive transfer gap and error in (0,1) required')
    bits=max(0,error.denominator.bit_length()-error.numerator.bit_length())
    while error*(1<<bits)<1:bits+=1
    while bits and error*(1<<(bits-1))>=1:bits-=1
    lower=F(2*seams,gap)*(vertices+bits)
    return (lower.numerator+lower.denominator-1)//lower.denominator


def limiting_weights(region, density):
    integer(region, 0, 'region'); fraction(density, 'density')
    if not 0 <= density <= 1 or region > 64: raise ValueError('density in [0,1] and at most 64 regional ports required')
    return tuple(comb(region,k)*density**k*(1-density)**(region-k)
                 for k in range(region+1))


def regional_error(vertices, occupation, region, density):
    dicke_weights(vertices,occupation,region); limiting_weights(region,density)
    return min(F(1),F(region*(region-1),2*vertices)
               +region*abs(F(occupation,vertices)-density))


def local_moments(support, degree=6):
    integer(support, 0, 'support'); integer(degree, 0, 'degree')
    return F(2*degree*support),F(4*degree**2*support*(support+1))


def history_budget(supports, times, seams, alpha=F(10,21), degree=6):
    if type(supports) is not tuple or type(times) is not tuple or len(supports)!=len(times) or len(supports)>64:
        raise ValueError('matching tuples of supports and evolution intervals required')
    integer(seams,1,'seams');fraction(alpha,'transfer floor')
    if not 0<alpha<=1:raise ValueError('transfer floor must lie in (0,1]')
    first=F();second=F();error=F();prefixes=[]
    for k,t in zip(supports,times,strict=True):
        fraction(t,'time');c,d=local_moments(k,degree)
        second+=2*c*first+d;first+=c
        error+=abs(t)*second/(4*seams*alpha)
        prefixes.append({'first_energy_moment_norm_upper':str(first),
                         'second_energy_moment_norm_upper':str(second)})
    return {'prefixes':prefixes,'history_vector_norm_error_upper':error,
            'complete_path_probability_error_upper':2*error}


def cooling_error(dimension, gap, seams, length):
    integer(dimension,1,'dimension');integer(seams,1,'seams');integer(length,0,'length')
    fraction(gap,'exchange gap')
    if not 0<gap<=2*seams:raise ValueError('gap outside the positive transfer domain')
    excited=(dimension-1)*(1-gap/(2*seams))**length
    return excited/(1+excited)


def thermal_witness():
    """Exact literal H66 witness; no eigensolver or energy fit."""
    import numpy as np
    from model.primitive_observer import transfer
    source=transfer.native_counts();transfer.exact_audit(source)
    counts=np.asarray(source.counts,dtype=object)
    occupations=np.asarray([state[0] for state in source.states],dtype=object)
    commutator=counts*occupations[None,:]-occupations[:,None]*counts
    moment=F(int(np.trace(counts@commutator.T@commutator)),60**3*56)
    if moment!=F(37,86400):raise ArithmeticError('the literal thermal witness changed')
    negative=F(18,5)*moment
    if negative!=F(37,24000):raise ArithmeticError('the spectral lower bound changed')
    return {'source_dimension':66,'transfer_trace':56,'transfer_floor':'2/3',
            'observable':'occupation n_0; norm one; native number preserving',
            'weighted_transfer_commutator_norm_squared':str(moment),
            'negative_energy_spectral_weight_lower':str(negative),
            'clock':'H=-60 log(P); counted t=n/60',
            'positive_underlying_H_contradicted':False,
            'thermal_QFT_excluded':False,'same_state_positive_energy_vacuum_excluded':True}
