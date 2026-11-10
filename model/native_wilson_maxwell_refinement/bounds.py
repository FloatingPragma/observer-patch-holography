"""Rational total-photon-cutoff bounds for the fixed native Wilson cube."""
from fractions import Fraction as F
from math import factorial, isqrt


def rational(value):
    if type(value) not in (F,int):raise TypeError('exact Fraction or integer required')
    return F(value)


def rounded(value, upper=True, digits=70):
    value=rational(value)
    scale=10**digits
    numerator=value.numerator*scale
    integer=-(-numerator//value.denominator) if upper else numerator//value.denominator
    return F(integer,scale)


def root(value, upper=True, digits=75):
    value=rational(value)
    if value<0:raise ValueError('nonnegative rational required')
    scale=10**digits
    integer=isqrt(value.numerator*scale*scale//value.denominator)
    if upper and integer*integer*value.denominator<value.numerator*scale*scale:integer+=1
    return F(integer,scale)


def exp_upper(value, digits=75):
    value=rational(value)
    if value<0:raise ValueError('nonnegative rational required')
    term=total=F(1);n=0;goal=F(1,10**digits)
    while True:
        following=term*value/(n+1);ratio=value/(n+2)
        if ratio<1:
            remainder=following/(1-ratio)
            if remainder<=goal:return rounded(total+remainder,True,digits)
        term=following;total+=term;n+=1


def displacement_radius(g=F(1,64)):
    """|g| sqrt(<vac,Q_e^2 vac>), rounded outward without floating roots."""
    g=rational(g)
    variance=F(1,16)+1/(12*root(6,False))
    return root(g*g*variance)


def coherent_word_error(cap,radius):
    """Finite product of literal truncated displacements versus full product."""
    radius=rational(radius)
    if type(cap) is not int or cap<0 or radius<0:raise ValueError('nonnegative cap/radius required')
    tail=radius**(cap+1)/root(factorial(cap+1),False)
    return min(F(2),tail*(1+root(cap+1)*radius))


def instrument_word_error(cap,radius):
    """Allows interleaved retained number-preserving contractions/isometries."""
    radius=rational(radius)
    if type(cap) is not int or cap<0 or radius<0:raise ValueError('nonnegative cap/radius required')
    tail=radius**(cap+1)/root(factorial(cap+1),False)
    if tail>=2:return F(2)
    return min(F(2),tail*(exp_upper(radius*radius)+root(cap+1)*radius*exp_upper(radius*radius/2)))


def dyson_tail(x,order):
    """Two times the omitted Dyson absolute series, with exact geometric tail."""
    x=rational(x)
    if x<0 or order<0 or x>=order+2:raise ValueError('require 0 <= x < order+2')
    return 2*x**(order+1)/factorial(order+1)/(1-x/F(order+2))


def cutoff_bound(cap,g=F(1,64),tau=F(1,64),order=24,apparatus_radius=F(0),histories=False):
    """Complete state/record comparison; tau sums actual autonomous waits.

    24 oriented matter coefficients each have norm <=2m. The photon-free
    baseline includes every native/complement/reference/Coulomb term.
    """
    g,tau,apparatus_radius=map(rational,(g,tau,apparatus_radius))
    if g<0 or tau<0 or apparatus_radius<0:raise ValueError('nonnegative parameters required')
    x=48*tau;r=displacement_radius(g)
    word=instrument_word_error if histories else coherent_word_error
    if apparatus_radius and not histories:raise ValueError('apparatus displacement needs history bound')
    weight=F(1);partial=F(0)
    for n in range(order+1):
        if n:weight*=x/n
        partial+=weight*word(cap,n*r+apparatus_radius)
    remainder=dyson_tail(x,order)
    error=rounded(min(F(2),partial+remainder))
    return {'cap':cap,'dimensionless_g':str(g),'total_dimensionless_autonomous_time':str(tau),
        'oriented_hopping_norm_sum_over_m':'48','dyson_weight_x':str(x),
        'single_link_displacement_radius_upper':str(r),'dyson_order':order,
        'retained_dyson_word_sum_upper':str(rounded(partial)),
        'omitted_dyson_words_vector_upper':str(rounded(remainder)),
        'apparatus_displacement_radius_upper':str(apparatus_radius),
        'retained_number_preserving_history_bound':histories,
        'state_vector_error_upper':str(error),'complete_record_TV_upper':str(min(F(1),error)),
        'photon_overflow_probability_upper':str(rounded(min(F(1),error*error))),
        'diagnostic_vector_error':float(error)}


def convergent_order(cap,g=F(1,64),apparatus_radius=F(0)):
    """An explicit diagonal schedule for fixed g and apparatus radius.

    For sufficiently large cap, J*r+s <= sqrt(cap+1)/4. Then the
    general-history word bound decays exponentially and the Dyson tail
    factorially. Small cap choices have no asymptotic significance.
    """
    r=displacement_radius(g)
    return max(0,int(root(cap+1,False)/(8*max(F(1),r))))
