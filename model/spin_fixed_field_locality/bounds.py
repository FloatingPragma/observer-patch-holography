"""Rational certificates for the positive-field interior marginal modulus."""
from fractions import Fraction as F


def positive(value,name):
    if type(value) not in (int,F) or value<=0:raise ValueError(name+' must be positive rational')
    return F(value)


def marginal(coupling,field,distance,sites,crossing_edges):
    j=positive(coupling,'coupling');b=positive(field,'field')
    if type(distance)is not int or distance<1:raise ValueError('positive integer distance required')
    if type(sites)is not int or sites<1:raise ValueError('positive observation site count required')
    if type(crossing_edges)is not int or crossing_edges<0:raise ValueError('nonnegative crossing-edge count required')
    g=b/8;t=F(distance)/(144*j);alpha=g/(2*t)
    rate=min(F(1,3),g/(288*j));exponent=(rate*distance).numerator//(rate*distance).denominator
    # pi>=3, sqrt(x)<=1+x and e>=2 give a rational upper certificate.
    prefactor=j*crossing_edges*(sites*t/2+1+2*t/(3*g)+1/g)
    upper=prefactor/F(2**exponent)
    return {'gap_lower':g,'time_cutoff':t,'gaussian_alpha':alpha,
        'common_decay_rate':rate,'binary_decay_exponent':exponent,
        'rational_prefactor':prefactor,'half_trace_upper':min(F(1),upper),
        'uncapped_half_trace_upper':upper,'exterior_volume_in_bound':False,
        'zero_field_modulus_claimed':False}


def box(coupling,field,inner_radius,outer_radius,observed_sites):
    if type(inner_radius)is not int or inner_radius<0:raise ValueError('nonnegative observation radius required')
    if type(outer_radius)is not int or outer_radius<=inner_radius:raise ValueError('larger positive source box required')
    if type(observed_sites)is not int or not 1<=observed_sites<=(2*inner_radius+1)**3:raise ValueError('observation support exceeds its box')
    value=marginal(coupling,field,outer_radius-inner_radius,observed_sites,6*(2*outer_radius+1)**2)
    return {'observation_radius':inner_radius,'source_radius':outer_radius,
        'source_sites':(2*outer_radius+1)**3,'crossing_edges_upper':6*(2*outer_radius+1)**2,**value}


def choose_box(coupling,field,inner_radius,observed_sites,accuracy):
    eps=positive(accuracy,'accuracy')
    if eps>=1:raise ValueError('accuracy must be less than one')
    radius=max(1,inner_radius+1);attempts=0
    while True:
        row=box(coupling,field,inner_radius,radius,observed_sites);attempts+=1
        if row['half_trace_upper']<=eps:return {**row,'target_accuracy':eps,'doubling_steps':attempts,'state_diagonalization_used':False}
        radius*=2
