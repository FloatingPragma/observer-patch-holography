"""Explicit Chebyshev band filter with complete finite-carrier error accounting."""
from functools import lru_cache
import sympy as s
from .isolation import data,DIMENSION


def smallest_power(base,target):
    """Smallest integer m with base^m >= target, using only integer arithmetic."""
    base,target=map(s.Rational,(base,target))
    assert base>1 and target>1
    def reached(m):return int(base.p)**m*int(target.q)>=int(base.q)**m*int(target.p)
    hi=1
    while not reached(hi):hi*=2
    lo=0
    while hi-lo>1:
        mid=(lo+hi)//2
        if reached(mid):hi=mid
        else:lo=mid
    assert reached(hi) and not reached(hi-1)
    return hi


@lru_cache(maxsize=1)
def design():
    d=data();rank=128;target=s.Rational(1,10**12)
    epsilon=s.Rational(1,10)
    while (4*epsilon+2*s.Rational(DIMENSION,rank)*epsilon**2)/(1-epsilon)**2>target:epsilon/=10
    a,b=d['band'];c=d['outside_lower'];upper=d['norm']
    za=(upper+c-2*a)/(upper-c);zb=(upper+c-2*b)/(upper-c)
    assert za>=zb>1
    root_lower=s.Integer(1)
    while root_lower**2>zb**2-1:root_lower/=2
    assert root_lower**2<=zb**2-1<4*root_lower**2
    for power in (1,2,4,8,16,32,64,128,256):
        stop=epsilon/(4*power)
        degree=smallest_power(1+root_lower,2/stop)
        pass_loss=degree*(za-zb)/root_lower
        if pass_loss<1 and pass_loss**power<=epsilon:break
    else:raise ArithmeticError('finite filter flattening budget exceeded')
    assert power*stop<1 and power*stop/(1-power*stop)<=epsilon
    trace_bound=(4*epsilon+2*s.Rational(DIMENSION,rank)*epsilon**2)/(1-epsilon)**2
    lower_success=s.Rational(rank,DIMENSION)*(1-epsilon)**2
    upper_success=s.Rational(rank,DIMENSION)+epsilon**2
    wrong_band=s.Rational(DIMENSION,rank)*epsilon**2/(1-epsilon)**2
    return {'a':a,'b':b,'c':c,'upper':upper,'za':za,'zb':zb,'root_lower':root_lower,
        'degree':degree,'flattening_power':power,'polynomial_degree':degree*power,
        'stop_amplitude':stop,'pass_loss':pass_loss,'epsilon':epsilon,
        'target_trace_error':target,'conditional_trace_error':trace_bound,
        'preparation_success_interval':(lower_success,upper_success),
        'outside_band_probability_upper':wrong_band,'rank':rank}


def certificate():
    d=design();p0=s.Rational(d['rank'],DIMENSION)
    return {'admissibility_projector':'joint physical Gauss, raw-L7, raw-B3=0 and periodic momentum zero; all other outcomes retained',
        'filter_parameters':{key:str(d[key]) for key in ('a','b','c','upper','za','zb','root_lower')},
        'Chebyshev_argument':'x(E)=(2E-upper-c)/(upper-c)',
        'normalized_polynomial':'f(E)=(-1)^m T_m(x(E))/T_m(za)',
        'flattened_polynomial':'M(E)=1-(1-f(E))^K inside the admissibility projector; zero outside',
        'Chebyshev_degree':d['degree'],'flattening_power':d['flattening_power'],
        'total_polynomial_degree':d['polynomial_degree'],
        'source_H_applications_for_one_polynomial_vector_evaluation':d['polynomial_degree'],
        'polynomial_evaluation_count_is_primitive_unitary_compilation_cost':False,
        'stop_band_amplitude_upper':str(d['stop_amplitude']),
        'unflattened_pass_band_error_upper':str(d['pass_loss']),
        'operator_error_to_known_band_projector_upper':str(d['epsilon']),
        'complete_Kraus_pair':['M','sqrt(I-M^2)'],
        'complete_instrument_reference_uniform_unhalved_trace_error_bound':'2(epsilon+sqrt(2epsilon))',
        'full_carrier_contraction_and_completeness_proved':True,
        'initial_full_source_density':'I/D',
        'ideal_zero_momentum_band_filter_success':str(p0),
        'polynomial_filter_original_input_success_interval':[str(x) for x in d['preparation_success_interval']],
        'conditional_trace_distance_to_uniform_known_band_upper':str(d['conditional_trace_error']),
        'conditional_outside_band_probability_upper':str(d['outside_band_probability_upper']),
        'accepted_density':'Pi_admissible M(H)^2 / Tr(Pi_admissible M(H)^2)',
        'accepted_density_exactly_stationary_under_actual_H':True,
        'failure_density_exactly_stationary_under_actual_H':True,
        'unknown_lowest_eigenprojector_required_by_polynomial':False,
        'full_interacting_source_matrix_materialized':False,
        'Stinespring_dilation_is_declared_finite_apparatus_not_primitive_compilation':True}
