"""Complete finite band filtering and polynomial polar correction."""
from functools import lru_cache
import sympy as s
from model.openwave_neutral_stationary.isolation import data
from model.openwave_neutral_stationary.filter import smallest_power


@lru_cache(maxsize=1)
def parameters():
    d=data();epsilon=s.Rational(1,10**16)
    a,b=d['band'];c=d['outside_lower'];upper=d['norm']
    za=(upper+c-2*a)/(upper-c);zb=(upper+c-2*b)/(upper-c)
    root=s.Integer(1)
    while root**2>zb**2-1:root/=2
    for power in (1,2,4,8):
        stop=epsilon/(4*power);degree=smallest_power(1+root,2/stop)
        loss=degree*(za-zb)/root
        if loss**power<=epsilon:break
    assert power*stop/(1-power*stop)<=epsilon
    eta2=d['eta']**2
    q0=(1-eta2)*(1-epsilon)**2
    alpha=1-q0
    # B1(A)=I+(I-A)/2 underestimates A^(-1/2), with all
    # omitted positive binomial coefficients bounded by one.
    tail=alpha**2/(1-alpha)
    success=q0*(1-tail)**2
    subspace=epsilon/q0 # 1/sqrt(q0)<=1/q0
    trace=4*tail/(1-tail)**2+2*subspace
    plain_trace=4*eta2/(1-eta2)+8*epsilon/(1-eta2)
    assert trace<s.Rational(1,10**14) and success>1-s.Rational(4,10**12)
    return {'epsilon':epsilon,'a':a,'b':b,'c':c,'upper':upper,'za':za,'zb':zb,
        'root_lower':root,'degree':degree,'flattening_power':power,'polynomial_degree':degree*power,
        'pass_loss':loss,'stop':stop,'q0':q0,'alpha':alpha,'polar_tail':tail,
        'success_lower':success,'trace_bound':trace,'plain_stationarity_trace_bound':plain_trace,
        'operator_evaluations':3*degree*power}


def certificate():
    p=parameters();d=data()
    return {'polynomial':'M=1-(1-f(H))^K, f=(-1)^m T_m(x(H))/T_m(za), on the physical L7 B0 sector and zero outside',
        'filter_parameters':{key:str(p[key]) for key in ('a','b','c','upper','za','zb','root_lower')},
        'Chebyshev_degree':p['degree'],'flattening_power':p['flattening_power'],
        'band_filter_degree':p['polynomial_degree'],'operator_error_upper':str(p['epsilon']),
        'finite_pass_loss':str(p['pass_loss']),'finite_stop_amplitude':str(p['stop']),
        'all_eight_momenta_retained':True,
        'complete_plain_filter_Kraus_pair':['M','sqrt(I-M^2)'],
        'plain_local_product_filter_success_lower':str(p['q0']),
        'plain_uniform_in_time_stationarity_defect_upper':str(p['plain_stationarity_trace_bound']),
        'atomic_compression':'A=P M^2 P on Ran(P), q0 I<=A<=I',
        'polar_polynomial':'B1(A)=(3I-A)/2',
        'complete_corrected_success_operator':'K=sqrt(q0) M P B1(A) P, extended by zero on I-P',
        'complete_corrected_failure_operator':'sqrt(I-K^dagger K)',
        'exact_success_contraction':'K^dagger K=q0 B1(A) A B1(A)<=q0 P<=I',
        'polar_tail_upper':str(p['polar_tail']),
        'corrected_original_input_success_interval':[str(p['success_lower']),str(p['q0'])],
        'corrected_original_input_success_decimal_lower':float(p['success_lower']),
        'conditional_reduced_density_trace_error_to_E_over_1024_upper':str(p['trace_bound']),
        'conditional_reduced_density_trace_error_decimal_upper':float(p['trace_bound']),
        'uniform_all_time_source_stationarity_defect_upper':str(2*p['trace_bound']),
        'unconditioned_reference_uniform_complete_Kraus_instrument':True,
        'stationarity_claim_about_arbitrary_joint_reference':False,
        'source_H_evaluations_per_success_vector':p['operator_evaluations'],
        'atomic_P_services_per_success_vector':3,
        'exact_sum_of_three_polynomial_terms_is_declared_operator_service':True,
        'Stinespring_gate_count_or_runtime_derived':False,
        'physical_band_leakage_probability_upper':str((p['epsilon']/p['q0'])**2),
        'eight_way_momentum_readout':'all eight outcomes retained; the reduced prepared density is translation invariant, and each has probability1/8 by the exact plane/link symmetries',
        'zero_momentum_original_input_success_lower':str(p['success_lower']/8),
        'source_band_isolation_and_ground_coefficients_unchanged':True}
