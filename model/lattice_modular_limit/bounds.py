"""Exact finite-residual error budgets, with no endpoint spectral gap."""
from fractions import Fraction as F
from math import isqrt


def rational(value,*,positive=False):
    if type(value) not in (int,F):raise ValueError('exact rational required')
    value=F(value)
    if value<0 or (positive and value==0):raise ValueError('nonnegative rational required')
    return value


def regular_log_budget(delta,terms,residual,probe_norm=1,map_norm=3):
    delta=rational(delta,positive=True);residual=rational(residual)
    probe_norm=rational(probe_norm);map_norm=rational(map_norm,positive=True)
    if type(terms) is not int or terms<1:raise ValueError('positive integer term count required')
    r=1/(1+2*delta);power=2*terms+1
    tail=2*r**power/(power*(1-r*r))
    return 2*map_norm*probe_norm*tail+4*r/(1-r*r)*residual


def core_flow_budget(degree,time_bound,residual,probe_norm=1,map_norm=3):
    if type(degree) is not int or degree<1:raise ValueError('positive integer degree required')
    time_bound=rational(time_bound);residual=rational(residual)
    probe_norm=rational(probe_norm);map_norm=rational(map_norm,positive=True)
    # 1/floor(sqrt N) is a rational upper bound for 1/sqrt N.
    return map_norm*(1+time_bound)*probe_norm/isqrt(degree)+F(degree*3**degree,4)*residual


def core_regularization_budget(delta,time_bound,probe_norm=1,map_norm=3):
    return rational(delta)*rational(time_bound)*rational(probe_norm)*rational(map_norm,positive=True)


def certificate():
    return {'uniform_eight_valley_map_norm_upper':'3',
            'regular_log_series':'2 sum_(j=0)^(N-1) r^(2j+1)*(1-2C)^(2j+1)/(2j+1), r=1/(1+2delta)',
            'regular_log_tail':'2*r^(2N+1)/((2N+1)*(1-r^2))',
            'regular_log_error':'2*M*||f||*tail+4*r/(1-r^2)*max_(j<2N-1)||R_a((1-2C)^j f)||',
            'core_flow_error':'M*(1+T)*||g||/sqrt(N)+(N*3^N/4)*max_(j<N)||R_a(C^j g)||',
            'raw_core_input_extra_error':'||R_a(g)||+||R_a(Cg)||',
            'core_joint_regularization_extra':'M*T*delta*||g||',
            'core_generator_regularization_extra':'M*delta*||g||',
            'regular_log_example_upper':str(regular_log_budget(F(1,4),8,F(1,10**8))),
            'core_flow_example_upper':str(core_flow_budget(144,1,F(1,4*144*3**144))),
            'core_flow_example_inputs':'N=144,T=1,||g||=1,M=3,residual=1/(4*N*3^N)',
            'endpoint_mass_core_bound':'||1_[0,eps]union[1-eps,1](C)*C*(1-C)*g|| <= eps*||g||',
            'fixed_region_required':True,
            'quantitative_spacing_rate_for_all_iterated_sharp_region_probes_claimed':False}
