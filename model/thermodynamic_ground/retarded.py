"""Upper-half-plane local response along the proved ground-state subsequence.

No real-axis continuation, isolated pole or computable state convergence is
inferred. Bounds concern bounded local operators in the original dynamics.
"""
from fractions import Fraction as F
from math import factorial

from .volume import averaged_dynamics_error


def _positive(value,name,zero=False):
    if type(value) not in (int,F) or value<0 or (not zero and value==0):
        raise ValueError('exact '+('nonnegative ' if zero else 'positive ')+name+' required')
    return F(value)


def exp_negative_upper(value,order=32):
    """1/positive Taylor partial sum >=exp(-value), with no asymptotic premise."""
    x=_positive(value,'exponent',True)
    if type(order) is not int or not 1<=order<=64:
        raise ValueError('integer Taylor order1..64 required')
    term=total=F(1)
    for k in range(1,order+1):
        term*=x/k;total+=term
    return 1/total


def tail_bound(eta,time,first_norm=F(1),second_norm=F(1)):
    eta=_positive(eta,'imaginary energy');time=_positive(time,'cutoff time',True)
    a=_positive(first_norm,'first norm',True);b=_positive(second_norm,'second norm',True)
    return 2*a*b*exp_negative_upper(eta*time)/eta


def derivative_bound(eta,order,first_norm=F(1),second_norm=F(1)):
    eta=_positive(eta,'imaginary energy')
    a=_positive(first_norm,'first norm',True);b=_positive(second_norm,'second norm',True)
    if type(order) is not int or not 0<=order<=64:raise ValueError('derivative order0..64 required')
    return 2*a*b*factorial(order)/eta**(order+1)


def comparison_bound(half_side,buffer_radius,support_radius,support_sites,time,eta,
                     *,state_discrepancy,speed=F(1,8),first_norm=F(1),second_norm=F(1)):
    """Conditional finite response error, retaining the unknown state term.

    state_discrepancy bounds sup_[0,time]|nu_L({tau_t(A),B*})-
    omega({tau_t(A),B*})| and must be separately supplied. Compactness proves
    it vanishes along a subsequence but supplies no finite numerical bound.
    """
    time=_positive(time,'cutoff time',True);eta=_positive(eta,'imaginary energy')
    a=_positive(first_norm,'first norm',True);b=_positive(second_norm,'second norm',True)
    state=_positive(state_discrepancy,'state discrepancy',True)
    dynamics=averaged_dynamics_error(half_side,buffer_radius,support_radius,support_sites,time,speed)
    return 2*a*b*time*dynamics+time*state+2*tail_bound(eta,time,a,b)


def certificate():
    return {'definition':'G_AB(z)=-i integral_0^infinity exp(i z t) omega({tau_t(A),B^*})dt,Imz>0',
            'finite_comparator':'spatial average of the same finite full-ground anticommutator correlations using finite-box dynamics',
            'same_ground_subsequence_as_existence_theorem':True,
            'convergence':'uniform on every compact subset of Imz>0',
            'individual_tail':'2||A||||B|| exp(-eta T)/eta forImz>=eta>0',
            'derivative_bound':'|d_z^n G_AB(z)|<=2||A||||B|| n!/eta^(n+1)',
            'comparison_bound':'2||A||||B||T[2B_R(T)+2bad_fraction]+T epsilon_L(T)+4||A||||B||exp(-eta T)/eta, capped dynamics as declared',
            'epsilon_L':'sup_0<=t<=T |nu_L({tau_t(A),B^*})-omega({tau_t(A),B^*})|; tends to0 along the common subsequence but no rate is computed',
            'rational_tail_examples':[{'eta':'1','time':str(t),'upper':str(tail_bound(1,t))} for t in (8,16,32)],
            'analytic_interacting_response_exists':True,
            'infinite_response_values_computed':False,
            'real_axis_limit_or_poles_or_gap_proved':False,
            'spatial_local_inverse_or_metric_proved':False}
