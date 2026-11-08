"""Rational complete-record lifetime bounds on one source clock."""
from fractions import Fraction as F
from model.protected_spin_limit.bounds import history


def _q(value, positive=False):
    if type(value) is not F or value < 0 or positive and not value:
        raise ValueError('nonnegative canonical Fraction required')
    return value


def contrast_lower(initial, time, coupling, degree, free_intervals,
                   tau_over_U, tau3_over_U2, per_protocol_half_error=F()):
    """Zero-field comparison on the same finite input; both outputs retained."""
    for value in (initial, time, tau_over_U, tau3_over_U2, per_protocol_half_error):
        _q(value)
    _q(coupling, True)
    if initial > 1 or type(degree) is not int or degree < 0:
        raise ValueError('probability contrast and nonnegative integer degree required')
    if tau_over_U > F(1, 4):
        raise ValueError('the source comparison requires U >= 4 tau')
    if type(free_intervals) is not int or free_intervals < 1:
        raise ValueError('positive count of free intervals required')
    repair_error = 8*free_intervals*tau_over_U+6*time*tau3_over_U2
    return max(F(), initial-degree*coupling*time-repair_error-2*per_protocol_half_error)


def trajectory_example():
    # 12 cube edges, ten complete repairs, all seams/laws and adaptive policies.
    edges = 12
    lower = contrast_lower(F(67, 100), F(1, 16), F(1), 3, 11,
                           F(1, edges**3), F(1, edges), F(1, 10**8))
    inherited = history(edges**4, edges**8, edges, F(1, 16), 11)
    if inherited != F(8*11, edges**3)+F(6, 16*edges):
        raise ArithmeticError('the retained source history estimate changed')
    if lower <= F(2, 5):
        raise ArithmeticError('the complete source lifetime lost its strict signal')
    return {'edges': edges, 'J': '1', 'c': str(edges**4), 'U': str(edges**8),
            'time_horizon': '1/16', 'maximum_complete_repairs': 10,
            'free_intervals': 11, 'all_seams_and_adaptive_C_E_choices': True,
            'per_protocol_half_trace_error_allowance': '1/100000000',
            'inherited_complete_record_trace_bound': str(inherited),
            'complete_X_pointer_contrast_lower': str(lower),
            'complete_X_pointer_contrast_strictly_greater_than': '2/5',
            'finite_signal_lifetime_is_a_particle_lifetime': False}
