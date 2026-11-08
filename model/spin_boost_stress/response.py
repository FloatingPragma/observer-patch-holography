"""Exact finite boost-response estimates from the actual source defect."""
from fractions import Fraction as F
from math import factorial


def hyperbolic_bounds(rapidity, order=48):
    if type(rapidity) is not F or abs(rapidity) > 1:
        raise ValueError('exact rapidity in [-1,1] required')
    if type(order) is not int or order < 4:
        raise ValueError('at least four series terms required')
    x = abs(rapidity)
    cosh = sum((x**(2*j)/factorial(2*j) for j in range(order)), F(0))
    sinh = sum((x**(2*j+1)/factorial(2*j+1) for j in range(order)), F(0))
    ct = x**(2*order)/factorial(2*order)/(1-F(1, (2*order+1)*(2*order+2)))
    st = x**(2*order+1)/factorial(2*order+1)/(1-F(1, (2*order+2)*(2*order+3)))
    return {'cosh_minus_one_lower': cosh-1, 'cosh_minus_one_upper': cosh-1+ct,
            'absolute_sinh_lower': sinh, 'absolute_sinh_upper': sinh+st}


def history_bound(rapidity, defect_norm, clock_durations, candidate_flux_durations):
    if type(defect_norm) is not F or defect_norm < 0:
        raise ValueError('exact nonnegative source defect bound required')
    clock_durations, candidate_flux_durations = tuple(clock_durations), tuple(candidate_flux_durations)
    if any(type(t) is not F for t in clock_durations+candidate_flux_durations):
        raise ValueError('exact declared operator-clock durations required')
    functions = hyperbolic_bounds(rapidity)
    clock = defect_norm*functions['cosh_minus_one_upper']
    flux = defect_norm*functions['absolute_sinh_upper']
    total = sum(abs(t) for t in clock_durations)*clock+sum(abs(t) for t in candidate_flux_durations)*flux
    return {'rapidity': rapidity, 'clock_generator_error_upper': clock,
            'candidate_flux_generator_error_upper': flux,
            'complete_history_halftrace_error_upper': min(F(1), total),
            'all_intervening_complete_channels_and_references_allowed': True,
            'candidate_flux_identified_as_physical_translation': False,
            'finite_defect_removed_by_small_rapidity_steps': False}
