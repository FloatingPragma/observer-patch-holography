"""Exact original-probability and common-clock bounds from source theorems."""
from fractions import Fraction as F
from .order import CERTIFIED_M0


def band(m0=CERTIFIED_M0, radius=F(1)):
    if any(type(x) not in (int, F) for x in (m0, radius)) or not 0 < m0 <= F(1, 2) or radius <= 0:
        raise ValueError('positive rational order lower bound and momentum radius required')
    m0, radius = F(m0), F(radius)
    squared_commutator = m0*m0/(4*13**3)
    lower = 3*squared_commutator/2
    upper = 2/squared_commutator
    tail = lower/12+1/(4*upper)
    if tail != squared_commutator/4:
        raise ArithmeticError('source conjugate-field tail budget differs')
    return {'order_lower_bound': str(m0), 'momentum_annulus': [str(radius/2), str(radius)],
        'squared_commutator_expectation_lower': str(squared_commutator),
        'energy_window_in_same_scaled_clock': [str(radius*lower), str(radius*upper)],
        'dimensionless_energy_over_radius_window': [str(lower), str(upper)],
        'conjugate_field_tail_norm_squared_upper': str(tail),
        'each_charge_original_joint_band_probability_lower': str(squared_commutator/2),
        'each_charge_total_original_probability_upper': '1/4',
        'each_creator_unnormalized_first_energy_moment_upper': str(3*radius),
        'three_preparation_records': ['charge_plus', 'charge_minus', 'complement'],
        'nonlinear_preparation_claimed_entirely_momentum_annular': False,
        'positive_band_weight_uniform_as_radius_tends_to_zero': True,
        'sharp_speed_or_mass_shell_derived': False}


def event_budget(m0=CERTIFIED_M0, complete_trace_error=F(1, 10**12)):
    row = band(m0)
    if type(complete_trace_error) not in (int, F) or not 0 <= complete_trace_error <= 2:
        raise ValueError('complete unhalved trace error must lie in [0,2]')
    original = F(row['each_charge_original_joint_band_probability_lower'])
    return {'ideal_original_probability_lower': str(original),
        'complete_trace_error': str(complete_trace_error),
        'observed_original_probability_lower': str(max(F(), original-F(complete_trace_error)/2)),
        'success_probability_renormalized': False}
