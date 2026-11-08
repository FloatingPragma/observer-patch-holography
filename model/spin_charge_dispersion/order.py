"""Stronger ordered-source bound from the existing exact infrared integral."""
from fractions import Fraction as F

from model.spin_goldstone_sector import infrared
from model.centered_particle_spectrum.certificate import sqrt_bounds, decimal_bounds

CERTIFIED_M0 = F(1, 10)


def biased_margin(integral_upper, bias):
    if type(integral_upper) not in (int, F) or type(bias) not in (int, F) or integral_upper < 0 or bias <= 0:
        raise ValueError('nonnegative rational integral bound and positive rational bias required')
    lower_sqrt, upper_sqrt = sqrt_bounds(F(24))
    if not lower_sqrt*lower_sqrt <= 24 <= upper_sqrt*upper_sqrt:
        raise ArithmeticError('outward sqrt24 enclosure required')
    bragg = F(1, 12)-F(integral_upper)/lower_sqrt
    if bragg <= 0:
        raise ValueError('this infrared bound does not certify positive order')
    magnetization = 2*bias*bragg/(1+bias*bias*bragg)
    return {'sqrt24_outward_interval': [str(lower_sqrt), str(upper_sqrt)],
        'staggered_Bragg_weight_lower': str(bragg), 'fixed_trial_bias': str(bias),
        'selected_magnetization_lower': str(magnetization),
        'magnetization_decimal_lower': decimal_bounds(magnetization, magnetization, 12)[0],
        'finite_trial_energy_upper_coefficient_in_J_over_N': str(3*bias*bias),
        'projected_excitation_operator_norm_upper': 'J/sqrt(N)'}


def certificate():
    rows = []
    for divisions, bias, admitted in ((64, F(14), F(7, 100)), (256, F(10), CERTIFIED_M0)):
        actual = infrared.certificate(divisions)
        row = biased_margin(F(actual['integral_upper']), bias)
        if F(row['selected_magnetization_lower']) <= admitted:
            raise ArithmeticError('actual infrared integral and biased ground do not prove the claimed order')
        rows.append({'existing_exact_infrared_producer': actual, 'biased_ground_bound': row,
            'simple_strict_magnetization_lower': str(admitted)})
    return {'controls': rows, 'selected_magnetization_strict_lower': str(CERTIFIED_M0),
        'source_twist_kappa_strict_lower': str(CERTIFIED_M0**2/2),
        'bias_changes_source_Hamiltonian': False,
        'trial_vector_claimed_exact_selected_ground': False,
        'same_maximal_field_selected_order_theorem_consumed': True,
        'empirical_magnetization_or_particle_data_used': False}
