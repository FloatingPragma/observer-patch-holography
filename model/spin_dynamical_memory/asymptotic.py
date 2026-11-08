"""Uniform actual-source filter bounds and a field-level implication control."""
from fractions import Fraction as F
from model.protected_defect_catalog.field import preparation_bounds


def source_filter_bounds(order_lower=F(1, 1000), regularizer=F(1), resolution=F(1)):
    m0, lam, z = order_lower, regularizer, resolution
    if any(type(x) not in (int, F) for x in (m0, lam, z)):
        raise TypeError('exact rational source bounds required')
    m0, lam, z = map(F, (m0, lam, z))
    if not 0 < m0 <= F(1, 2) or min(lam, z) <= 0:
        raise ValueError('positive filter scales and admissible order lower bound required')
    inherited = preparation_bounds(m0, lam, z)
    n0 = F(inherited['liminf_original_odd_preparation_probability'])
    h0 = F(inherited['limsup_unnormalized_scaled_energy'])
    asymptotic_floor = F(inherited['liminf_odd_filtered_original_probability'])/lam**2
    eventual_filter = F(inherited['sufficiently_small_scale_filtered_probability_lower'])/lam**2
    return {'order_lower': m0, 'Q_Weyl_amplitude': 4/m0,
        'comparison_P_Weyl_amplitude': m0/4,
        'odd_sine_amplitude': 2/m0,
        'bounded_probe': 'sin((2/m0)Q)',
        'source_parity': 'odd',
        'chosen_reference_component_removed': True,
        'all_zero_energy_vectors_removed': False,
        'unfiltered_norm_liminf_lower': n0,
        'unfiltered_scaled_energy_limsup_upper': h0,
        'regularizer': lam, 'memory_resolution': z,
        'filtered_norm_liminf_lower': asymptotic_floor,
        'eventual_filtered_norm_lower': eventual_filter,
        'eventual_original_filter_success_lower': lam**2*eventual_filter,
        'eventual_filtered_memory_upper': 1/(z*eventual_filter),
        'preparation_and_parity_producer': 'model.protected_defect_catalog.field.preparation_bounds',
        'effective_small_p_threshold_supplied': False,
        'vanishing_memory_inferred': False}


def two_speed_control():
    """Exact Fock-field control of theorem consequences, not an OPH source."""
    import sympy as s
    m, kappa = s.Rational(1, 4), s.Rational(1, 4)
    a = 12*kappa/m
    b1, b2 = kappa/(2*m), 3*kappa/(2*m)
    c1, c2 = s.sqrt(a*b1), s.sqrt(a*b2)
    q_norm = s.simplify(m*a/4*(1/c1+1/c2))
    p_norm = s.simplify(m/(4*a)*(c1+c2))
    q_first = m*a/2
    p_first = m*(b1+b2)/4
    q_inverse = s.simplify(m*a/4*(1/c1**2+1/c2**2))
    p_inverse = m/(2*a)
    ratio_defect = s.simplify(2*(c1*c1+c2*c2)/(c1+c2)**2-1)
    assert q_first == 6*kappa and p_first == kappa/2
    assert q_inverse == s.Rational(1, 6) and p_inverse == s.Rational(1, 96)
    assert ratio_defect == 7-4*s.sqrt(3)
    assert m*m/2 <= kappa <= s.Rational(1, 2)
    assert q_norm < 3 and p_norm < s.Rational(1, 4)
    assert q_inverse <= 4 and p_inverse <= s.Rational(1, 12)
    assert c2*s.sqrt(3) < 72
    return {'order': str(m), 'bond_coefficient': str(kappa),
        'kinetic_coefficient': str(a), 'gradient_coefficients': [str(b1), str(b2)],
        'two_speeds_squared': [str(c1*c1), str(c2*c2)],
        'Q_norm_coefficient_over_inverse_momentum': str(q_norm),
        'P_norm_coefficient_over_momentum': str(p_norm),
        'Q_first_energy_moment': str(q_first),
        'P_first_energy_coefficient_over_momentum_squared': str(p_first),
        'Q_inverse_energy_coefficient_over_inverse_momentum_squared': str(q_inverse),
        'P_inverse_energy_moment': str(p_inverse),
        'P_energy_over_momentum_ratio_defect': str(ratio_defect),
        'canonical_commutator_coefficient': str(m),
        'positive_energy_spatial_rotations_and_finite_cone': True,
        'source_leading_first_moments_and_stated_bounds_satisfied': True,
        'one_common_Lorentz_speed': False,
        'is_a_state_or_scaling_limit_of_the_OPH_spin_source': False,
        'nonzero_scattering': False}
