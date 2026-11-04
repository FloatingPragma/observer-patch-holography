"""Exact rational budgets for the interior coarea proof and common clock."""
from fractions import Fraction as F

from ..contact_continuum.continuum import rational
from ..moving_bound_pair.limit import moving_pole_certificate


def coarea_certificate():
    k2, h = F(1, 16), F(1, 128)
    f_min = 1-k2/2
    derivative_min = f_min-k2/4
    derivative_slope = F(5, 2)*k2
    t_lower = F(5, 3)*(1-k2/4-h/6)
    assert t_lower > F(3, 2)
    assert f_min == F(31, 32) and derivative_min == F(61, 64)
    assert derivative_min > F(7, 8) and derivative_slope < F(1, 4)
    assert F(9, 4)*(F(32, 31)**5) < 4
    threshold_error = (F(5, 2)+F(1, 9))/f_min
    assert threshold_error < 3
    angular_derivative = 2*(2*1*2*2+2*1*2*2+2*1*1*2+2*1*2*1)
    assert angular_derivative == 48
    interior_real = F(4, 9)*(F(64, 7)+F(96, 7))
    inner_and_tail = 512+interior_real+4
    assert inner_and_tail < 528
    momentum_error, energy_error = F(21, 8)+F(1, 48), F(3173, 12)+F(1, 96)
    assert momentum_error < 3 and energy_error < 265
    return {'maximum_lattice_momentum': F(1, 4),
            'maximum_scaled_interior_energy': h, 'transformed_ball_radius': F(1, 4),
            'radial_shape_lower': f_min, 'radial_squared_map_derivative_lower': derivative_min,
            'radial_squared_map_second_derivative_upper': derivative_slope,
            'angular_amplitude_derivative_upper': angular_derivative,
            'angular_threshold_deviation_coefficient': 2,
            'threshold_W_momentum_error_coefficient': threshold_error,
            'interior_W_energy_remainder_coefficient': inner_and_tail,
            'rounded_W_energy_remainder_coefficient': 528,
            'rounded_W_momentum_remainder_coefficient': 1,
            'scaled_green_momentum_error_coefficient': momentum_error,
            'scaled_green_energy_error_coefficient': energy_error,
            'rounded_scaled_green_momentum_coefficient': 3,
            'rounded_scaled_green_energy_coefficient': 265,
            'upper_boundary_value_sign': 'positive E regulator gives positive denominator sign',
            'all_edge_valleys_retained': 8}


def green_certificate(spacing, momentum_cap, energy_cap):
    a = rational(spacing, 'spacing', positive=True)
    p = rational(momentum_cap, 'componentwise model momentum cap')
    cap = rational(energy_cap, 'relative energy cap', positive=True)
    if a*p > F(1, 4) or a*a*cap > F(1, 128):
        raise ValueError('requires a*Pmax<=1/4 and a^2*S<=1/128')
    return {'spacing': a, 'momentum_cap': p, 'relative_energy_cap': cap,
            'scaled_green_error_upper': a*(3*p*p+265*cap),
            'comparison': '(g_K(B_K-s+i0)-a^2 W0/2)/a^3 -> -i sqrt(s)/pi',
            'uniform_domain': '|P_j|<=Pmax, K=aP, 0<=s<=S',
            'PV_bound_proved_separately_from_exterior_bound': True}


def scattering_certificate(spacing, gamma, momentum_cap, energy_cap):
    result = green_certificate(spacing, momentum_cap, energy_cap)
    a, g = result['spacing'], rational(gamma, 'inverse length', positive=True)
    error = result['scaled_green_error_upper']
    if a*g > 1 or 44*error > 7*g:
        raise ValueError('requires a*gamma<=1 and 44*GreenError<=7*gamma')
    return {**result, 'gamma': g, 'coupling_positive': True,
            'denominator_modulus_lower_after_a3': 7*g/44,
            'uniform_scattering_operator_error_upper': min(F(2), F(88, 7)*error/g),
            'source_scattering_limit': '(gamma+i sqrt(s))/(gamma-i sqrt(s))',
            'reflected_scattering_limit': '(gamma-i sqrt(s))/(gamma+i sqrt(s))',
            'comparison_space': 'direct-integral contact cyclic spectral channel L2(dP ds)',
            'orthogonal_noncontact_channels': 'identity, retained in the source carrier',
            'bound_poles_discarded': False}


def joined_certificate(spacing, gamma, momentum_cap, energy_cap, duration):
    result = scattering_certificate(spacing, gamma, momentum_cap, energy_cap)
    pole = moving_pole_certificate(spacing, gamma, momentum_cap)
    time = rational(duration, 'absolute model time')
    a, p = result['spacing'], result['momentum_cap']
    return {'scattering': result, 'moving_bound_band': pole,
            'same_coupling': '1/U_a=a^2 W0/2-a^3 gamma/pi',
            'same_global_generator': 'J_a=6/a^2-H_a',
            'relative_scattering_energy': 'J_a free energy = 6/a^2-B_a(aP)+s',
            'limiting_scattering_energy': '|P|^2/4+s',
            'limiting_bound_energy': '|P|^2/4-gamma^2',
            'scattering_free_clock_operator_error_upper': min(F(2), time*a*a*p**4/4),
            'bound_clock_operator_error_upper': min(F(2), time*pole['reflected_band_energy_error_upper']),
            'relativistic_or_many_body_QFT_constructed': False,
            'source_selected_gamma': False}


def spacing_certificate(gamma, momentum_cap, energy_cap, tolerance, duration=1):
    g = rational(gamma, 'inverse length', positive=True)
    p = rational(momentum_cap, 'momentum cap')
    cap = rational(energy_cap, 'relative energy cap', positive=True)
    eps = rational(tolerance, 'operator tolerance', positive=True)
    if eps > 2:
        raise ValueError('operator tolerance is at most two')
    a = min(1/g, 1/(4*(1+p)), 1/(8*(1+cap)),
            g/(256*(g*g+p*p)), 7*eps*g/(88*(3*p*p+265*cap)))
    return {'spacing_upper': a, 'scattering_tolerance': eps,
            **joined_certificate(a, g, p, cap, duration)}
