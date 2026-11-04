"""Rational spectral, localization and contact-profile bounds."""
from fractions import Fraction as F

from ..native_mass_channels.algebra import rational


def spectral_window(spacing=1, velocity=1, mass=4, hopping_bound=2):
    a, v, m, k = (rational(x) for x in (spacing, velocity, mass, hopping_bound))
    if min(a, v, m, k) <= 0:
        raise ValueError('positive exact spacing, velocity, mass and hopping bound required')
    speed = v/a
    if k*k < 3*speed*speed:
        raise ValueError('hopping bound must dominate sqrt(3)*v/a')
    if m <= k:
        raise ValueError('mass must strictly exceed the hopping bound')
    gap = m-k
    coupling_squared = 3*speed*speed/2
    origin = 1/(1+coupling_squared/(gap*gap))
    return {'spacing': a, 'velocity': v, 'mass': m, 'hopping_norm_upper': k,
            'one_body_gap_lower': gap, 'kernel_dimension_exact': 2,
            'origin_coupling_norm_squared': coupling_squared,
            'hopping_ratio_upper': k/m,
            'origin_weight_lower': origin, 'origin_weight_upper': F(1),
            'ipr_lower': origin*origin, 'ipr_upper': F(1),
            'origin_particle_number_deviation_contrast_lower': origin*origin,
            'kernel_is_inversion_even': True,
            'local_zero_frame_Gram_is_scalar': True,
            'bounds_uniform_in_infinite_lattice_and_symmetric_finite_cubes': True,
            'fixed_mass_velocity_continuum_refinement_stays_in_domain': False}


def _integer(value, minimum, name):
    if type(value) is not int or value < minimum:
        raise ValueError(f'integer {name} at least {minimum} required')
    return value


def neumann_error_squared(terms, **parameters):
    terms = _integer(terms, 0, 'term count')
    value = spectral_window(**parameters)
    return (value['origin_coupling_norm_squared']/value['one_body_gap_lower']**2
            *value['hopping_ratio_upper']**(2*terms))


def tail_probability_upper(distance, **parameters):
    distance = _integer(distance, 1, 'distance from the origin')
    return min(F(1), neumann_error_squared(distance-1, **parameters))


def finite_volume_projector_error_squared(radius, **parameters):
    radius = _integer(radius, 1, 'cube radius')
    return min(F(1), 4*neumann_error_squared(radius, **parameters))


def encode(value):
    if type(value) is F:
        return str(value)
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(item) for item in value]
    return value


def certificate():
    examples = []
    for mass in (4, 64):
        row = encode(spectral_window(mass=mass))
        row['local_neumann_graph_error_squared'] = [
            {'terms': terms, 'upper': str(neumann_error_squared(terms, mass=mass))}
            for terms in (0, 2, 4, 8)]
        row['finite_volume_projector_error_squared'] = [
            {'radius': radius, 'upper': str(finite_volume_projector_error_squared(radius, mass=mass))}
            for radius in (1, 2, 4, 8)]
        examples.append(row)
    return {'sharp_strong_mass_domain': 'm>sqrt(3)*v/a',
            'rational_admission': 'k^2>=3*(v/a)^2 and 0<k<m; d=m-k, r=k/m',
            'exact_kernel_dimension': 2,
            'nonzero_spectrum_excluded_from': '(-d,d)',
            'origin_coupling_norm_squared': 'b^2=3*v^2/(2*a^2)',
            'origin_weight_lower': '1/(1+b^2/d^2)',
            'IPR_lower': '(1/(1+b^2/d^2))^2',
            'zero_mode_tail_probability': 'sum_{|n|_infinity>=L}|psi(n)|^2 <= min(1,(b^2/d^2)*r^(2*(L-1))), L>=1',
            'finite_volume_projector_bound': '||Z_R-Z|| <= min(1,(2*b/d)*r^R), R>=1',
            'finite_volume_convergence_is_quantitative': True,
            'all_nonzero_bound_levels_outside_central_gap_classified': False,
            'fixed_physical_mass_continuum_limit_proved': False,
            'examples': examples}
