"""Rational localization constants for the exact bounded-source zero modes."""
from fractions import Fraction as F

from ..native_mass_channels.algebra import rational


def parameters(spacing=1, velocity=1, mass=1):
    values = tuple(rational(x) for x in (spacing, velocity, mass))
    if any(x <= 0 for x in values):
        raise ValueError('positive exact spacing, velocity and mass required')
    return values


def localization(spacing=1, velocity=1, mass=1, core_support_radius=0):
    a, v, m = parameters(spacing, velocity, mass)
    if type(core_support_radius) is not int or core_support_radius < 0:
        raise ValueError('nonnegative integer perturbation support radius required')
    scale = 12*v/(a*m)
    radius = max(1, -(-scale.numerator//scale.denominator), core_support_radius)
    q = 1+a*m/(12*v)
    anticommutator = 6*m*v/(a*(radius+1))
    squared_gap = m*m-anticommutator
    weight_error = 3*v*(q-1)/a
    cutoff_error = 3*v/(a*radius)
    weighted_inverse = 4/m
    if not (squared_gap >= m*m/2 and weight_error == m/4
            and cutoff_error <= m/4 and weighted_inverse*cutoff_error <= 1):
        raise ArithmeticError('the exterior gap or weighted inverse estimate failed')
    return {'spacing': a, 'velocity': v, 'mass': m,
            'core_support_radius': core_support_radius, 'reference_radius_sites': radius,
            'physical_reference_radius': a*radius, 'weight_per_site': q,
            'exterior_anticommutator_norm_upper': anticommutator,
            'exterior_squared_gap_lower': squared_gap,
            'exterior_gap_conservative_lower': m/2,
            'weight_conjugation_error_upper': weight_error,
            'weighted_exterior_inverse_norm_upper': weighted_inverse,
            'cutoff_commutator_norm_upper': cutoff_error,
            'weighted_tail_prefactor_after_radius_2R': weighted_inverse*cutoff_error,
            'small_spacing_physical_decay_domain': a*m <= 12*v,
            'physical_probability_decay_length_upper_in_small_spacing_domain': 12*v/m,
            'unperturbed_radius_2aR_upper': 24*v/m+2*a if core_support_radius == 0 else None}


def tail_bound(steps, spacing=1, velocity=1, mass=1):
    if type(steps) is not int or steps < 0:
        raise ValueError('nonnegative integer tail distance required')
    value = localization(spacing, velocity, mass)
    return value['weight_per_site']**(-2*steps)


def encode(value):
    if type(value) is F:
        return str(value)
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, list):
        return [encode(item) for item in value]
    return value


def certificate():
    return {'localization_statement': 'For every normalized exact zero mode, sum_{|n|_infinity>=2R+L}|psi(n)|^2 <= q^(-2L), integer L>=0.',
            'reference_radius': 'R=max(1,ceil(12*v/(a*m)),support_radius)',
            'weight': 'q=1+a*m/(12*v)',
            'localization_applies_to_original_nearest_neighbor_profile': True,
            'finite_core_perturbations_require_support_inside_reference_cube': True,
            'localization_extended_to_arbitrary_nonlocal_compact_perturbations': False,
            'small_spacing_domain': 'a*m<=12*v',
            'small_spacing_probability_bound': 'q^(-2L)<=exp(-m*a*L/(12*v))',
            'small_spacing_unperturbed_physical_radius': '2*a*R<=24*v/m+2*a',
            'examples': [encode(localization(a)) for a in (F(1), F(1, 2), F(1, 8))],
            'tail_probability_examples_at_unit_parameters': [
                {'additional_lattice_steps': steps, 'upper': str(tail_bound(steps))}
                for steps in (0, 12, 60, 120)]}
