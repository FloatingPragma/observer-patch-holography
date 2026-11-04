"""Rational bounds derived from the complete compact-source gap."""
from fractions import Fraction as F

from model.compact_defect_pair.topology import spectral_bounds
from .source import parameters


def certificate(mass=64, half_separation=1, radius=3, eta=F(1, 512), rounds=1):
    base = spectral_bounds(mass, half_separation, radius)
    eta = parameters(eta)
    if type(rounds) is not int or rounds < 0:
        raise ValueError('nonnegative integer round count required')
    r = F(base['band_abs_energy_upper'])
    g = F(base['complement_abs_energy_lower'])
    contour = (g+r)/2
    distance = (g-r)/2
    epsilon = 4*F(mass)*eta
    beta = contour*epsilon/distance**2
    if beta >= 1:
        raise ValueError('relative budget does not certify unique polar transport')
    return {'mass': str(mass), 'half_separation': half_separation,
            'relative_budget': str(eta), 'fixed_core_zero_coordinates': 4,
            'isolated_band_dimension': 4,
            'band_abs_energy_upper': str(r), 'complement_abs_energy_lower': str(g),
            'core_probability_lower': base['core_probability_lower'],
            'any_two_admitted_source_difference_norm_upper': str(epsilon),
            'resolvent_contour_radius': str(contour),
            'resolvent_distance_to_each_spectrum_lower': str(distance),
            'spectral_projection_difference_norm_upper': str(beta),
            'single_quench_leakage_probability_upper': str(beta**2),
            'uncontrolled_rounds': rounds,
            'uncontrolled_history_leakage_probability_upper': str(min(F(1), (rounds*beta)**2)),
            'controlled_record_endpoint_band_leakage_probability': '0',
            'controlled_one_body_band_work_abs_upper': str(2*r),
            'bounds_uniform_in_exterior_volume': True,
            'local_topology_guards_needed_for_band_count': False}
