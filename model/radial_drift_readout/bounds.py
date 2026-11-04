"""Exact local drift and classifier bounds, independent of source volume."""
from fractions import Fraction as F
from math import isqrt

from model.radial_geometry_join.bounds import joined_bounds


def rational(value, *, positive=False):
    if type(value) not in (int, F):
        raise ValueError('an exact rational input is required')
    value = F(value)
    if value < 0 or positive and value == 0:
        raise ValueError('a positive input is required' if positive else 'a nonnegative input is required')
    return value


def sqrt_interval(square, bits=80):
    square = rational(square)
    if type(bits) is not int or not 1 <= bits <= 4096:
        raise ValueError('integer precision between 1 and 4096 required')
    scale = 1 << bits
    lower = F(isqrt(square.numerator*scale*scale//square.denominator), scale)
    upper = lower if lower*lower == square else lower+F(1, scale)
    if not lower*lower <= square <= upper*upper:
        raise ArithmeticError('integer square-root enclosure failed')
    return lower, upper


def local_drift(total_time, *, velocity=F(1), spacing=F(1), coupling=F(1, 65536), dimensions=3):
    total_time = rational(total_time)
    velocity, spacing = rational(velocity, positive=True), rational(spacing, positive=True)
    coupling = rational(coupling)
    if type(dimensions) is not int or dimensions not in (1, 3):
        raise ValueError('one or three spatial directions required')
    star_squared = 2*dimensions*(velocity/spacing)**2
    local_squared = star_squared+coupling*coupling/16
    star, local = sqrt_interval(star_squared), sqrt_interval(local_squared)
    error = min(F(1), total_time*local[1])
    return {
        'total_elapsed_Hamiltonian_time': total_time,
        'dimensions': dimensions,
        'velocity': velocity, 'spacing': spacing, 'centered_contact_coupling': coupling,
        'star_Fock_norm_squared_exact': star_squared,
        'origin_parity_projector_commutator_norm_squared_exact': star_squared,
        'origin_sign_parity_commutator_norm_squared_exact': 4*star_squared,
        'star_plus_center_contact_norm_squared_exact': local_squared,
        'star_Fock_norm_interval': star,
        'star_plus_center_contact_norm_interval': local,
        'interround_only_complete_half_trace_error_upper': min(F(1), total_time*star[1]),
        'timed_control_complete_half_trace_error_upper': error,
        'timed_control_complete_diamond_error_upper': 2*error,
        'source_volume_factor': 1,
        'mass_texture_factor': 0,
        'field_reference_and_full_archive_included': True,
    }


def classifier_bound(total_time=F(1, 10**6), *, rounds=32):
    if type(rounds) is not int or rounds != 32:
        raise ValueError('the retained classifier has exactly 32 scored rounds')
    source = joined_bounds()
    drift = local_drift(total_time)
    spectral_error = source['ground_D0_expectation_interval'][1]
    tail, controls = F(1, 10**6), rounds*F(1, 10**9)
    error = spectral_error+tail+controls+drift['timed_control_complete_half_trace_error_upper']
    return {
        'promised_source': 'actual mass64, U=1/65536, radius1 full-Fock ground mixture or either charged singlet',
        'rounds': rounds, 'tie_rule': 'charged at count >=16; every failure and halted slot scored zero',
        'actual_spectral_wrong_parity_upper': spectral_error,
        'conditional_binomial_tail_strict_upper': tail,
        'complete_native_control_error_strict_upper': controls,
        'drift': drift,
        'misclassification_strict_upper': error,
        'below_one_in_2500': error < F(1, 2500),
        'record_law_is_exact_frozen_parity_mixture': False,
        'complete_native_load_mean_coordinate_error_upper': 2*min(F(1), controls+drift['timed_control_complete_half_trace_error_upper']),
        'input_spectral_preparation_is_supplied': True,
        'postmeasurement_ground_or_charged_sector_preserved': False,
    }


def certificate():
    example = classifier_bound()
    if not example['below_one_in_2500']:
        raise ArithmeticError('positive-time complete source classifier is not certified')
    square = example['drift']['star_plus_center_contact_norm_squared_exact']
    if not square < F(25, 4):
        raise ArithmeticError('simple local norm bound 5/2 failed')
    return {
        'positive_time_witness': example,
        'simple_norm_strict_upper': F(5, 2),
        'simple_drift_strict_upper_at_witness': F(1, 400000),
        'generic_norm_has_no_volume_or_mass_factor': True,
        'interacting_spectral_preparation_bound_is_volume_uniform': False,
        'elapsed_time_must_include_every_control_transfer_and_record_delay': True,
        'clock_rate_or_physical_control_bandwidth_derived': False,
    }
