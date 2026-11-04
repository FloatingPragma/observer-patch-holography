"""Exact full-carrier pole-residue and complete-archive bounds."""
from fractions import Fraction as F

from ..native_mass_channels.algebra import rational


def residue_bounds(*, sites, coupling, gap_lower, inverse_participation,
                   detector_probability):
    """Consequence of the source hypotheses proved in PROOF.md.

    These numbers alone do not establish a ground or its symmetries.
    The certificate binds separately proved sources.
    """
    if type(sites) is not int or sites < 1 or sites % 2 != 1:
        raise ValueError('a positive odd site count is required')
    u, delta, ipr, density = map(rational, (
        coupling, gap_lower, inverse_participation, detector_probability))
    if not (u > 0 and delta > 0 and 0 < ipr <= 1 and 0 <= density <= 1):
        raise ValueError('positive exact coupling/gap and normalized weights required')
    norm = sites*u/4
    if 4*norm >= delta:
        raise ValueError('the isolated four-level hypothesis is not satisfied')
    remainder = norm**2/(delta-2*norm)
    charged = u*ipr/2
    if charged <= 2*remainder:
        raise ValueError('the singly occupied doublet is not certified as ground')
    angle = norm/(delta-2*norm)
    error = 3*angle
    if error >= F(1, 2):
        raise ValueError('the charged pole rank is not resolved by this estimate')
    return {
        'sites': sites, 'full_Fock_dimension': 2**(2*sites),
        'coupling': u, 'full_contact_norm': norm,
        'one_body_and_Fock_gap_lower': delta,
        'zero_orbital_inverse_participation': ipr,
        'projection_angle_bound': angle,
        'normalized_ground_trace_norm_error': 2*angle,
        'pole_matrix_operator_norm_error': error,
        'charged_gap_lower': charged-2*remainder,
        'charged_gap_upper': charged+2*remainder,
        'other_charged_excitations_threshold_lower': delta-2*norm,
        'each_addition_removal_pole_matrix_rank': 2,
        'each_zero_mode_low_pole_weight_lower': F(1, 2)-error,
        'each_zero_mode_low_pole_weight_upper': F(1, 2),
        'each_zero_mode_high_energy_weight_upper': error,
        'detector_zero_orbital_probability': density,
        'each_detector_low_pole_weight_lower': density/2-error,
        'each_detector_low_pole_weight_upper': density/2+error,
        'detector_pole_certified_visible': density/2 > error,
        'zero_mode_anticommutator_nonpole_remainder_bound': 2*error,
    }


def poisson_archive(*, edges, rate, time, ground_angle=F(0)):
    """Exact log probability; tiny probabilities are never rounded to zero."""
    if type(edges) is not int or edges < 1:
        raise ValueError('a positive exact seam count is required')
    rate, time, angle = map(rational, (rate, time, ground_angle))
    if rate < 0 or time < 0 or not 0 <= angle <= 1:
        raise ValueError('nonnegative exact rate/time and a normalized angle required')
    exponent = edges*rate*time
    return {
        'native_seams': edges, 'rate_per_seam': rate, 'duration': time,
        'trace_distance_convention': 'one_half_trace_norm',
        'negative_log_no_event_probability': exponent,
        'no_event_probability': 'exp(-'+str(exponent)+')',
        'exact_complete_archive_distance_for_same_input': '1-exp(-'+str(exponent)+')',
        'trial_ground_input_distance_upper': angle,
        'combined_distance_upper_formula': '1-exp(-x)+exp(-x)*ground_angle',
        'every_event_and_abort_output_retained': True,
        'conditional_no_event_carrier_law_is_exact_Hamiltonian_law': True,
        'autonomous_spectral_protection_inferred': False,
    }


def encode(value):
    if type(value) is F:
        return str(value)
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, list):
        return [encode(item) for item in value]
    return value
