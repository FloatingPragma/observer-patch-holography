"""A certified complete-clock event on the unchanged selected eight-site source."""
from fractions import Fraction as F
from pathlib import Path
import numpy as np

from model.native_ordered_ground import certificate, source
from model.native_joint_spin_spectroscopy import algebra, spectral
from model.native_spin_observer.precision import pi_interval
from model.centered_particle_spectrum.certificate import decimal_bounds
from .bounds import bin_ranges, tube_bounds

ROOT = Path(__file__).resolve().parents[2]


def energy_chart(energy_upper, tau):
    if type(energy_upper) is not F or type(tau) is not F or energy_upper < 0 or tau <= 0:
        raise ValueError('exact nonnegative full excitation bound and positive clock step required')
    pi_lo, pi_hi = pi_interval(110)
    if tau*energy_upper >= 2*pi_lo:
        raise ValueError('the full source excitation spectrum aliases across the clock period')
    return {'full_excitation_turn_upper': tau*energy_upper/(2*pi_lo),
            'clock_period_lower': 2*pi_lo/tau,
            'clock_period_upper': 2*pi_hi/tau}


def kernel_interval(energy, bin_index, size, tau):
    """Exact Fejer interval using the quotient away from its removable zeros."""
    if (type(energy) is not F or type(tau) is not F or tau <= 0
            or type(size) is not int or size < 2 or type(bin_index) is not int
            or not 0 <= bin_index < size):
        raise ValueError('exact energy, clock and recorded bin required')
    a, b = spectral.cosine_interval(tau*energy, F(-2*bin_index, size))
    denominator = (1-b, 1-a)
    if denominator[0] <= F(1, 10**12):
        return spectral.fejer_interval(energy, bin_index, size, tau)
    # M times the bin phase is an integer multiple of 2pi.
    a, b = spectral.cosine_interval(size*tau*energy)
    numerator = (max(F(), 1-b), max(F(), 1-a))
    lo = max(F(), numerator[0]/(size**2*denominator[1]))
    hi = min(F(1), numerator[1]/(size**2*denominator[0]))
    # Outward dyadic rounding keeps sums over the full witness practical;
    # arbitrary quotient denominators would multiply across all256 levels.
    scale = 1 << 80
    lower, upper = lo*scale, hi*scale
    return (F(lower.numerator//lower.denominator, scale),
            F(-(-upper.numerator//upper.denominator), scale))


def certificate_report():
    with np.load(ROOT/'model/native_ordered_ground/witness.npz', allow_pickle=False) as archive:
        witness = {name: archive[name] for name in archive.files}
    checked = certificate.verify(witness)
    h, scale = source.matrix()
    group = algebra.audit(h, scale)
    vectors = witness['vectors'].astype(object)
    ground = vectors[:, 0]
    ground_norm = sum(int(x)**2 for x in ground)
    divisor = 1 << certificate.BITS
    energies = [F(int(x), divisor) for x in witness['energies']]
    origin = energies[0]
    e0lo, e0hi = map(F, checked['ground_energy_interval'])
    origin_error = max(origin-e0lo, e0hi-origin)
    radius = F(checked['complete_polar_Weyl_radius'])
    state_error = F(checked['normalized_ground_vector_error'])
    size, tau, character = 128, F(1, 2), 7
    pi_lo, pi_hi = pi_interval(110)
    excitation_upper = energies[-1]+radius-e0lo
    chart = energy_chart(excitation_upper, tau)
    # Finite source-selected test band [4/5,7/5], centered at11/10.
    # Conservative rational turn radii include pi uncertainty explicitly.
    center = (tau*F(11, 10)/(2*pi_hi), tau*F(11, 10)/(2*pi_lo))
    width = tau*F(3, 10)/(2*pi_hi)
    margin = tau*F(1, 5)/(2*pi_hi)
    ranges = bin_ranges(size, center, width, margin)
    probe = ground*np.array([1-2*(mask & 1) for mask in range(256)], dtype=object)
    projected = algebra.character_projection(probe, character)
    coefficients = vectors.T @ projected
    weights = [F(int(x)**2, divisor**2*64*ground_norm) for x in coefficients]
    success = F(sum(int(x)**2 for x in projected), 64*ground_norm)
    success_interval = (max(F(), success-2*state_error), min(F(1), success+2*state_error))
    one_bin_error = spectral.operator_error(size, tau, radius, origin_error,
                                            F(checked['Gram_defect']), state_error)
    events = {}
    for name in ('inner', 'outer'):
        first, last = ranges[name]
        lo = hi = F()
        for j in range(first, last+1):
            for weight, energy in zip(weights, energies, strict=True):
                a, b = kernel_interval(energy-origin, j, size, tau)
                lo += weight*a
                hi += weight*b
        error = max(0, last-first+1)*one_bin_error
        events[name] = (max(F(), lo-error), min(F(1), hi+error))
    bound = tube_bounds(success_interval, events['inner'], events['outer'],
                        ranges['complete_clock_tail_upper'])
    if bound['original_tube_weight_lower'] <= F(1, 5):
        raise ArithmeticError('complete finite tube readout lost its positive source weight')
    # Every true eigenvalue lies within radius of its listed proposal, and
    # the ground origin has its independent certified error. Certify that
    # every level in [4/5,7/5] belongs to the first excitation cluster.
    cluster = []
    for i, energy in enumerate(energies):
        lo, hi = energy-radius-e0hi, energy+radius-e0lo
        phase_lo, phase_hi = tau*lo/(2*pi_hi), tau*hi/(2*pi_lo)
        if hi < F(4, 5) or lo > F(7, 5):
            if not (phase_hi < center[0]-width or phase_lo > center[1]+width):
                raise ArithmeticError('a noncluster level intersects the certified phase tube')
            continue
        if not F(1) < lo <= hi < F(6, 5):
            raise ArithmeticError('the source test band contains an additional spectral cluster')
        if not center[1]-width < phase_lo <= phase_hi < center[0]+width:
            raise ArithmeticError('the first source cluster is outside the inner phase tube')
        cluster.append(i)
    if not cluster:
        raise ArithmeticError('the certified source gap cluster disappeared')
    return {'source': 'selected open eight-site spin cube, J=1, B=1/4',
            'source_group': group, 'character_bits': '111',
            'character_is_translation_momentum': False,
            'complete_clock_size': size, 'tau': str(tau),
            'all_8_times_128_joint_bins_retained': True,
            'all_source_energies_in_one_clock_period': True,
            'full_excitation_turn_upper': str(chart['full_excitation_turn_upper']),
            'source_selected_energy_band': ['4/5', '7/5'],
            'center_turn_interval': list(map(str, center)),
            'turn_width': str(width), 'turn_margin': str(margin),
            'inner_bins_inclusive': list(ranges['inner']),
            'outer_bins_inclusive': list(ranges['outer']),
            'clock_tail_upper': str(ranges['complete_clock_tail_upper']),
            'success_probability_interval': decimal_bounds(*success_interval, 11),
            'inner_event_probability_interval': decimal_bounds(*events['inner'], 11),
            'outer_event_probability_interval': decimal_bounds(*events['outer'], 11),
            'original_band_weight_interval': decimal_bounds(
                bound['original_tube_weight_lower'], bound['original_tube_weight_upper'], 10),
            'band_weight_strict_lower': '1/5',
            'finite_gap_interval': checked['finite_field_gap_interval'],
            'source_spectral_cluster_indices': cluster,
            'full256_vector_integer_certificate_rechecked': True,
            'clock_probability_certificate_uses_no_eigenvector_gap': True,
            'additional_native_128_bin_tape_executed': False,
            'rounded_source_apparatus_error_included_in_spectral_table': False,
            'infinite_ordered_state_identified_with_finite_selected_ground': False,
            'sharp_scaling_mass_shell_or_physical_mass_proved': False}
