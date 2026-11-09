"""Exact outward bounds for the actual CAR response near a magnetic node."""
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt

import sympy as s

from .local_matter import annihilator
from .local_spectral import dimer_source
from .spectral_window_access import _radical_interval


FOURIER_TAIL_CONSTANT = 6709248
LN2_LOWER = F(2, 3)


def _rational(value, *, positive=False):
    if type(value) not in (int, F):
        raise ValueError('exact rational input required')
    value = F(value)
    if value < 0 or (positive and not value):
        raise ValueError('nonnegative rational in the admitted domain required')
    return value


def _integer(value, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError('integer in the admitted domain required')
    return value


def _ceil(value):
    return -((-value.numerator)//value.denominator)


def _pow2(exponent):
    return F(2**exponent) if exponent >= 0 else F(1, 2**(-exponent))


def ceil_log2(value):
    """Exact smallest integer j with value<=2**j; no floating logarithm."""
    value = _rational(value, positive=True)
    j = value.numerator.bit_length()-value.denominator.bit_length()
    if _pow2(j) < value:
        j += 1
    if not _pow2(j-1) < value <= _pow2(j):
        raise ArithmeticError('integer logarithm enclosure failed')
    return j


def fourier_tail(radius, duration, coupling=1):
    """Uniform 16-by-16 cell tail from the full odd-field cut estimate."""
    radius = _integer(radius, 1)
    time = _rational(duration); c = _rational(coupling, positive=True)
    if not time:
        return F(0)
    # exp(49ct)<=2**ceil(49ct/ln(2))<=2**ceil(147ct/2).
    exponent = _ceil(F(147, 2)*c*time)-radius
    return FOURIER_TAIL_CONSTANT*_pow2(exponent)


def fourier_cutoff(duration, accuracy, coupling=1):
    time = _rational(duration); epsilon = _rational(accuracy, positive=True)
    c = _rational(coupling, positive=True)
    if epsilon > 1:
        raise ValueError('accuracy must be at most one')
    radius = ceil_log2(F(FOURIER_TAIL_CONSTANT)/epsilon)+_ceil(F(147, 2)*c*time)
    tail = fourier_tail(radius, time, c)
    if radius < 1 or tail > epsilon:
        raise ArithmeticError('uniform source Fourier tail exceeds accuracy')
    return {'cell_radius': radius, 'original_response_tail_upper': tail,
            'requested_accuracy': epsilon, 'source_duration': time,
            'source_coupling': c, 'matrix_components': 16,
            'ground_clustering_assumed': False}


def lorentzian_difference(distance_l1, eta, coupling=1, accuracy=None):
    """A rational upper bound for ||Q_eta(k)-Q_eta(k_*)|| in one Bloch lift."""
    d = _rational(distance_l1); eta = _rational(eta, positive=True)
    c = _rational(coupling, positive=True)
    epsilon = min(d, F(1)) if accuracy is None else _rational(accuracy, positive=True)
    if epsilon > 1:
        raise ValueError('Fourier accuracy must be at most one')
    if not d:
        return F(0)
    degree_offset = ceil_log2(F(FOURIER_TAIL_CONSTANT)/epsilon)+1
    # Integrate the integer cutoff bound; retain the physical-cell frame term.
    bound = 2*epsilon+2*d*((1+epsilon)*(degree_offset+F(147, 2)*c/eta)+1)
    return min(F(1), bound)


def dyadic_log_remainder(power):
    """Uniform remainder for every 0<distance<=2**(-power), including jumps."""
    power = _integer(power)
    return _pow2(-power)*(12+4*(ceil_log2(FOURIER_TAIL_CONSTANT)+power))


def conditional_window(atom_lower, packet_radius_l1, velocity=1):
    """Conditional finite physical-window budget; atom_lower is NOT derived."""
    z = _rational(atom_lower, positive=True)
    k = _rational(packet_radius_l1, positive=True); v = _rational(velocity, positive=True)
    # Exact comparison with chi=(31+sqrt(577))/(95+sqrt(577)).
    # Rearrangement is squared only when its left side is positive.
    if z >= 1 or (95*z-31 > 0 and (95*z-31)**2 > 577*(1-z)**2):
        raise ValueError('hypothesized source atom must respect the paired-star upper bound')
    tau = max(F(1), 2352*v*k/z)
    power = 0
    while dyadic_log_remainder(power) > z/8:
        power += 1
    spacing = _pow2(-power)/k
    n = max(1, isqrt(_ceil(4/z)))
    while F(1, 1+n*n) > z/4:
        n += 1
    outside = F(1, 1+n*n)
    dynamic = 294*v*k/tau
    total = dynamic+dyadic_log_remainder(power)
    true_score = z-total
    service = z/8
    measured_score = true_score-service
    inferred_score = measured_score-service
    original_floor = (true_score-outside)/(1-outside)
    inferred_floor = (inferred_score-outside)/(1-outside)
    if total > z/4 or original_floor < z/2 or inferred_floor < z/4:
        raise ArithmeticError('conditional original-weight margin failed')
    return {'assumed_nodal_atom_lower': z,
            'actual_nodal_atom_lower_proved': False,
            'admitted_atom_ceiling_exact': '(31+sqrt(577))/(95+sqrt(577))',
            'fixed_normalized_source_c': F(1), 'fixed_normalized_source_U': F(8),
            'same_normalized_reference_at_every_spacing_required': True,
            'packet_radius_l1_coefficient': k, 'physical_velocity_scale': v,
            'physical_Lorentzian_scale': tau, 'spacing_upper': spacing,
            'dyadic_distance_power': power,
            'joint_frequency_momentum_error_upper': total,
            'original_Lorentzian_score_lower': true_score,
            'window_multiplier': n, 'physical_energy_window': [F(0), n*tau],
            'original_window_weight_lower': original_floor,
            'complete_score_error_allowance': service,
            'measured_score_lower_if_allowance_met': measured_score,
            'inferred_original_window_weight_lower': inferred_floor,
            'zero_energy_included': True,
            'native_growing_readout_error_met': False,
            'positive_energy_annulus_proved': False,
            'dispersion_or_particle_identification_proved': False}


def _simple(matrix):
    return matrix.applyfunc(s.simplify)


@lru_cache(maxsize=1)
def _finite_control():
    """Complete actual interacting dimer spectral Gram, without a pole fit."""
    h, rho, _, projections = dimer_source()
    modes = tuple(annihilator(2, j) for j in range(4))
    size = len(modes)
    signed_moments = []
    for degree in range(3):
        kernel = s.zeros(size)
        for energy, projector in projections:
            for i, ai in enumerate(modes):
                for j, aj in enumerate(modes):
                    kernel[i, j] += energy**degree*s.trace(rho*(
                        ai*projector*aj.H+(-1)**degree*aj.H*projector*ai))
        signed_moments.append(_simple(kernel))
    m0, m1, m2 = signed_moments
    expected_h = s.zeros(4)
    expected_h[0, 2] = expected_h[2, 0] = s.Rational(1, 16)
    expected_h[1, 3] = expected_h[3, 1] = s.Rational(1, 16)
    if m0 != s.eye(4) or m1 != expected_h or m2-m1*m1 != s.eye(4)/4:
        raise ArithmeticError('full source retarded contraction/memory moments failed')
    q = s.zeros(4)
    for energy, projector in projections:
        for i, ai in enumerate(modes):
            for j, aj in enumerate(modes):
                q[i, j] += s.trace(rho*(ai*projector*aj.H+aj.H*projector*ai))/(1+energy**2)
    q = _simple(q)
    scalar = s.simplify(q[0, 0])
    if q != scalar*s.eye(4) or not 0 < scalar < 1:
        raise ArithmeticError('complete source Lorentzian is not the expected effect')
    # The individual charge Gram has total local probability 1/2. Its omission
    # cannot be hidden by renormalizing the surviving record.
    add = _simple(s.Matrix(4, 4, lambda i, j: s.trace(rho*modes[i]*modes[j].H)))
    remove = _simple(s.Matrix(4, 4, lambda i, j: s.trace(rho*modes[j].H*modes[i])))
    if add+remove != s.eye(4) or any(add[j, j] != s.Rational(1, 2) for j in range(4)):
        raise ArithmeticError('complete addition/removal normalization failed')
    return {'actual_Fock_dimension': h.rows, 'actual_CAR_modes': size,
            'source_c': F(1, 8), 'source_U': F(1),
            'both_signed_charge_measures_retained': True,
            'retarded_total_matrix_mass_identity': True,
            'retarded_second_moment_minus_first_square': F(1, 4),
            'retarded_response_is_not_assumed_unitary': True,
            'eta': F(1), 'original_Lorentzian_interval': _radical_interval(scalar),
            'each_local_addition_and_removal_weight': F(1, 2),
            'finite_dimer_replaces_thermodynamic_reference': False,
            'native_tape_executed': False}


def _detach(value):
    if type(value) is F:
        return F(value.numerator, value.denominator)
    if type(value) is dict:
        return {k: _detach(v) for k, v in value.items()}
    if type(value) is list:
        return [_detach(v) for v in value]
    if type(value) is tuple:
        return tuple(_detach(v) for v in value)
    return value


def finite_control():
    return _detach(_finite_control())


def certificate():
    """Source theorem controls and explicitly hypothetical atom-dependent plan."""
    return {'source': 'same centered full CAR source, fixed c=1,U=8 selected reference',
            'periodic_cell_matrix_components': 16,
            'Fourier_tail_constant': FOURIER_TAIL_CONSTANT,
            'natural_log_2_lower_bound': F(LN2_LOWER.numerator, LN2_LOWER.denominator),
            'cutoff_example': fourier_cutoff(F(1, 4), F(1, 1000)),
            'joint_modulus_example': lorentzian_difference(F(1, 1000000), F(1, 10)),
            'conditional_example': conditional_window(F(1, 100), 1),
            'actual_source_control': finite_control(),
            'atom_positivity_equivalent_to_bounded_window_survival_on_fixed_reference': True,
            'actual_nodal_atom_positivity_proved': False,
            'physical_mass_assigned': None,
            'sharp_dispersion_proved': False,
            'growing_native_instrument_executed': False}
