"""Exact native small-gap family and complete coordinate-symbol controls.

Small spectral weight does not bound inverse-gap electric response. The
finite native matrices here are coordinate symbols, not point-localized
normalizable gauge preparations or measured particle mass spectra.
"""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import json

import sympy as s

from model.particle_local_wilson_density import density


ROOT = Path(__file__).resolve().parents[2]
SOURCE_INPUTS = (
    'model/particle_local_wilson_density/density.py',
    'model/particle_local_wilson_density/PROOF.md',
    'model/particle_local_gauge_state/PROOF.md',
    'model/particle_observer_outputs/field_route.py',
    'model/particle_observer_outputs/FIELD_ROUTE.md',
    'model/interacting_extension/chiral_lattice/PROOF.md',
    'model/interacting_extension/chiral_lattice/source_data.py',
    'model/common_source/overlap_regulator.py',
    'model/observer_particles/results/representation_identification_operators.json',
)


def positive_gap(value):
    if type(value) not in (int, str, F, s.Rational):
        raise ValueError('the gap requires exact rational data')
    d = s.Rational(str(value))
    if not 0 < d < s.Rational(1, 2):
        raise ValueError('the small source gap must lie strictly between zero and one half')
    return d


@lru_cache(maxsize=1)
def native_blocks():
    """All sixty source columns, with the actual four weak doublets retained."""
    clean = density.clean
    raw = json.loads((ROOT/SOURCE_INPUTS[-1]).read_text())
    weak = s.zeros(15)
    for i, j, value in raw['even_identification_operators']['weak_C2']['entries']:
        weak[i, j] = s.Rational(4, 3)*s.Rational(value)
    weak = clean(weak)
    singlet = s.eye(15)-weak
    tau = []
    for axis in (0, 1):
        generator = clean(density.native_weak_link(axis, 0, 1)-singlet)
        if generator.H != -generator or clean(generator*generator) != -weak:
            raise ArithmeticError('the actual native weak generator differs')
        if -s.trace(generator*generator) != 8:
            raise ArithmeticError('the source trace metric of the weak direction differs')
        tau.append(clean(-s.I*generator))
    gamma = density.gamma_data()
    beta = clean(s.kronecker_product(gamma[3], s.eye(15)))
    pw = clean(s.kronecker_product(s.eye(4), weak))
    ps = s.eye(60)-pw
    axes = tuple(clean(s.kronecker_product(s.I*gamma[3]*gamma[mu], tau[mu]))
                 for mu in (0, 1))
    psmall = clean((pw-axes[0]*axes[1])/2)
    plarge = clean((pw+axes[0]*axes[1])/2)
    if density.comm(*axes) != s.zeros(60):
        raise ArithmeticError('the two spin/internal involutions do not commute')
    for p, rank in ((psmall, 16), (plarge, 16), (ps, 28)):
        if clean(p*p) != p or p != p.H or p.trace() != rank:
            raise ArithmeticError('the complete native block projector differs')
    if clean(psmall+plarge+ps) != s.eye(60):
        raise ArithmeticError('the native block menu is incomplete')
    for axis in axes:
        if clean(axis*axis) != pw or clean(beta*axis+axis*beta) != s.zeros(60):
            raise ArithmeticError('the source Clifford/internal block relation differs')
    return {'beta': beta, 'weak': pw, 'singlet': ps,
            'small': psmall, 'large': plarge, 'axes': axes,
            'internal_generator_trace_metric_norm_squared': 8,
            'opposite_two_link_direction_norm_squared': 16}


def symbol_family(gap):
    """Exact source H_W and electric derivative at the stated native links."""
    d = positive_gap(gap)
    c = (1-d)/2
    sine = s.sqrt(1-c*c)
    b = native_blocks()
    beta, a1, a2 = b['beta'], *b['axes']
    h = density.clean(d*beta*b['weak']-beta*b['singlet']+sine*(a1+a2))
    radius = s.sqrt(3+2*d)
    sign = density.clean(beta*b['small']-beta*b['singlet']
                         +(d*beta+sine*(a1+a2))*b['large']/radius)
    dh = density.clean(c*(a1-a2))
    derivative = density.clean((1-d)*a1*b['small']/d)
    if density.clean(h*h-(d*d*b['small']+b['singlet']+(3+2*d)*b['large'])) != s.zeros(60):
        raise ArithmeticError('the exact native spectral decomposition differs')
    if density.clean(sign*sign) != s.eye(60) or density.clean(density.comm(h, sign)) != s.zeros(60):
        raise ArithmeticError('the native sign formula differs')
    if density.clean(sign*derivative+derivative*sign) != s.zeros(60):
        raise ArithmeticError('the sign derivative fails its involution tangent')
    if density.clean(density.comm(h, derivative)-density.comm(sign, dh)) != s.zeros(60):
        raise ArithmeticError('the sign derivative fails its commutation tangent')
    # The curve changes two independent links with native trace-metric
    # squared speeds8+8. The unit electric direction is therefore (1/4)d/du.
    electric_derivative = derivative/4
    hs = s.simplify(s.trace(electric_derivative.H*electric_derivative))
    if hs != (1-d)**2/d**2:
        raise ArithmeticError('the unit-electric derivative norm differs')
    return {'gap': d, 'cos_theta': c, 'sin_theta': sine,
            'kernel': h, 'sign': sign, 'kernel_u_derivative': dh,
            'sign_u_derivative': derivative,
            'unit_electric_sign_derivative': electric_derivative,
            'squared_kernel_spectral_ranks': {d*d: 16, s.S.One: 28, 3+2*d: 16},
            'low_mode_rank': 16,
            'unit_electric_sign_derivative_norm': (1-d)/(4*d),
            'unit_electric_sign_derivative_HS_squared': hs}


def pulse_records(gap):
    """Complete native coordinate-symbol record law for two sign readouts.

    Each of three invariant blocks and all four sign pairs is retained.
    The auxiliary input is I60/60. The source links are at a specified
    coordinate, not asserted to be normalizable quantum field states.
    """
    d = positive_gap(gap)
    c, sine = (1-d)/2, s.sqrt(1-((1-d)/2)**2)
    cu, su = s.sqrt(1-d*d), d  # u=arcsin(delta), with 0<u<2 delta.
    B = 1-2*c*cu
    difference, total = 2*c*su, 2*sine*cu
    small_corr = B/s.sqrt(B*B+difference*difference)
    large_corr = (d*B+2*sine*total)/(s.sqrt(3+2*d)*s.sqrt(B*B+total*total))
    correlations = [('small', 16, small_corr), ('large', 16, large_corr), ('singlet', 28, s.S.One)]
    records = []
    for label, rank, corr in correlations:
        for first in (-1, 1):
            for last in (-1, 1):
                records.append({'block': label, 'first_sign': first, 'last_sign': last,
                                'original_probability': s.Rational(rank, 240)*(1+first*last*corr)})
    if s.simplify(sum(r['original_probability'] for r in records)) != 1:
        raise ArithmeticError('the original native symbol records are incomplete')
    return {'gap': d, 'pulse_angle': 'arcsin(delta)',
            'unit_electric_length': '4 arcsin(delta)',
            'all_native_columns': 60, 'original_records': records,
            'small_block_sign_correlation': small_corr,
            'small_block_original_probability': s.Rational(4, 15),
            'small_block_conditional_flip_probability': (1-small_corr)/2,
            'complete_original_probability_sum': s.S.One,
            'pulse_and_fixed_block_menu_Gauss_preservation_proved': False,
            'coordinate_symbol_not_normalizable_preparation': True}


def rare_packet_counter(gap):
    d = positive_gap(gap)
    if d > s.Rational(1, 4):
        raise ValueError('the explicit finite-pulse lower bound uses delta<=1/4')
    return {'gap': d, 'rare_branch_original_weight': d,
            'rare_symbol_low_mode_trace': 16*d,
            'rare_symbol_normalized_low_effect_probability': 4*d/15,
            'rare_symbol_unit_electric_derivative_HS_squared': (1-d)**2/d,
            'complete_probe_original_rare_flip_lower': d/60,
            'mean_square_record_change_per_unit_electric_length_squared_lower': 1/(3840*d),
            'rare_branch_conditioning_used': False,
            'normalizable_packet_realization': 'smooth neighborhoods and finite-word approximation with all preparation flags retained',
            'packet_word_degree_or_uniform_electric_energy_bound_supplied': False,
            'counterfamily_is_the_existing_prepared_density': False,
            'fixed_duration_interacting_response_impossible_claimed': False}


@lru_cache(maxsize=1)
def exact_native_control():
    d = s.Rational(1, 4)
    c = (1-d)/2
    sine = s.sqrt(1-c*c)
    shifts = (density.native_weak_link(0, c, sine),
              density.native_weak_link(1, c, sine), s.eye(15))
    actual = density.decomposition(shifts)
    closed = symbol_family(d)
    if actual['kernel'] != closed['kernel']:
        raise ArithmeticError('the explicit source links do not realize the claimed family')
    # Direct differentiation of the two normalized2x2 sign blocks.
    u = s.symbols('u', real=True)
    B, D = 1-2*c*s.cos(u), 2*c*s.sin(u)
    R = s.sqrt(B*B+D*D)
    small_beta = s.simplify(s.diff(B/R, u).subs(u, 0))
    small_axis = s.simplify(s.diff(D/R, u).subs(u, 0))
    if small_beta != 0 or small_axis != (1-d)/d:
        raise ArithmeticError('the analytic small-block derivative differs')
    return {'exact_source_link_kernel_match': True,
            'native_spectral_control': {str(k): v for k, v in closed['squared_kernel_spectral_ranks'].items()},
            'small_block_beta_derivative': small_beta,
            'small_block_axis_derivative': small_axis,
            'unit_electric_HS_squared': closed['unit_electric_sign_derivative_HS_squared'],
            'all_sixty_columns_retained': True,
            'two_time_complete_coordinate_symbol_records': pulse_records(d)}


def certificate():
    return density.serial({'schema': 'oph.sign_response_boundary.v1',
        'source_inputs': list(SOURCE_INPUTS),
        'exact_native_control': exact_native_control(),
        'rare_packet_sequence': [rare_packet_counter(s.Rational(1, n)) for n in (4, 16, 256, 4096)],
        'positive_sufficient_input': 'inverse-gap-weighted electric matrix-element control, or a proved smooth-regulator local comparison',
        'low_mode_density_alone_controls_electric_derivative': False,
        'specific_prepared_state_failure_proved': False,
        'physical_mass_or_species_assignments': [],
        'complete_interacting_local_Duhamel_bound_proved': False})
