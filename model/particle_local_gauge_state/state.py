"""Exact density planner for a gauge-twirled product of finite link words.

The full mixed source, canonical charged reference and action are declared.
The planner proves local rough-plaquette bounds; it does not calculate a
local Wilson low-mode bound, physical mass gap or relativistic limit.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path

import sympy as s

from model.particle_field_preparation import preparation
from model.particle_observer_outputs import field_route


SOURCE_INPUTS = (
    'model/particle_field_preparation/preparation.py',
    'model/particle_field_preparation/PROOF.md',
    'model/particle_field_preparation/REFINEMENT.md',
    'model/particle_observer_outputs/field_route.py',
    'model/particle_observer_outputs/FIELD_ROUTE.md',
    'model/particle_observer_outputs/READOUT_CONTRACT.md',
    'model/interacting_extension/chiral_lattice/source_data.py',
    'model/native_repair_spin_limit/PROOF.md',
    'model/spin_reference_limit/PROOF.md',
)


def integer(value, name, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError(name+' requires an integer at least '+str(minimum))
    return value


def rational(value, name, positive=False):
    if type(value) not in (int, str, F):
        raise ValueError(name+' requires exact rational data')
    result = F(value)
    if result < 0 or (positive and not result):
        raise ValueError(name+' has the wrong sign')
    return result


def link_bound(root_degree, threshold_bits):
    """One normalized link, no volume-dependent normalization denominator.

    The threshold delta is for a plaquette; the bad-link threshold is
    delta/4. Frobenius distance bounds the represented operator distance.
    """
    n = integer(root_degree, 'root degree', 1)
    r = integer(threshold_bits, 'threshold bits', 1)
    L = n*n
    delta = F(1, 2**r)
    cover_bits = (8*n).bit_length()
    exponent = 2+450*cover_bits-L//(480*4**r)
    return {
        'root_degree': n, 'per_link_word_degree': L,
        'plaquette_operator_threshold': delta,
        'bad_link_operator_threshold': delta/4,
        'normalization_squared_lower': {
            'rational_factor': F(59, 60)**2,
            'integer_base': 1+8*n, 'inverse_power': 450},
        'bad_link_probability_upper': {'base': 2, 'exponent': min(0, exponent)},
        'uncapped_binary_exponent': exponent,
        'electric_energy_per_link_upper': F(29, 10)*L*L,
        'normalization_contains_total_volume': False,
    }


@lru_cache(maxsize=16)
def density_plan(target_bits=20, degree=128, spacing='1'):
    """Fix L and kappa once for every periodic volume and every model time."""
    bits = integer(target_bits, 'target roughness bits', 1)
    m = integer(degree, 'Legendre degree')
    a = rational(spacing, 'spacing', True)
    eta = F(1, 2**bits)
    r = 1
    while F(3, 2*4**r) > eta/(4*48000):
        r += 1
    qbits = 1
    while F(24, 2**qbits) > eta/(4*48000):
        qbits += 1
    n = 1
    while True:
        link = link_bound(n, r)
        if link['uncapped_binary_exponent'] <= -qbits:
            break
        n *= 2
    delta, q = F(1, 2**r), F(1, 2**qbits)
    magnetic_density = 3*delta**2/2+24*q
    electric_density = 3*link['electric_energy_per_link_upper']
    poly = field_route.polynomial_bounds(m, a)
    kinetic_density = 60*poly['global_one_body_H_norm_upper']
    kappa = 96000*(electric_density+2*kinetic_density)/eta
    local_bound = 48000*(magnetic_density+(electric_density+2*kinetic_density)/kappa)
    if magnetic_density > eta/96000 or local_bound > eta:
        raise ArithmeticError('the chosen densities exceed the roughness allowance')
    return {
        'target_probability_upper': eta, 'target_bits': bits,
        'Legendre_degree': m, 'spacing': a,
        'single_link': link,
        'chosen_bad_link_probability_upper': q,
        'initial_rough_plaquette_probability_at_delta_upper': min(F(1), 4*q),
        'initial_cost_per_plaquette_upper': delta**2/2+8*q,
        'initial_magnetic_cost_per_cell_upper': magnetic_density,
        'initial_electric_energy_per_cell_upper': electric_density,
        'full_Fock_kinetic_norm_per_cell_upper': kinetic_density,
        'declared_plaquette_coefficient': kappa,
        'all_time_rough_plaquette_probability_upper': local_bound,
        'roughness_operator_threshold': F(1, 40),
        'state': 'site-gauge twirl of independently normalized link-polynomial product density tensor maximally mixed full Fock density',
        'fermion_modes_per_cell': 60,
        'fermion_density': 'I_Fock / 2^(60 V), with canonical conjugate reference',
        'translation_invariant_source_marginal': True,
        'no_orientation_factor_three_division': True,
        'local_word_degree_independent_of_volume': True,
        'plaquette_coefficient_independent_of_volume': True,
        'bound_holds_for_every_model_time': True,
        'field_only_Gauss_singlet_support_claimed': False,
        'global_Wilson_gap_probability_bound_claimed': False,
        'local_Wilson_low_mode_bound_claimed': False,
        'physical_mass_or_species_claimed': False,
        'interacting_relativistic_refinement_claimed': False,
        'measured_particle_inputs': [],
    }


def volume_instance(plan, side, extra_degree=0):
    """Full finite word/Fock carrier and retained canonical purification.

    The source marginal is exactly defined by Haar matrix integrals on a
    finite code. This certificate does not allocate that huge density.
    """
    side = integer(side, 'periodic side', 1)
    extra = integer(extra_degree, 'extra total word degree')
    volume = side**3
    links = 3*volume
    L = plan['single_link']['per_link_word_degree']
    initial = links*L
    K = initial+extra
    return {
        'periodic_side': side, 'cells': volume, 'links_and_positive_plaquettes': links,
        'per_link_initial_word_degree': L,
        'initial_total_word_degree_upper': initial,
        'evolution_total_word_cutoff': K,
        'gauge_code_dimension_upper': {'base': 1+450*links, 'exponent': K},
        'full_Fock_dimension': {'base': 2, 'exponent': 60*volume},
        'purifier_dimension_upper': {
            'factors': [{'base': 451, 'exponent': initial},
                        {'base': 2, 'exponent': 60*volume}]},
        'same_kappa': plan['declared_plaquette_coefficient'],
        'same_all_time_local_probability_bound': plan['all_time_rough_plaquette_probability_upper'],
        'source_H': 'P_K [H_E+dGamma(h_m)-kappa sum_p Re Tr_M15(U_p)/15] P_K',
        'detector_effects': ['P_K 1_{||R(U_p)-I||>1/40} P_K', 'I_code minus rough effect'],
        'all_Fock_occupations_retained': True,
        'complete_detector_probability_sum': F(1),
        'purifier_gauge_representation': 'conjugate of the source support representation',
        'reference_Hamiltonian_during_observation': 0,
        'joint_source_reference_Gauss_invariant': True,
        'source_only_Gauss_projection_performed': False,
        'autonomous_cooling_or_repair_preparation_claimed': False,
        'large_density_or_rough_detector_matrix_materialized': False,
        'proper_clock_or_GeV_units_assigned': False,
        'operator_error_for_compiled_rough_detector_included': False,
    }


@lru_cache(maxsize=1)
def exact_one_link_control():
    """L=1 source character norm and complete independent-endpoint twirl.

    This is a single open link with two independent endpoint gauges. A
    periodic lattice's site twirl correlates adjacent links and is not
    replaced by independent link twirls in the state family above.
    """
    rows = preparation.native_casimirs()
    norm = s.Rational(1, 4)+10*s.Rational(1, 60)**2
    if norm != s.Rational(91, 360):
        raise ArithmeticError('the native character normalization differs')
    outcomes = [{'irrep': 'trivial', 'dimension': 1, 'matrix_coefficient_states': 1,
                 'original_probability': s.Rational(90, 91),
                 'density_eigenvalue': s.Rational(90, 91), 'electric_energy': s.S.Zero}]
    for row in rows:
        d = row['dimension']
        for conjugate in (False, True):
            outcomes.append({
                'irrep': {'q6Y': -row['q6Y'] if conjugate else row['q6Y'],
                          'native_block_dimension': d, 'conjugate': conjugate},
                'dimension': d, 'matrix_coefficient_states': d*d,
                'original_probability': s.Rational(1, 910),
                'density_eigenvalue': s.Rational(1, 910*d*d),
                'electric_energy': row['trace_metric_C2']/2})
    probability = sum(row['original_probability'] for row in outcomes)
    energy = sum(row['original_probability']*row['electric_energy'] for row in outcomes)
    rank = sum(row['matrix_coefficient_states'] for row in outcomes)
    if probability != 1 or energy != s.Rational(33, 9100) or rank != 119:
        raise ArithmeticError('complete endpoint-twirled native character law differs')
    return {'single_link_normalization_squared': norm,
            'positive_state_support_dimension': rank,
            'complete_irrep_records': outcomes,
            'outside_support_probability': s.S.Zero,
            'all_original_probabilities_sum': probability,
            'electric_energy_expectation': energy,
            'nontrivial_reference_charge_probability': s.Rational(1, 91),
            'source_only_Gauss_singlet': False,
            'periodic_joint_twirl_replaced_by_independent_link_twirl': False,
            'local_roughness_target_certified_for_L1': False}


def serial(value):
    if isinstance(value, (F, s.Basic)):
        return str(value)
    if isinstance(value, dict):
        return {str(k): serial(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [serial(v) for v in value]
    return value


def certificate():
    plan = density_plan()
    instances = [volume_instance(plan, side, extra_degree=4) for side in (1, 2, 4, 16)]
    if len({instance['same_kappa'] for instance in instances}) != 1:
        raise ArithmeticError('the finite volume family changed its local action')
    return serial({
        'schema': 'oph.local_gauge_state.v1', 'source_inputs': list(SOURCE_INPUTS),
        'uniform_density_plan': plan, 'finite_volume_instances': instances,
        'exact_single_open_link_control': exact_one_link_control(),
        'source_density_construction': 'finite exact Haar twirl, no smooth-outcome conditioning',
        'preparation_service': 'canonical conjugate purification via a complete controlled unitary',
        'comparison_to_occupied_origin_preparation_claimed': False,
        'remaining_local_operator_premise': 'volume-uniform low-Wilson-mode and local sign-response control on this same quantum state and clock',
        'numerical_exponential_used_as_proof': False,
    })
