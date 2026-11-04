"""Exact finite identities; no continuum or physical particle identification."""
from fractions import Fraction as F

import sympy as s

from model.common_source.chiral_candidate import action
from model.radial_defect_gap import bounds


def annihilators(modes=4):
    if type(modes) is not int or not 1 <= modes <= 8:
        raise ValueError('explicit CAR matrices require one to eight modes')
    result = []
    for mode in range(modes):
        matrix = s.zeros(2**modes)
        for state in range(2**modes):
            if state & (1 << mode):
                matrix[state ^ (1 << mode), state] = (-1)**((state & ((1 << mode)-1)).bit_count())
        result.append(matrix)
    return result


def block(mass, coupling):
    c = annihilators()
    number = sum((a.H*a for a in c), s.zeros(16))
    pairing = sum((c[i].H*c[j].H+c[j]*c[i] for i,j in ((0,1),(2,3))), s.zeros(16))
    perturbation = (number-2*s.eye(16))**2
    return mass*pairing+coupling*perturbation, number, pairing, perturbation


def exact_certificate():
    labels = action.ordering()['all_left_components']
    q = {'Q': 1, 'u^c': -1, 'd^c': -1, 'L': -3, 'e^c': 3}
    charge15 = s.diag(*(q[row['internal_type']] for row in labels))
    assert all(charge15*g == g*charge15 for g in action.generators(families=False))
    channels = action._representation()[4]
    for channel in channels.values():
        for tensor in channel['components']:
            assert charge15*tensor+tensor*charge15.T == s.zeros(15)
    neutrino = next(i for i,row in enumerate(labels)
                   if row['internal_type'] == 'L' and row['electric_Q'] == '0')
    majorana = s.zeros(15)
    majorana[neutrino,neutrino] = 1
    assert charge15*majorana+majorana*charge15.T == -6*majorana
    charges90 = [int(charge15[i,i]) for _ in range(6) for i in range(15)]
    m, u = s.symbols('m u', real=True)
    z = s.Symbol('z')
    h, n, pairing, v = block(m,u)
    expected = z**5*(z-4*u)*((z-2*u)**2-4*(u*u+m*m))*((z-u)**2-m*m)**4
    assert s.expand(h.charpoly(z).as_expr()-expected) == 0
    parity = s.diag(*((-1)**i.bit_count() for i in range(16)))
    charge = s.diag(*(sum(((-3,3,-3,3)[j] if i & (1<<j) else 0) for j in range(4))
                      for i in range(16)))
    assert h*parity == parity*h and h*charge == charge*h
    assert v*pairing != pairing*v and h*n != n*h
    c = annihilators()
    # A commutator with a creation operator has occupation-dependent cubic terms.
    commutator = v*c[0].H-c[0].H*v
    assert commutator == c[0].H*(2*n-3*s.eye(16))
    return {
        'native_currents_commuting_with_core_number_interaction': 12,
        'charge_3B_minus_3L_by_native_type': q,
        'charged_Higgs_tensors_preserve_this_charge': True,
        'charge_on_core_Fock_space_minimum': sum(q for q in charges90 if q < 0),
        'charge_on_core_Fock_space_maximum': sum(q for q in charges90 if q > 0),
        'Majorana_neutrino_pair_changes_this_charge': -6,
        'core_annihilators': 90,
        'core_number_interaction': 'V=lambda*(N_W-45)^2',
        'normal_ordered_quartic_coefficient': '2*lambda*sum_{i<j} n_i*n_j',
        'core_number_conserved_by_pairing': False,
        'fermion_parity_conserved': True,
        'four_mode_characteristic_polynomial': 'z^5*(z-4u)*((z-2u)^2-4*(u^2+m^2))*((z-u)^2-m^2)^4',
        'four_mode_ground_energy': '2u-2*sqrt(u^2+m^2), m>0,u>=0',
        'four_mode_lowest_odd_gap': '2*sqrt(u^2+m^2)-u-m',
        'four_mode_single_fermion_poles': ['2*sqrt(u^2+m^2)-u-m', '2*sqrt(u^2+m^2)-u+m'],
        'addition_plus_removal_residues': ['(1+m/sqrt(u^2+m^2))/2', '(1-m/sqrt(u^2+m^2))/2'],
        'ground_connected_cross_pair_density': '-u/(4*sqrt(u^2+m^2))',
        'four_mode_block_is_full_core_ground_sector': False,
        'quartic_dynamics_nontrivial_within_fixed_charge_sector': True}


def radius_uniform_bound(gain, coupling, p, time, steps):
    values = (gain,coupling,p,time)
    if any(isinstance(x, (float,bool)) for x in values):
        raise ValueError('exact rational bounds required')
    g,u,p,t = map(F, values)
    if g < 0 or u < 0 or not 0 <= p < 1 or t < 0 or type(steps) is not int or steps < 1:
        raise ValueError('invalid pulse-bound domain')
    norm = 45*(1+g)+2025*u
    transient = 4*norm*t/(steps*(1-p))
    stationary = (8/(1-p)+4)*norm*norm*t*t/steps
    return {'core_Fock_norm_upper': str(norm), 'arbitrary_initial_CQ_error_upper': str(transient+stationary),
            'stationary_product_memory_error_upper': str(stationary),
            'full_history_error_squared_upper': str(2*(transient+stationary)),
            'stationary_product_full_history_error_squared_upper': str(2*stationary),
            'full_history_rate': 'O(n^-1/2), capped by trace norm2',
            'initially_correlated_memory_ideal': 'retain p(m,history) and apply the mean unitary separately to each initial conditional core state rho_m',
            'uniform_in_cube_radius': True, 'uniform_in_lattice_spacing': False,
            'comparison_algebra': 'the same90 localized annihilators and an untouched finite reference'}


def localization_bound(distance):
    tail = bounds.tail_probability_upper(distance, mass=64)
    # Projection onto |x|_infinity < distance is denoted A_L.
    # ||W W* - A_L W W* A_L|| <= 2 sqrt(tail).
    return {'distance': distance, 'one_particle_tail_probability_upper': str(tail),
            'projector_truncation_norm_squared_upper': str(4*tail),
            'core_number_truncation_norm_squared_upper': str(180**2*4*tail),
            'interaction_truncation_error': '||V-V_L|| <= lambda*||N_W-N_L||*(90+||N_W-N_L||)',
            'finite_range': False, 'exponentially_localized_about_prescribed_origin': True}
