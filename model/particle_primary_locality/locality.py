"""Instantiate the existing all-state cone on the primary finite source.

The product link cutoff preserves each interaction's support and norm
bound. Gauss support is imposed on states and invariant instruments;
no globally projected Hamiltonian is declared to be a local interaction.
"""
from __future__ import annotations

from fractions import Fraction as F
from itertools import product

from model.particle_common_contract import contract
from model.particle_physical_boundary_state import boundary
from model.particle_physical_causality import cone
from model.particle_spatial_ground_state import ground

SOURCE_INPUTS = (
    'model/particle_common_contract/contract.py',
    'model/particle_common_contract/PROOF.md',
    'model/particle_physical_boundary_state/boundary.py',
    'model/particle_physical_boundary_state/PROOF.md',
    'model/particle_physical_causality/cone.py',
    'model/particle_physical_causality/PROOF.md',
    'model/particle_field_causality/bounds.py',
    'model/particle_field_causality/PROOF.md',
    'model/particle_spatial_ground_state/ground.py',
    'model/particle_spatial_ground_state/PROOF.md',
    'model/particle_gauss_ground_state/gauss.py',
    'model/particle_gauss_ground_state/PROOF.md',
    'model/particle_local_wilson_density/density.py',
    'model/particle_local_wilson_density/PROOF.md',
)


def _cell(value, side):
    if (not isinstance(value, (list, tuple)) or len(value) != 3
            or any(type(x) is not int or not 0 <= x < side for x in value)):
        raise ValueError('an actual interior cell of the selected box is required')
    return tuple(value)


def _step(cell, axis):
    return tuple(x+int(i == axis) for i, x in enumerate(cell))


def plaquette(side, anchor, axes=(1, 2)):
    """Complete oriented plaquette and actual outgoing-link owner support."""
    side = ground._i(side, 'box side', 2)
    anchor = _cell(anchor, side)
    if (not isinstance(axes, (list, tuple)) or len(axes) != 2
            or any(type(i) is not int or i not in range(3) for i in axes)
            or axes[0] >= axes[1]):
        raise ValueError('two increasing coordinate axes are required')
    first, second = axes
    a = _cell(_step(anchor, first), side)
    b = _cell(_step(anchor, second), side)
    far = _cell(_step(a, second), side)
    return {'anchor': anchor, 'axes': tuple(axes),
            'vertices': (anchor, a, far, b),
            'oriented_links': ((anchor, first, 1), (a, second, 1),
                               (b, first, -1), (anchor, second, -1)),
            'owner_support': tuple(sorted((anchor, a, b))),
            'normalized_ReTr_operator_norm_upper': F(1)}


def selected_source_bound(member):
    """Authenticate the exact member and instantiate its actual coefficients."""
    contract.verify_member(member)
    source = member['source']
    if source['word_projection'] != 'tensor_product_per_link_Haar_orthogonal':
        raise ValueError('a local product link cutoff is required')
    a = F(source['spacing'])
    kappa = F(source['magnetic_coefficient'])
    electric = F(source['electric_coefficient'])/(2*a)
    config = cone.family('wilson', magnetic_scale=kappa, weight='2')
    actual = cone.spacing_member(config, a, magnetic=kappa/a,
                                  electric_coefficient=electric)
    return {'family': config, 'actual_member': actual,
            'product_link_cutoff_preserves_term_support': True,
            'compression_norm_increases': False,
            'global_Gauss_projection_used_as_local_Hamiltonian': False,
            'ambient_bound_restricted_to_invariant_instruments': True,
            'all_sixty_CAR_modes_per_cell_retained': True,
            'electric_part_removed_by_onsite_interaction_picture': True,
            'per_link_word_cutoff_enters_propagation_rate': False,
            'source_clock_changed': False,
            'reference_idle_H': member['reference']['idle_H']}


def exponential_increment_upper(value, order=12):
    """Positive Taylor sum plus a geometric bound for exp(x)-1, 0<=x<=2."""
    x = ground._q(value, 'exponential argument')
    order = ground._i(order, 'Taylor order', 1)
    if x > 2 or x >= order+2:
        raise ValueError('the finite short-time exponential control requires x<=2')
    term = F(1); total = F(0)
    for j in range(1, order+1):
        term *= x/j
        total += term
    tail = term*x/F(order+1)/(1-x/F(order+2))
    return {'argument': x, 'Taylor_order': order,
            'positive_sum': total, 'tail_upper': tail,
            'exp_minus_one_upper': total+tail}


def complete_boundary_experiment(member, duration='1/10000', remote_anchor=None):
    """Actual same-input experiment: remote plaquette phase versus identity.

    Observer output is the complete origin-cell charge pointer and local
    output, with the untouched canonical reference. The remote field and
    any remote private apparatus are traced, not directly communicated.
    """
    source_bound = selected_source_bound(member)
    state_readout = boundary.common_contract_payload(member)
    source = member['source']; side = source['side']; a = F(source['spacing'])
    time = ground._q(duration, 'nonnegative source duration')
    anchor = (side-1, 0, 0) if remote_anchor is None else remote_anchor
    remote = plaquette(side, anchor, (1, 2))
    distance = min(sum(abs(x) for x in cell) for cell in remote['owner_support'])
    if distance <= 0:
        raise ArithmeticError('the remote plaquette overlaps the observation cell')
    family = source_bound['family']
    W = family['weighted_rate_times_spacing_upper']
    mu = family['log_weight_rational_lower']
    separation = a*distance
    margin = mu*separation-2*W*time
    # Exact one-cell observer support; no physical-volume growth prefactor.
    # min(1,2 exp(-margin/a)) <= min(1,2^(1-floor(margin/a))) if margin>=0.
    bits = max(0, (margin/a).__floor__()-1) if margin >= 0 else 0
    short = None
    argument = 2*W*time/a
    if 0 < argument <= 2:
        exponential = exponential_increment_upper(argument)
        prefactor = 2*exponential['exp_minus_one_upper']
        bits = max(bits, distance-cone.ceil_log2(prefactor), 0)
        short = {'exponential_control': exponential, 'coefficient': prefactor,
                 'base': 2, 'negative_exponent': distance}
    if time == 0:
        error = {'exact_zero': True, 'base': 2, 'negative_exponent': None}
    else:
        error = {'exact_zero': False, 'base': 2, 'negative_exponent': bits}
    return {'duration': time, 'source_unitary': member['clock']['Schrodinger'],
            'source_bound': source_bound,
            'observer_owner_support': ((0, 0, 0),),
            'observer_support_cell_count': 1,
            'observer_instrument': state_readout['complete_central_charge_measurement'],
            'intermediate_local_intervention_count': 0,
            'remote': remote,
            'remote_intervention': 'exp(i B_p), B_p=compressed ReTr_M15(U_p)/15',
            'comparison_intervention': 'identity',
            'remote_operation_even_and_Gauss_invariant': True,
            'graph_support_distance': distance,
            'assigned_physical_support_separation': separation,
            'rational_exponential_margin': margin,
            'finite_chain_short_time_upper': short,
            'complete_output_half_trace_upper': error,
            'both_inputs_same_selected_canonical_purification': True,
            'all_observer_records_and_quantum_reference_retained': True,
            'remote_field_and_private_record_traced': True,
            'outcome_postselection': False,
            'number_of_record_labels_multiplies_error_bound': False,
            'preparation_or_boundary_approximation_error_added': F(0),
            'full_spatial_experiment_executed': False}


def incidence_control(side=3, spacing='1', magnetic='1/7', weight='2'):
    """Enumerate actual open-box Wilson term supports, not continuum cells."""
    side = ground._i(side, 'control side', 2)
    if side > 5:
        raise ValueError('the independent finite incidence control uses side at most five')
    a = ground._q(spacing, 'spacing', True)
    k = ground._q(magnetic, 'magnetic coefficient')
    lam = ground._q(weight, 'spatial weight', True)
    if lam <= 1:
        raise ValueError('spatial weight must exceed one')
    cells = tuple(product(range(side), repeat=3))
    values = {cell: F(0) for cell in cells}
    bond_count = 0; plaquette_count = 0
    for cell in cells:
        for axis in range(3):
            if cell[axis]+1 < side:
                support = (cell, _step(cell, axis))
                for owner in support:
                    values[owner] += 2*30*lam/a
                bond_count += 1
        for first in range(3):
            for second in range(first+1, 3):
                if cell[first]+1 < side and cell[second]+1 < side:
                    support = plaquette(side, cell, (first, second))['owner_support']
                    diameter = max(sum(abs(x-y) for x, y in zip(left, right))
                                   for left in support for right in support)
                    for owner in support:
                        values[owner] += len(support)*k*lam**diameter/a
                    plaquette_count += 1
    envelope = (360*lam+48*k*lam*lam)/a
    actual = max(values.values())
    if actual > envelope:
        raise ArithmeticError('the complete finite interaction exceeds the reused cone')
    return {'side': side, 'bond_count': bond_count, 'plaquette_count': plaquette_count,
            'actual_support_weighted_norm_bound_max': actual,
            'reused_weighted_norm_envelope': envelope,
            'all_open_box_interactions_enumerated': True,
            'onsite_terms_excluded_by_interaction_picture': True}


def common_contract_payload(member):
    source = selected_source_bound(member)
    return contract.serial({
        'selected_source_locality': source,
        'complete_boundary_record_comparison': complete_boundary_experiment(member),
        'local_cutoff_and_Gauss_bound_proved_for_every_selected_member': True,
        'uniform_conservative_cone_in_assigned_coordinates': True,
        'joint_state_observable_refinement_proved': False,
        'Lorentz_covariance_or_measured_light_speed_proved': False,
        'particle_mass_dispersion_or_scattering_identified': False,
        'conditional_limit_statement': 'if compatible record laws converge, remote influence vanishes outside the conservative cone',
    })


def bind_member(member):
    payload = common_contract_payload(member)
    return {'schema': 'oph.primary_locality.v1', 'member_id': member['member_id'],
            **member['ids'], 'result': payload, 'result_sha256': contract.digest(payload)}


def verify_bound(member, receipt):
    if receipt != bind_member(member):
        raise ValueError('locality receipt differs from the selected source or result')


def common_contract_bound(member):
    return contract.bind_result(member,
        'complete local boundary records on the selected physical source',
        common_contract_payload(member))


def certificate():
    fixed_distance_member = contract.member(6)
    return contract.serial({
        'schema': 'oph.primary_locality_certificate.v1',
        'selected_member': bind_member(contract.member(0)),
        'finite_members_with_separation': [bind_member(contract.member(n)) for n in (1, 2, 3)],
        'fixed_assigned_distance_one_comparison': {
            'member_id': fixed_distance_member['member_id'], **fixed_distance_member['ids'],
            'result': complete_boundary_experiment(fixed_distance_member, remote_anchor=(64, 0, 0))},
        'independent_incidence_controls': [incidence_control(n) for n in (2, 3, 4)],
        'actual_native_bond_spectra': cone.native_wilson_bonds(),
        'full_spatial_field_run_executed': False,
        'same_source_locality_premise_proved': True,
        'interacting_relativistic_continuum_proved': False,
    })
