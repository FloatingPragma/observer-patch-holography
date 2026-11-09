"""Bind the finite energy Ward identity to the selected Wilson product cutoff.

Cell energy includes its allocated interaction terms. The source coordinate
partition is explicit; no proper metric or Lorentz stress is inferred.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product

import sympy as s

from model.particle_common_contract import contract
from model.particle_gravity_response import response
from model.particle_primary_locality import locality
from model.particle_spatial_ground_state import ground

SOURCE_INPUTS = (
    'model/particle_common_contract/contract.py',
    'model/particle_common_contract/PROOF.md',
    'model/particle_spatial_ground_state/ground.py',
    'model/particle_spatial_ground_state/PROOF.md',
    'model/particle_gauss_ground_state/PROOF.md',
    'model/particle_primary_locality/locality.py',
    'model/particle_primary_locality/PROOF.md',
    'model/particle_gravity_response/response.py',
    'model/particle_gravity_response/PROOF.md',
    'model/particle_geometric_response/PROOF.md',
    'model/particle_spatial_modular_response/PROOF.md',
    'model/particle_spatial_modular_response/regional.py',
    'model/native_gravity_identification/entropy.py',
    'model/prepared_regional_entropy/PROOF.md',
    'model/matched_einstein_limit/PROOF.md',
    'Lean/QFT/CoarseEnergyWard.lean',
)
CLAIM = 'local source-energy continuity and complete energy-current records'


def _step(cell, axis, amount=1):
    return tuple(value+amount*(axis == j) for j, value in enumerate(cell))


def _inside(cell, side):
    return all(0 <= value < side for value in cell)


def _term(kind, anchor, axes, support, norm, share):
    return {'id': (kind, tuple(anchor), tuple(axes)), 'kind': kind,
            'anchor': tuple(anchor), 'axes': tuple(axes),
            'owner_support': tuple(sorted(support)),
            'unallocated_norm_upper': norm, 'cell_share': share}


def _cell_terms(source, cell):
    """The literal local assignment, with no enumeration of the whole box."""
    side = source['side']; cell = locality._cell(cell, side)
    a = F(source['spacing']); degree = source['per_link_word_degree']
    electric = F(source['electric_coefficient'])*F(29, 10)*degree**2/a
    magnetic = F(source['magnetic_coefficient'])/a
    terms = [_term('onsite_matter', cell, (), (cell,), F(90)/a, F(1))]
    terms += [_term('electric', cell, (axis,), (cell,), electric, F(1))
              for axis in range(3)]
    for axis in range(3):
        for anchor in (cell, _step(cell, axis, -1)):
            endpoint = _step(anchor, axis)
            if _inside(anchor, side) and _inside(endpoint, side):
                terms.append(_term('bond', anchor, (axis,), (anchor, endpoint), F(30)/a, F(1, 2)))
    for axes in combinations(range(3), 2):
        for anchor in (cell, _step(cell, axes[0], -1), _step(cell, axes[1], -1)):
            if (not _inside(anchor, side)
                    or not all(_inside(_step(anchor, axis), side) for axis in axes)):
                continue
            plaquette = locality.plaquette(side, anchor, axes)
            if cell not in plaquette['owner_support']:
                raise ArithmeticError('the cell is not an owner of this plaquette term')
            terms.append(_term('plaquette', anchor, axes, plaquette['owner_support'], magnetic, F(1, 3)))
    if len({row['id'] for row in terms}) != len(terms):
        raise ArithmeticError('a cell energy contains duplicate terms')
    return terms


def _cell_data(source, cell):
    cell = locality._cell(cell, source['side'])
    terms = _cell_terms(source, cell)
    support = tuple(sorted({owner for row in terms for owner in row['owner_support']}))
    norm = sum((row['cell_share']*row['unallocated_norm_upper'] for row in terms), F(0))
    if any(sum(abs(u-v) for u, v in zip(cell, owner)) > 2 for owner in support):
        raise ArithmeticError('the actual cell energy exceeds its support radius')
    return {'cell': cell, 'terms': terms, 'owner_support': support, 'norm_upper': norm}


def cell_energy(member, cell=(0, 0, 0)):
    contract.verify_member(member)
    return _cell_data(member['source'], cell)


def _nearby(cell, side, radius=4):
    for shift in product(range(-radius, radius+1), repeat=3):
        if sum(abs(value) for value in shift) <= radius:
            candidate = tuple(x+d for x, d in zip(cell, shift))
            if _inside(candidate, side):
                yield candidate


def regional_bounds(member, region=((0, 0, 0),)):
    """Complete potentially nonzero crossing pairs for a finite cell region."""
    contract.verify_member(member)
    source = member['source']; side = source['side']
    if not isinstance(region, (list, tuple)) or not region:
        raise ValueError('a nonempty list of distinct source cells is required')
    cells = tuple(locality._cell(cell, side) for cell in region)
    if len(set(cells)) != len(cells):
        raise ValueError('a regional cell cannot be counted twice')
    inside = set(cells); cached = {}
    def data(cell):
        if cell not in cached: cached[cell] = _cell_data(source, cell)
        return cached[cell]
    pairs = []
    for x in sorted(cells):
        dx = data(x)
        for y in _nearby(x, side):
            if y in inside: continue
            dy = data(y)
            if set(dx['owner_support']).isdisjoint(dy['owner_support']): continue
            pairs.append({'from': x, 'to': y,
                          'current_norm_upper': 2*dx['norm_upper']*dy['norm_upper']})
    return {'region': tuple(sorted(cells)),
            'regional_energy_norm_upper': sum((data(x)['norm_upper'] for x in cells), F(0)),
            'crossing_pairs': pairs,
            'outward_flux_norm_upper': sum((row['current_norm_upper'] for row in pairs), F(0)),
            'pair_list_is_support_upper_envelope_not_nonzero_current_claim': True,
            'all_crossing_pairs_retained': True,
            'separately_project_region_to_Gauss_kernel': False}


def complete_record_budget(member, region=((0, 0, 0),), duration='1/10000',
                           half_trace_error='1/1000000', quadrature_slices=8):
    """Bound a measured continuity residual with left-endpoint flux sampling.

    The error is a uniform complete-output half-trace comparison at each
    readout time, including preparation, evolution and apparatus errors.
    """
    bounds = regional_bounds(member, region)
    time = ground._q(duration, 'source duration')
    error = ground._q(half_trace_error, 'complete-output half-trace error')
    if error > 1: raise ValueError('half-trace error must not exceed one')
    count = ground._i(quadrature_slices, 'positive quadrature count', 1)
    cap = F(member['state']['spectral_upper'])
    cap = max(cap, abs(F(member['state']['spectral_lower'])))
    energy = bounds['regional_energy_norm_upper']; flux = bounds['outward_flux_norm_upper']
    return {'duration': time, 'complete_output_half_trace_upper': error,
            'quadrature_slices': count, 'source_H_norm_upper': cap,
            'regional_energy_norm_upper': energy, 'outward_flux_norm_upper': flux,
            'two_endpoint_mean_error_upper': 4*energy*error,
            'integrated_flux_record_error_upper': 2*time*flux*error,
            'left_rule_quadrature_error_upper': time*time*cap*flux/count,
            'measured_continuity_residual_upper':
                4*energy*error+2*time*flux*error+time*time*cap*flux/count,
            'exact_integral_replaced_by_endpoint_current_without_error': False,
            'readouts_use_separate_copies_of_the_same_initial_preparation': True,
            'sequential_measurement_disturbance_ignored': False,
            'error_is_uniform_over_every_sampled_complete_output': True,
            'all_quantum_references_and_record_outcomes_retained': True}


def regional_first_law_bound(member):
    """Reuse the regional entropy theorem on this same state's actual support.

    The entropy functional concerns the source regional marginal. Its
    preparation coin and conjugate reference remain in the observer output.
    """
    from model.particle_spatial_modular_response import regional
    contract.verify_member(member)
    support = _cell_data(member['source'], (0, 0, 0))['owner_support']
    inherited = regional.common_contract_payload(member, region_cells=len(support))
    source_floor_bits = inherited['physical_density_floor_factored']['negative_exponent']
    regional_bits = (60+27*member['source']['per_link_word_degree'])*len(support)
    exponent = source_floor_bits+regional_bits
    return {
        'region_owner_cells': support,
        'inherited_selected_support_consumer': inherited,
        'regional_dimension_bits_upper': regional_bits,
        'supported_regional_density_floor': {'base': 2, 'negative_exponent': exponent,
                                            'operator': 'P_R=support Tr_out(P_G)'},
        'floor_argument': 'supported n_lambda/d_lambda >= 1/d_lambda >= 1/D_R',
        'regional_modular_norm_upper': exponent,
        'unitary': 'U_R=exp(i theta h_origin) on its actual owner support, theta=1/1000',
        'complete_preparation_instrument': 'sqrt(1-s) I, sqrt(s) U_R; both coin outcomes retained',
        'path': 'sigma_s=(1-s)sigma+s U_R sigma U_R_dagger',
        'tangent': 'X=U_R sigma U_R_dagger-sigma',
        'first_law': 'dS(sigma_s)/ds at zero = Tr(X K_R), K_R=-log(sigma)|P_R',
        'relative_entropy_remainder_upper': 's^2 ||X||_HS^2/(2m_R) <= s^2/m_R',
        'small_coin_probability_example': {'base': 2, 'negative_exponent': exponent+4},
        'small_coin_relative_entropy_remainder_upper': {'base': 2, 'negative_exponent': exponent+8},
        'entropy_functional': 'source regional marginal, after mathematical trace over the coin and reference',
        'complete_coin_and_quantum_reference_retained': True,
        'coin_plus_region_entropy': 'H2(s)+S(sigma), with a singular derivative at s=0',
        'marginal_first_law_is_full_observer_entropy': False,
        'ambient_region_full_rank_claimed': False,
        'regional_modular_operator_is_beta_times_cell_energy': False,
        'geometric_area_or_modular_diamond_identified': False,
        'full_spatial_density_or_logarithm_computed': False,
    }


def common_contract_payload(member):
    contract.verify_member(member)
    source = member['source']; a = F(source['spacing']); ell = source['per_link_word_degree']
    if source['word_projection'] != 'tensor_product_per_link_Haar_orthogonal':
        raise ValueError('this local current construction requires the product link cutoff')
    cap = (F(180)+F(1827, 2)*ell**2+F(3, 7))/a
    origin = _cell_data(source, (0, 0, 0))
    regional = regional_bounds(member)
    if origin['norm_upper'] > cap: raise ArithmeticError('cell norm exceeds the all-cell cap')
    return contract.serial({
        'program': contract.PROGRAM, 'member_id': member['member_id'], **member['ids'],
        'claim': CLAIM, 'cell_partition':
            'onsite matter + owned electric + half incident bonds + one-third owner plaquettes',
        'source_level': source['level'], 'all_sixty_CAR_modes_and_link_words_retained': True,
        'energy_zero_and_clock_unchanged': True, 'reference_idle_H': member['reference']['idle_H'],
        'local_energy_support_radius_graph': 2, 'current_zero_beyond_graph_distance': 4,
        'all_cell_energy_norm_upper': cap, 'all_pair_current_norm_upper': 2*cap*cap,
        'norm_bound_uniform_in_refinement_claimed': False,
        'origin_energy': origin, 'origin_region': regional,
        'current_definition': 'J_xy=i[h_x,h_y]',
        'operator_identity': 'd alpha_t(h_x)/dt + sum_y alpha_t(J_xy)=0',
        'stationary_selected_source_regional_net_flux_mean': F(0),
        'each_individual_current_mean_claimed_zero': False,
        'idle_reference_purification_claimed_stationary': False,
        'operator_identity_applies_after_any_physical_intervention': True,
        'complete_current_instrument': {
            'effects': '(I+J/M)/2, (I-J/M)/2 for M>=||J||; M>0',
            'record_labels': ('+M', '-M'), 'zero_current_case': 'one identity outcome labeled zero',
            'mean': 'Tr rho J', 'gauge_invariant': True, 'outcomes_discarded': False,
            'full_spatial_instrument_executed': False},
        'source_coefficient_response': {
            'family': 'H_f(t)=sum_x f_x(t) h_x',
            'fixed_state_derivative': 'partial H_f/partial f_x=h_x',
            'local_work_identity': 'd< f_x h_x >/dt + sum_y f_x f_y <J_xy>=dot(f_x)<h_x>',
            'base_member_at_all_f_equal_one': True,
            'perturbed_ground_state_reprepared': False,
            'coefficients_identified_with_observer_proper_lapse': False},
        'readout_error_example': complete_record_budget(member),
        'same_state_regional_first_law': regional_first_law_bound(member),
        'reused_Lean_theorems': [
            'OPH.QFT.CoarseEnergyWard.coarse_continuity',
            'OPH.QFT.CoarseEnergyWard.coarse_compose',
            'OPH.QFT.CoarseEnergyWard.time_increments_telescope'],
        'full_spatial_matrices_allocated': False,
        'proper_metric_or_symmetric_Lorentz_stress_identified': False,
        'regional_modular_diamond_kernel_or_area_identified': False,
        'Einstein_small_ball_premises_instantiated': False,
    })


def common_contract_bound(member):
    """Self-contained consumer binding; parent registration is separate."""
    payload = common_contract_payload(member)
    return {'schema': 'oph.primary_stress_bound.v1', 'program': contract.PROGRAM,
            'member_id': member['member_id'], **member['ids'], 'claim': CLAIM,
            'result': payload, 'result_sha256': contract.digest(payload)}


def verify_bound(member, receipt):
    if receipt != common_contract_bound(member):
        raise ValueError('the local energy consumer differs from this primary member')


def incidence_control(side=3):
    """Count each term's full allocation in an independently assembled box."""
    side = ground._i(side, 'small open box side', 2)
    if side > 4: raise ValueError('the finite incidence replay uses side at most four')
    source = {'side': side, 'spacing': F(1), 'per_link_word_degree': 1,
              'electric_coefficient': F(105), 'magnetic_coefficient': F(1, 7)}
    allocations = {}; local = {}
    for cell in product(range(side), repeat=3):
        local[cell] = _cell_data(source, cell)
        for term in local[cell]['terms']:
            allocations[term['id']] = allocations.get(term['id'], F(0))+term['cell_share']
    expected = set()
    for cell in product(range(side), repeat=3):
        expected.add(('onsite_matter', cell, ()))
        for axis in range(3):
            expected.add(('electric', cell, (axis,)))
            if cell[axis]+1 < side: expected.add(('bond', cell, (axis,)))
        for axes in combinations(range(3), 2):
            if all(cell[axis]+1 < side for axis in axes): expected.add(('plaquette', cell, axes))
    if set(allocations) != expected or any(value != 1 for value in allocations.values()):
        raise ArithmeticError('cell allocations do not recover every source term exactly')
    overlapping = 0; max_distance = 0
    for x, y in combinations(local, 2):
        if not set(local[x]['owner_support']).isdisjoint(local[y]['owner_support']):
            overlapping += 1
            max_distance = max(max_distance, sum(abs(u-v) for u, v in zip(x, y)))
    if max_distance > 4: raise ArithmeticError('support overlap exceeds current range')
    return {'side': side, 'cells': side**3, 'allocated_terms': len(expected),
            'each_source_term_allocated_once': True,
            'onsite_count': side**3, 'electric_count': 3*side**3,
            'bond_count': 3*side**2*(side-1), 'plaquette_count': 3*side*(side-1)**2,
            'maximum_plaquettes_per_owner_cell': max(sum(t['kind'] == 'plaquette' for t in d['terms']) for d in local.values()),
            'maximum_cell_energy_norm_upper': max(d['norm_upper'] for d in local.values()),
            'potentially_nonzero_unordered_currents': overlapping,
            'maximum_overlap_graph_distance': max_distance}


def _annihilation(mode, modes=4):
    entries = {}
    for mask in range(1 << modes):
        if mask & (1 << mode):
            entries[mask ^ (1 << mode), mask] = (-1)**((mask & ((1 << mode)-1)).bit_count())
    return s.ImmutableSparseMatrix(1 << modes, 1 << modes, entries)


@lru_cache(maxsize=1)
def quantum_control():
    """Separate full48-state neutral two-cell CAR/U(1) quantum-link control."""
    annihilation = [_annihilation(i) for i in range(4)]
    number = [c.H*c for c in annihilation]
    identity_f = s.eye(16); identity_l = s.eye(3)
    electric = s.diag(-1, 0, 1)
    shift = s.Matrix([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    lift = lambda op: s.kronecker_product(op, identity_l)
    n0 = lift(number[0]+number[1]); n1 = lift(number[2]+number[3])
    e = s.kronecker_product(identity_f, electric)
    forward = (s.kronecker_product(annihilation[0].H*annihilation[2], shift)/3
               +s.kronecker_product(annihilation[1].H*annihilation[3], shift.H)/5)
    bond = forward+forward.H
    h0 = n0+s.Rational(3, 2)*e*e+bond/2
    h1 = 2*n1+bond/2
    total = h0+h1
    gauss = (lift(number[0]-number[1])-e, lift(number[2]-number[3])+e)
    physical = s.diag(*(int(gauss[0][i, i] == gauss[1][i, i] == 0) for i in range(48)))
    amplitude = physical*(10*s.eye(48)-total)
    normalization = s.trace(amplitude*amplitude)
    rho = amplitude*amplitude/normalization
    seed = s.zeros(48, 1); seed[3*3+1, 0] = 1
    generic, currents = response.partition_certificate(total, (h0, h1), seed, symmetry=gauss[0])
    current = currents[0, 1]
    if any(response.commutator(op, g) != s.zeros(48) for op in (h0, h1, current) for g in gauss):
        raise ArithmeticError('the literal control term violates a Gauss constraint')
    if total*rho != rho*total or physical*rho != rho or s.trace(rho) != 1:
        raise ArithmeticError('the control polynomial is not a stationary physical density')
    phase = s.diag(*(s.I if (mask >> 0) & 1 else 1 for mask in range(16)))
    phase = s.kronecker_product(phase, identity_l)
    perturbed = phase*rho*phase.H
    cap = max(sum(abs(current[i, j]) for j in range(48)) for i in range(48))
    if cap <= 0: raise ArithmeticError('the full control current is zero')
    effects = tuple((s.eye(48)+sign*current/cap)/2 for sign in (1, -1))
    probabilities = tuple(s.simplify(s.trace(perturbed*effect)) for effect in effects)
    mean = s.simplify(s.trace(perturbed*current))
    rate = s.simplify(s.trace(perturbed*s.I*(total*h0-h0*total)))
    if sum(probabilities) != 1 or any(p < 0 for p in probabilities) or cap*(probabilities[0]-probabilities[1]) != mean:
        raise ArithmeticError('the complete current records do not recover their mean')
    if mean == 0 or rate+mean != 0:
        raise ArithmeticError('the intervention must exhibit nonzero exact boundary transport')
    return {'H': total, 'h': (h0, h1), 'J': current, 'Gauss': gauss,
            'physical_projector': physical, 'rho': rho, 'phase': phase,
            'intervened_rho': perturbed, 'effects': effects, 'record_cap': cap,
            'certificate': generic,
            'record_probabilities': probabilities, 'intervened_outward_flux': mean,
            'intervened_cell0_energy_rate': rate,
            'selected_control_mean_current': s.simplify(s.trace(rho*current)),
            'idle_purifier_energy_variance': s.simplify(s.trace(rho*total*total)-s.trace(rho*total)**2),
            'physical_dimension': s.trace(physical),
            'same_as_selected_spatial_primary': False}


def _exact(value):
    if isinstance(value, s.Rational): return str(value)
    if isinstance(value, dict): return {key: _exact(val) for key, val in value.items()}
    if isinstance(value, (list, tuple)): return [_exact(val) for val in value]
    return contract.serial(value)


def certificate():
    control = quantum_control()
    return _exact({'schema': 'oph.primary_stress.v1',
        'selected_member_bounds': [common_contract_bound(contract.member(level)) for level in (0, 1, 3)],
        'incidence_control': incidence_control(),
        'independent_quantum_control': {
            'ambient_dimension': 48, 'physical_dimension': control['physical_dimension'],
            'scope': 'two neutral CAR species per cell and one three-state U(1) link; separate control',
            'same_as_selected_spatial_primary': False,
            'finite_partition_check': control['certificate'],
            'complete_current_record_labels': (control['record_cap'], -control['record_cap']),
            'complete_current_probabilities_after_phase': control['record_probabilities'],
            'intervened_outward_flux': control['intervened_outward_flux'],
            'intervened_cell0_energy_rate': control['intervened_cell0_energy_rate'],
            'stationary_source_net_current': control['selected_control_mean_current'],
            'idle_purifier_energy_variance': control['idle_purifier_energy_variance'],
            'all_ambient_states_and_both_record_outcomes_retained': True},
        'same_source_finite_energy_conservation_constructed': True,
        'continuum_stress_or_Einstein_dynamics_proved': False})
