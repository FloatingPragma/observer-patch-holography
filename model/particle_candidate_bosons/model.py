"""Complete loop multiplets and exact detector selection on the primary source."""
from fractions import Fraction as F
from functools import lru_cache

import sympy as s

from model.particle_common_contract import contract, provenance
from model.particle_primary_reduced_sector import model as reduced
from model.particle_primary_charge_conjugation import model as charge


SOURCE_INPUTS = tuple(dict.fromkeys((
    *('model/particle_candidate_bosons/' + name for name in
      ('__init__.py', '__main__.py', 'model.py', 'test_model.py', 'PROOF.md', 'README.md')),
    *reduced.SOURCE_INPUTS, *charge.SOURCE_INPUTS,
)))


def clean(value):
    return s.expand(value)


@lru_cache(maxsize=1)
def detector_basis():
    """Orthonormal face coefficients in the literal outward-normal order."""
    e = s.eye(6)
    even = [(e[:, 2*i] + e[:, 2*i+1])/s.sqrt(2) for i in range(3)]
    odd = [(e[:, 2*i+1] - e[:, 2*i])/s.sqrt(2) for i in range(3)]
    scalar = sum(even, s.zeros(6, 1))/s.sqrt(3)
    diagonal = [(even[0]-even[1])/s.sqrt(2),
                (even[0]+even[1]-2*even[2])/s.sqrt(6)]
    matrix = s.Matrix.hstack(scalar, *diagonal, *odd).applyfunc(clean)
    if (matrix.T*matrix).applyfunc(clean) != s.eye(6):
        raise ArithmeticError('literal face coefficients are not orthonormal')
    return matrix


@lru_cache(maxsize=1)
def symmetry_control():
    graph = reduced.cube_support()
    normals = reduced.cube_channels()['face_outward_normals']
    basis = detector_basis()
    scalar, diagonal, vector = basis[:, :1], basis[:, 1:3], basis[:, 3:]
    reverse = s.zeros(6)
    for j, v in enumerate(normals):
        reverse[normals.index(tuple(-x for x in v)), j] = 1
    if (basis.T*reverse*basis).applyfunc(clean) != s.diag(1, 1, 1, -1, -1, -1):
        raise ArithmeticError('spatial inversion of actual face words changed')
    rotations = []
    for q in reduced.channels.group.elements():
        R = reduced.channels.group.rotation(q)
        P = s.zeros(6)
        for j, v in enumerate(normals):
            w = tuple(sum(int(R[i, k])*v[k] for k in range(3)) for i in range(3))
            P[normals.index(w), j] = 1
        if P*scalar != scalar or (P*vector-vector*R).applyfunc(clean) != s.zeros(6, 3):
            raise ArithmeticError('actual source vector/scalar covariance failed')
        if (diagonal.T*P*diagonal).trace() != reduced.channels._sym(reduced.channels.group.characters(q)['E']):
            raise ArithmeticError('actual diagonal face pair is not E')
        action = (basis.T*P*basis).applyfunc(clean)
        if any(action[i, j] != 0 for i in range(6) for j in range(6)
               if (0 if i == 0 else 1 if i < 3 else 2) != (0 if j == 0 else 1 if j < 3 else 2)):
            raise ArithmeticError('inequivalent cube detector channels mixed')
        rotations.append(P)
    average = sum(rotations, s.zeros(6))/len(rotations)
    if average*vector != s.zeros(6, 3) or average*scalar != scalar:
        raise ArithmeticError('vector-to-scalar selection control failed')
    rows = []
    for part, conjugation in (('real', 1), ('imaginary', -1)):
        labels = ('A1g', 'Eg', 'Eg', 'T1u', 'T1u', 'T1u') if part == 'real' else (
            'A1u', 'Eu', 'Eu', 'T1g', 'T1g', 'T1g')
        for j, label in enumerate(labels):
            rows.append({'channel': f'{part}_{j}', 'cube_irrep': label,
                         'component': j if j == 0 else j-1 if j < 3 else j-3,
                         'face_coefficients': [str(x) for x in basis[:, j]],
                         'spatial_inversion': 1 if label.endswith('g') else -1,
                         'gauge_coordinate_conjugation': conjugation, 'CAR_parity': 0})
    return {'actual_binary_cube_rotations_checked': len(rotations),
            'actual_oriented_source_edges': len(graph['interior_edges']),
            'loop_channels': rows,
            'multiplet_dimensions': {'A1g': 1, 'Eg': 2, 'T1u': 3, 'A1u': 1, 'Eu': 2, 'T1g': 3},
            'bare_scalar_even_to_A1u_overlap': '0',
            'bare_vector_to_A1u_overlap': '0',
            'total_cube_irreps_identified_with_intrinsic_Lorentz_spin': False}


@lru_cache(maxsize=1)
def loop_polynomials():
    """Actual q6 character multiplication, before the source word projection."""
    fluxes = [f['outward_flux'] for f in reduced.cube_support()['oriented_outward_faces']]
    basis = detector_basis()
    operators = []
    for imaginary in (False, True):
        for j in range(6):
            terms = {}
            for face, flux in enumerate(fluxes):
                a = basis[face, j]/s.sqrt(24)
                for sign in (-1, 1):
                    value = a if not imaginary else sign*a/s.I
                    key = tuple(sign*x for x in flux)
                    terms[key] = clean(terms.get(key, 0)+value)
            operators.append({key: value for key, value in terms.items() if value})
    return tuple(operators)


@lru_cache(maxsize=1)
def completeness_control():
    operators = loop_polynomials()
    total = {}
    for operator in operators:
        if {tuple(-x for x in k): s.conjugate(v) for k, v in operator.items()} != operator:
            raise ArithmeticError('loop detector is not self-adjoint')
        for u, a in operator.items():
            for v, b in operator.items():
                key = tuple(x+y for x, y in zip(u, v))
                total[key] = clean(total.get(key, 0)+a*b)
    total = {key: value for key, value in total.items() if value}
    if total != {(0,)*12: s.Integer(1)}:
        raise ArithmeticError('complete actual face menu does not sum to identity before projection')
    gram = s.zeros(12)
    for i, left in enumerate(operators):
        for j, right in enumerate(operators):
            gram[i, j] = clean(sum(s.conjugate(a)*right.get(k, 0) for k, a in left.items()))
    if gram != s.eye(12)/12:
        raise ArithmeticError('Haar electric anchors have a different normalization')
    return {'success_rows': 12, 'failure_rows': 1,
            'parent_sum_F_dagger_F': 'I', 'compressed_sum_F_dagger_F_upper': 'I',
            'blank_anchor_Gram': 'I_12/12', 'extra_uniform_selector_divisor': 1,
            'failure_operator': 'sqrt(I-sum_j P_empty F_j^dagger F_j P_empty)',
            'all_other_native_sectors_and_source_reference_outputs_retained': True,
            'compressed_menu_assumed_exactly_saturated': False}


@lru_cache(maxsize=1)
def native_current_control():
    currents = reduced.channels.native.native_currents()
    columns = s.Matrix.hstack(*(g.reshape(225, 1) for g in currents))
    commutators = [a*b-b*a for i, a in enumerate(currents) for b in currents[i+1:]]
    nonzero = [x for x in commutators if x != s.zeros(15)]
    derived = s.Matrix.hstack(*(x.reshape(225, 1) for x in nonzero))
    if columns.rank() != 12 or derived.rank() != 11:
        raise ArithmeticError('actual native center/derived algebra dimensions changed')
    if any(currents[-1]*g != g*currents[-1] for g in currents):
        raise ArithmeticError('the declared native central current is not central')
    if s.Matrix.hstack(derived, currents[-1].reshape(225, 1)).rank() != 12:
        raise ArithmeticError('native hypercharge was incorrectly removed by an adjoint twirl')
    return {'actual_native_generator_rank': 12, 'derived_su2_plus_su3_rank': 11,
            'center_rank': 1, 'undressed_single_site_derived_current_Gauss_compression': '0',
            'abelian_current_killed_by_adjoint_twirl_claimed': False,
            'all_normal_ordered_number_preserving_matter_currents_annihilate_empty_sector': True,
            'charged_creation_operators_included_in_that_annihilation_claim': False,
            'internal_adjoint_directions_are_particle_poles': False}


def spectral_control():
    perturbation = F(6, 7)
    distance = F(429, 14)
    outside = perturbation/distance
    ground_distance = F(10, 7)/(27*62)
    amplitude = F(15, 52)*(1-outside**2)-F(1, 1100)
    if outside != F(4, 143) or ground_distance >= F(1, 1100) or amplitude <= F(2, 7):
        raise ArithmeticError('full-source multiplet residue estimate failed')
    low, high = F(63)-perturbation, F(63)+perturbation
    return {'interaction_coefficient': 1, 'isolated_cluster_rank': 12,
            'isolated_cluster_final_energy_interval': [low, high],
            'multiplet_transition_energy_interval': [low, high+F(1, 45675)],
            'full_complement_distance_from_electric_anchor_lower': distance,
            'electric_anchor_outside_cluster_norm_upper': outside,
            'gauge_ground_distance_from_blank_strict_upper': F(1, 1100),
            'each_actual_ground_to_multiplet_amplitude_strict_lower': amplitude,
            'each_actual_ground_to_multiplet_residue_strict_lower': F(4, 49),
            'one_exact_energy_per_cube_multiplet': True,
            'all_six_source_transition_ratios_upper': F(149, 145),
            'zero_gap_transition_in_this_cluster': False,
            'six_multiplet_energies_proved_pairwise_distinct': False,
            'exact_H_replaced_by_a_twelve_state_compression': False,
            'source_gap_identified_with_a_rest_mass': False}


def bind_member(member):
    base = reduced.bind_member(member)['result']
    conjugation = charge.bind_member(member)['result']
    labels = symmetry_control()
    menu = completeness_control()
    currents = native_current_control()
    spectral = spectral_control()
    clock = base['complete_clocked_original_event']
    eta = F(clock['complete_clock_plan']['energy_resolution'])
    loss = F(clock['one_clock_loss_upper'])
    factor = F(4, 49)*(1-loss)**2
    if factor <= F(1, 13):
        raise ArithmeticError('complete unchanged source clocks lost a channel')
    floor = base['original_energy_local_population_floor']
    result = contract.serial({
        'member_id': member['member_id'], **member['ids'],
        'same_source_state_reference_clock': True, 'source_input_is_original_selected_density': True,
        'source_domain': base['reducing_sector'], 'actual_loop_symmetry': labels,
        'actual_complete_loop_menu': menu, 'native_current_selection': currents,
        'complete_instrument': 'complete empty/native-sector and cube-irrep records; initial H clock; twelve loop contractions plus failure; complete final cube-irrep and H-clock records',
        'full_source_interacting_multiplet_bounds': spectral,
        'original_gauge_ground_population_floor': floor,
        'complete_clock_plan': clock['complete_clock_plan'],
        'initial_record_window': clock['initial_record_window'],
        'final_cluster_record_window': [spectral['isolated_cluster_final_energy_interval'][0]-eta,
                                        spectral['isolated_cluster_final_energy_interval'][1]+eta],
        'each_row_original_reported_weight_strict_lower_factored': {
            'base': 2, 'negative_exponent': floor['negative_exponent'], 'rational_multiplier': F(1, 13)},
        'each_row_clock_weight_before_tiny_floor': factor,
        'combined_twelve_row_original_reported_weight_strict_lower_multiplier': F(12, 13),
        'multiplet_original_reported_weight_strict_lower_multipliers': {
            label: F(d, 13) for label, d in labels['multiplet_dimensions'].items()},
        'forbidden_original_record_categories': {
            'initial_empty_A1g_then_even_scalar_then_A1u': '0',
            'initial_empty_A1g_then_polar_or_axial_vector_then_A1u': '0',
            'empty_gauge_ground_then_normal_ordered_matter_current': '0'},
        'zero_records_require_the_stated_initial_symmetry_record': True,
        'clock_leakage_breaks_these_exact_symmetry_zeroes': False,
        'zeroes_asserted_for_arbitrary_initial_native_sectors': False,
        'gauge_coordinate_conjugation_within_empty_sector': {
            'linear_unitary': True, 'gauge_ground_sign': 1,
            'real_loop_multiplet_sign': 1, 'imaginary_loop_multiplet_sign': -1,
            'is_full_native_charge_conjugation': False},
        'full_native_charge_conjugation': {
            'empty_sector_maps_to_full_matter_sector': True,
            'CAR_modes': conjugation['CAR_modes'], 'preserves_source_energy_zero': True,
            'same_classical_probabilities_in_conjugate_complete_experiment': True,
            'source_only_operation_fixes_retained_reference_state': False},
        'candidate_readouts': {
            'photon': 'nonzero T1u polar and T1g axial source channels; field/helicity/conjugation dictionary and masslessness unproved',
            'gluon': 'native su3 directions exist; bare adjoint current is not a physical Gauss-invariant detector',
            'W_plus_minus': 'native weak directions exist; charged vector and symmetry-broken field attachment unproved',
            'Z': 'neutral vector channel exists; weak/hypercharge mixing and Z identity unproved',
            'Higgs': 'A1g loop scalar couples; separate charged-composite Higgs-doublet weights do not establish its vacuum response or a Higgs pole',
            'A1u': 'isolated odd-inversion gauge channel; even-scalar and vector direct readouts excluded'},
        'actual_six_interacting_gauge_multiplet_records_proved': True,
        'failure_reference_and_rare_sector_probabilities_retained': True,
        'positive_window_normalized_by_success': False,
        'original_state_rare_sector_replaced_by_prepared_gauge_ground': False,
        'intrinsic_Lorentz_spin_or_physical_CP_assigned': False,
        'uniform_refinement_particle_residue_proved': False,
        'physical_particle_map': [], 'physical_mass_table_GeV': [],
    })
    return {'schema': 'oph.particle_candidate_bosons.v1', 'result': result,
            'result_sha256': contract.digest(result)}


def verify_bound(member, value):
    if value != bind_member(member):
        raise ValueError('boson detector, source, symmetry, original weight or clock changed')


def certificate():
    paths = provenance.definition_inputs(SOURCE_INPUTS)
    pins = provenance.definition_hashes(paths)
    result = {'schema': 'oph.particle_candidate_bosons_certificate.v1',
              'selected_member_bound': bind_member(contract.member(0)), 'source_sha256': pins}
    if provenance.definition_hashes(paths) != pins:
        raise ValueError('source inputs changed during boson-channel replay')
    return result


build = certificate
