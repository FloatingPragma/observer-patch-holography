"""Full-complement finite loop fusion on the unchanged primary source."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, permutations

from model.particle_common_contract import contract, provenance
from model.particle_primary_reduced_sector import model as reduced


SOURCE_INPUTS = tuple(dict.fromkeys((
    *('model/particle_primary_loop_fusion/'+name for name in
      ('__init__.py', '__main__.py', 'model.py', 'test_model.py', 'PROOF.md', 'README.md')),
    *reduced.SOURCE_INPUTS,
    'model/particle_primary_reduced_sector/test_model.py',
)))


def cycle_control():
    graph = reduced.cube_support()
    edges = graph['interior_edges']
    faces = graph['oriented_outward_faces']
    hexagons = []
    for selected in combinations(range(12), 6):
        degree = Counter(v for edge in selected for v in edges[edge])
        if any(d != 2 for d in degree.values()):
            continue
        neighbors = {v: [] for v in degree}
        for i in selected:
            x, y = edges[i]
            neighbors[x].append(y)
            neighbors[y].append(x)
        visited = {min(neighbors)}
        frontier = list(visited)
        while frontier:
            for y in neighbors[frontier.pop()]:
                if y not in visited:
                    visited.add(y)
                    frontier.append(y)
        if len(visited) != 6:
            raise ArithmeticError('a six-edge support has an additional component')
        hexagons.append(selected)
    if len(hexagons) != 16 or graph['supports_with_no_degree_one_vertex_counts'][6] != len(hexagons):
        raise ArithmeticError('actual cube hexagon count changed')
    pairs = []
    for first, second in permutations(range(6), 2):
        left, right = faces[first], faces[second]
        shared = set(left['edges']) & set(right['edges'])
        flux = tuple(a+b for a, b in zip(left['outward_flux'], right['outward_flux']))
        divergence = [0]*8
        for (x, y), charge in zip(edges, flux):
            divergence[x] -= charge
            divergence[y] += charge
        if any(divergence) or max(map(abs, flux)) != 1:
            raise ArithmeticError('the sequential pair leaves actual Gauss or link cutoff')
        active = tuple(i for i, q in enumerate(flux) if q)
        kind = 'adjacent_fused_hexagon' if shared else 'opposite_disjoint_squares'
        if shared:
            if len(shared) != 1 or any(flux[i] for i in shared) or active not in hexagons:
                raise ArithmeticError('the common edge did not cancel exactly')
        elif len(active) != 8:
            raise ArithmeticError('opposite source faces do not give disjoint loops')
        pairs.append({'first_face': first, 'second_face': second, 'kind': kind,
                      'resulting_q6_flux': flux, 'nonconstant_links': len(active),
                      'electric_energy': F(63, 4)*len(active)})
    counts = Counter(row['kind'] for row in pairs)
    if counts != {'adjacent_fused_hexagon': 24, 'opposite_disjoint_squares': 6}:
        raise ArithmeticError('complete ordered face-pair menu changed')
    return {'oriented_native_charge': 6, 'all_ordered_distinct_face_pairs': pairs,
            'adjacent_ordered_pairs': 24, 'opposite_ordered_pairs': 6,
            'six_edge_simple_cycles': len(hexagons),
            'q6_hexagon_electric_eigenspace_dimension': 2*len(hexagons),
            'weak_square_electric_eigenspace_dimension': 12,
            'combined_electric_cluster_dimension': 44,
            'each_intermediate_and_final_word_admitted_at_degree_one': True,
            'pair_intervention_uses_unprojected_parent_dynamics': False}


def fusion_bounds():
    eps = F(6, 7)
    center = F(189, 2)
    next_floor = F(441, 4)
    previous = F(63)
    outside_distance = min(center-previous-eps, next_floor-center-eps)
    electric_residual_outside = (eps/outside_distance)**2
    ground_distance = F(10, 7)/(27*62)
    if ground_distance >= F(1, 1000):
        raise ArithmeticError('gauge ground comparison exceeds its exact angle bound')
    if 1-electric_residual_outside <= F(998, 1000)**2:
        raise ArithmeticError('full-complement fused window bound failed')
    fused_weight = (F(998, 1000)-F(1, 1000))**2
    disjoint_weight = (F(9, 10)-F(1, 1000))**2
    if F(40, 49) <= F(9, 10)**2 or fused_weight <= F(99, 100) or disjoint_weight <= F(4, 5):
        raise ArithmeticError('the original-ground probe comparison failed')
    fundamental = reduced.spectral_bounds()
    e0 = fundamental['interacting_gauge_ground_energy_interval']
    e1 = fundamental['interacting_A1u_energy_interval']
    cluster = [center-eps, center+eps]
    deficit = [2*e1[0]-cluster[1]-e0[1], 2*e1[1]-cluster[0]-e0[0]]
    return {'same_source_interaction_coefficient': 1,
            'full_magnetic_norm_upper': eps,
            'electric_fusion_cluster_energy': center,
            'electric_next_level_lower': next_floor,
            'interacting_fusion_cluster_rank': 44,
            'interacting_fusion_cluster_interval': cluster,
            'distance_of_electric_anchor_energy_to_full_cluster_complement_lower': outside_distance,
            'fused_anchor_weight_outside_cluster_upper': electric_residual_outside,
            'actual_ground_to_fused_cluster_weight_strict_lower': F(99, 100),
            'actual_ground_to_fused_cluster_weight_comparison_lower': fused_weight,
            'actual_ground_to_disjoint_window_weight_strict_lower': F(4, 5),
            'actual_ground_to_disjoint_window_weight_comparison_lower': disjoint_weight,
            'disjoint_final_energy_window': [F(124), F(128)],
            'fusion_transfer_energy_interval': [cluster[0]-e0[1], cluster[1]-e0[0]],
            'deficit_relative_to_twice_isolated_A1u_transition': deficit,
            'deficit_decimal_display_only': [format(float(x), '.12f') for x in deficit],
            'small_subspace_assumed_invariant': False,
            'bound_state_below_a_two_particle_scattering_threshold_proved': False,
            'moving_two_packet_collision_executed': False}


def bind_member(member):
    base = reduced.bind_member(member)['result']
    graph = cycle_control()
    bounds = fusion_bounds()
    clock = base['complete_clocked_original_event']
    delta = F(clock['one_clock_loss_upper'])
    eta = F(clock['complete_clock_plan']['energy_resolution'])
    floor = base['original_energy_local_population_floor']
    categories = []
    for kind, interval, ideal, simple in (
            ('adjacent_fused_hexagon', bounds['interacting_fusion_cluster_interval'], F(99, 100), F(49, 50)),
            ('opposite_disjoint_squares', bounds['disjoint_final_energy_window'], F(4, 5), F(3, 4))):
        factor = ideal*(1-delta)**2
        if factor <= simple:
            raise ArithmeticError('complete two-clock loss exceeds the original weight budget')
        categories.append({'kind': kind, 'ideal_final_energy_window': interval,
                           'reported_final_energy_window': [interval[0]-eta, interval[1]+eta],
                           'original_weight_lower_factored': {'base': 2, 'negative_exponent': floor['negative_exponent'],
                                                             'rational_multiplier': factor},
                           'simple_original_weight_strict_lower_multiplier': simple,
                           'complete_thirty_pair_menu_each_row_multiplier': factor/30})
    result = contract.serial({
        'member_id': member['member_id'], **member['ids'],
        'same_source_state_reference_clock': True,
        'original_selected_density_retained': True,
        'reducing_sector': base['reducing_sector'],
        'source_loop_geometry': graph, 'full_complement_fusion_bounds': bounds,
        'original_gauge_ground_population_floor': floor,
        'complete_clock_plan': clock['complete_clock_plan'],
        'initial_record_window': clock['initial_record_window'],
        'one_clock_loss_upper': delta,
        'three_sequential_probe_outcomes': ['C_g C_f P_empty', 'sqrt(I-C_g* C_g) C_f P_empty', 'sqrt(I-P_empty C_f* C_f P_empty)'],
        'compressed_single_face_contraction': 'C_f=P_word chi6(outward_face_f) P_word on physical Gauss support',
        'complete_ordered_pair_selector_probability': F(1, 30),
        'complete_pair_and_probe_outcomes': 90,
        'complete_original_categories': categories,
        'each_complete_pair_row_original_weight_strict_lower_factored': {'base': 2, 'negative_exponent': floor['negative_exponent'], 'rational_multiplier': F(1, 40)},
        'reference_outputs_and_all_failures_retained': True,
        'initial_extreme_sector_weight_conditioned_out': False,
        'retained_pure_reference_Loewner_bound_asserted': False,
        'actual_static_loop_fusion_energy_separation_proved': True,
        'identified_physical_particle_binding_or_mass_proved': False,
        'relativistic_dispersion_or_conserved_spatial_momentum_proved': False,
        'extreme_sector_supported_weight_vanishes_on_selected_refinement': True,
        'physical_mass_table_GeV': [],
    })
    return {'schema': 'oph.primary_loop_fusion.v1', 'result': result,
            'result_sha256': contract.digest(result)}


def verify_bound(member, receipt):
    if receipt != bind_member(member):
        raise ValueError('loop fusion source, spectrum, original weight or apparatus changed')


def certificate():
    paths = provenance.definition_inputs(SOURCE_INPUTS)
    pins = provenance.definition_hashes(paths)
    value = {'schema': 'oph.primary_loop_fusion_certificate.v1',
             'selected_member_bound': bind_member(contract.member(0)), 'source_sha256': pins}
    if provenance.definition_hashes(paths) != pins:
        raise ValueError('loop fusion source changed during replay')
    return value
