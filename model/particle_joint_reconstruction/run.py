"""Replay the unified field/spectrum dependencies on their actual domains."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from model.particle_common_contract import contract, provenance
from model.particle_joint_reconstruction import candidate, positive_energy
from model.particle_primary_field_algebra import algebra
from model.particle_primary_stress import stress
from model.particle_positive_window import window
from model.particle_filtered_probe import filter as filtered_probe
from model.particle_primary_survival import survival
from model.particle_primary_mass_readout import readout as mass_readout
from model.particle_primary_spatial_algebra import spatial
from model.particle_primary_capacity import capacity
from model.particle_compact_mass_probe import compact
from model.particle_primary_momentum import momentum
from model.particle_primary_channel_labels import channels
from model.particle_primary_spectral_scale import scale
from model.particle_primary_endpoint_bounds import endpoint
from model.particle_primary_variational_readout import variational
from model.particle_primary_joint_spectroscopy import spectroscopy
from model.particle_primary_occupation_bound import model as occupation
from model.particle_primary_residual_refinement import model as local_residual
from model.particle_primary_response_moments import model as pair_response
from model.particle_primary_pair_channels import model as pair_channels
from model.particle_primary_refinement_response import model as refinement_response
from model.particle_primary_gravity_premises import model as gravity_premises
from model.particle_primary_branch_resolution import model as branch_resolution
from model.particle_primary_origin_momentum import model as origin_momentum
from model.particle_primary_energy_tightness import model as energy_tightness
from model.particle_primary_open_coframe import model as open_coframe
from model.particle_primary_charged_string import model as charged_string
from model.particle_primary_odd_residual import model as odd_residual
from model.particle_primary_momentum_transfer import model as momentum_transfer
from model.particle_primary_local_fluctuation import model as local_fluctuation
from model.particle_primary_spectral_resolution import model as spectral_resolution
from model.particle_primary_charged_readout import model as charged_readout
from model.particle_primary_staggered_response import model as staggered_response
from model.particle_primary_current_tail import model as current_tail
from model.particle_primary_sector_threshold import model as sector_threshold
from model.particle_primary_populated_momentum import model as populated_momentum
from model.particle_primary_weighted_survival import model as weighted_survival
from model.particle_primary_charge_conjugation import model as charge_conjugation
from model.particle_primary_composite_transport import model as composite_transport
from model.particle_primary_composite_spin import model as composite_spin
from model.particle_primary_selected_odd_response import model as selected_odd_response
from model.particle_primary_covariant_stress import model as covariant_stress
from model.particle_primary_selected_spin_response import model as selected_spin
from model.particle_primary_odd_resolvent import model as odd_resolvent
from model.particle_primary_odd_low_energy import model as odd_filter
from model.particle_primary_geometric_cut import model as geometric_cut
from model.particle_primary_coherent_odd_response import model as coherent_odd
from model.particle_primary_original_shell_moments import model as original_shells
from model.particle_primary_original_branch import model as original_branch
from model.particle_primary_critical_refinement import model as critical_refinement
from model.particle_primary_boundary_information import model as boundary_information
from model.particle_primary_branch_certificate import model as branch_certificate
from model.particle_primary_odd_sum_rules import model as odd_sum_rules
from model.particle_primary_work_spectrum import model as work_spectrum
from model.particle_primary_boundary_correlation import model as boundary_correlation
from model.particle_primary_reduced_sector import model as reduced_sector
from model.particle_primary_regional_fermion import model as regional_fermion
from model.particle_primary_weighted_spectral_tail import model as weighted_spectral_tail
from model.particle_primary_energy_response import model as energy_response
from model.particle_primary_loop_fusion import model as loop_fusion
from model.particle_primary_equal_work_control import model as equal_work_control
from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE / 'results.json'
NODE_IDS = (
    'fixed_program', 'physical_preparation', 'local_records', 'joint_refinement',
    'relativistic_field', 'chiral_species_attachment', 'particle_observables',
    'independent_comparison', 'pure_gauge_gap', 'geometric_gravity',
)
FINITE_CONSUMERS = {
    'fixed_program': ('source_bounds', 'common_quadratic_clock'),
    'physical_preparation': ('physical_preparation', 'physical_regional_support',
                             'physical_boundary_preparation'),
    'local_records': ('complete_local_records',),
}
FINITE_DOMAINS = {
    'fixed_program': 'The selected open Wilson family, with sixty CAR modes per cell, product link-word cutoff and the original energy zero.',
    'physical_preparation': 'The exact Gauss-polynomial state of each selected finite member and its retained canonical conjugate reference.',
    'local_records': 'The actual selected Wilson action and prepared state, the complete origin-cell boundary-charge instrument, and a remote physical plaquette phase on the same input and clock.',
}
REQUIRED_EDGES = {
    node: ([] if index == 0 else ['relativistic_field'] if node in ('pure_gauge_gap', 'geometric_gravity')
           else [NODE_IDS[index-1]]) for index, node in enumerate(NODE_IDS)
}
COMPONENT_REQUIREMENTS = {
    'field_common_theory_ten': {str(n) for n in range(1, 11)},
    'native_full_spectrum_ten': {str(n) for n in range(1, 11)},
    'defect_identification_ten': {str(n) for n in range(1, 11)},
    'yang_mills_Y1_Y6': {f'Y{n}' for n in range(1, 7)},
}


def validate_dependencies(graph):
    if graph.get('schema') != 'oph.joint_reconstruction_dependencies.v1':
        raise ValueError('unknown dependency schema')
    if graph.get('primary_program') != contract.PROGRAM:
        raise ValueError('dependency graph belongs to another program')
    rows = graph['nodes']
    if [row['id'] for row in rows] != list(NODE_IDS):
        raise ValueError('the unified graph must contain its ten distinct nodes')
    if [row['number'] for row in rows] != list(range(1, 11)):
        raise ValueError('dependency numbering differs')
    known = set(NODE_IDS)
    assumptions = {row['id'] for row in candidate.assumption_ledger()}
    covered = set()
    completed = set()
    for row in rows:
        dependencies = row['depends_on']
        if len(dependencies) != len(set(dependencies)) or not set(dependencies) <= completed:
            raise ValueError('dependency graph is cyclic or not in acceptance order')
        if dependencies != REQUIRED_EDGES[row['id']]:
            raise ValueError('required acceptance dependency was removed or changed')
        completed.add(row['id'])
        expected = ('finite_constructed' if row['id'] in FINITE_CONSUMERS else
                    'conditional_criterion' if row['id'] == 'pure_gauge_gap' else
                    'physical_proof_required')
        if row['status'] != expected:
            raise ValueError('a physical premise cannot be promoted by editing its status')
        if row['id'] in FINITE_DOMAINS and row['domain'] != FINITE_DOMAINS[row['id']]:
            raise ValueError('finite consumer domain differs from the actual experiment')
        if not row['domain'] or not row['acceptance'] or not row['evidence']:
            raise ValueError('every requirement needs a domain, acceptance conditions and evidence')
        for item in row['evidence']:
            path = (ROOT / item['path']).resolve()
            if not path.is_relative_to(ROOT) or not path.is_file() or not item['scope']:
                raise ValueError('evidence needs an existing repository file and its scope')
        listed = set(row['conditional_assumptions'])
        if not listed <= assumptions:
            raise ValueError('unregistered conditional assumption')
        covered |= listed
    if covered != assumptions:
        raise ValueError('conditional assumptions are omitted from the unified graph')
    extension = graph['conditional_candidate']
    if (extension['id'] != candidate.EXTENSION_ID
            or extension['discharges_primary_physical_dependencies']
            or extension['calibrated_tree_coefficients_are_independent_predictions']
            or extension['full_theory_positive_global_gap_claimed']):
        raise ValueError('conditional extension cannot discharge primary physics or a global photon gap')
    pure = rows[8]
    if (pure['required_for_full_particle_spectrum']
            or pure['implies_full_theory_positive_global_gap']
            or pure['required_shared_nodes'] != list(NODE_IDS[:5])):
        raise ValueError('pure Yang-Mills must instantiate its own source through reconstruction')
    gravity = rows[9]
    if gravity['required_for_full_particle_spectrum'] or gravity['proper_geometry_or_Einstein_limit_proved']:
        raise ValueError('finite stress cannot be promoted to a physical Einstein limit')
    maps = graph['requirement_maps']
    if set(maps) != set(COMPONENT_REQUIREMENTS):
        raise ValueError('component requirement coverage is incomplete')
    for name, mapping in maps.items():
        if set(mapping) != COMPONENT_REQUIREMENTS[name]:
            raise ValueError('component requirement coverage is incomplete')
        for targets in mapping.values():
            if not targets or any(target not in known and not target.startswith('outside_target:')
                                  for target in targets):
                raise ValueError('a component requirement has no unified destination')


def build():
    stored_candidate = strict_json.load(candidate.OUTPUT)
    graph = strict_json.load(HERE / 'dependencies.json')
    validate_dependencies(graph)
    paths = (set(stored_candidate['source_sha256']) | set(contract.SOURCE_INPUTS)
             | set(positive_energy.SOURCE_INPUTS))
    paths.update(item['path'] for node in graph['nodes'] for item in node['evidence'])
    for module in (algebra, stress, window, filtered_probe, survival,
                   mass_readout, spatial, capacity, compact, momentum, channels, scale,
                   endpoint, variational, spectroscopy, occupation, local_residual,
                   pair_response, pair_channels, refinement_response, gravity_premises,
                   branch_resolution, origin_momentum, energy_tightness, open_coframe,
                   charged_string, odd_residual, momentum_transfer, local_fluctuation,
                   spectral_resolution, charged_readout, staggered_response, current_tail,
                   sector_threshold, populated_momentum, weighted_survival, charge_conjugation,
                   composite_transport, composite_spin, selected_odd_response, covariant_stress,
                   selected_spin, odd_resolvent, odd_filter, geometric_cut, coherent_odd,
                   original_shells, original_branch, critical_refinement, boundary_information,
                   branch_certificate, odd_sum_rules, work_spectrum, boundary_correlation,
                   reduced_sector, regional_fermion, weighted_spectral_tail, energy_response,
                   loop_fusion, equal_work_control):
        paths.update(module.SOURCE_INPUTS)
        paths.update(str(path.relative_to(ROOT)) for path in Path(module.__file__).parent.iterdir()
                     if path.is_file() and path.suffix in ('.py', '.md'))
    paths.update(str(path.relative_to(ROOT)) for path in HERE.iterdir()
                 if path.suffix in ('.py', '.md') or path.name in
                 ('dependencies.json', 'candidate_results.json'))
    pins = provenance.definition_hashes(provenance.definition_inputs(paths))
    common = contract.certificate()
    member = common['fixed_primary_member']
    conditional = candidate.build_candidate()
    if conditional['primary_member'] != member:
        raise ValueError('candidate ancestry differs from the selected finite member')
    # Replaying the candidate also verifies its frozen calibration ancestry.
    if candidate.OUTPUT.read_text() != candidate.canonical(conditional):
        raise ValueError('candidate receipt differs; regenerate the actual candidate first')
    admissions = []
    for row in graph['nodes']:
        consumers = []
        for name in FINITE_CONSUMERS.get(row['id'], ()):
            receipt = common['bound_consequences'][name]
            contract.verify_result(member, receipt)
            consumers.append({'name': name, 'claim': receipt['claim'],
                              'result_sha256': receipt['result_sha256']})
        admissions.append({
            'id': row['id'], 'number': row['number'], 'status': row['status'],
            'depends_on': row['depends_on'], 'domain': row['domain'],
            'selected_primary_finite_acceptance': bool(consumers),
            'verified_consumers': consumers,
            'transfers_to_chiral_extension_or_pure_gauge': False,
            'physical_continuum_acceptance': False,
        })
    positive = positive_energy.certificate()
    positive_energy.verify_bound(member, positive['selected_member_bounds'][0])
    field_algebra = algebra.certificate()
    for level, receipt in zip((1, 6), field_algebra['selected_members'], strict=True):
        algebra.verify_bound(contract.member(level), receipt)
    source_stress = stress.certificate()
    for level, receipt in zip((0, 1, 3), source_stress['selected_member_bounds'], strict=True):
        stress.verify_bound(contract.member(level), receipt)
    filtered = filtered_probe.certificate()
    filtered_probe.verify_bound(member, filtered['selected_member_plan'],
                                **filtered_probe.precision_parameters(member))
    menu_member = contract.member(1)
    menu_parameters = filtered_probe.precision_parameters(menu_member, half_width='1/100000')
    for kind, receipt in zip(('charge', 'plaquette'), filtered['authenticated_field_menu_plans'], strict=True):
        filtered_probe.verify_bound(menu_member, receipt,
            probe_descriptor=algebra.probe(menu_member, kind=kind), **menu_parameters)
    spectral_window = window.certificate()
    if spectral_window['selected_requirements_without_a_positive_contrast_witness'] != window.selected_requirement(contract.member(2), 0, 0):
        raise ValueError('the spectral diagnostic substituted its primary source or contrast')
    finite_survival = survival.certificate()
    for level, receipt in zip((1, 2), finite_survival['finite_primary_witnesses'], strict=True):
        survival.verify_bound(contract.member(level), receipt)
    spatial_algebra = spatial.certificate()
    for level, receipt in zip((2, 4, 6), spatial_algebra['selected_members'], strict=True):
        spatial.verify_bound(contract.member(level), receipt)
    mass = mass_readout.certificate()
    for level, receipt in zip((0, 2, 4), mass['selected_members'], strict=True):
        mass_readout.verify_bound(contract.member(level), receipt)
    capacities = capacity.certificate()
    compact_probe = compact.certificate()
    for level, receipt in zip((1, 2), compact_probe['selected_members'], strict=True):
        compact.verify_bound(contract.member(level), receipt)
    spatial_readout = momentum.certificate()
    momentum.verify_bound(contract.member(2), spatial_readout['selected_aperture'])
    channel_labels = channels.certificate()
    for level, receipt in zip((1, 2), channel_labels['finite_members'], strict=True):
        channels.verify_bound(contract.member(level), receipt)
    spectral_scale = scale.certificate()
    endpoint_bounds = endpoint.certificate()
    for level, receipt in zip((0, 1), endpoint_bounds['source_ground_bounds'], strict=True):
        endpoint.verify_bound(contract.member(level), receipt)
    variational_bounds = variational.certificate()
    for level, receipt in zip((0, 1, 2), variational_bounds['members'], strict=True):
        variational.verify_bound(contract.member(level), receipt)
    joint_spectroscopy = spectroscopy.certificate()
    spectroscopy.verify_bound(contract.member(2), joint_spectroscopy['actual_source_instrument'])
    occupation_bounds = occupation.certificate()
    occupation.verify_bound(contract.member(0), occupation_bounds['selected_member_bound'])
    exact_residual = local_residual.certificate()
    for level, receipt in zip((0, 1, 2), exact_residual['members'], strict=True):
        local_residual.verify_bound(contract.member(level), receipt)
    neutral_response = pair_response.certificate()
    pair_response.verify_bound(contract.member(0), neutral_response['selected_member_bound'])
    neutral_channels = pair_channels.certificate()
    pair_channels.verify_bound(contract.member(0), neutral_channels['selected_member'])
    scalable_response = refinement_response.certificate()
    for level, receipt in zip((0, 1, 2), scalable_response['members'], strict=True):
        refinement_response.verify_bound(contract.member(level), receipt)
    modular_work = gravity_premises.certificate()
    gravity_premises.verify_bound(contract.member(0), modular_work['selected_member_bound'])
    concentrated_response = branch_resolution.certificate()
    branch_resolution.verify_bound(member, concentrated_response['selected_member_bound'])
    relative_momentum = origin_momentum.certificate()
    origin_momentum.verify_bound(member, relative_momentum['selected_member_bound'])
    for level, receipt in zip((6, 13, 16), relative_momentum['refinement_members'], strict=True):
        origin_momentum.verify_bound(contract.member(level), receipt)
    sharper_response = energy_tightness.certificate()
    for level, receipt in zip((0, 1, 2), sharper_response['members'], strict=True):
        energy_tightness.verify_bound(contract.member(level), receipt)
    energy_tightness.verify_bound(contract.member(7), sharper_response['retained_origin_member'])
    matched_strain = open_coframe.certificate()
    for level, receipt in zip((0, 1, 2), matched_strain['members'], strict=True):
        open_coframe.verify_bound(contract.member(level), receipt)
    charged_pairs = charged_string.certificate()
    charged_string.verify_bound(member, charged_pairs['selected_member_bound'])
    odd_excitation = odd_residual.certificate()
    odd_residual.verify_bound(member, odd_excitation['selected_member_bound'])
    momentum_detector = momentum_transfer.certificate()
    momentum_transfer.verify_bound(contract.member(8), momentum_detector['selected_member_bound'])
    for level, receipt in zip((0, 4, 12), momentum_detector['other_members'], strict=True):
        momentum_transfer.verify_bound(contract.member(level), receipt)
    local_energy = local_fluctuation.certificate()
    for level, receipt in zip((0, 1, 2), local_energy['members'], strict=True):
        local_fluctuation.verify_bound(contract.member(level), receipt)
    local_fluctuation.verify_bound(contract.member(32), local_energy['retained_origin_member'])
    prepared_spectrum = spectral_resolution.certificate()
    spectral_resolution.verify_bound(member, prepared_spectrum['selected_member_bound'])
    joint_charge = charged_readout.certificate()
    charged_readout.verify_bound(member, joint_charge['selected_member_bound'])
    checkerboard = staggered_response.certificate()
    staggered_response.verify_bound(member, checkerboard['selected_member_bound'])
    phase_histories = current_tail.certificate()
    for level, receipt in zip((0, 1, 2, 8), phase_histories['members'], strict=True):
        current_tail.verify_bound(contract.member(level), receipt)
    physical_sectors = sector_threshold.certificate()
    sector_threshold.verify_bound(member, physical_sectors['selected_member_bound'])
    scanned_records = populated_momentum.certificate()
    populated_momentum.verify_bound(member, scanned_records['selected_member_bound'])
    strength_boundary = weighted_survival.certificate()
    for level, receipt in zip((0, 1, 2, 8), strength_boundary['members'], strict=True):
        weighted_survival.verify_bound(contract.member(level), receipt)
    conjugate_sectors = charge_conjugation.certificate()
    for level, receipt in zip((0, 1, 2, 8), conjugate_sectors['members'], strict=True):
        charge_conjugation.verify_bound(contract.member(level), receipt)
    composite_motion = composite_transport.certificate()
    composite_transport.verify_bound(member, composite_motion['selected_member_bound'])
    composite_labels = composite_spin.certificate()
    composite_spin.verify_bound(member, composite_labels['selected_member_bound'])
    original_odd_records = selected_odd_response.certificate()
    for level, receipt in zip((0, 1, 2, 8), original_odd_records['members'], strict=True):
        selected_odd_response.verify_bound(contract.member(level), receipt)
    cubic_energy = covariant_stress.certificate()
    for level, receipt in zip((0, 1, 2, 8), cubic_energy['members'], strict=True):
        covariant_stress.verify_bound(contract.member(level), receipt)
    original_spin_records = selected_spin.certificate()
    for level, receipt in zip((0, 1, 2, 8), original_spin_records['members'], strict=True):
        selected_spin.verify_bound(contract.member(level), receipt)
    full_odd_resolvent = odd_resolvent.certificate()
    odd_resolvent.verify_bound(member, full_odd_resolvent['selected_member_bound'])
    complete_odd_filters = odd_filter.certificate()
    for level, receipt in zip((0, 1, 2, 8), complete_odd_filters['members'], strict=True):
        odd_filter.verify_bound(contract.member(level), receipt)
    regional_cuts = geometric_cut.certificate()
    for level, receipt in zip((0, 1, 2, 8), regional_cuts['members'], strict=True):
        geometric_cut.verify_bound(contract.member(level), receipt)
    spatial_odd_records = coherent_odd.certificate()
    coherent_odd.verify_bound(member, spatial_odd_records['members'][0])
    density_shells = original_shells.certificate()
    original_shells.verify_bound(member, density_shells['selected_member_bound'])
    original_branches = original_branch.certificate()
    original_branch.verify_bound(member, original_branches['members'][0])
    critical_profiles = critical_refinement.certificate()
    for level, receipt in zip((0, 1, 2, 8), critical_profiles['members'], strict=True):
        critical_refinement.verify_bound(contract.member(level), receipt)
    boundary_records = boundary_information.certificate()
    for level, receipt in zip((0, 1, 2, 8), boundary_records['members'], strict=True):
        boundary_information.verify_bound(contract.member(level), receipt)
    admitted_branches = branch_certificate.certificate()
    branch_certificate.verify_bound(member, admitted_branches['selected_member_bound'])
    odd_sum_records = odd_sum_rules.certificate()
    odd_sum_rules.verify_bound(member, odd_sum_records['members'][0])
    for level, receipt in zip((0, 1, 2, 8), odd_sum_records['chirality_members'], strict=True):
        odd_sum_rules.verify_chirality_bound(contract.member(level), receipt)
    phase_spectra = work_spectrum.certificate()
    for level, receipt in zip((0, 1, 2, 8), phase_spectra['members'], strict=True):
        work_spectrum.verify_bound(contract.member(level), receipt)
    boundary_correlations = boundary_correlation.certificate()
    for level, receipt in zip((0, 1, 2, 8), boundary_correlations['members'], strict=True):
        boundary_correlation.verify_bound(contract.member(level), receipt)
    reduced_transitions = reduced_sector.certificate()
    reduced_sector.verify_bound(member, reduced_transitions['selected_member_bound'])
    regional_fermions = regional_fermion.certificate()
    regional_fermion.verify_bound(member, regional_fermions['members'][0])
    weighted_tails = weighted_spectral_tail.certificate()
    for level, receipt in zip((0, 1, 2, 8), weighted_tails['members'], strict=True):
        weighted_spectral_tail.verify_bound(contract.member(level), receipt)
    controlled_energy = energy_response.certificate()
    for level, receipt in zip((0, 1, 2, 8), controlled_energy['members'], strict=True):
        energy_response.verify_bound(contract.member(level), receipt)
    fused_loops = loop_fusion.certificate()
    loop_fusion.verify_bound(member, fused_loops['selected_member_bound'])
    equal_work_records = equal_work_control.certificate()
    for level, receipt in zip((0, 1, 2, 8), equal_work_records['members'], strict=True):
        equal_work_control.verify_bound(contract.member(level), receipt)
    if any(hashlib.sha256((ROOT/path).read_bytes()).hexdigest() != value for path, value in pins.items()):
        raise ValueError('reconstruction inputs changed during replay')
    return {
        'schema': 'oph.joint_reconstruction.v1', 'program': contract.PROGRAM,
        'member_id': member['member_id'], **member['ids'],
        'nodes': admissions, 'dependency_count': len(admissions),
        'finite_primary_acceptance_count': sum(row['selected_primary_finite_acceptance'] for row in admissions),
        'candidate_id': candidate.EXTENSION_ID,
        'conditional_candidate_scope': conditional['scope'],
        'conditional_positive_energy_transfer': positive,
        'primary_common_filtered_algebra': field_algebra,
        'primary_complete_filtered_apparatus': filtered,
        'primary_finite_energy_and_regional_response': source_stress,
        'original_weight_positive_window_diagnostic': spectral_window,
        'primary_finite_excitation_and_extreme_sector_weight': finite_survival,
        'primary_regular_spatial_algebra': spatial_algebra,
        'primary_excitation_to_mass_readout': mass,
        'primary_capacity_family_in_nats': capacities,
        'primary_compact_four_outcome_probe': compact_probe,
        'primary_complete_spatial_Fourier_readout': spatial_readout,
        'primary_conserved_occupations_and_cubic_channel': channel_labels,
        'primary_decoupled_expansion_and_spectral_scale': spectral_scale,
        'primary_interacting_endpoint_ground_bounds': endpoint_bounds,
        'primary_interacting_variational_bounds': variational_bounds,
        'primary_complete_joint_energy_wavevector_readout': joint_spectroscopy,
        'primary_occupation_sensitive_ground_and_state_bound': occupation_bounds,
        'primary_exact_local_moment_residual': exact_residual,
        'primary_original_weight_neutral_pair_response': neutral_response,
        'primary_populated_cubic_parity_pair_channels': neutral_channels,
        'primary_scalable_collective_response': scalable_response,
        'primary_modular_work_and_geometric_premises': modular_work,
        'primary_concentrated_pair_response_and_odd_trial': concentrated_response,
        'primary_retained_origin_relative_momentum': relative_momentum,
        'primary_sharper_capacity_dependent_energy_windows': sharper_response,
        'primary_boundary_completed_spatial_strain': matched_strain,
        'primary_charged_endpoint_pair_response': charged_pairs,
        'primary_scalar_odd_excitation_and_coherent_response': odd_excitation,
        'primary_relative_momentum_detector_attachment': momentum_detector,
        'primary_local_Weyl_original_energy_tightness': local_energy,
        'primary_prepared_scalar_odd_spectral_window': prepared_spectrum,
        'primary_joint_charge_transfer_and_energy_records': joint_charge,
        'primary_checkerboard_coherent_response': checkerboard,
        'primary_common_unfiltered_phase_histories': phase_histories,
        'primary_physical_occupation_sector_thresholds': physical_sectors,
        'primary_scanned_joint_spatial_energy_records': scanned_records,
        'primary_complete_phase_strength_boundary': strength_boundary,
        'primary_exact_conjugate_sector_energy_matching': conjugate_sectors,
        'primary_physical_composite_transport_moments': composite_motion,
        'primary_physical_composite_spin_and_exchange': composite_labels,
        'primary_selected_odd_creation_energy_records': original_odd_records,
        'primary_cubic_covariant_energy_and_boundary_improvement': cubic_energy,
        'primary_selected_spin_resolved_creation_energy_records': original_spin_records,
        'primary_actual_coupling_odd_full_complement_resolvent': full_odd_resolvent,
        'primary_complete_selected_odd_energy_filters': complete_odd_filters,
        'primary_covariant_regional_cuts_and_boundary_entropy': regional_cuts,
        'primary_original_coherent_spatial_odd_energy_records': spatial_odd_records,
        'primary_original_density_shell_energy_moments': density_shells,
        'primary_original_total_cube_and_number_branches': original_branches,
        'primary_critical_profile_work_refinement': critical_profiles,
        'primary_boundary_center_information_and_record_work': boundary_records,
        'primary_admitted_two_link_branch_and_symmetry_reduction': admitted_branches,
        'primary_original_odd_sum_rules_and_native_chirality_pairing': odd_sum_records,
        'primary_critical_work_spectrum_and_fluctuation_records': phase_spectra,
        'primary_physical_boundary_correlation_upper_bound': boundary_correlations,
        'primary_extreme_gauge_isolated_transition': reduced_transitions,
        'primary_regional_charged_spin_and_pair_energy_records': regional_fermions,
        'primary_selected_weighted_spectral_tail': weighted_tails,
        'primary_original_energy_current_response': controlled_energy,
        'primary_extreme_gauge_isolated_transition_proved': True,
        'primary_regional_charged_spin_records_proved': True,
        'primary_selected_weighted_spectral_tail_proved': True,
        'primary_original_energy_current_response_proved': True,
        'primary_static_loop_fusion_energy_records': fused_loops,
        'primary_same_source_equal_work_elastic_control': equal_work_records,
        'primary_static_loop_fusion_energy_records_proved': True,
        'primary_same_source_equal_work_control_proved': True,
        'primary_original_finite_graded_neighbor_contrast_proved': True,
        'primary_selected_native_parity_chirality_pairing_proved': True,
        'primary_odd_sector_exact_multiplicity_space_reduction_proved': True,
        'primary_critical_phase_finite_fluctuation_premises_instantiated': True,
        'primary_boundary_correlation_bound_with_polynomial_remainder_proved': True,
        'primary_common_subsequential_algebra_and_positive_GNS_proved': True,
        'primary_finite_positive_energy_weight_proved': True,
        'primary_uniform_fixed_band_positive_energy_weight_proved': False,
        'primary_spatial_action_on_regular_algebra_proved': True,
        'primary_selected_state_spatial_invariance_proved': False,
        'relative_observation_state_spatial_invariance_proved': True,
        'relative_observation_commuting_positive_energy_and_momentum_proved': True,
        'relative_observation_nonzero_momentum_response_proved': False,
        'relative_observation_Fourier_detector_spectral_kernel_proved': True,
        'primary_scalar_odd_excited_eigenvalue_enclosed': True,
        'primary_charged_endpoint_original_energy_band_proved': True,
        'primary_local_Weyl_energy_probability_tightness_proved': True,
        'primary_local_Weyl_limit_proved_inelastic': False,
        'primary_unfiltered_Weyl_and_filtered_probes_common_time_representation_proved': True,
        'primary_unfiltered_Weyl_spatial_momentum_attachment_proved': False,
        'primary_prepared_excited_window_with_retained_displaced_state_proved': True,
        'primary_undisturbed_joint_charge_transfer_energy_history_proved': True,
        'primary_member_zero_positive_odd_CAR_sector_thresholds_proved': True,
        'primary_exact_charge_conjugation_on_same_source_proved': True,
        'primary_member_zero_odd_composite_transport_proved': True,
        'primary_local_composite_gauge_spin_kernels_and_graded_exchange_proved': True,
        'primary_selected_odd_creation_original_energy_band_proved': True,
        'primary_local_odd_probe_has_uniform_energy_escape_under_refinement': True,
        'primary_cubic_energy_allocation_and_finite_Ward_premises_proved': True,
        'primary_selected_G1_H_probe_original_energy_records_proved': True,
        'primary_actual_coupling_prepared_odd_full_complement_resolvent_bounded': True,
        'primary_complete_odd_filter_with_original_band_weight_proved': True,
        'primary_covariant_regional_entropy_and_boundary_representation_bound_proved': True,
        'primary_member_zero_each_coherent_spatial_odd_row_populated': True,
        'primary_original_coherent_spatial_odd_dispersion_proved': False,
        'primary_original_density_spatial_energy_mean_contrast_proved': True,
        'primary_native_block_hopping_expectations_strictly_negative': True,
        'primary_original_total_cube_and_native_number_joint_records_proved': True,
        'primary_critical_profile_Dirichlet_work_limit_proved': True,
        'primary_boundary_center_record_source_work_and_coherence_proved': True,
        'primary_selected_odd_filter_proves_fixed_energy_continuum_residue': False,
        'primary_composite_transport_is_an_isolated_particle_band': False,
        'primary_scanned_analyzer_labels_are_conserved_momenta': False,
        'primary_subcritical_normalized_phase_histories_lose_fixed_positive_energy': True,
        'primary_finite_aperture_Fourier_instrument_proved': True,
        'primary_conserved_particle_momentum_proved': False,
        'primary_interacting_endpoint_gap_proved': False,
        'primary_interacting_ground_energy_enclosure_proved': True,
        'primary_isolated_excitation_branch_proved': False,
        'primary_isolated_excitation_branch_domain': 'A physical-vacuum particle branch with controlled dispersion and original spectral weight; the separately certified extreme-gauge finite transition does not establish this branch.',
        'declared_SI_and_GeV_comparison_conversion_available': True,
        'source_derived_particle_mass_table_GeV': [],
        'primary_proper_geometry_or_Einstein_limit_proved': False,
        'conditional_candidate_calibration': 'inherited common-scale action coefficients; no physical pole matching supplied',
        'full_particle_spectrum_derived': False,
        'interacting_relativistic_primary_limit_proved': False,
        'pure_Yang_Mills_gap_proved': False,
        'all_reference_targets_numerically_matched': conditional['conditional_physical_completion_arm']['all_reference_targets_numerically_matched'],
        'next_shared_acceptance': 'joint_refinement',
        'input_sha256': pins,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = candidate.canonical(build())
    if args.write:
        OUTPUT.write_text(result)
    elif not OUTPUT.is_file() or OUTPUT.read_text() != result:
        raise ValueError('unified reconstruction receipt differs from actual producers')
    print('PASS: ten dependencies; three finite primary acceptances; physical interfaces explicit')


if __name__ == '__main__':
    main()
