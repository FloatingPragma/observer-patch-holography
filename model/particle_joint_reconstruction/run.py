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
    paths = (set(stored_candidate['source_sha256']) | set(contract.SOURCE_INPUTS)
             | set(positive_energy.SOURCE_INPUTS))
    for module in (algebra, stress, window, filtered_probe, survival,
                   mass_readout, spatial, capacity, compact, momentum, channels, scale,
                   endpoint, variational, spectroscopy, occupation, local_residual,
                   pair_response):
        paths.update(module.SOURCE_INPUTS)
        paths.update(str(path.relative_to(ROOT)) for path in Path(module.__file__).parent.iterdir()
                     if path.is_file() and path.suffix in ('.py', '.md'))
    paths.update(str(path.relative_to(ROOT)) for path in HERE.iterdir()
                 if path.suffix in ('.py', '.md') or path.name in
                 ('dependencies.json', 'candidate_results.json'))
    pins = provenance.definition_hashes(provenance.definition_inputs(paths))
    graph = strict_json.load(HERE / 'dependencies.json')
    validate_dependencies(graph)
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
        'primary_common_subsequential_algebra_and_positive_GNS_proved': True,
        'primary_finite_positive_energy_weight_proved': True,
        'primary_uniform_fixed_band_positive_energy_weight_proved': False,
        'primary_spatial_action_on_regular_algebra_proved': True,
        'primary_selected_state_spatial_invariance_proved': False,
        'primary_finite_aperture_Fourier_instrument_proved': True,
        'primary_conserved_particle_momentum_proved': False,
        'primary_interacting_endpoint_gap_proved': False,
        'primary_interacting_ground_energy_enclosure_proved': True,
        'primary_isolated_excitation_branch_proved': False,
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
