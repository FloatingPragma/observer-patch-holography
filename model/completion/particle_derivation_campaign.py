"""Validate the finite source-particle campaign and its precise domains.

The campaign binds the six mathematical results and the source-particle
map; it does not replace their proofs, scientific replay or physical matching.
"""
import argparse
import hashlib
import json
from pathlib import Path

from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
PATH = Path('model/completion/particle_derivation_campaign.json')
SOURCE_MAP_PATH = 'model/completion/SOURCE_PARTICLE_MAP.md'
COMPARISON_PATH = 'model/completion/source_flavor_comparison.json'
SCHEMA = 'oph.particle_derivation_campaign.v1'
CAMPAIGN = 'OPH-SOURCE-PARTICLE-DERIVATION-2026-11-04'
ISSUES = [1139, 1140, 1141, 1142, 1143]
COMPONENTS = {
    'certified_winding_core': 'FA-137',
    'radial_current_map': 'FA-138',
    'native_flavor_orbit': 'FA-139',
    'lattice_texture_selection': 'FA-140',
    'mcpr_native_memory_bridge': 'FA-141',
    'defect_current_bridge': 'FA-142',
}
BOUNDARY = {
    'simulator_observer_scheduler_repair_memory_specs_admitted': True,
    'finite_winding_spectra_enclosed_with_multiplicity': True,
    'full_Fock_contact_bounds_retain_all_virtual_sectors': True,
    'finite_spectral_certification_is_interacting_continuum': False,
    'radial_eta_identified_as_Standard_Model_or_Lorentz_current': False,
    'all_nonlocal_or_broken_gauge_currents_classified': False,
    'native_same_current_basis_changes_preserve_flavor_invariants': True,
    'source_flavor_architecture_and_duty_policy_are_supplied': True,
    'empirical_Yukawa_coefficients_consumed_by_flavor_calculation': False,
    'static_mean_real_projector_feedback_has_zero_CP_invariant': True,
    'static_mean_zero_CP_applies_to_all_correlated_feedback': False,
    'finite_unread_feedback_equals_unitary_mean_action': False,
    'P_dependent_memory_configuration_is_a_declared_design': True,
    'target_informed_mass_architecture_is_an_independent_prediction': False,
    'localized_native_current_mass_join_uses_extended_carrier': True,
    'same_localized_CAR_current_and_P_mass_operator_constructed': True,
    'record_driven_quadratic_limit_controlled_at_fixed_volume': True,
    'local_Gauss_constraint_or_dynamical_gauge_fields_derived': False,
    'localized_fermion_occupation_creates_texture_winding': False,
    'finite_energy_topology_bound_uses_selected_ordered_rotation_orbit': True,
    'finite_energy_derives_eventual_shell_admissibility_and_zero_degree': True,
    'topology_energy_crossing_floor_is_particle_lifetime': False,
    'bare_protected_repair_identity_equals_full_Hamiltonian_evolution': False,
    'source_texture_vacuum_and_absolute_mass_scale_selected': False,
    'same_source_current_map_to_observed_particles_complete': False,
    'common_interacting_continuum_and_Einstein_model_constructed': False,
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def _bound_file(root, name, digest, role):
    require(type(name) is str and type(digest) is str, 'invalid '+role+' binding')
    path = (root/name).resolve()
    require(path.is_relative_to(root) and path.is_file(), 'missing or external '+role)
    require(hashlib.sha256(path.read_bytes()).hexdigest() == digest, 'stale '+role)


def validate(packet, root=ROOT):
    root = Path(root).resolve()
    require(type(packet) is dict and packet.get('schema') == SCHEMA
            and packet.get('campaign_id') == CAMPAIGN, 'wrong component campaign')
    require(packet.get('selected_issues') == ISSUES, 'five-issue selection changed')
    require(packet.get('scope') == 'FINITE_COMPONENT_FAMILY'
            and packet.get('whole_program_frozen') is False
            and packet.get('integrated_gate_accepted') is False
            and packet.get('canonical_premises_discharged') == []
            and type(packet.get('physical_comparison_count_added')) is int
            and packet['physical_comparison_count_added'] == 0, 'physical scope promoted')
    boundary = packet.get('composition_boundary')
    require(type(boundary) is dict and set(boundary) == set(BOUNDARY)
            and all(boundary[k] is v for k, v in BOUNDARY.items()), 'component domain changed')
    source_map = packet.get('source_map')
    require(type(source_map) is dict and set(source_map) == {'path', 'sha256'}
            and source_map['path'] == SOURCE_MAP_PATH, 'source-particle map binding changed')
    _bound_file(root, source_map['path'], source_map['sha256'], 'source-particle map')
    comparison = packet.get('downstream_comparison')
    require(type(comparison) is dict and set(comparison) == {'path', 'sha256'}
            and comparison['path'] == COMPARISON_PATH, 'downstream comparison binding changed')
    _bound_file(root, comparison['path'], comparison['sha256'], 'downstream comparison')
    from .source_flavor_comparison import validate as validate_comparison
    validate_comparison(strict_json.load(root/COMPARISON_PATH), root)
    finite = strict_json.load(root/'model/completion/finite_attachment_map.json')
    attachments = {row['id']: row for row in finite['attachments']}
    require(set(COMPONENTS.values()) <= set(attachments), 'missing finite attachment')
    program = strict_json.load(root/'claims/completion_program.json')
    stages = {row['id']: row for row in program['stages']}
    rows = packet.get('component_results')
    require(type(rows) is list and all(type(row) is dict for row in rows)
            and [row.get('id') for row in rows] == list(COMPONENTS), 'missing or duplicate component')
    issue_set = set()
    for row in rows:
        name = row['id']
        attachment = attachments[COMPONENTS[name]]
        require(set(attachment['owners']) <= set(stages), 'missing component owner')
        issues = [stages[owner]['issue'] for owner in attachment['owners']]
        issue_set.update(issues)
        require(row.get('attachment') == COMPONENTS[name] and row.get('issues') == issues,
                'component ownership changed')
        require(attachment['classification'] == 'realized_on_claimed_domain'
                and attachment['physical_identification'] is False, 'attachment scope promoted')
        for key, source in [('result', 'conclusion'), ('interface', 'interface'), ('remaining', 'remaining')]:
            require(row.get(key) == attachment[source], 'component result or domain changed')
        evidence = row.get('evidence')
        require(type(evidence) is dict and set(evidence) == set(attachment['witnesses']),
                'incomplete component source custody')
        for name, digest in evidence.items():
            _bound_file(root, name, digest, 'component source')
        receipt = f'model/{row["id"]}/results.json'
        require(receipt in evidence, 'missing component receipt')
        for owner in attachment['owners']:
            progress = [p for p in stages[owner].get('progress_evidence', [])
                        if p.get('campaign_id') == CAMPAIGN and p.get('id') == row['id']]
            require(len(progress) == 1 and progress[0]['path'] == receipt
                    and progress[0]['sha256'] == evidence[receipt], 'unbound completion progress')
            require(not any(p.get('campaign_id') == CAMPAIGN
                            for p in stages[owner].get('accepted_evidence', [])),
                    'component results entered physical acceptance')
    require(issue_set == set(ISSUES), 'component coverage differs from the five selected issues')
    return {'selected_issues': ISSUES, 'reviewed_components': len(rows),
            'source_map_bound': SOURCE_MAP_PATH, 'canonical_premises_discharged': 0,
            'source_selection_solved': False, 'physical_particle_map_complete': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.parse_args()
    print(json.dumps(validate(strict_json.load(ROOT/PATH)), indent=2))


if __name__ == '__main__':
    main()
