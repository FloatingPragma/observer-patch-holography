"""Check the five-issue interacting-core, exact-index and responsive-geometry campaign and its explicit non-composition.

This verifies custody and scope; mathematical proofs remain in the source
packages. Separate finite algebras and references are never joined by name.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from tools import strict_json


ROOT = Path(__file__).resolve().parents[2]
PATH = Path('model/completion/defect_geometry_campaign.json')
CAMPAIGN = 'OPH-DEFECT-GEOMETRY-CLOSURE-2026-11-04'
OWNERS = {
    'defect_core_response': ('FA-123', ['CP-04', 'CP-06', 'CP-07'], [1140, 1142, 1143]),
    'bounded_defect_index': ('FA-124', ['CP-03', 'CP-06', 'CP-07'], [1139, 1142, 1143]),
    'native_odd_box_core': ('FA-125', ['CP-04', 'CP-06', 'CP-07'], [1140, 1142, 1143]),
    'matter_geometry_response': ('FA-126', ['CP-08'], [1144]),
}
ISSUES = [1139, 1140, 1142, 1143, 1144]
INTERFACES = {
    'actual_interacting_ground_spectral_measures_have_positive_pole_weight': True,
    'unmodified_odd_cube_uses_no_supplied_mass_texture': True,
    'all_Fock_virtual_states_retained_in_core_bounds': True,
    'odd_cube_zero_mode_is_boundary_standing_wave': True,
    'bounded_radial_source_has_exact_lattice_zeros': True,
    'bounded_index_grading_is_nonlocal_spatial_inversion': True,
    'bounded_index_changes_sign_with_texture_orientation': False,
    'bounded_kernel_exhaustion_proved': False,
    'native_repair_protects_Hamiltonian_zero_modes_or_poles': False,
    'continuing_small_positive_repair_rate_is_a_specified_parameter': True,
    'continuing_comparison_conditions_on_no_repair_event': False,
    'source_energy_changes_Gram_under_one_fixed_instrument': True,
    'fresh_probe_independence_is_geometry_iteration_input': True,
    'arbitrary_correlated_probe_reuse_has_same_mean_law': False,
    'ideal_code_completeness_certifies_physical_service_failures': False,
    'bounded_radial_index_inherits_finite_linear_texture_interaction_proof': False,
    'finite_core_source_mapped_to_responsive_native_geometry': False,
    'regional_Einstein_premises_instantiated': False,
    'massive_interacting_relativistic_continuum_proved': False,
    'physical_particle_identification': False,
    'one_common_interacting_QFT_and_responsive_Einstein_source': False,
}
JOIN = {
    'result_component': 'defect_core_response',
    'source_results': ['FA-120 linear finite core', 'FA-125 unmodified odd open core'],
    'state': 'Each parent separately supplies its actual N27 interacting Kramers ground mixture, charged N26/N28 projectors and full Fock complement.',
    'clock': 'The corresponding centered Hamiltonian time and its supplied coefficients; no cross-source energy calibration.',
    'conclusion': 'Positive addition/removal spectral atoms and controlled one-body residue matrices on each parent source, with native number and energy-current identities and complete repair-record comparison.',
    'positive_rate_witness': 'Each of 54 seams has chosen Gamma=g_lower/54000. At T=1/g_lower the complete archive distance is below1/1000 for the exact ground, with every repair event retained; the trial mixture adds t.',
    'cross_parent_states_or_energy_levels_identified': False,
    'native_interacting_ground_preparation_executed': False,
    'native_repair_invariance_inferred_from_Hamiltonian_poles': False,
    'join_to_bounded_radial_index_or_responsive_Gram': False,
}



def require(condition, message):
    if not condition:
        raise ValueError(message)


def validate(packet, root=ROOT):
    root = Path(root).resolve()
    require(packet.get('schema') == 'oph.defect_geometry_campaign.v1'
            and packet.get('campaign_id') == CAMPAIGN, 'wrong defect-geometry contract')
    require(packet.get('scope') == 'FINITE_COMPONENT_FAMILY'
            and packet.get('whole_program_frozen') is False,
            'finite campaign cannot freeze the whole program')
    require(packet.get('selected_issues') == ISSUES, 'five-issue ownership changed')
    require(packet.get('physical_comparison_count_added') == 0
            and type(packet.get('physical_comparison_count_added')) is int
            and packet.get('canonical_premises_discharged') == []
            and packet.get('integrated_gate_accepted') is False
            and packet.get('common_interacting_QFT_and_Einstein_model') is False,
            'component evidence cannot accept physical completion')
    require(packet.get('specification_readout_allowed') is True,
            'the declared simulator and observer specifications are admitted inputs')
    rows = packet.get('component_results', [])
    require([r.get('id') for r in rows] == list(OWNERS), 'missing or duplicate component')
    finite = strict_json.load(root/'model/completion/finite_attachment_map.json')
    attachments = {a['id']: a for a in finite['attachments']}
    program = strict_json.load(root/'claims/completion_program.json')
    stages = {stage['id']: stage for stage in program['stages']}
    for row in rows:
        identifier = row['id']
        attachment, owners, issues = OWNERS[identifier]
        require(row.get('attachment') == attachment and row.get('issues') == issues,
                f'component owner changed: {identifier}')
        source = attachments[attachment]
        require(source['owners'] == owners
                and source['classification'] == 'realized_on_claimed_domain'
                and source['physical_identification'] is False,
                f'finite map scope changed: {identifier}')
        require(row.get('result') == source['conclusion']
                and row.get('remaining') == source['remaining'],
                f'component result or boundary differs from its reviewed attachment: {identifier}')
        require(row.get('interface') == source['interface'],
                f'component state/algebra/operation/clock/refinement changed: {identifier}')
        evidence = row.get('evidence', {})
        require(isinstance(evidence, dict) and bool(evidence), f'missing source custody: {identifier}')
        require(set(evidence) == set(source['witnesses']),
                f'incomplete component artifact custody: {identifier}')
        for name, digest in evidence.items():
            path = (root/name).resolve()
            require(path.is_relative_to(root) and path.is_file(), f'missing source: {name}')
            require(hashlib.sha256(path.read_bytes()).hexdigest() == digest,
                    f'stale component source: {name}')
        receipt = f'model/{identifier}/results.json'
        require(receipt in evidence, f'unbound receipt: {identifier}')
        for owner in owners:
            progress = [p for p in stages[owner].get('progress_evidence', [])
                        if p.get('campaign_id') == CAMPAIGN and p.get('id') == identifier]
            require(len(progress) == 1 and progress[0]['path'] == receipt
                    and progress[0]['sha256'] == evidence[receipt],
                    f'completion progress does not bind this component: {identifier}/{owner}')
            require(not any(p.get('campaign_id') == CAMPAIGN
                            for p in stages[owner].get('accepted_evidence', [])),
                    f'component campaign entered accepted evidence: {owner}')
    boundary = packet.get('composition_boundary', {})
    require(boundary.get('automatic_cross_domain_join') is False
            and boundary.get('common_state_reference_constructed') is False
            and boundary.get('common_clock_and_refinement_constructed') is False
            and boundary.get('distinct_domains') == list(OWNERS)
            and bool(boundary.get('reason')),
            'distinct source states cannot be silently composed')
    interface = packet.get('conditional_interfaces', {})
    require(all(interface.get(k) is v for k, v in INTERFACES.items()),
            'source, interaction, index or geometry domain promoted or erased')
    require(all(bool(interface.get(k)) for k in ('particle_domain', 'geometry_domain')),
            'missing explicit particle or geometry source domain')
    require(packet.get('supported_source_join') == JOIN,
            'finite pole source join broadened or erased')
    return {'selected_issues': ISSUES, 'reviewed_components': len(rows),
            'canonical_premises_discharged': 0, 'integrated_gate_accepted': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.parse_args()
    print(json.dumps(validate(strict_json.load(ROOT/PATH)), indent=2))


if __name__ == '__main__':
    main()
