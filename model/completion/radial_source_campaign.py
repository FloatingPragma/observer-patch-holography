"""Validate the same radial source kernel, interaction, entropy and native geometry join.

Custody and interface checks preserve the distinction between a finite
controlled source experiment and an interacting continuum/Einstein theory.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
from tools import strict_json
ROOT=Path(__file__).resolve().parents[2]
PATH=Path('model/completion/radial_source_campaign.json')
CAMPAIGN = 'OPH-RADIAL-SOURCE-CLOSURE-2026-11-04'
OWNERS = {'radial_defect_gap': ('FA-127', ['CP-03', 'CP-06', 'CP-07'], [1139, 1142, 1143]), 'radial_interacting_core': ('FA-128', ['CP-03', 'CP-06', 'CP-07'], [1139, 1142, 1143]), 'parity_entropy_response': ('FA-129', ['CP-08', 'CP-09'], [1144, 1145]), 'radial_geometry_join': ('FA-130', ['CP-03', 'CP-06', 'CP-07', 'CP-08'], [1139, 1142, 1143, 1144]), 'radial_reused_readout': ('FA-131', ['CP-07', 'CP-08'], [1143, 1144])}
ISSUES = [1139, 1142, 1143, 1144, 1145]
INTERFACES = {'same_bounded_radial_source_for_kernel_and_full_Fock_interaction': True, 'strong_mass_kernel_exhaustion_and_positive_complement_gap_proved': True, 'supplied_radial_mass_texture_is_selected_by_native_law': False, 'full_Fock_virtual_complement_retained': True, 'actual_interacting_ground_and_charged_poles_attached_to_local_parity': True, 'actual_ground_central_M4_reduction_faithful': True, 'actual_ground_full_Fock_state_faithful': False, 'same_full_H_faithful_affine_family_has_regional_entropy_energy_relation': True, 'source_native_regional_modular_word_retains_complement_failures': True, 'regional_modular_flow_has_stated_local_contact_time_factor': True, 'regional_modular_clock_is_full_Hamiltonian_flow': False, 'same_actual_parity_drives_rank_three_native_Gram': True, 'controlled_readout_work_and_complete_entropy_records_accounted': True, 'field_preparation_and_auxiliary_services_are_supplied': True, 'fresh_field_copies_required_inside_readout_only_QND_block': False, 'reused_field_block_mean_is_latent_parity_mixture': True, 'reused_field_block_mean_equals_power_of_averaged_mean': False, 'repeated_readout_estimates_ensemble_parity_from_one_field': False, 'Hamiltonian_drift_absent_or_compensated_during_readout_block': True, 'arbitrary_native_repair_protects_radial_kernel_or_charged_poles': False, 'signed_kinetic_bias_identified_with_geometric_area': False, 'regional_Einstein_fixed_volume_premises_instantiated': False, 'autonomous_reciprocal_Einstein_dynamics_derived': False, 'massive_interacting_relativistic_continuum_proved': False, 'physical_particle_identification': False, 'one_common_interacting_QFT_and_responsive_Einstein_source': False}
JOIN = {'components': ['radial_defect_gap', 'radial_interacting_core', 'parity_entropy_response', 'radial_geometry_join', 'radial_reused_readout'], 'source': 'The bounded radial radius-one 27-site Hamiltonian at a=v=1,m=64,k=2,U=1/65536, with the complete 54-mode Fock algebra and central D0=(n_origin-1)^2.', 'kernel_to_interaction': 'FA127 supplies the exact two-mode kernel and positive complement gap of the same bounded radial source. FA128 retains the full Fock complement and derives its actual N27 Kramers ground doublet, unique PH-related neighboring charged states and positive spectral residues.', 'interaction_to_geometry': 'FA130 copies the actual central occupation parity into the FA126 charge-probe input, retains every native quantum failure and all scheduler/current records, and derives the same-source charged-versus-ground Gram contrast with measurement work.', 'interaction_to_entropy': 'FA129 proves the actual ground central M4 reduction faithful, constructs a stationary faithful affine family of this same H, and relates its regional modular generator, entropy, full energy and native response bias. The modular/contact time factor differs from the full Gibbs clock, vanishes at q=1/2 and changes sign across that point.', 'finite_reused_readout': 'FA131 retains the complete history of a finite readout-only block on one prepared field. The native mean is (1-q)M0^t+q M1^t; fresh auxiliaries and compensated Hamiltonian drift are inputs. The 32-round classifier assumes the promised ground or charged input; its postmeasurement field need not remain in that spectral sector. No ensemble-mean power law or ensemble-q estimate from a single field is inferred.', 'state_preparation': 'Actual spectral-state selection and preparation, the radial texture, blank work memory and addressed control accuracy are supplied inputs. Mathematical spectral projectors are not executed preparations.', 'clock': 'The source Hamiltonian defines the pole/energy clock. Native modular/contact time, pulse duration and idle-inclusive scheduler attempt counts have their stated conversions and control conditions. No physical Lorentzian clock is inferred.', 'common_finite_matter_reference_and_controlled_response_constructed': True, 'common_continuum_QFT_Einstein_reference_and_refinement_constructed': False, 'regional_entropy_is_geometric_area': False, 'full_interacting_Hamiltonian_runs_uncompensated_during_QND_block': False, 'all_canonical_physical_premises_discharged': False, 'interleaved_field_repair_admitted_during_QND_block': False, 'readout_preserves_input_spectral_sector': False, 'classification_bound_requires_promised_ground_or_charged_input': True, 'ensemble_Gram_inferred_from_single_reused_field': False}
BOUNDARY = {'automatic_cross_domain_join': False, 'common_finite_matter_reference_and_controlled_response_constructed': True, 'common_continuum_QFT_Einstein_reference_and_refinement_constructed': False, 'components': ['radial_defect_gap', 'radial_interacting_core', 'parity_entropy_response', 'radial_geometry_join', 'radial_reused_readout'], 'reason': 'The same bounded radial source is attached through its actual central parity to regional entropy and controlled native geometry. Its finite spectral preparation and readout schedule are supplied operations. Neither a volume-uniform interacting continuum nor regional area/stress/fixed-volume Einstein premises follow from this finite join.'}


DOMAINS = {'particle_domain': 'The same bounded radial source has an exact two-dimensional zero kernel and controlled complement gap. On its finite27-site full Fock space, the centered interaction gives actual ground/charged states and positive spectral residues; the interaction estimate depends on volume. The strong-mass condition supplies no fixed-mass spatial continuum limit.', 'geometry_domain': 'The actual central parity drives the complete native measurement and all-seam scheduler. The same regional state has a native modular word and an entropy/energy/bias relation. A finite readout-only block retains one field and every record, with fresh apparatus and compensated Hamiltonian drift. The rank-three Gram and signed kinetic bias are not an area, Lorentzian stress or Einstein premise instance.'}

def require(condition, message):
    if not condition:
        raise ValueError(message)


def validate(packet, root=ROOT):
    root = Path(root).resolve()
    require(packet.get('schema') == 'oph.radial_source_campaign.v1'
            and packet.get('campaign_id') == CAMPAIGN, 'wrong radial-source contract')
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
    require(packet.get('composition_boundary') == BOUNDARY,
            'finite same-source join or continuum boundary changed')
    interface = packet.get('conditional_interfaces', {})
    require(all(interface.get(k) is v for k, v in INTERFACES.items()),
            'source, entropy, readout or Einstein domain promoted or erased')
    require(set(interface) == set(INTERFACES)|set(DOMAINS)
            and all(interface.get(k) == value for k, value in DOMAINS.items()),
            'explicit particle or geometry source domain changed')
    require(packet.get('supported_source_join') == JOIN,
            'radial source join broadened or erased')
    return {'selected_issues': ISSUES, 'reviewed_components': len(rows),
            'canonical_premises_discharged': 0, 'integrated_gate_accepted': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.parse_args()
    print(json.dumps(validate(strict_json.load(ROOT/PATH)), indent=2))


if __name__ == '__main__':
    main()
