"""Validate source custody and the domains of the four defect/field results.

This registry check does not replace the component mathematics or replay.
"""
import argparse
import hashlib
import json
from pathlib import Path

from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
PATH = Path('model/completion/defect_selection_campaign.json')
CAMPAIGN = 'OPH-DEFECT-SELECTION-CONSTRAINTS-2026-11-04'
ISSUES = [1139, 1141, 1142, 1143, 1145]
COMPONENTS = {
    'winding_sector_search': 'FA-133',
    'radial_drift_readout': 'FA-134',
    'majorana_defect_bridge': 'FA-135',
    'higgs_curvature_response': 'FA-136',
}
BOUNDARY = {
    'simulator_observer_scheduler_repair_memory_specs_admitted': True,
    'winding_index_is_measured_generation_count': False,
    'drift_bound_requires_no_compensation': True,
    'timed_readout_rate_and_port_embedding_are_supplied': True,
    'Majorana_coefficient_selected_by_native_current_constraints': False,
    'conserved_additive_charges_must_respect_actual_mixing_support': True,
    'varying_Higgs_curvature_coefficient_uses_generalized_volume': True,
    'physical_entropy_stationarity_instantiated': False,
    'same_current_defect_particle_map_constructed': False,
    'nature_couplings_selected': False,
    'common_interacting_continuum_and_Einstein_model_constructed': False,
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def validate(packet, root=ROOT):
    root = Path(root).resolve()
    require(packet.get('schema') == 'oph.defect_selection_campaign.v1'
            and packet.get('campaign_id') == CAMPAIGN, 'wrong component campaign')
    require(packet.get('selected_issues') == ISSUES, 'five-issue selection changed')
    require(packet.get('scope') == 'FINITE_COMPONENT_FAMILY'
            and packet.get('whole_program_frozen') is False
            and packet.get('integrated_gate_accepted') is False
            and packet.get('canonical_premises_discharged') == []
            and type(packet.get('physical_comparison_count_added')) is int
            and packet['physical_comparison_count_added'] == 0, 'physical scope promoted')
    boundary = packet.get('composition_boundary', {})
    require(set(boundary) == set(BOUNDARY)
            and all(boundary[k] is v for k, v in BOUNDARY.items()), 'component domain changed')
    finite = strict_json.load(root/'model/completion/finite_attachment_map.json')
    attachments = {row['id']: row for row in finite['attachments']}
    program = strict_json.load(root/'claims/completion_program.json')
    stages = {row['id']: row for row in program['stages']}
    rows = packet.get('component_results', [])
    require([row.get('id') for row in rows] == list(COMPONENTS), 'missing or duplicate component')
    issue_set = set()
    for row in rows:
        name = row['id']
        attachment = attachments[COMPONENTS[name]]
        issues = [stages[owner]['issue'] for owner in attachment['owners']]
        issue_set.update(issues)
        require(row.get('attachment') == COMPONENTS[name] and row.get('issues') == issues,
                'component ownership changed')
        require(attachment['classification'] == 'realized_on_claimed_domain'
                and attachment['physical_identification'] is False, 'attachment scope promoted')
        for key, source in [('result', 'conclusion'), ('interface', 'interface'), ('remaining', 'remaining')]:
            require(row.get(key) == attachment[source], 'component result or domain changed')
        evidence = row.get('evidence', {})
        require(set(evidence) == set(attachment['witnesses']), 'incomplete component source custody')
        for name, digest in evidence.items():
            path = (root/name).resolve()
            require(path.is_relative_to(root) and path.is_file(), 'missing component source')
            require(hashlib.sha256(path.read_bytes()).hexdigest() == digest, 'stale component source')
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
            'canonical_premises_discharged': 0, 'source_selection_solved': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.parse_args()
    print(json.dumps(validate(strict_json.load(ROOT/PATH)), indent=2))


if __name__ == '__main__':
    main()
