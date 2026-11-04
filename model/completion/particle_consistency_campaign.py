"""Check the finite particle consistency results and their composition limits."""
import argparse
import json
from pathlib import Path

from tools import strict_json
from .particle_derivation_campaign import ROOT, _bound_file, require

PATH = Path('model/completion/particle_consistency_campaign.json')
CAMPAIGN = 'OPH-PARTICLE-CONSISTENCY-COMPONENTS-2026-11-04'
ISSUES = [1139, 1140, 1142, 1143, 1145]
COMPONENTS = {
    'pixel_mass_closure': 'FA-143',
    'compact_defect_pair': 'FA-144',
    'interacting_defect_response': 'FA-145',
}
BOUNDARY = {
    'simulator_design_choices_admitted': True,
    'pixel_closure_uses_declared_response_model': True,
    'pixel_closure_is_native_electromagnetic_response': False,
    'interacting_spectral_response_reinserted_into_pixel_loop': False,
    'compact_pair_has_finite_total_texture_energy': True,
    'compact_pair_is_an_autonomous_stable_particle': False,
    'projected_quartic_is_a_spatially_local_interaction': False,
    'core_memory_limit_is_a_local_continuum_limit': False,
    'final_classical_history_and_core_comparison_controlled': True,
    'finite_Gibbs_first_law_establishes_geometric_stationarity': False,
    'three_results_share_one_proved_physical_limit': False,
    'observed_particle_map_complete': False,
    'interacting_QFT_and_Einstein_dynamics_established': False,
}


def validate(packet, root=ROOT):
    root = Path(root).resolve()
    require(packet.get('schema') == 'oph.particle_consistency_campaign.v1'
            and packet.get('campaign_id') == CAMPAIGN, 'wrong consistency campaign')
    require(packet.get('selected_issues') == ISSUES, 'five-issue ownership changed')
    require(packet.get('scope') == 'FINITE_COMPONENT_FAMILY'
            and packet.get('whole_program_frozen') is False
            and packet.get('integrated_gate_accepted') is False
            and packet.get('canonical_premises_discharged') == []
            and type(packet.get('physical_comparison_count_added')) is int
            and packet['physical_comparison_count_added'] == 0, 'scope promotion')
    boundary = packet.get('composition_boundary')
    require(type(boundary) is dict and set(boundary) == set(BOUNDARY)
            and all(boundary[k] is v for k, v in BOUNDARY.items()),
            'unsupported physical or composition claim')
    bindings = packet.get('source_bindings')
    require(type(bindings) is dict and bindings, 'missing source bindings')
    for path, digest in bindings.items():
        _bound_file(root, path, digest, 'particle consistency source')
    finite = strict_json.load(root/'model/completion/finite_attachment_map.json')
    attachments = {row['id']: row for row in finite['attachments']}
    program = strict_json.load(root/'claims/completion_program.json')
    stages = {row['id']: row for row in program['stages']}
    rows = packet.get('component_results')
    require(type(rows) is list and all(type(row) is dict for row in rows)
            and [row.get('id') for row in rows] == list(COMPONENTS),
            'missing or duplicated result')
    issues = set()
    for row in rows:
        name = row['id']
        require(COMPONENTS[name] in attachments, 'missing finite attachment')
        a = attachments[COMPONENTS[name]]
        require(a['classification'] == 'realized_on_claimed_domain'
                and a['physical_identification'] is False, 'unsupported attachment status')
        owned = [stages[g]['issue'] for g in a['owners']]
        issues.update(owned)
        require(row.get('attachment') == a['id'] and row.get('issues') == owned,
                'result ownership mismatch')
        for key, source in [('result', 'conclusion'), ('interface', 'interface'),
                            ('remaining', 'remaining')]:
            require(row.get(key) == a[source], 'result differs from reviewed attachment')
        require(set(a['witnesses'] + a['countercontrols']) <= set(bindings),
                'missing reviewed source or control')
        receipt = f'model/{name}/results.json'
        require(receipt in bindings, 'missing owning receipt')
        for g in a['owners']:
            progress = [p for p in stages[g].get('progress_evidence', [])
                        if p.get('campaign_id') == CAMPAIGN and p.get('id') == name]
            require(len(progress) == 1 and progress[0]['path'] == receipt
                    and progress[0]['sha256'] == bindings[receipt],
                    'unbound component progress')
            require(not any(p.get('campaign_id') == CAMPAIGN
                            for p in stages[g].get('accepted_evidence', [])),
                    'finite result entered physical acceptance')
    require(issues == set(ISSUES), 'result coverage differs from selected issues')
    return {'selected_issues': ISSUES, 'reviewed_components': len(rows),
            'canonical_premises_discharged': 0, 'physical_particle_map_complete': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.parse_args()
    print(json.dumps(validate(strict_json.load(ROOT/PATH)), indent=2))


if __name__ == '__main__':
    main()
