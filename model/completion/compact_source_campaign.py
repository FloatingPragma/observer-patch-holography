"""Validate the compact-source component family and its physical boundaries.

This checks source custody and exact correspondence with the reviewed finite
attachments. Each component's independent mathematical replay remains required.
"""
import argparse
import json
from pathlib import Path

from tools import strict_json
from .particle_derivation_campaign import ROOT, _bound_file, require

PATH = Path('model/completion/compact_source_campaign.json')
SCHEMA = 'oph.compact_source_campaign.v1'
CAMPAIGN = 'OPH-COMPACT-SOURCE-CONSISTENCY-2026-11-04'
ISSUES = [1139, 1140, 1141, 1142, 1145]
COMPONENTS = {
    'compact_defect_transport': ('FA-146', ['CP-03', 'CP-05']),
    'compact_interacting_response': ('FA-147', ['CP-04', 'CP-06']),
    'compact_band_response': ('FA-148', ['CP-04', 'CP-09']),
    'native_quantum_gauss': ('FA-149', ['CP-05', 'CP-06']),
}
OWNERS = {'CP-03': 1139, 'CP-04': 1140, 'CP-05': 1141, 'CP-06': 1142, 'CP-09': 1145}
PACKET_FIELDS = {'schema', 'campaign_id', 'selected_issues', 'scope', 'whole_program_frozen',
                 'integrated_gate_accepted', 'canonical_premises_discharged',
                 'physical_comparison_count_added', 'component_results', 'source_bindings'}
RESULT_FIELDS = {'id', 'attachment', 'issues', 'result', 'interface', 'remaining'}
INTERFACE_FIELDS = {'state', 'algebra', 'operation', 'clock', 'refinement'}
FORBIDDEN_PROMOTIONS = {
    'compact_defect_transport': (
        'autonomous_texture_formation_or_moving_core_dynamics',
        'local_Gauss_dressing_or_physical_species_identification',
        'relativistic_dispersion_or_interacting_continuum', 'conserved_gravitational_stress',
        'unchanged_native_CE_pulse_compilation'),
    'compact_interacting_response': (
        'physical_rest_masses_or_observed_species_derived',
        'interacting_chiral_relativistic_continuum_derived',
        'conserved_relativistic_stress_or_Einstein_response_derived',
        'four_mode_slice_is_full_core_vacuum', 'four_mode_slice_reduces_each_individual_record',
        'same_interaction_is_strictly_site_local',
        'interacting_electromagnetic_pixel_feedback_closed'),
    'compact_band_response': (
        'conserved_relativistic_stress_constructed',
        'hopping_drive_identified_with_metric_variation',
        'regional_area_entropy_or_Einstein_stationarity_derived',
        'Kato_transporter_is_autonomous_native_evolution'),
    'native_quantum_gauss': (
        'Wilson_unitary_links', 'same_compact_defect_band_gauge_dressing',
        'continuum_chirality_or_regulator_anomaly_cancellation',
        'physical_mass_or_coupling_calibration', 'bare_CE_compilation'),
}


def _unique_rows(rows, keys, role):
    require(type(rows) is list and all(type(row) is dict for row in rows), 'invalid '+role)
    result = {}
    for key in keys:
        matching = [row for row in rows if row.get('id') == key]
        require(len(matching) == 1, 'missing or duplicated '+role)
        result[key] = matching[0]
    return result


def validate(packet, root=ROOT):
    root = Path(root).resolve()
    require(type(packet) is dict and set(packet) == PACKET_FIELDS
            and packet.get('schema') == SCHEMA and packet.get('campaign_id') == CAMPAIGN,
            'wrong compact campaign schema')
    require(packet.get('selected_issues') == ISSUES, 'five-issue selection changed')
    require(packet.get('scope') == 'FINITE_COMPONENT_FAMILY'
            and packet.get('whole_program_frozen') is False
            and packet.get('integrated_gate_accepted') is False
            and packet.get('canonical_premises_discharged') == []
            and type(packet.get('physical_comparison_count_added')) is int
            and packet['physical_comparison_count_added'] == 0, 'physical or canonical scope promoted')
    rows = packet.get('component_results')
    require(type(rows) is list and all(type(row) is dict and set(row) == RESULT_FIELDS for row in rows)
            and [row['id'] for row in rows] == list(COMPONENTS), 'missing, duplicated or altered component')
    bindings = packet.get('source_bindings')
    require(type(bindings) is dict and bindings, 'missing source bindings')
    checked = {}

    def bound(path, digest, role):
        require(type(path) is str and type(digest) is str, 'invalid '+role+' binding')
        require(not path.startswith('docs/'), 'research note used as compact source evidence')
        if path in checked:
            require(checked[path] == digest, 'conflicting source bindings')
        else:
            _bound_file(root, path, digest, role)
            checked[path] = digest

    for path, digest in bindings.items():
        bound(path, digest, 'compact campaign source')
    finite = strict_json.load(root/'model/completion/finite_attachment_map.json')
    attachments = _unique_rows(finite.get('attachments'), [value[0] for value in COMPONENTS.values()],
                               'finite attachment')
    program = strict_json.load(root/'claims/completion_program.json')
    stages = _unique_rows(program.get('stages'), OWNERS, 'completion owner')
    require(all(stages[key].get('issue') == issue for key, issue in OWNERS.items()),
            'completion owner issue changed')
    expected_bindings, covered_issues = set(), set()
    receipts = {f'model/{name}/results.json' for name in COMPONENTS}
    for row in rows:
        name = row['id']
        attachment_id, owners = COMPONENTS[name]
        attachment = attachments[attachment_id]
        require(attachment.get('classification') == 'realized_on_claimed_domain'
                and attachment.get('physical_identification') is False
                and attachment.get('new_lean_inhabitant') is False, 'unsupported attachment status')
        require(attachment.get('owners') == owners and attachment.get('reviewed_sources') == [name],
                'attachment ownership changed')
        issues = [OWNERS[owner] for owner in owners]
        require(row.get('attachment') == attachment_id and row.get('issues') == issues,
                'component ownership mismatch')
        covered_issues.update(issues)
        require(type(row['interface']) is dict and set(row['interface']) == INTERFACE_FIELDS
                and all(type(value) is str and value for value in row['interface'].values()),
                'missing component interface')
        require(type(row['remaining']) is list and row['remaining']
                and all(type(value) is str and value for value in row['remaining']),
                'missing component limits')
        for key, source in [('result', 'conclusion'), ('interface', 'interface'), ('remaining', 'remaining')]:
            require(row[key] == attachment.get(source), 'result or domain differs from reviewed attachment')
        witnesses, controls = attachment.get('witnesses'), attachment.get('countercontrols')
        require(type(witnesses) is list and witnesses and type(controls) is list and controls
                and all(type(path) is str for path in witnesses+controls), 'missing reviewed evidence')
        expected_bindings.update(witnesses+controls)
        require(set(witnesses+controls) <= set(bindings), 'missing reviewed source or control')
        finite_bindings = finite.get('source_bindings', {})
        require(all(finite_bindings.get(path) == bindings[path] for path in witnesses+controls),
                'finite attachment source custody differs')
        receipt_path = f'model/{name}/results.json'
        require(receipt_path in bindings and f'model/{name}/PROOF.md' in witnesses,
                'missing owning proof or receipt')
        receipt = strict_json.load(root/receipt_path)
        require(type(receipt) is dict and receipt.get('schema') == f'oph.{name}.v1',
                'wrong component receipt schema')
        receipt_scope = receipt.get('scope')
        require(type(receipt_scope) is dict and receipt_scope.get('integrated_gate_accepted') is False
                and receipt.get('canonical_physical_premises_discharged') == []
                and all(receipt_scope.get(key) is False for key in FORBIDDEN_PROMOTIONS[name]),
                'component receipt scope promoted')
        if 'physical_comparison_count_added' in receipt:
            require(type(receipt['physical_comparison_count_added']) is int
                    and receipt['physical_comparison_count_added'] == 0,
                    'component receipt added a physical comparison')
        dependency_bindings = receipt.get('source_bindings')
        require(type(dependency_bindings) is dict and dependency_bindings, 'missing component dependencies')
        for path, digest in dependency_bindings.items():
            bound(path, digest, 'compact component dependency')
        for owner in owners:
            progress = [item for item in stages[owner].get('progress_evidence', [])
                        if item.get('campaign_id') == CAMPAIGN and item.get('id') == name]
            require(len(progress) == 1 and progress[0].get('path') == receipt_path
                    and progress[0].get('sha256') == bindings[receipt_path]
                    and progress[0].get('independent_check') in controls
                    and progress[0].get('scope') == row['result']+' '+row['remaining'][0],
                    'unbound component progress')
    require(set(bindings) == expected_bindings, 'campaign source set differs from reviewed evidence')
    require(covered_issues == set(ISSUES), 'component coverage differs from selected issues')
    for stage in program['stages']:
        require(not any(item.get('campaign_id') == CAMPAIGN or item.get('path') in receipts
                        or item.get('id') in COMPONENTS for item in stage.get('accepted_evidence', [])),
                'finite component entered physical acceptance')
        for item in stage.get('progress_evidence', []):
            if item.get('campaign_id') == CAMPAIGN:
                require(item.get('id') in COMPONENTS
                        and stage.get('id') in COMPONENTS[item['id']][1],
                        'component progress assigned to an unrelated owner')
    campaign_rows = _unique_rows(program.get('component_campaigns'), [CAMPAIGN], 'program campaign')
    registered = campaign_rows[CAMPAIGN]
    require(registered.get('scope') == 'FINITE_COMPONENT_FAMILY'
            and registered.get('whole_program_frozen') is False
            and registered.get('supports') == list(OWNERS), 'program campaign scope or coverage changed')
    contract = registered.get('contract')
    require(type(contract) is dict and set(contract) == {'path', 'sha256', 'independent_check', 'scope'}
            and contract.get('path') == str(PATH)
            and contract.get('independent_check') == 'model/completion/test_compact_source_campaign.py'
            and type(contract.get('scope')) is str and contract['scope'], 'invalid program campaign contract')
    _bound_file(root, contract['path'], contract['sha256'], 'compact campaign contract')
    return {'selected_issues': ISSUES, 'reviewed_components': len(rows),
            'canonical_premises_discharged': 0, 'physical_particle_map_complete': False,
            'interacting_continuum_and_Einstein_response_established': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.parse_args()
    print(json.dumps(validate(strict_json.load(ROOT/PATH)), indent=2))


if __name__ == '__main__':
    main()
