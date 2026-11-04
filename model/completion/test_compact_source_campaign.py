"""Guard custody, exact component interfaces and finite-only acceptance."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import shutil

import pytest

from tools import strict_json
from .compact_source_campaign import CAMPAIGN, COMPONENTS, ISSUES, PATH, ROOT, validate


def packet():
    return strict_json.load(ROOT/PATH)


def write(path, value):
    path.write_text(json.dumps(value, sort_keys=True, indent=2)+'\n')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture
def repository(tmp_path):
    value = packet()
    paths = set(value['source_bindings'])
    for name in COMPONENTS:
        paths.update(strict_json.load(ROOT/f'model/{name}/results.json')['source_bindings'])
    paths.update((str(PATH), 'model/completion/finite_attachment_map.json', 'claims/completion_program.json'))
    for name in paths:
        target = tmp_path/name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT/name, target)
    return tmp_path


def test_four_exact_components_cover_five_issues_without_physical_promotion():
    result = validate(packet())
    assert result['selected_issues'] == ISSUES
    assert result['reviewed_components'] == 4
    assert result['canonical_premises_discharged'] == 0
    assert result['physical_particle_map_complete'] is False
    assert result['interacting_continuum_and_Einstein_response_established'] is False


@pytest.mark.parametrize('damage', [
    'stale_source', 'missing_receipt', 'missing_control', 'reassigned_issue',
    'altered_result', 'altered_clock', 'erased_remaining', 'missing_component',
    'duplicated_component', 'physical_acceptance', 'canonical_discharge',
    'frozen_whole_program', 'added_comparison', 'boolean_comparison', 'unknown_physical_claim',
    'external_source',
])
def test_campaign_custody_and_scientific_boundary_mutations_fail(damage):
    value = deepcopy(packet())
    if damage == 'stale_source':
        value['source_bindings'][next(iter(value['source_bindings']))] = '0'*64
    elif damage == 'missing_receipt':
        del value['source_bindings']['model/native_quantum_gauss/results.json']
    elif damage == 'missing_control':
        del value['source_bindings']['model/compact_interacting_response/test_packet.py']
    elif damage == 'reassigned_issue':
        value['component_results'][0]['issues'] = [1147]
    elif damage == 'altered_result':
        value['component_results'][1]['result'] = 'The complete measured particle spectrum is derived.'
    elif damage == 'altered_clock':
        value['component_results'][2]['interface']['clock'] = 'The physical metric and calibrated clock.'
    elif damage == 'erased_remaining':
        value['component_results'][3]['remaining'] = []
    elif damage == 'missing_component':
        value['component_results'].pop()
    elif damage == 'duplicated_component':
        value['component_results'].append(deepcopy(value['component_results'][0]))
    elif damage == 'physical_acceptance':
        value['integrated_gate_accepted'] = True
    elif damage == 'canonical_discharge':
        value['canonical_premises_discharged'] = ['PR-23']
    elif damage == 'frozen_whole_program':
        value['whole_program_frozen'] = True
    elif damage == 'added_comparison':
        value['physical_comparison_count_added'] = 1
    elif damage == 'boolean_comparison':
        value['physical_comparison_count_added'] = False
    elif damage == 'unknown_physical_claim':
        value['physical_particle_map_complete'] = True
    elif damage == 'external_source':
        value['source_bindings']['../outside-repository'] = '0'*64
    with pytest.raises(ValueError):
        validate(value)


@pytest.mark.parametrize('damage', [
    'attachment_physical', 'attachment_lean', 'attachment_duplicate', 'attachment_owner',
    'progress_digest', 'progress_duplicate', 'progress_reassigned', 'unlabelled_acceptance',
    'contract_digest', 'contract_coverage',
])
def test_reviewed_ledger_and_program_mutations_fail(repository, damage):
    finite_path = repository/'model/completion/finite_attachment_map.json'
    program_path = repository/'claims/completion_program.json'
    finite, program = strict_json.load(finite_path), strict_json.load(program_path)
    attachment = next(row for row in finite['attachments'] if row['id'] == 'FA-146')
    owner = next(row for row in program['stages'] if row['id'] == 'CP-03')
    progress = next(row for row in owner['progress_evidence'] if row.get('campaign_id') == CAMPAIGN)
    registered = next(row for row in program['component_campaigns'] if row['id'] == CAMPAIGN)
    if damage == 'attachment_physical':
        attachment['physical_identification'] = True
    elif damage == 'attachment_lean':
        attachment['new_lean_inhabitant'] = True
    elif damage == 'attachment_duplicate':
        finite['attachments'].append(deepcopy(attachment))
    elif damage == 'attachment_owner':
        attachment['owners'] = ['CP-07']
    elif damage == 'progress_digest':
        progress['sha256'] = '0'*64
    elif damage == 'progress_duplicate':
        owner['progress_evidence'].append(deepcopy(progress))
    elif damage == 'progress_reassigned':
        wrong_owner = next(row for row in program['stages'] if row['id'] == 'CP-07')
        wrong_owner.setdefault('progress_evidence', []).append(deepcopy(progress))
    elif damage == 'unlabelled_acceptance':
        # Omitting the campaign label cannot turn the same finite receipt into
        # accepted physical evidence in an unrelated stage.
        wrong_owner = next(row for row in program['stages'] if row['id'] == 'CP-07')
        wrong_owner.setdefault('accepted_evidence', []).append({'path': progress['path']})
    elif damage == 'contract_digest':
        registered['contract']['sha256'] = '0'*64
    elif damage == 'contract_coverage':
        registered['supports'].remove('CP-09')
    write(finite_path, finite)
    write(program_path, program)
    with pytest.raises(ValueError):
        validate(strict_json.load(repository/PATH), repository)


def test_transitive_source_edit_invalidates_an_unchanged_component_receipt(repository):
    value = strict_json.load(repository/PATH)
    receipt = strict_json.load(repository/'model/compact_interacting_response/results.json')
    dependency = next(path for path in receipt['source_bindings'] if path not in value['source_bindings'])
    path = repository/dependency
    path.write_bytes(path.read_bytes()+b'\n')
    with pytest.raises(ValueError, match='stale compact component dependency'):
        validate(value, repository)


def test_rehashing_a_promoted_receipt_does_not_create_a_physical_join(repository):
    receipt_name = 'model/native_quantum_gauss/results.json'
    receipt_path = repository/receipt_name
    receipt = strict_json.load(receipt_path)
    receipt['scope']['same_compact_defect_band_gauge_dressing'] = True
    write(receipt_path, receipt)
    changed_digest = digest(receipt_path)
    value = strict_json.load(repository/PATH)
    value['source_bindings'][receipt_name] = changed_digest
    write(repository/PATH, value)
    finite_path = repository/'model/completion/finite_attachment_map.json'
    finite = strict_json.load(finite_path)
    finite['source_bindings'][receipt_name] = changed_digest
    write(finite_path, finite)
    program_path = repository/'claims/completion_program.json'
    program = strict_json.load(program_path)
    for stage in program['stages']:
        for row in stage.get('progress_evidence', []):
            if row.get('campaign_id') == CAMPAIGN and row.get('path') == receipt_name:
                row['sha256'] = changed_digest
    registered = next(row for row in program['component_campaigns'] if row['id'] == CAMPAIGN)
    registered['contract']['sha256'] = digest(repository/PATH)
    write(program_path, program)
    with pytest.raises(ValueError, match='component receipt scope promoted'):
        validate(value, repository)
