"""Reject source substitution and physical promotion in the finite join."""
from copy import deepcopy
import json
import shutil

import pytest

from tools import strict_json
from .finite_observable_campaign import CAMPAIGN, COMPONENTS, ISSUES, PATH, ROOT, validate


def packet():
    return strict_json.load(ROOT/PATH)


def test_three_mathematical_results_cover_five_issues():
    result = validate(packet())
    assert result == {'selected_issues': ISSUES, 'finite_components': 3,
                      'canonical_premises_discharged': 0, 'integrated_gate_accepted': False}


@pytest.mark.parametrize('damage', [
    'custody', 'missing', 'duplicate', 'clock', 'boundary', 'owner',
    'acceptance', 'premise', 'comparison_boolean', 'specification', 'unknown_claim'])
def test_altered_contract_is_rejected(damage):
    value = deepcopy(packet())
    if damage == 'custody':
        value['source_bindings'][next(iter(value['source_bindings']))] = '0'*64
    elif damage == 'missing':
        value['component_results'].pop()
    elif damage == 'duplicate':
        value['component_results'].append(deepcopy(value['component_results'][0]))
    elif damage == 'clock':
        value['component_results'][0]['interface']['clock'] = 'SI calibrated rest masses'
    elif damage == 'boundary':
        value['component_results'][0]['remaining'] = []
    elif damage == 'owner':
        value['component_results'][1]['issues'] = [1147]
    elif damage == 'acceptance':
        value['integrated_gate_accepted'] = True
    elif damage == 'premise':
        value['canonical_premises_discharged'] = ['PR-23']
    elif damage == 'comparison_boolean':
        value['physical_comparison_count_added'] = False
    elif damage == 'specification':
        value['specification_readout_allowed'] = False
    else:
        value['physical_particles_discovered'] = True
    with pytest.raises(ValueError):
        validate(value)


@pytest.fixture
def repository(tmp_path):
    value = packet()
    paths = set(value['source_bindings'])
    for name in COMPONENTS:
        receipt = strict_json.load(ROOT/f'model/{name}/results.json')
        paths.update(receipt.get('source_bindings', receipt.get('source_sha256', {})))
    paths.update([str(PATH), 'claims/completion_program.json',
                  'model/completion/finite_attachment_map.json'])
    for name in paths:
        destination = tmp_path/name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT/name, destination)
    return tmp_path


@pytest.mark.parametrize('damage', ['promote', 'duplicate', 'receipt', 'wrong_stage', 'dependency'])
def test_independently_changed_ledger_or_source_is_rejected(repository, damage):
    path = repository/'claims/completion_program.json'
    program = strict_json.load(path)
    stage = next(s for s in program['stages'] if s['id'] == 'CP-03')
    progress = next(r for r in stage['progress_evidence'] if r.get('campaign_id') == CAMPAIGN)
    if damage == 'promote':
        stage['accepted_evidence'].append(deepcopy(progress))
    elif damage == 'duplicate':
        stage['progress_evidence'].append(deepcopy(progress))
    elif damage == 'receipt':
        progress['sha256'] = '0'*64
    elif damage == 'wrong_stage':
        other = next(s for s in program['stages'] if s['id'] == 'CP-10')
        other.setdefault('progress_evidence', []).append(deepcopy(progress))
    else:
        receipt = strict_json.load(repository/'model/charged_composite_search/results.json')
        dependency = next(iter(receipt.get('source_bindings', receipt.get('source_sha256', {}))))
        with (repository/dependency).open('a') as handle:
            handle.write('\n# altered dependency\n')
    path.write_text(json.dumps(program))
    with pytest.raises(ValueError):
        validate(packet(), repository)
