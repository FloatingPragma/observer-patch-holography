"""The new service inventory has its own source-bound, limited contract."""
import json

import pytest

from . import joins_contract as j


def test_successor_preserves_the_original_diagnostic_and_names_new_services():
    value = j.load_contract()
    assert value['predecessor_campaign'] == j.load_diagnostic()['campaign_id']
    assert value['service_inventory_extended'] is True
    assert not value['native_family_changed']
    assert not value['integrated_gate_accepted']


@pytest.mark.parametrize('edit', ('omit_source', 'boolean_count', 'hidden_services', 'promote'))
def test_changed_join_scope_or_source_is_rejected(tmp_path, monkeypatch, edit):
    value = j.load_contract()
    if edit == 'omit_source':
        value['source_sha256'].pop(next(iter(value['source_sha256'])))
    elif edit == 'wrong_hash':
        value['source_sha256'][next(iter(value['source_sha256']))] = '0'*64
    elif edit == 'boolean_count':
        value['physical_comparison_count_added'] = False
    elif edit == 'hidden_services':
        value['service_inventory_extended'] = False
    else:
        value['integrated_gate_accepted'] = True
    (tmp_path/'joins_campaign.json').write_text(json.dumps(value))
    monkeypatch.setattr(j, 'HERE', tmp_path)
    with pytest.raises(ValueError):
        j.load_contract()


def test_duplicate_contract_keys_are_rejected(tmp_path, monkeypatch):
    original = (j.HERE/'joins_campaign.json').read_text()
    changed = original.replace('"service_inventory_extended": true',
                               '"service_inventory_extended": false, "service_inventory_extended": true')
    assert changed != original
    (tmp_path/'joins_campaign.json').write_text(changed)
    monkeypatch.setattr(j, 'HERE', tmp_path)
    with pytest.raises(ValueError):
        j.load_contract()
