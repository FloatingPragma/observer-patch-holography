"""Only the named copied child and its four actual source scopes are admitted."""
import copy

import pytest

from tools import strict_json
from tools import renew_canonical_live_receipts as live
from . import model


def candidate():
    approved = model.authority()
    old = strict_json.load(model.ROOT / approved['retained_receipts'][model.COMPLETION]['retained_path'])
    value = copy.deepcopy(old)
    spec = approved['exact_added_subtree']
    parent = live.at_path(value, spec['target_json_path'][:-1])
    parent[spec['target_json_path'][-1]] = live.at_path(
        strict_json.load(model.ROOT / spec['source_receipt']), spec['source_json_path'])
    maps = {where: pins for _, pins, where in live.source_maps(value)}
    for where, additions in approved['exact_source_additions_by_location'].items():
        maps[where].update(additions)
    return old, value, approved


def test_exact_four_locations_and_one_child_keep_all_older_values():
    old, value, approved = candidate()
    assert len(approved['exact_source_additions_by_location']) == 4
    assert sum(map(len, approved['exact_source_additions_by_location'].values())) == 16
    assert model.compare_completion(old, value, approved) == value
    restored = copy.deepcopy(value)
    spec = approved['exact_added_subtree']
    del live.at_path(restored, spec['target_json_path'][:-1])[spec['target_json_path'][-1]]
    maps = {where: pins for _, pins, where in live.source_maps(restored)}
    for where, additions in approved['exact_source_additions_by_location'].items():
        for name in additions:
            del maps[where][name]
    assert restored == old


@pytest.mark.parametrize('mutation', ['old_science', 'child_science', 'physical_scope',
                                     'extra_map', 'extra_path', 'wrong_location'])
def test_unreviewed_science_or_scope_cannot_be_admitted(mutation):
    old, value, approved = candidate()
    if mutation == 'old_science':
        value['unreviewed_scientific_clause'] = True
    elif mutation in {'child_science', 'physical_scope'}:
        spec = approved['exact_added_subtree']
        child = live.at_path(value, spec['target_json_path'])
        if mutation == 'child_science':
            child['Gram_calibration_events'] += 1
        else:
            child['scope']['physical_clock_or_particle_energy_identified'] = True
    elif mutation == 'extra_map':
        value['unreviewed_source'] = {'source_sha256': {}}
    elif mutation == 'extra_path':
        value['source_sha256']['model/unreviewed.py'] = '0' * 64
    else:
        maps = {where: pins for _, pins, where in live.source_maps(value)}
        where, additions = next(iter(approved['exact_source_additions_by_location'].items()))
        name = next(iter(additions))
        del maps[where][name]
        unrelated = next(p for w, p in maps.items()
                         if w not in approved['exact_source_additions_by_location'] and name not in p)
        unrelated[name] = additions[name]
    with pytest.raises(ValueError):
        model.compare_completion(old, value, approved)


def test_existing_source_hashes_remain_checked():
    value = {'source_sha256': {model.CLASS: '0' * 64}}
    with pytest.raises(ValueError, match='live source changed'):
        model.check_bindings(model.COMPLETION, value)


def test_nested_additions_match_their_executed_owning_receipts():
    old, _, approved = candidate()
    owners = {
        'native_clock_comparison': 'model/native_clock_comparison/results.json',
        'native_controlled_transfer_experiment': 'model/native_transfer_control/results.json',
        'weighted_native_complete_record_response': 'model/native_record_response/results.json',
    }
    for section, receipt in owners.items():
        fresh = strict_json.load(model.ROOT / receipt)
        if section == 'weighted_native_complete_record_response':
            del fresh['short_record_population_realization']
        additions = next(v for k, v in approved['exact_source_additions_by_location'].items()
                         if k.startswith('root.' + section + '.'))
        carried = live.carry_metadata(receipt, old[section], fresh, scope_append=additions)
        assert live.canonical(carried) == live.canonical(fresh)


def test_final_current_record_and_exact_scope():
    value = model.check()
    assert value['old_scientific_fields_and_statuses_exactly_preserved'] is True
    assert value['complete_scientific_payload_unchanged'] is False
    assert value['added_child_result_is_metadata_only'] is False
    assert value['class_modified'] is False
    assert value['physical_statuses_changed'] is False
    assert value['class_does_not_bind_consumers_back'] is True


@pytest.mark.parametrize('field', ['class_sha256', 'producer_source_sha256', 'actual_unchanged_owner_executed'])
def test_execution_record_must_bind_the_actual_class_and_owner(field):
    approved = model.authority()
    execution = strict_json.load(model.ROOT / model.EXECUTIONS)
    execution[model.COMPLETION][field] = False if field.startswith('actual_') else '0' * 64
    with pytest.raises(ValueError, match='actual late-owner execution'):
        model.validate_execution(execution, {p: model.sha(model.ROOT / p)
                                            for p in (model.COMPLETION, model.REGIONAL)}, approved)
