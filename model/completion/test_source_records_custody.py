"""The current field successor cannot alter retained science or review scope."""
from copy import deepcopy
import pytest
from tools import strict_json
from . import source_records_custody as custody


def documents():
    return {name: strict_json.load(custody.ROOT / name) for name in custody.FILES}


def test_current_successor_keeps_science_status_and_old_review_attribution():
    report = custody.validate()
    assert report['executed_live_receipts'] == 27
    assert report['program_pointer_schema']['locations'] == 813
    assert report['program_pointer_schema']['targets'] == 342
    assert report['navigation']['claim_count'] == 392
    assert report['navigation']['lean_module_count'] == 681
    assert report['scientific_statuses_changed'] is False
    assert report['old_review_signed_current_bytes'] is False
    assert report['historical_evidence_rewritten'] is False


@pytest.mark.parametrize('damage', ['gate_status', 'acceptance', 'pointer_target', 'extra_pointer',
                                    'class_parameter', 'class_clock', 'review_claim',
                                    'map_classification', 'source_append', 'source_delete'])
def test_scientific_contract_or_source_scope_cannot_change(damage):
    values = documents()
    if damage == 'gate_status':
        stage = next(row for row in values[custody.PROGRAM]['stages'] if row['status'] != 'ACCEPTED')
        stage['status'] = 'ACCEPTED'
    elif damage == 'acceptance':
        values[custody.PROGRAM]['stages'][0]['acceptance'] = 'No conditions.'
    elif damage in ('pointer_target', 'extra_pointer'):
        path, _ = custody.pointers(values[custody.PROGRAM])[0]
        row = custody.live.at_path(values[custody.PROGRAM], path)
        if damage == 'pointer_target':
            row['path'] = custody.CLASS
        else:
            values[custody.PROGRAM]['extra_evidence'] = deepcopy(row)
    elif damage == 'class_parameter':
        values[custody.CLASS]['laws']['local_alternative_family']['port_weights']['max'] = 64
    elif damage == 'class_clock':
        values[custody.CLASS]['resource_envelope']['host_limits_are_physical_clock'] = True
    elif damage == 'review_claim':
        values[custody.REVIEW]['verdict'] = 'PHYSICS_ESTABLISHED'
    elif damage == 'map_classification':
        values[custody.MAP]['attachments'][0]['classification'] = 'new_physical_theorem'
    elif damage == 'source_append':
        values[custody.INTERFACES]['source_bindings']['unreviewed/source.py'] = 'f' * 64
    else:
        values[custody.INTERFACES]['source_bindings'].pop(next(iter(values[custody.INTERFACES]['source_bindings'])))
    with pytest.raises(ValueError, match='scientific contract or binding scope'):
        custody.conserve_documents(values, custody.load_authority())


@pytest.mark.parametrize('key,value', [('new_mathematical_acceptance', True), ('gate_statuses_changed', True),
                                      ('canonical_physical_premises_discharged', ['P1']),
                                      ('historical_evidence_rewritten', True),
                                      ('old_review_signed_current_bytes', True),
                                      ('current_bindings', {}), ('classification_counts', {'realized': 156})])
def test_record_cannot_promote_science_or_relabel_prior_review(key, value):
    packet = strict_json.load(custody.ROOT / custody.RECORD)
    packet[key] = value
    with pytest.raises(ValueError, match='record changed'):
        custody.validate(packet)


def test_nonrenewable_source_hash_is_not_silently_repaired():
    values = documents()
    allowed = set(custody.load_authority()['renewable_targets'])
    name = next(name for name in values[custody.INTERFACES]['source_bindings'] if name not in allowed)
    values[custody.INTERFACES]['source_bindings'][name] = 'f' * 64
    with pytest.raises(ValueError, match='unreviewed source'):
        custody.check_binding_sources(values, custody.load_authority(), require_current=False)


def test_interpreter_change_is_explicit_and_old_signed_digest_stays_retained():
    values = documents()
    row = values[custody.REVIEW][custody.SUCCESSOR_KEY]
    assert row == custody.transition(values)
    assert row['old_review_signed_current_bytes'] is False
    old = strict_json.load(custody.ROOT / row['retained_review']['path'])
    assert old['reviewed_contract_digest'] == row['retained_reviewed_contract_digest']
    assert row['retained_interpreter_receipt_sha256'] != row['current_interpreter_receipt_sha256']
    assert row['retained_reviewed_contract_digest'] != row['current_reviewed_contract_digest']


def test_old_checker_and_authority_remain_exact_historical_bytes():
    authority = custody.load_authority()
    for name in ('model/completion/field_custody.py',
                 'model/completion/test_field_custody.py',
                 'model/completion/field_custody_authority.json',
                 'model/completion/field_custody_current.json'):
        assert custody.sha(custody.ROOT, name) == authority['historical_sha256'][name]
