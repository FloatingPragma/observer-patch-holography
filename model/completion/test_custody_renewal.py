"""Protect the historical/current distinction in the dated custody record."""
import copy
import json
import shutil
import pytest

from model.completion import custody_renewal as renewal
from tools import strict_json


def test_current_custody_preserves_scientific_status():
    result = renewal.validate()
    assert result['classification_counts'] == {'domains': 156, 'realized': 151}
    assert result['new_mathematical_acceptance'] is False


def test_integration_authority_binds_every_existing_live_pointer():
    program = strict_json.load(renewal.ROOT / renewal.PROGRAM)
    pointers = renewal.program_pointers(program)
    assert len(pointers) == 736
    assert len({target for _, target in pointers}) == 267
    assert renewal.INTEGRATION_BASELINE['commit'] == '85a0cc208f52b385cfa674b288bcb079156abfec'
    assert renewal.INTEGRATION_CLASS_DIGEST == '54fe2116bef3a7cb2eab1678ff20281d80f01032d21d5024009424fe5ff8d877'
    assert renewal.SOURCE_REFRESH_BASELINE['commit'] == '12047f556a3f4d6bae5b48aa5b85980a0c5eda40'
    assert renewal.CLASS_DIGEST == 'b6bcea47eeef27e3bd09102d0ea94a9f3222690f76face79d833ff56cf29c77b'


@pytest.mark.parametrize('damage', ('add', 'retarget', 'move'))
def test_live_pointer_schema_cannot_expand_or_retarget(damage):
    program = strict_json.load(renewal.ROOT / renewal.PROGRAM)
    path, _ = renewal.program_pointers(program)[0]
    if damage == 'add':
        program['extra_evidence'] = copy.deepcopy(renewal.at_path(program, path))
    elif damage == 'retarget':
        renewal.at_path(program, path)['path'] = renewal.MAP
    else:
        program['moved_evidence'] = renewal.at_path(program, path[:-1]).pop(path[-1])
    with pytest.raises(ValueError, match='immutable program pointer schema'):
        renewal.program_semantics(program, renewal.HISTORICAL_SCOPE)


@pytest.mark.parametrize('damage, message', [
    ('count', 'classification counts differ'),
    ('history', 'immutable historical restoration changed'),
    ('physical_acceptance', 'cannot promote scientific acceptance'),
    ('binding', 'current custody binding changed'),
    ('coverage', 'current custody binding coverage changed'),
    ('integration_baseline', 'immutable integration baseline changed'),
    ('integration_source', 'immutable integration baseline changed'),
    ('pointer_authority', 'immutable program pointer authority changed'),
    ('original_authority', 'immutable original conservation authority changed'),
    ('refresh_baseline', 'immutable accepted-main refresh authority changed'),
    ('refresh_source', 'immutable accepted-main refresh authority changed'),
])
def test_new_renewal_cannot_relabel_history_or_scientific_acceptance(damage, message):
    packet = copy.deepcopy(strict_json.load(renewal.ROOT / renewal.RECORD))
    if damage == 'count':
        packet['current_classification_counts']['domains'] = 145
    elif damage == 'history':
        packet['historical_snapshot_restorations'][0]['restored_sha256'] = '0' * 64
    elif damage == 'physical_acceptance':
        packet['new_mathematical_acceptance'] = True
    elif damage == 'coverage':
        del packet['current_bindings'][renewal.INTERFACES]
    elif damage == 'integration_baseline':
        packet['integration_baseline']['commit'] = '0' * 40
    elif damage == 'integration_source':
        packet['integration_baseline']['source_sha256'][renewal.REVIEW] = '0' * 64
    elif damage == 'pointer_authority':
        packet['program_pointer_schema']['locations'] = 69
    elif damage == 'original_authority':
        packet['original_conserved']['review_history'] = '0' * 64
    elif damage == 'refresh_baseline':
        packet['source_refresh_baseline']['commit'] = '0' * 40
    elif damage == 'refresh_source':
        packet['source_refresh_baseline']['source_sha256'][renewal.CLASS] = '0' * 64
    else:
        first = next(iter(packet['current_bindings']))
        packet['current_bindings'][first] = '0' * 64
    with pytest.raises(ValueError, match=message):
        renewal.validate(packet)


def _read(root, name):
    return strict_json.load(root / name)


def _write(root, name, value):
    (root / name).write_text(json.dumps(value, indent=2, sort_keys=True) + '\n',
                            encoding='utf-8', newline='\n')


@pytest.fixture
def isolated_custody_root(tmp_path):
    names = renewal.REQUIRED_BINDINGS | {
        renewal.RECORD, renewal.REVIEW, renewal.CLASS, renewal.PROGRAM, renewal.NAVIGATION_SOURCE}
    names |= {target for _, target in renewal.program_pointers(
        strict_json.load(renewal.ROOT / renewal.PROGRAM))}
    for name in names:
        target = tmp_path / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(renewal.ROOT / name, target)
    renewal.validate(root=tmp_path)
    return tmp_path


def _rebind_outer_hashes(root):
    """An attacker may re-sign data pointers, but not trusted code anchors."""
    packet = _read(root, renewal.RECORD)
    packet['current_bindings'] = {
        name: renewal.file_digest(root / name) for name in renewal.REQUIRED_BINDINGS}
    _write(root, renewal.RECORD, packet)
    review = _read(root, renewal.REVIEW)
    for name in (renewal.MAP, renewal.INTERFACES, renewal.RECORD):
        review['source_bindings'][name] = renewal.file_digest(root / name)
    review['current_source_custody_renewal']['sha256'] = renewal.file_digest(root / renewal.RECORD)
    _write(root, renewal.REVIEW, review)
    klass = _read(root, renewal.CLASS)
    klass['adoption']['independent_review']['sha256'] = renewal.file_digest(root / renewal.REVIEW)
    _write(root, renewal.CLASS, klass)
    program = _read(root, renewal.PROGRAM)
    for path, target in renewal.program_pointers(program):
        renewal.at_path(program, path)['sha256'] = renewal.file_digest(root / target)
    _write(root, renewal.PROGRAM, program)


@pytest.mark.parametrize('damage', [
    'historical_snapshot', 'other_history', 'canonical_owner', 'map_cleanup',
    'program_promotion', 'adoption_proposal', 'adoption_registry', 'unapproved_hash',
    'binding_omission',
])
def test_coordinated_data_rewrite_cannot_redefine_preservation(isolated_custody_root, damage):
    root = isolated_custody_root
    packet = _read(root, renewal.RECORD)
    if damage in ('historical_snapshot', 'other_history'):
        review = _read(root, renewal.REVIEW)
        if damage == 'historical_snapshot':
            row = packet['historical_snapshot_restorations'][0]
            review['review_refreshes'][row['refresh_index']]['source_receipt_snapshots'][
                row['receipt']] = row['restored_sha256'] = 'f' * 64
        else:
            review['review_refreshes'][0]['scope'] += ' forged extension'
        packet['conserved']['review_history'] = renewal.digest(renewal.review_history(review))
        _write(root, renewal.REVIEW, review)
    elif damage in ('canonical_owner', 'map_cleanup', 'binding_omission'):
        finite = _read(root, renewal.MAP)
        if damage == 'canonical_owner':
            finite['canonical_ownership'][0]['owners'] = ['CP-02']
        elif damage == 'map_cleanup':
            finite['cleanup'] = {'forged': 'outside the old four selected projections'}
        else:
            del finite['source_bindings'][next(iter(finite['source_bindings']))]
        packet['conserved']['finite_map'] = renewal.digest(renewal.without_source_bindings(finite))
        _write(root, renewal.MAP, finite)
    elif damage in ('program_promotion', 'unapproved_hash'):
        program = _read(root, renewal.PROGRAM)
        if damage == 'program_promotion':
            next(row for row in program['stages'] if row['id'] == 'CP-02')['status'] = 'ACCEPTED'
        else:
            program['unapproved_receipt_sha256'] = 'f' * 64
        packet['conserved']['completion_program'] = renewal.digest(
            renewal.program_semantics(program, packet['historical_scope_text']))
        _write(root, renewal.PROGRAM, program)
    else:
        klass = _read(root, renewal.CLASS)
        if damage == 'adoption_proposal':
            klass['adoption']['status'] = 'PROPOSED_FOR_INDEPENDENT_REVIEW'
        else:
            klass['adoption']['changes_current_registry'] = False
        _write(root, renewal.CLASS, klass)
    _write(root, renewal.RECORD, packet)
    _rebind_outer_hashes(root)
    with pytest.raises(ValueError, match='immutable|adoption|source binding path coverage'):
        renewal.validate(root=root)


@pytest.mark.parametrize('index', (0, 368, 735))
def test_live_pointer_hashes_must_match_target_bytes(isolated_custody_root, index):
    root = isolated_custody_root
    program = _read(root, renewal.PROGRAM)
    path, _ = renewal.program_pointers(program)[index]
    renewal.at_path(program, path)['sha256'] = 'f' * 64
    _write(root, renewal.PROGRAM, program)
    with pytest.raises(ValueError, match='program custody target changed'):
        renewal.validate(root=root)


@pytest.mark.parametrize('damage, message', [
    ('class_clause', 'adopted scientific class clauses changed'),
    ('class_other_source', 'adopted scientific class clauses changed'),
    ('class_live_source', 'adopted scientific class clauses changed'),
    ('navigation_count', 'current finite-map navigation changed'),
    ('navigation_hash', 'current finite-map navigation changed'),
    ('navigation_owner', 'immutable current navigation authority changed'),
    ('review_class', 'current review class digest changed'),
])
def test_named_metadata_refresh_cannot_expand_its_scope(isolated_custody_root, damage, message):
    root = isolated_custody_root
    if damage.startswith('class_'):
        klass = _read(root, renewal.CLASS)
        if damage == 'class_clause':
            klass['interpreter_binding']['fresh_trillion_pulse_replay_claimed'] = True
        elif damage == 'class_other_source':
            klass['interpreter_binding']['source_hashes']['model/native_program/interpreter.py'] = 'f'*64
        else:
            klass['interpreter_binding']['source_hashes'][renewal.CLASS_LIVE_RECEIPT] = 'f'*64
        _write(root, renewal.CLASS, klass)
    elif damage.startswith('navigation_'):
        if damage == 'navigation_owner':
            register = _read(root, renewal.NAVIGATION_SOURCE)
            register['model_supply']['contract_coverage']['claim_ids'].pop()
            _write(root, renewal.NAVIGATION_SOURCE, register)
        else:
            finite = _read(root, renewal.MAP)
            navigation = finite['existing_navigation_coverage']
            if damage == 'navigation_count': navigation['claim_count'] = 389
            else: navigation['contract_coverage_sha256'] = 'f'*64
            _write(root, renewal.MAP, finite)
    else:
        review = _read(root, renewal.REVIEW)
        review['reviewed_contract_digest'] = renewal.INTEGRATION_CLASS_DIGEST
        _write(root, renewal.REVIEW, review)
    _rebind_outer_hashes(root)
    with pytest.raises(ValueError, match=message):
        renewal.validate(root=root)
