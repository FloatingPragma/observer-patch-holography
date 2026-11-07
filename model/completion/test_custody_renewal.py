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


@pytest.mark.parametrize('damage, message', [
    ('count', 'classification counts differ'),
    ('history', 'immutable historical restoration changed'),
    ('physical_acceptance', 'cannot promote scientific acceptance'),
    ('binding', 'current custody binding changed'),
    ('coverage', 'current custody binding coverage changed'),
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
        renewal.RECORD, renewal.REVIEW, renewal.CLASS, renewal.PROGRAM}
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
    for path, target in renewal.PROGRAM_POINTERS:
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
            allowed = {path for path, _ in renewal.PROGRAM_POINTERS}
            candidates = []
            def walk(value, path=()):
                if isinstance(value, dict):
                    if 'path' in value and 'sha256' in value and path not in allowed:
                        candidates.append(value)
                    for key, child in value.items():
                        walk(child, path + (key,))
                elif isinstance(value, list):
                    for index, child in enumerate(value):
                        walk(child, path + (index,))
            walk(program)
            assert candidates
            candidates[0]['sha256'] = 'f' * 64
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
