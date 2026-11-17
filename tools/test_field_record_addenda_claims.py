"""The sealed supplement admits custody renewal, never re-signed new science."""
from copy import deepcopy
import json
from pathlib import Path
import shutil

import pytest

from tools import check_field_record_addenda_claims as claims


def load(root, name):
    return claims.strict_json.load(root/name)


def write(root, name, value):
    (root/name).write_text(json.dumps(value, indent=2, sort_keys=True)+'\n')


@pytest.fixture
def data():
    return load(claims.ROOT, 'claims/addenda/field_record_addenda.json')


@pytest.fixture
def repository(tmp_path):
    klass = load(claims.ROOT, claims.program_class_v2.PATH)
    review_path = klass['adoption']['independent_review']['path']
    review = load(claims.ROOT, review_path)
    catalog = load(claims.ROOT, claims.renewal.CATALOG)
    spec = catalog['owners']['native_program']
    member = load(claims.ROOT, spec['receipt'])
    paths = {str(claims.program_class_v2.PATH), review_path, claims.renewal.CATALOG,
             spec['receipt'], spec['retained_path'], 'claims/executable_model_hypothesis.json',
             'tracking/observation_ledger.json', klass['laws']['local_alternative_family']['contract']}
    paths.update(klass['interpreter_binding']['source_hashes'])
    paths.update(klass['review_implementation_bindings'])
    paths.update(claims.program_class_v2.EXPOSURE_SNAPSHOTS.get(p, p)
                 for p in klass['exposure']['source_hashes'])
    paths.update(review['source_bindings'])
    paths.update(catalog['reviewed_executable_source_sha256'])
    paths.update(member['source_sha256'])
    for name in paths:
        target = tmp_path/name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(claims.ROOT/name, target)
    try:
        yield tmp_path
    finally:
        shutil.rmtree(tmp_path)


def resign(root):
    """A coherent forgery updates both class and review, not merely one hash."""
    klass = load(root, claims.program_class_v2.PATH)
    member = claims.custody_renewal.CLASS_LIVE_RECEIPT
    klass['interpreter_binding']['source_hashes'][member] = claims.sha(root/member)
    review_path = klass['adoption']['independent_review']['path']
    review = load(root, review_path)
    review['source_bindings'][member] = claims.sha(root/member)
    digest = claims.program_class_v2.reviewed_contract_digest(klass)
    review['reviewed_contract_digest'] = digest
    write(root, review_path, review)
    klass['adoption']['independent_review'].update(
        sha256=claims.sha(root/review_path), reviewed_contract_digest=digest)
    write(root, claims.program_class_v2.PATH, klass)
    return klass


def test_current_class_is_equivalent_without_rewriting_the_seal(data):
    before = deepcopy(data)
    claims.check_parent_class(data)
    assert data == before
    assert data['sealed_program_class_sha256'] == claims.SEALED_CLASS_SHA256


@pytest.mark.parametrize('field', ('sealed_parent_commit', 'sealed_program_class_sha256'))
def test_replacing_the_sealed_authority_is_rejected(data, field):
    data[field] = '0'*len(data[field])
    with pytest.raises(ValueError, match='Sealed parent authority'):
        claims.check_parent_class(data)


@pytest.mark.parametrize('damage', ('science', 'scope'))
def test_coherently_resigned_member_cannot_change_science_or_scope(data, repository, damage):
    name = claims.custody_renewal.CLASS_LIVE_RECEIPT
    member = load(repository, name)
    if damage == 'science':
        member['scope']['physical_clock_calibration_claimed'] = True
    else:
        del member['source_sha256'][next(iter(member['source_sha256']))]
    write(repository, name, member)
    klass = resign(repository)
    # Coherent review signatures cannot authorize changed live member science.
    claims.validate_live_review_signatures(klass, repository)
    with pytest.raises(ValueError, match='scientific payload|source map scope'):
        claims.check_parent_class(data, repository)


@pytest.mark.parametrize('damage', ('clause', 'adoption', 'adoption_type', 'other_pin', 'member_pin'))
def test_class_meaning_adoption_and_actual_hashes_are_required(data, repository, damage):
    path = claims.program_class_v2.PATH
    klass = load(repository, path)
    if damage == 'clause':
        klass['resource_envelope']['source_word_literals_max'] += 1
    elif damage == 'adoption':
        klass['adoption']['changes_current_registry'] = False
    elif damage == 'adoption_type':
        klass['adoption']['changes_current_registry'] = 1
    else:
        pins = klass['interpreter_binding']['source_hashes']
        target = (claims.custody_renewal.CLASS_LIVE_RECEIPT if damage == 'member_pin'
                  else 'model/native_program/interpreter.py')
        pins[target] = '0'*64
    write(repository, path, klass)
    if damage != 'member_pin':
        resign(repository)
    with pytest.raises(ValueError, match='class semantics|adoption semantics|owning bytes'):
        claims.check_parent_class(data, repository)


@pytest.mark.parametrize('damage', ('source', 'archive', 'catalog'))
def test_actual_sources_and_original_member_authorities_remain_exact(data, repository, damage):
    catalog = load(repository, claims.renewal.CATALOG)
    spec = catalog['owners']['native_program']
    target = {'source': 'model/native_program/interpreter.py',
              'archive': spec['retained_path'], 'catalog': claims.renewal.CATALOG}[damage]
    with (repository/target).open('a') as output:
        output.write('\n')
    with pytest.raises(ValueError, match='stale reviewed source|reviewed transitive executable source changed|Retained member|owner authority'):
        claims.check_parent_class(data, repository)


def test_sealed_component_does_not_adopt_current_hypothesis(data, repository):
    write(repository, 'claims/executable_model_hypothesis.json', {'unrelated_current_law': True})
    claims.check_parent_class(data, repository)
    result = claims.check_sealed_class(claims.ROOT)
    assert result == {'sealed_class_valid': True, 'official_CCF_LT2_adoption_claimed': False}


def test_original_sealed_class_bytes_remain_exact(monkeypatch):
    original = claims.subprocess.check_output
    def changed(command, **kwargs):
        raw = original(command, **kwargs)
        if command[:2] == ['git', 'show'] and command[2].endswith(':' + str(claims.program_class_v2.PATH)):
            return raw + b' '
        return raw
    monkeypatch.setattr(claims.subprocess, 'check_output', changed)
    with pytest.raises(ValueError, match='Sealed parent class bytes'):
        claims.check_sealed_class.__wrapped__(claims.ROOT)


def test_inventory_helper_change_needs_exact_consumed_semantics(repository, monkeypatch):
    catalog = load(repository, claims.renewal.CATALOG)
    target = repository / claims.INVENTORY_HELPER
    text = target.read_text().replace('if not condition:', 'if condition:', 1)
    target.write_text(text)
    monkeypatch.setattr(claims, 'CURRENT_INVENTORY_SHA256', claims.sha(target))
    with pytest.raises(ValueError, match='Consumed inventory helper semantics'):
        claims.verify_live_executable_sources(catalog['reviewed_executable_source_sha256'], repository)
