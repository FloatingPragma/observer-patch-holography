"""Successor custody protects history, scientific semantics and every live pin."""
import copy
import json
import os
import shutil
import pytest
from model.completion import r2063_custody_renewal as renewal
from tools import strict_json


def test_current_custody_preserves_scientific_status_and_historical_bytes():
    result = renewal.validate()
    assert result['classification_counts'] == {'domains': 156, 'realized': 151}
    assert result['new_mathematical_acceptance'] is False
    assert result['historical_evidence_rewritten'] is False
    assert result['navigation']['lean_module_count'] == 681
    assert result['navigation']['claim_count'] == 392


def test_successor_binds_accepted_main_without_relabelling_november_seven():
    authority = renewal.load_authority()
    assert authority['immutable_base_commit'] == 'bd5e6ba80bb3e24c09c6e2289e9b72a627a039fa'
    pointers = renewal.program_pointers(strict_json.load(renewal.ROOT / renewal.PROGRAM))
    assert len(pointers) == 813 and len({p for _, p in pointers}) == 342
    old = strict_json.load(renewal.ROOT / 'model/completion/custody_renewal_20261107.json')
    assert old['date'] == '2026-11-07'
    assert old['program_pointer_schema']['locations'] == 736
    assert old['program_pointer_schema']['targets'] == 267


@pytest.mark.parametrize('field,value', [
    ('date', '2026-11-07'), ('new_mathematical_acceptance', True),
    ('gate_statuses_changed', True), ('canonical_physical_premises_discharged', ['A1']),
    ('historical_evidence_rewritten', True), ('authority', {}),
    ('program_pointer_schema', {'locations': 736}), ('classification_counts', {'domains': 999}),
    ('conserved_semantics', {}), ('navigation', {}), ('current_bindings', {}),
])
def test_record_cannot_promote_science_or_forge_its_authority(field, value):
    packet = strict_json.load(renewal.ROOT / renewal.RECORD)
    packet[field] = value
    with pytest.raises(ValueError):
        renewal.validate(packet)


@pytest.mark.parametrize('damage', ('add', 'retarget', 'move'))
def test_live_pointer_schema_cannot_expand_or_retarget(damage):
    packet = strict_json.load(renewal.ROOT / renewal.PROGRAM)
    path, _ = renewal.program_pointers(packet)[0]
    if damage == 'add':
        packet['extra_evidence'] = copy.deepcopy(renewal.at_path(packet, path))
    elif damage == 'retarget':
        renewal.at_path(packet, path)['path'] = renewal.MAP
    else:
        packet['moved_evidence'] = renewal.at_path(packet, path[:-1]).pop(path[-1])
    with pytest.raises(ValueError, match='immutable live pointer schema'):
        renewal.program_pointers(packet)


@pytest.fixture(scope='module')
def custody_template(tmp_path_factory):
    """Hardlink read-only inputs once; every mutation atomically replaces."""
    root = tmp_path_factory.mktemp('custody-template')
    names = {renewal.AUTHORITY, renewal.RECORD, renewal.PROGRAM, renewal.CLASS,
             renewal.MAP, renewal.REVIEW, renewal.INTERFACES, renewal.TRACKING}
    authority = renewal.load_authority()
    names.update(authority['historical_file_sha256'])
    names.update(authority['reviewed_live_interface_sha256'])
    for name, row in authority['permitted_integration_configuration_renewals'].items():
        names.update((name, row['retained_path']))
        names.update(strict_json.load(renewal.ROOT / name)['source_bindings'])
    from tools import renew_canonical_live_receipts as live
    for name, row in authority['permitted_canonical_live_renewals'].items():
        names.update((name, row['retained_path']))
        value = strict_json.load(renewal.ROOT / name)
        if name != 'model/completion/campaign.json':
            names.update(str(target.relative_to(renewal.ROOT))
                         for _, _, _, target in live.bindings(renewal.ROOT, name, value))
        for pointer in row.get('reviewed_receipt_pointers', ()):
            names.update((pointer['target'], pointer['retained_path']))
    names.update(renewal.CURRENT_BINDINGS)
    for name in (renewal.MAP, renewal.REVIEW, renewal.INTERFACES):
        names.update(strict_json.load(renewal.ROOT / name)['source_bindings'])
    names.update(p for _, p in renewal.program_pointers(strict_json.load(renewal.ROOT / renewal.PROGRAM)))
    for name in sorted(names):
        target = root / name
        target.parent.mkdir(parents=True, exist_ok=True)
        os.link(renewal.ROOT / name, target)
    return root


@pytest.fixture
def isolated_root(custody_template, tmp_path):
    root = tmp_path / 'root'
    shutil.copytree(custody_template, root, copy_function=os.link)
    return root


def read(root, name):
    return strict_json.load(root / name)


def write(root, name, value):
    # Never write through a hardlink to the repository or another test.
    path = root / name
    path.unlink()
    path.write_text(json.dumps(value, indent=2) + '\n')


@pytest.mark.parametrize('damage', (
    'gate_status', 'attachment', 'canonical_contract', 'interface_scope',
    'class_clause', 'review_history', 'assumption_union', 'extra_assumption',
    'canonical_discharge', 'unreviewed_group_text', 'navigation_count',
    'historical_record', 'authority', 'owner_science',
    'r2063_assumption', 'r2063_remaining', 'r2063_source',
))
def test_rebuilt_outer_hashes_do_not_authorize_scientific_or_historical_changes(isolated_root, damage):
    root = isolated_root
    if damage == 'gate_status':
        name = renewal.PROGRAM
        value = read(root, name)
        next(r for r in value['stages'] if r['id'] == 'CP-02')['status'] = 'ACCEPTED'
    elif damage in ('attachment', 'canonical_contract'):
        name = renewal.MAP
        value = read(root, name)
        if damage == 'attachment':
            value['attachments'][0]['classification'] = 'needs_correction'
        else:
            value['canonical_ownership'][0]['contract_sha256'] = 'f' * 64
    elif damage in ('interface_scope', 'class_clause', 'review_history'):
        name = {'interface_scope': renewal.INTERFACES, 'class_clause': renewal.CLASS,
                'review_history': renewal.REVIEW}[damage]
        value = read(root, name)
        value['unreviewed_scientific_change'] = True
    elif damage in ('assumption_union', 'extra_assumption', 'canonical_discharge',
                    'unreviewed_group_text', 'navigation_count', 'r2063_assumption',
                    'r2063_remaining', 'r2063_source'):
        name = renewal.TRACKING
        value = read(root, name)
        if damage == 'r2063_assumption':
            value['model_supply']['contract_coverage']['bindings'][1]['assumptions'].remove(
                'declared_collective_native_block_architecture')
        elif damage == 'r2063_remaining':
            value['model_supply']['groups'][18]['remaining'] += ' Unreviewed physical acceptance.'
        elif damage == 'r2063_source':
            value['model_supply']['groups'][8]['sources'].pop()
        elif damage == 'canonical_discharge':
            value['rows'][0]['discharge']['status'] = 'discharged'
        elif damage == 'unreviewed_group_text':
            value['model_supply']['groups'][8]['constructed'] += ' New unreviewed result.'
        elif damage == 'navigation_count':
            value['model_supply']['contract_coverage']['lean_files'].pop()
        else:
            assumptions = value['model_supply']['contract_coverage']['bindings'][1]['assumptions']
            if damage == 'assumption_union':
                assumptions.remove('declared_family_collision_clock')
            else:
                assumptions.append('unreviewed_extra_assumption')
    elif damage == 'historical_record':
        name = 'model/completion/custody_renewal_20261107.json'
        value = read(root, name)
        value['date'] = '2026-11-11'
    elif damage == 'authority':
        name = renewal.AUTHORITY
        value = read(root, name)
        value['baseline_semantics']['program'] = 'f' * 64
    else:
        name = renewal.OWNER
        value = read(root, name)
        value['scope']['physical_gate_accepted'] = True
    write(root, name, value)
    packet = read(root, renewal.RECORD)
    if name in packet['current_bindings']:
        packet['current_bindings'][name] = renewal.file_digest(root / name)
    with pytest.raises(ValueError):
        renewal.validate(packet, root)


@pytest.mark.parametrize('index', (0, 406, 812))
def test_every_live_pointer_is_checked_against_target_bytes(isolated_root, index):
    value = read(isolated_root, renewal.PROGRAM)
    path, _ = renewal.program_pointers(value)[index]
    renewal.at_path(value, path)['sha256'] = 'f' * 64
    write(isolated_root, renewal.PROGRAM, value)
    with pytest.raises(ValueError, match='live program custody target changed'):
        renewal.validate(root=isolated_root)


@pytest.mark.parametrize('damage', ('omit', 'add', 'hash'))
def test_current_transitive_source_custody_remains_strict(isolated_root, damage):
    value = read(isolated_root, renewal.INTERFACES)
    if damage == 'omit':
        del value['source_bindings'][next(iter(value['source_bindings']))]
    elif damage == 'add':
        value['source_bindings']['invented/source.py'] = 'f' * 64
    else:
        value['source_bindings']['Lean/lakefile.lean'] = 'f' * 64
    write(isolated_root, renewal.INTERFACES, value)
    with pytest.raises(ValueError, match='live (binding path coverage|source binding) changed'):
        renewal.validate(root=isolated_root)


def test_renewal_rejects_unreviewed_scientific_source_before_writing(isolated_root):
    path = isolated_root / 'model/matched_einstein_limit/frame.py'
    original = path.read_text()
    path.unlink()
    path.write_text(original + '\n# unreviewed source edit\n')
    before = (isolated_root / renewal.MAP).read_bytes()
    with pytest.raises(ValueError, match='unreviewed source change cannot be renewed'):
        renewal.renew(isolated_root)
    assert (isolated_root / renewal.MAP).read_bytes() == before


def test_self_consistent_live_repin_cannot_replace_immutable_scientific_source(isolated_root):
    name = 'model/matched_einstein_limit/frame.py'
    path = isolated_root / name
    original = path.read_text()
    path.unlink()
    path.write_text(original + '\n# unreviewed scientific implementation\n')
    document = read(isolated_root, renewal.INTERFACES)
    document['source_bindings'][name] = renewal.file_digest(path)
    with pytest.raises(ValueError, match='immutable nonrenewable source custody changed'):
        renewal.check_bindings(document, renewal.INTERFACES, isolated_root, renewal.load_authority())


@pytest.mark.parametrize('damage', ('acceptance', 'extra_source', 'deleted_source', 'scope'))
def test_configuration_transition_retains_science_even_if_outer_digest_is_forged(isolated_root, damage):
    authority = copy.deepcopy(renewal.load_authority())
    name = 'model/pair_cyclic_control/campaign.json'
    value = read(isolated_root, name)
    if damage == 'acceptance':
        value['integrated_gate_accepted'] = True
    elif damage == 'extra_source':
        value['source_bindings']['tools/unreviewed.py'] = 'f' * 64
    elif damage == 'deleted_source':
        value['source_bindings'].pop(next(iter(value['source_bindings'])))
    else:
        value['scope'] = 'PHYSICAL_COMPLETION'
    write(isolated_root, name, value)
    authority['permitted_integration_configuration_renewals'][name]['replay_verified_sha256'] = renewal.file_digest(isolated_root / name)
    with pytest.raises(ValueError, match='scientific configuration changed'):
        renewal.validate_integration_configurations(isolated_root, authority)


@pytest.mark.parametrize('damage', ('science', 'source_inventory'))
def test_live_canonical_science_survives_forged_outer_receipt_pin(isolated_root, damage):
    authority = copy.deepcopy(renewal.load_authority())
    name = 'model/native_program/results.json'
    value = read(isolated_root, name)
    if damage == 'science':
        value['unreviewed_physical_acceptance'] = True
    else:
        value['source_sha256'].pop(next(iter(value['source_sha256'])))
    write(isolated_root, name, value)
    authority['permitted_canonical_live_renewals'][name]['replay_verified_sha256'] = renewal.file_digest(isolated_root / name)
    with pytest.raises(ValueError, match='scientific payload|source map scope'):
        renewal.validate_canonical_live_inputs(isolated_root, authority)


@pytest.mark.parametrize('damage', ('old_attribution', 'receipt', 'class_fingerprint',
                                   'review_fingerprint', 'class_clause', 'status'))
def test_interpreter_custody_successor_cannot_rewrite_review_or_class(damage):
    authority = renewal.load_authority()
    documents = {p: read(renewal.ROOT, p) for p in (
        renewal.PROGRAM, renewal.MAP, renewal.CLASS, renewal.REVIEW, renewal.INTERFACES)}
    klass, review = documents[renewal.CLASS], documents[renewal.REVIEW]
    if damage == 'old_attribution':
        review['current_interpreter_custody_successor']['old_review_signed_current_bytes'] = True
    elif damage == 'receipt':
        klass['interpreter_binding']['source_hashes']['model/native_program/results.json'] = 'f' * 64
    elif damage == 'class_fingerprint':
        klass['adoption']['independent_review']['reviewed_contract_digest'] = 'f' * 64
    elif damage == 'review_fingerprint':
        review['reviewed_contract_digest'] = 'f' * 64
    elif damage == 'class_clause':
        klass['services']['coherent_services_derived_from_passive_repair'] = True
    else:
        klass['adoption']['status'] = 'PROPOSED_FOR_INDEPENDENT_REVIEW'
    if damage == 'status':
        bodies = renewal.semantic_bodies(documents, authority)
        assert renewal.digest(bodies['class']) != authority['baseline_semantics']['class']
    else:
        with pytest.raises(ValueError):
            renewal.semantic_bodies(documents, authority)
