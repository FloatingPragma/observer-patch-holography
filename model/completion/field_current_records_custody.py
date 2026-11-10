"""Current custody successor for the reviewed field/current package field/current-record navigation.

The sealed a036 source-record class, review, registry projection and scientific receipt bytes
remain explicit predecessors. This module authenticates executed live outputs,
the exact ledger delta and existing source-map keys. It does not reinterpret a
retained review as a signature of new bytes or accept any physical gate.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import re
from pathlib import Path

from tools import strict_json
from tools import renew_canonical_live_receipts as live
from tools.build_completion_program import validate_integration_sources
from . import field_current_records_navigation as field_navigation
from . import program_class_v2

ROOT = Path(__file__).resolve().parents[2]
AUTHORITY = 'model/completion/field_current_records_custody_authority.json'
AUTHORITY_SHA256 = 'e9fe3efea315ffcc1e877df916eb411a28959eff122f8b6bd57addf91c2b63c0'
RECORD = 'model/completion/field_current_records_custody_current.json'
PROGRAM = 'claims/completion_program.json'
CLASS = 'model/completion/program_class_v2.json'
REVIEW = 'model/completion/program_class_v2_review.json'
MAP = 'model/completion/finite_attachment_map.json'
INTERFACES = 'model/completion/component_interfaces.json'
FILES = (PROGRAM, CLASS, REVIEW, MAP, INTERFACES)
ARCHIVE = 'model/completion/historical/field_current_records_inputs/'
INTERPRETER = 'model/native_program/results.json'
SUCCESSOR_KEY = 'field_current_records_custody_successor'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(root, name):
    path = (Path(root) / name).resolve()
    require(path.is_relative_to(Path(root).resolve()) and path.is_file(),
            'missing or escaping source: ' + name)
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_authority(root=ROOT):
    require(sha(root, AUTHORITY) == AUTHORITY_SHA256, 'field custody authority changed')
    result = strict_json.load(Path(root) / AUTHORITY)
    require(result['immutable_checkpoint'] == 'a0367dc26b07ad959c915054dbfcb3f9828af295',
            'retained checkpoint changed')
    for name, expected in result['historical_sha256'].items():
        require(sha(root, name) == expected, 'historical custody changed: ' + name)
    for name, expected in result['unchanged_executable_sha256'].items():
        require(sha(root, name) == expected, 'scientific executable changed: ' + name)
    for name, expected in result['maintenance_source_sha256'].items():
        require(sha(root, name) == expected, 'reviewed maintenance source changed: ' + name)
    source = (Path(root) / 'model/completion/field_current_records_custody.py').read_text()
    normalized, count = re.subn(r"AUTHORITY_SHA256 = '[a-f0-9]{64}'",
                               "AUTHORITY_SHA256 = '<AUTHORITY_DIGEST>'", source)
    require(count == 1 and hashlib.sha256(normalized.encode()).hexdigest()
            == result['validator_template_sha256'], 'reviewed custody validator changed')
    require(result['field_navigation_authority_sha256'] == field_navigation.AUTHORITY_SHA256,
            'field navigation authority differs')
    return result


def pointers(program):
    result = []
    def visit(value, path=()):
        if isinstance(value, dict):
            if 'path' in value and 'sha256' in value:
                result.append((path, value['path']))
            for key in sorted(value):
                visit(value[key], path + (key,))
        elif isinstance(value, list):
            for index, child in enumerate(value):
                visit(child, path + (index,))
    visit(program)
    return result


def science_view(name, document, additions=()):
    """Remove only the enumerated custody fields, retaining their path scopes."""
    value = copy.deepcopy(document)
    if name == PROGRAM:
        for path, _ in pointers(value):
            live.at_path(value, path)['sha256'] = '<CURRENT_CUSTODY>'
    if name in (INTERFACES, MAP, REVIEW):
        if name == INTERFACES:
            for target in additions:
                value['source_bindings'].pop(target, None)
        value['source_bindings'] = {p: '<CURRENT_CUSTODY>' for p in value['source_bindings']}
    if name == MAP:
        value['existing_navigation_coverage']['contract_coverage_sha256'] = '<REVIEWED_ASSUMPTION_UNION>'
    if name == CLASS:
        value['interpreter_binding']['source_hashes'][INTERPRETER] = '<CURRENT_INTERPRETER_RECEIPT>'
        value['adoption']['independent_review']['sha256'] = '<CURRENT_REVIEW>'
        value['adoption']['independent_review']['reviewed_contract_digest'] = '<CURRENT_CLASS_DIGEST>'
    if name == REVIEW:
        value.pop(SUCCESSOR_KEY, None)
        value['reviewed_contract_digest'] = '<CURRENT_CLASS_DIGEST>'
        value['current_source_custody_renewal'] = '<EXPLICIT_CURRENT_CUSTODY_RECORD>'
    return value


def conserve_documents(documents, approved, root=ROOT):
    for name in FILES:
        original = strict_json.load(Path(root) / (ARCHIVE + name))
        additions = approved['maintenance_input_scope_append'] if name == INTERFACES else ()
        require(live.canonical(science_view(name, original)) == live.canonical(science_view(name, documents[name], additions)),
                'retained scientific contract or binding scope changed: ' + name)
    scope = [[list(path), target] for path, target in pointers(documents[PROGRAM])]
    require(field_navigation.digest(scope) == approved['program_pointer_schema']['sha256'],
            'completion pointer scope changed')


def validate_live_receipts(approved, root=ROOT):
    for name, row in approved['executed_live_receipts'].items():
        require(sha(root, name) == row['current_sha256'],
                'canonical receipt differs from executed output: ' + name)
        old = strict_json.load(Path(root) / row['retained_path'])
        current = strict_json.load(Path(root) / name)
        if name == 'model/completion/campaign.json':
            pointer = row['reviewed_receipt_pointers'][0]
            require(live.at_path(old, pointer['json_path']) == pointer['retained_sha256'],
                    'campaign retained control changed')
            live.set_path(old, pointer['json_path'], sha(root, pointer['target']))
            require(old == current, 'campaign changed outside its one admitted control pointer')
            continue
        projected = live.carry_metadata(name, old, current,
            pointers=row.get('reviewed_receipt_pointers', ()),
            scope_append=row.get('reviewed_source_scope_append'), root=root)
        require(live.canonical(projected) == live.canonical(current),
                'canonical scientific payload changed: ' + name)
        for _, relative, expected, path in live.bindings(Path(root), name, current):
            require(hashlib.sha256(path.read_bytes()).hexdigest() == expected,
                    'executed live receipt source changed: ' + name + ':' + relative)
    for name in approved['integration_configurations']:
        current = strict_json.load(Path(root) / name)
        old = strict_json.load(Path(root) / (ARCHIVE + name))
        for target in ('model/completion/campaign.json', INTERPRETER):
            old['source_bindings'][target] = sha(root, target)
        require(current == old, 'integration configuration changed outside two named pins')
        validate_integration_sources(Path(root), current)


def check_binding_sources(documents, approved, root=ROOT, *, require_current):
    allowed = set(approved['renewable_targets'])
    original_program = strict_json.load(Path(root) / (ARCHIVE + PROGRAM))
    for path, target in pointers(original_program):
        original_hash = live.at_path(original_program, path)['sha256']
        if target not in allowed:
            require(sha(root, target) == original_hash,
                    'unreviewed completion target changed: ' + target)
    for name in (INTERFACES, MAP, REVIEW):
        original = strict_json.load(Path(root) / (ARCHIVE + name))['source_bindings']
        current = documents[name]['source_bindings']
        additions = set(approved['maintenance_input_scope_append']) if name == INTERFACES else set()
        scopes = (set(original) | additions,) if require_current else (set(original), set(original) | additions)
        require(set(current) in scopes, 'source binding inventory changed: ' + name)
        for target, expected in current.items():
            actual = sha(root, target)
            if target in additions:
                if require_current:
                    require(expected == actual, 'maintenance source changed: ' + target)
            elif target not in allowed:
                require(expected == original[target] == actual,
                        'unreviewed source changed: ' + target)
            elif require_current:
                require(expected == actual, 'stale current source: ' + target)


def transition(documents, root=ROOT):
    original = strict_json.load(Path(root) / (ARCHIVE + CLASS))
    return {'role': 'CUSTODY_ONLY_SUCCESSOR_OF_SEALED_A036_SOURCE_RECORD_REVIEW',
            'retained_class': {'path': ARCHIVE + CLASS, 'sha256': sha(root, ARCHIVE + CLASS)},
            'retained_review': {'path': ARCHIVE + REVIEW, 'sha256': sha(root, ARCHIVE + REVIEW)},
            'retained_reviewed_contract_digest': program_class_v2.reviewed_contract_digest(original),
            'current_reviewed_contract_digest': program_class_v2.reviewed_contract_digest(documents[CLASS]),
            'interpreter_receipt': INTERPRETER,
            'retained_interpreter_receipt_sha256': original['interpreter_binding']['source_hashes'][INTERPRETER],
            'current_interpreter_receipt_sha256': sha(root, INTERPRETER),
            'old_review_signed_current_bytes': False,
            'scientific_clauses_and_gate_statuses_exactly_preserved': True}


def record_value(documents, approved, root=ROOT):
    register = strict_json.load(Path(root) / field_navigation.TRACKING)
    coverage = register['model_supply']['contract_coverage']
    return {'schema': 'oph.field_current_records_custody_successor.v1',
            'authority': {'path': AUTHORITY, 'sha256': AUTHORITY_SHA256},
            'immutable_checkpoint': approved['immutable_checkpoint'],
            'scope': 'Exact field/current package field/current-record ledger projection and replay-verified current custody; retained physical identities are unchanged.',
            'new_mathematical_acceptance': False, 'gate_statuses_changed': False,
            'canonical_physical_premises_discharged': [], 'historical_evidence_rewritten': False,
            'old_review_signed_current_bytes': False,
            'classification_counts': {'domains': 156, 'realized': 151},
            'program_pointer_schema': approved['program_pointer_schema'],
            'navigation': {'claim_count': len(coverage['claim_ids']), 'lean_module_count': len(coverage['lean_files']),
                           'contract_coverage_sha256': field_navigation.digest(coverage)},
            'class_custody_successor': transition(documents, root),
            'executed_live_receipts': approved['executed_live_receipts'],
            'maintenance_input_scope_append': {name: sha(root, name)
                for name in approved['maintenance_input_scope_append']},
            'current_bindings': {name: sha(root, name) for name in approved['current_binding_paths']}}


def validate(packet=None, root=ROOT):
    root = Path(root).resolve()
    approved = load_authority(root)
    field_navigation.validate(root)
    validate_live_receipts(approved, root)
    documents = {name: strict_json.load(root / name) for name in FILES}
    conserve_documents(documents, approved, root)
    check_binding_sources(documents, approved, root, require_current=True)
    validate_integration_sources(root, documents[INTERFACES])
    packet = strict_json.load(root / RECORD) if packet is None else packet
    expected = record_value(documents, approved, root)
    require(packet == expected, 'current field custody record changed')
    require(all(documents[MAP]['existing_navigation_coverage'][k] == v
                for k, v in expected['navigation'].items()),
            'finite map navigation differs from reviewed unions')
    klass, review = documents[CLASS], documents[REVIEW]
    require(klass['interpreter_binding']['source_hashes'][INTERPRETER] == sha(root, INTERPRETER),
            'class interpreter receipt changed')
    fingerprint = program_class_v2.reviewed_contract_digest(klass)
    require(klass['adoption']['independent_review']['reviewed_contract_digest'] == fingerprint
            and review['reviewed_contract_digest'] == fingerprint
            and review[SUCCESSOR_KEY] == transition(documents, root),
            'current class review fingerprint or attribution changed')
    require(review['current_source_custody_renewal'] == {
        'date': '2026-11-11', 'path': RECORD, 'role': 'CUSTODY_ONLY_NO_NEW_MATHEMATICAL_ACCEPTANCE',
        'sha256': sha(root, RECORD)}, 'review does not bind the current custody successor')
    require(klass['adoption']['independent_review']['sha256'] == sha(root, REVIEW),
            'class does not bind current review bytes')
    for path, target in pointers(documents[PROGRAM]):
        require(live.at_path(documents[PROGRAM], path)['sha256'] == sha(root, target),
                'completion pointer stale: ' + target)
    return {'retained_checkpoint': approved['immutable_checkpoint'],
            'executed_live_receipts': len(approved['executed_live_receipts']),
            'program_pointer_schema': expected['program_pointer_schema'],
            'navigation': expected['navigation'], 'scientific_statuses_changed': False,
            'historical_evidence_rewritten': False, 'old_review_signed_current_bytes': False}


def write(root, name, value):
    (root / name).write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n')


def renew(root=ROOT):
    root = Path(root).resolve()
    approved = load_authority(root)
    field_navigation.validate(root)
    validate_live_receipts(approved, root)
    documents = {name: strict_json.load(root / name) for name in FILES}
    conserve_documents(documents, approved, root)
    # A reviewed draft may contain an earlier subset of the same enumerated
    # maintenance closure. Conservation above rejects every unlisted key.
    # Materialize the complete approved closure before checking its inventory.
    for target in approved['maintenance_input_scope_append']:
        documents[INTERFACES]['source_bindings'][target] = sha(root, target)
    check_binding_sources(documents, approved, root, require_current=False)
    documents[CLASS]['interpreter_binding']['source_hashes'][INTERPRETER] = sha(root, INTERPRETER)
    fingerprint = program_class_v2.reviewed_contract_digest(documents[CLASS])
    documents[CLASS]['adoption']['independent_review']['reviewed_contract_digest'] = fingerprint
    documents[REVIEW]['reviewed_contract_digest'] = fingerprint
    documents[REVIEW][SUCCESSOR_KEY] = transition(documents, root)
    expected = record_value(documents, approved, root)
    documents[MAP]['existing_navigation_coverage'].update(expected['navigation'])
    for name in (INTERFACES, MAP):
        targets = list(documents[name]['source_bindings'])
        if name == INTERFACES:
            targets += sorted(set(approved['maintenance_input_scope_append']) - set(targets))
        documents[name]['source_bindings'] = {target: sha(root, target) for target in targets}
        if name == INTERFACES:
            validate_integration_sources(root, documents[name])
        write(root, name, documents[name])
    write(root, RECORD, record_value(documents, approved, root))
    documents[REVIEW]['source_bindings'] = {target: sha(root, target) for target in documents[REVIEW]['source_bindings']}
    documents[REVIEW]['current_source_custody_renewal'] = {
        'date': '2026-11-11', 'path': RECORD, 'role': 'CUSTODY_ONLY_NO_NEW_MATHEMATICAL_ACCEPTANCE',
        'sha256': sha(root, RECORD)}
    write(root, REVIEW, documents[REVIEW])
    documents[CLASS]['adoption']['independent_review']['sha256'] = sha(root, REVIEW)
    write(root, CLASS, documents[CLASS])
    for path, target in pointers(documents[PROGRAM]):
        live.at_path(documents[PROGRAM], path)['sha256'] = sha(root, target)
    write(root, PROGRAM, documents[PROGRAM])
    return validate(root=root)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--renew', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    print(json.dumps(renew() if args.renew else validate(), indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
