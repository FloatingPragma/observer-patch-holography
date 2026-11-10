"""Renew explicit live custody while preserving the accepted-main semantics.

The November 7 record and checker remain historical bytes. This successor
authenticates the larger already-accepted program at bd5e6ba80; it neither
backdates that program nor promotes any scientific gate.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path

from tools import strict_json
from tools.build_completion_program import integration_source_paths, validate_integration_sources

ROOT = Path(__file__).resolve().parents[2]
AUTHORITY = 'model/completion/custody_authority_20261111.json'
AUTHORITY_SHA256 = '9537d11b1c352da10f0531f9bee002cc9695a307a86f1aefe831a386b45b19d4'
RECORD = 'model/completion/custody_renewal_20261111.json'
REVIEW = 'model/completion/program_class_v2_review.json'
CLASS = 'model/completion/program_class_v2.json'
MAP = 'model/completion/finite_attachment_map.json'
PROGRAM = 'claims/completion_program.json'
INTERFACES = 'model/completion/component_interfaces.json'
TRACKING = 'tracking/premise_register.json'
OWNER = 'model/matched_einstein_limit/results.json'
CHILD = 'model/completion/spectrum_einstein_campaign.json'
SELF = 'model/completion/current_custody_renewal.py'
TEST = 'model/completion/test_custody_renewal.py'
CURRENT_TOOL = 'tools/build_completion_program_current.py'
DATE = '2026-11-11'
BASE = 'bd5e6ba80bb3e24c09c6e2289e9b72a627a039fa'
ADDITIONS = {INTERFACES: {SELF, AUTHORITY, CURRENT_TOOL}, MAP: set(), REVIEW: {SELF, AUTHORITY, RECORD, CURRENT_TOOL}}
CURRENT_BINDINGS = {INTERFACES, MAP, OWNER, CHILD, SELF, TEST, AUTHORITY, CURRENT_TOOL,
                    'code/a5_closure/manifests/source_current_order_sensitive_inventory.json'}
RENEWABLE_TARGETS = {INTERFACES, MAP, OWNER, CHILD, SELF, TEST, AUTHORITY, CURRENT_TOOL, RECORD,
                     'Lean/lakefile.lean', 'tools/build_completion_program.py', 'tools/test_completion_program.py',
                     'model/completion/README.md'}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def file_digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_authority(root=ROOT):
    require(file_digest(Path(root) / AUTHORITY) == AUTHORITY_SHA256, 'immutable custody authority changed')
    value = strict_json.load(Path(root) / AUTHORITY)
    require(value['immutable_base_commit'] == BASE, 'immutable accepted-main baseline changed')
    require(value['renewable_targets'] == sorted(RENEWABLE_TARGETS), 'reviewed renewal target boundary changed')
    return value


def at_path(value, path):
    for key in path:
        value = value[key]
    return value


def program_pointers(program, authority=None):
    authority = load_authority() if authority is None else authority
    rows = []
    def walk(value, path=()):
        if isinstance(value, dict):
            if 'path' in value and 'sha256' in value:
                rows.append((path, value['path']))
            for key in sorted(value):
                walk(value[key], path + (key,))
        elif isinstance(value, list):
            for index, child in enumerate(value):
                walk(child, path + (index,))
    walk(program)
    schema = {'locations': len(rows), 'targets': len({p for _, p in rows}), 'sha256': digest(rows)}
    require(schema == authority['program_pointer_schema'], 'immutable live pointer schema changed')
    return rows


def semantic_bodies(documents, authority):
    program = copy.deepcopy(documents[PROGRAM])
    for path, _ in program_pointers(program, authority):
        del at_path(program, path)['sha256']
    finite = {k: copy.deepcopy(v) for k, v in documents[MAP].items() if k != 'source_bindings'}
    for key in ('claim_count', 'lean_module_count', 'contract_coverage_sha256'):
        finite['existing_navigation_coverage'][key] = '<current navigation>'
    klass = copy.deepcopy(documents[CLASS])
    klass['adoption']['independent_review']['sha256'] = '<current review>'
    return {
        'program': program,
        'finite': finite,
        'class': klass,
        'interfaces': {k: v for k, v in documents[INTERFACES].items() if k != 'source_bindings'},
        'review': {k: v for k, v in documents[REVIEW].items()
                   if k not in ('source_bindings', 'current_source_custody_renewal')},
    }


def tracking_transition(register, authority):
    """Reverse only the owner's eight exact, typed append/union changes."""
    result = copy.deepcopy(register)
    transition = authority['permitted_tracking_transition']
    for change in transition['groups']:
        group = result['model_supply']['groups'][change['group_index']]
        require(group['id'] == change['group_id'], 'navigation group identity changed')
        suffix = change['constructed_append']
        require(bool(suffix) and group['constructed'].endswith(suffix), 'unreviewed navigation description change')
        group['constructed'] = group['constructed'][:-len(suffix)]
        for addition in change['source_additions']:
            require(group['sources'].count(addition) == 1, 'navigation source append changed')
            group['sources'].remove(addition)
    additions = transition['assumptions_added']
    for index in transition['binding_indices']:
        assumptions = result['model_supply']['contract_coverage']['bindings'][index]['assumptions']
        require(all(assumptions.count(a) == 1 for a in additions), 'navigation assumption union changed')
        result['model_supply']['contract_coverage']['bindings'][index]['assumptions'] = [a for a in assumptions if a not in additions]
    require(digest(result) == authority['baseline_semantics']['tracking'],
            'tracking changed outside the exact reviewed append/union transition')
    coverage = register['model_supply']['contract_coverage']
    require(len(coverage['claim_ids']) == authority['baseline_navigation']['claim_count']
            and len(coverage['lean_files']) == authority['baseline_navigation']['lean_module_count'],
            'current navigation inventory changed')
    return {'claim_count': len(coverage['claim_ids']), 'lean_module_count': len(coverage['lean_files']),
            'contract_coverage_sha256': digest(coverage)}


def authenticated_documents(root, authority):
    documents = {p: strict_json.load(root / p) for p in (PROGRAM, MAP, CLASS, REVIEW, INTERFACES)}
    for name, value in semantic_bodies(documents, authority).items():
        require(digest(value) == authority['baseline_semantics'][name],
                'accepted-main scientific semantics changed: ' + name)
    for name, expected in authority['historical_file_sha256'].items():
        require(file_digest(root / name) == expected, 'historical evidence changed: ' + name)
    tracking = strict_json.load(root / TRACKING)
    navigation = tracking_transition(tracking, authority)
    owner = strict_json.load(root / OWNER)
    require(digest({k: v for k, v in owner.items() if k != 'source_bindings'})
            == authority['matched_einstein_scientific_sha256'], 'owning Einstein scientific payload changed')
    child = strict_json.load(root / CHILD)
    row = next(r for r in child['component_results'] if r['id'] == 'matched_einstein_limit')
    row['evidence'][OWNER] = '<live owning receipt>'
    require(digest(child) == authority['live_child_evidence']['semantic_sha256'],
            'spectrum-Einstein campaign changed outside its one live owning receipt')
    return documents, navigation


def check_binding_scope(document, path, authority):
    bindings = document['source_bindings']
    added = ADDITIONS[path]
    require(digest(sorted(set(bindings) - added)) == authority['baseline_binding_path_sha256'][path],
            'live binding path coverage changed: ' + path)
    static = {p: v for p, v in bindings.items() if p not in RENEWABLE_TARGETS}
    require(digest(static) == authority['baseline_static_binding_sha256'][path],
            'immutable nonrenewable source custody changed: ' + path)


def check_bindings(document, path, root, authority):
    check_binding_scope(document, path, authority)
    bindings = document['source_bindings']
    require(ADDITIONS[path] <= set(bindings), 'live binding path coverage changed: ' + path)
    for name, expected in bindings.items():
        source = (root / name).resolve()
        require(source.is_relative_to(root) and source.is_file(), 'missing or escaping live source: ' + name)
        require(file_digest(source) == expected, 'live source binding changed: ' + name)


def validate(packet=None, root=ROOT):
    root = Path(root).resolve()
    authority = load_authority(root)
    packet = strict_json.load(root / RECORD) if packet is None else packet
    require(packet.get('schema') == 'oph.current_custody_renewal.v2' and packet.get('date') == DATE,
            'wrong successor custody record')
    require(packet.get('authority') == {'path': AUTHORITY, 'sha256': AUTHORITY_SHA256, 'base_commit': BASE},
            'immutable successor authority changed')
    require(packet.get('new_mathematical_acceptance') is False and packet.get('gate_statuses_changed') is False
            and packet.get('canonical_physical_premises_discharged') == [], 'custody cannot promote scientific acceptance')
    require(packet.get('historical_evidence_rewritten') is False,
            'custody cannot rewrite historical evidence')
    require(packet.get('program_pointer_schema') == authority['program_pointer_schema'],
            'immutable program pointer authority changed')
    require(packet.get('classification_counts') == authority['classification_counts'], 'classification authority changed')
    documents, navigation = authenticated_documents(root, authority)
    require(packet.get('navigation') == navigation, 'successor navigation certificate changed')
    require(all(documents[MAP]['existing_navigation_coverage'][k] == v for k, v in navigation.items()),
            'finite-map live navigation changed')
    require(packet.get('conserved_semantics') == authority['baseline_semantics'], 'immutable conservation authority changed')
    require(packet.get('live_child_evidence') == authority['live_child_evidence'], 'live child evidence authority changed')
    child = strict_json.load(root / CHILD)
    require(next(r for r in child['component_results'] if r['id'] == 'matched_einstein_limit')['evidence'][OWNER]
            == file_digest(root / OWNER), 'live campaign owning receipt changed')
    for name in (INTERFACES, MAP, REVIEW):
        check_bindings(documents[name], name, root, authority)
    validate_integration_sources(root, documents[INTERFACES])
    for path, target in program_pointers(documents[PROGRAM], authority):
        require(at_path(documents[PROGRAM], path)['sha256'] == file_digest(root / target),
                'live program custody target changed: ' + target)
    require(set(packet.get('current_bindings', {})) == CURRENT_BINDINGS, 'current binding coverage changed')
    for path, expected in packet['current_bindings'].items():
        require(file_digest(root / path) == expected, 'current custody binding changed: ' + path)
    expected_pointer = {'date': DATE, 'path': RECORD, 'role': 'CUSTODY_ONLY_NO_NEW_MATHEMATICAL_ACCEPTANCE',
                        'sha256': file_digest(root / RECORD)}
    require(documents[REVIEW]['current_source_custody_renewal'] == expected_pointer,
            'review does not bind successor custody')
    require(documents[CLASS]['adoption']['independent_review']['sha256'] == file_digest(root / REVIEW),
            'adoption custody does not bind current review')
    return {'date': DATE, 'immutable_base_commit': BASE, 'classification_counts': authority['classification_counts'],
            'program_pointer_schema': authority['program_pointer_schema'], 'navigation': navigation,
            'new_mathematical_acceptance': False, 'scientific_status_changed': False,
            'historical_evidence_rewritten': False}


def write_json(root, name, value):
    (root / name).write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n')


def renew(root=ROOT):
    """Explicit maintenance operation; validation alone never repairs custody."""
    root = Path(root).resolve()
    authority = load_authority(root)
    documents, navigation = authenticated_documents(root, authority)
    # A live renewal is not permission to repin arbitrary altered science.
    # The unchanged dependency set is checked before any output is written.
    for name in (INTERFACES, MAP, REVIEW):
        check_binding_scope(documents[name], name, authority)
        changed = {p for p, expected in documents[name]['source_bindings'].items()
                   if file_digest(root / p) != expected}
        require(changed <= RENEWABLE_TARGETS,
                'unreviewed source change cannot be renewed: ' + ', '.join(sorted(changed - RENEWABLE_TARGETS)))
    # Execute the owning exact producer, not just its hash calculator.
    from model.matched_einstein_limit import packet as producer
    require(root == ROOT, 'producer-backed renewal requires the actual repository root')
    fresh = producer.report()
    require(digest({k: v for k, v in fresh.items() if k != 'source_bindings'})
            == authority['matched_einstein_scientific_sha256'], 'fresh owning science changed')
    (root / OWNER).write_text(producer.canonical(fresh))
    child = strict_json.load(root / CHILD)
    next(r for r in child['component_results'] if r['id'] == 'matched_einstein_limit')['evidence'][OWNER] = file_digest(root / OWNER)
    write_json(root, CHILD, child)
    for key, value in navigation.items():
        documents[MAP]['existing_navigation_coverage'][key] = value
    for name in (INTERFACES, MAP):
        keys = list(documents[name]['source_bindings'])
        keys += sorted(ADDITIONS[name] - set(keys))
        if name == INTERFACES:
            require(integration_source_paths(root, documents[name]) <= set(keys),
                    'new transitive dependency requires explicit custody review')
        documents[name]['source_bindings'] = {p: file_digest(root / p) for p in keys}
        write_json(root, name, documents[name])
    record = {
        'schema': 'oph.current_custody_renewal.v2', 'date': DATE,
        'authority': {'path': AUTHORITY, 'sha256': AUTHORITY_SHA256, 'base_commit': BASE},
        'scope': 'Live source custody and the exact owner-approved navigation append/union; no historical proof replay or new gate acceptance.',
        'new_mathematical_acceptance': False, 'gate_statuses_changed': False,
        'canonical_physical_premises_discharged': [], 'historical_evidence_rewritten': False,
        'conserved_semantics': authority['baseline_semantics'],
        'live_child_evidence': authority['live_child_evidence'],
        'classification_counts': authority['classification_counts'],
        'program_pointer_schema': authority['program_pointer_schema'], 'navigation': navigation,
        'current_bindings': {p: file_digest(root / p) for p in sorted(CURRENT_BINDINGS)},
        'owning_producer': {'module': 'model.matched_einstein_limit.packet',
                            'scientific_payload_sha256': authority['matched_einstein_scientific_sha256'],
                            'scientific_payload_exactly_preserved': True,
                            'fresh_Lean_kernel_execution_claimed_by_this_record': False},
    }
    write_json(root, RECORD, record)
    review = documents[REVIEW]
    keys = list(review['source_bindings'])
    keys += sorted(ADDITIONS[REVIEW] - set(keys))
    review['source_bindings'] = {p: file_digest(root / p) for p in keys}
    review['current_source_custody_renewal'] = {
        'date': DATE, 'path': RECORD, 'role': 'CUSTODY_ONLY_NO_NEW_MATHEMATICAL_ACCEPTANCE',
        'sha256': file_digest(root / RECORD)}
    write_json(root, REVIEW, review)
    documents[CLASS]['adoption']['independent_review']['sha256'] = file_digest(root / REVIEW)
    write_json(root, CLASS, documents[CLASS])
    for path, target in program_pointers(documents[PROGRAM], authority):
        at_path(documents[PROGRAM], path)['sha256'] = file_digest(root / target)
    write_json(root, PROGRAM, documents[PROGRAM])
    return validate(root=root)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--check', action='store_true')
    group.add_argument('--renew', action='store_true', help='explicit producer-backed live custody maintenance')
    args = parser.parse_args()
    print(json.dumps(renew() if args.renew else validate(), indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
