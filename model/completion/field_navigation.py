"""Exact reviewed ledger projection from the sealed r2063 checkpoint.

Four prose values are replaced, evidence is appended, and four declared
assumptions enter the two existing binding unions. Every other semantic value
is conserved. Baseline file/blob hashes authenticate the retained checkpoint;
the complete old narratives need not be duplicated in another receipt.
This module does not discharge a physical gate or replay scientific producers.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path

from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
AUTHORITY = 'model/completion/field_navigation_authority.json'
AUTHORITY_SHA256 = '0b8c5169a4396a62fd324c6109a4a677331b9a6ac4fec943321d77afbd47987a'
REGISTRY = 'claims/claim_registry.yaml'
TRACKING = 'tracking/premise_register.json'
DICTIONARY = 'claims/assumption_dictionary.md'
JSON_FILES = (REGISTRY, TRACKING)


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def digest(value):
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def bytes_digest(value):
    return hashlib.sha256(value).hexdigest()


def at(value, path):
    for key in path:
        value = value[key]
    return value


def require(test, message):
    if not test:
        raise ValueError(message)


def authority(root=ROOT):
    root = Path(root)
    raw = (root / AUTHORITY).read_bytes()
    require(bytes_digest(raw) == AUTHORITY_SHA256, 'reviewed field navigation authority changed')
    value = strict_json.loads(raw.decode())
    predecessor = value['predecessor_custody_authority']
    require(bytes_digest((root / predecessor['path']).read_bytes()) == predecessor['sha256'],
            'sealed predecessor authority changed')
    for name, expected in value['new_evidence_sha256'].items():
        path = (root / name).resolve()
        require(path.is_relative_to(root.resolve()) and path.is_file()
                and bytes_digest(path.read_bytes()) == expected,
                'reviewed field evidence changed: ' + name)
    return value


def normalize(documents, approved):
    """Invert only enumerated additions; compare all nonprose predecessor data."""
    result = copy.deepcopy(documents)
    for row in approved['prose_replacements']:
        container = at(result[row['file']], row['json_path'][:-1])
        key = row['json_path'][-1]
        require(container[key] == row['after'] and digest(container[key]) == row['after_value_sha256'],
                'unreviewed narrative replacement')
        container[key] = '<REVIEWED_PROSE_FIELD>'
    for row in approved['ordered_appends']:
        values = at(result[row['file']], row['json_path'])
        additions = row['items']
        require(bool(additions) and values[-len(additions):] == additions,
                'reviewed evidence or assumption append changed')
        require(all(values.count(v) == 1 for v in additions),
                'reviewed addition duplicated or previously present')
        del values[-len(additions):]
    for row in approved['assumption_unions']:
        container = at(result[row['file']], row['json_path'][:-1])
        key = row['json_path'][-1]
        require(container[key] == row['after'], 'reviewed binding assumption union changed')
        container[key] = copy.deepcopy(row['before'])
    for name, expected in approved['normalized_predecessor_semantics_sha256'].items():
        require(digest(result[name]) == expected,
                'nonlisted ledger content, evidence, status or scope changed: ' + name)
    return result


def project(baseline, approved):
    result = copy.deepcopy(baseline)
    for row in approved['prose_replacements']:
        container = at(result[row['file']], row['json_path'][:-1])
        key = row['json_path'][-1]
        require(digest(container[key]) == row['before_value_sha256'], 'retained narrative differs')
        container[key] = row['after']
    for row in approved['ordered_appends']:
        values = at(result[row['file']], row['json_path'])
        require(all(v not in values for v in row['items']), 'listed append is not new')
        values.extend(copy.deepcopy(row['items']))
    for row in approved['assumption_unions']:
        container = at(result[row['file']], row['json_path'][:-1])
        key = row['json_path'][-1]
        require(container[key] == row['before'], 'retained binding assumptions differ')
        container[key] = copy.deepcopy(row['after'])
    normalize(result, approved)
    return result


def validate(root=ROOT, approved=None):
    root = Path(root)
    approved = authority(root) if approved is None else approved
    documents = {name: strict_json.load(root / name) for name in JSON_FILES}
    normalize(documents, approved)
    raw = (root / DICTIONARY).read_bytes()
    suffix = approved['dictionary_append'].encode()
    require(bool(suffix) and raw.endswith(suffix), 'reviewed assumption definitions changed')
    require(bytes_digest(raw[:-len(suffix)]) == approved['baseline_files'][DICTIONARY]['sha256'],
            'retained assumption dictionary changed')
    return {'schema': 'oph.field_navigation_validation.v1',
            'retained_checkpoint': approved['immutable_checkpoint'],
            'prose_fields_replaced': len(approved['prose_replacements']),
            'added_evidence_files': len(approved['new_evidence_sha256']),
            'added_assumptions': list(approved['new_assumption_definitions']),
            'claim_count': len(documents[REGISTRY]['claims']),
            'physical_statuses_and_nonlisted_semantics_exactly_preserved': True,
            'scientific_producer_replay_claimed_by_this_projection': False,
            'current_files': {name: bytes_digest((root / name).read_bytes())
                              for name in (*JSON_FILES, DICTIONARY)}}


def apply(root=ROOT):
    root = Path(root)
    approved = authority(root)
    raw = {name: (root / name).read_bytes() for name in (*JSON_FILES, DICTIONARY)}
    if all(bytes_digest(v) == approved['baseline_files'][name]['sha256'] for name, v in raw.items()):
        result = project({name: strict_json.loads(raw[name].decode()) for name in JSON_FILES}, approved)
        outputs = {name: (json.dumps(result[name], indent=2, ensure_ascii=False) + '\n').encode()
                   for name in JSON_FILES}
        outputs[DICTIONARY] = raw[DICTIONARY] + approved['dictionary_append'].encode()
        for name, value in outputs.items():
            require((root / name).read_bytes() == raw[name], 'ledger changed during projection')
        for name, value in outputs.items():
            (root / name).write_bytes(value)
    return validate(root, approved)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--apply', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    print(json.dumps(apply() if args.apply else validate(), indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
