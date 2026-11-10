"""Admit two actual late-owner executions without changing the sealed class."""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
import subprocess
import sys

from tools import strict_json
from tools import renew_canonical_live_receipts as live

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
AUTHORITY = 'model/field_current_late_custody/authority.json'
AUTHORITY_SHA256 = '803db6ccf068a02841c834ecc62d731a119487a3a940137ee3a7e1b124b52282'
OUTPUT = 'model/completion/field_current_class_consumer_replays.json'
COMPLETION = 'model/particle_identification/completion_results.json'
REGIONAL = 'model/regional_clock/results.json'
CLASS = 'model/completion/program_class_v2.json'
EXECUTIONS = 'model/field_current_late_custody/executions.json'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def authority(root=ROOT):
    root = Path(root)
    if sha(root / AUTHORITY) != AUTHORITY_SHA256:
        raise ValueError('late-owner authority changed')
    value = strict_json.load(root / AUTHORITY)
    if set(value['retained_receipts']) != {COMPLETION, REGIONAL}:
        raise ValueError('late-owner inventory changed')
    for name, expected in value['input_sha256'].items():
        if sha(root / name) != expected:
            raise ValueError('sealed late-owner input changed: ' + name)
    for row in value['retained_receipts'].values():
        if sha(root / row['retained_path']) != row['retained_sha256']:
            raise ValueError('retained late-owner receipt changed')
    return value


def compare_completion(old, fresh, approved, *, root=ROOT):
    """Allow one exact child subtree and four exact source-map expansions."""
    expanded = copy.deepcopy(old)
    spec = approved['exact_added_subtree']
    subtree = live.at_path(strict_json.load(Path(root) / spec['source_receipt']), spec['source_json_path'])
    if hashlib.sha256(live.canonical(subtree).encode()).hexdigest() != spec['canonical_sha256']:
        raise ValueError('authorized population child result changed')
    parent = live.at_path(expanded, spec['target_json_path'][:-1])
    if spec['target_json_path'][-1] in parent:
        raise ValueError('additive result already existed in retained input')
    parent[spec['target_json_path'][-1]] = copy.deepcopy(subtree)
    maps = {where: pins for _, pins, where in live.source_maps(expanded)}
    fresh_maps = {where: pins for _, pins, where in live.source_maps(fresh)}
    if maps.keys() != fresh_maps.keys():
        raise ValueError('late source-map inventory changed')
    additions = approved['exact_source_additions_by_location']
    if not set(additions) <= maps.keys():
        raise ValueError('authorized source-map location missing')
    for where, pins in maps.items():
        extra = additions.get(where, {})
        if set(extra) & set(pins):
            raise ValueError('source addition overwrites a retained path')
        if set(fresh_maps[where]) != set(pins) | set(extra):
            raise ValueError('unreviewed source-map scope change: ' + where)
        for name, expected in extra.items():
            if sha(Path(root) / name) != expected or fresh_maps[where][name] != expected:
                raise ValueError('reviewed source addition differs: ' + where + ':' + name)
        pins.update(extra)
    admitted = live.carry_metadata(COMPLETION, expanded, fresh, root=root)
    if live.canonical(admitted) != live.canonical(fresh):
        raise ValueError('older completion science changed or child addition differs')
    return admitted


def check_bindings(name, value, root=ROOT):
    for _, relative, expected, path in live.bindings(Path(root), name, value):
        if sha(path) != expected:
            raise ValueError('late-owner live source changed: ' + name + ':' + relative)


def compare_receipts(completion, regional, *, root=ROOT):
    approved = authority(root)
    old = {name: strict_json.load(Path(root) / row['retained_path'])
           for name, row in approved['retained_receipts'].items()}
    compare_completion(old[COMPLETION], completion, approved, root=root)
    admitted = live.carry_metadata(REGIONAL, old[REGIONAL], regional, root=root)
    if live.canonical(admitted) != live.canonical(regional):
        raise ValueError('regional canonical science changed')
    check_bindings(COMPLETION, completion, root)
    check_bindings(REGIONAL, regional, root)
    return approved


def own_pins(root=ROOT):
    folder = Path(root) / HERE.relative_to(ROOT)
    return {str(p.relative_to(root)): sha(p) for p in sorted(folder.iterdir())
            if p.suffix in {'.py', '.md', '.json'}}


def validate_execution(execution, receipt_hashes, approved, *, root=ROOT):
    root = Path(root)
    expected_producers = {
        COMPLETION: ('model.particle_identification.completion_certificate:report',
                     'model/particle_identification/completion_certificate.py'),
        REGIONAL: ('model.regional_clock.verify:report', 'model/regional_clock/verify.py'),
    }
    if set(execution) != set(expected_producers):
        raise ValueError('late execution record inventory changed')
    for name, (producer, source) in expected_producers.items():
        row = execution[name]
        if (row['producer'] != producer or row['producer_source_sha256'] != sha(root / source)
                or row['class_sha256'] != approved['sealed_class_sha256']
                or row['admitted_sha256'] != receipt_hashes[name]
                or row['actual_unchanged_owner_executed'] is not True):
            raise ValueError('actual late-owner execution does not bind this result')
    if execution[COMPLETION]['fresh_sha256'] != receipt_hashes[COMPLETION]:
        raise ValueError('completion receipt differs from full executed report')
    if execution[REGIONAL]['original_complete_checker_executed'] is not True:
        raise ValueError('original regional numerical checker was not executed')


def record_value(*, root=ROOT):
    root = Path(root)
    approved = compare_receipts(strict_json.load(root / COMPLETION), strict_json.load(root / REGIONAL), root=root)
    execution = strict_json.load(root / EXECUTIONS)
    validate_execution(execution, {p: sha(root / p) for p in (COMPLETION, REGIONAL)}, approved, root=root)
    return {'schema': 'oph.field_current_late_class_consumer_admission.v1',
            'authority_sha256': AUTHORITY_SHA256, 'current_class_sha256': sha(root / CLASS),
            'source_sha256': own_pins(root), 'executed': execution,
            'internal_parent_execution_sha256': sha(root / 'model/completion/field_current_interface_parent_replays.json'),
            'exact_source_additions_by_location': approved['exact_source_additions_by_location'],
            'exact_added_scientific_subtree': approved['exact_added_subtree'],
            'old_scientific_fields_and_statuses_exactly_preserved': True,
            'complete_scientific_payload_unchanged': False,
            'added_child_result_is_metadata_only': False,
            'physical_statuses_changed': False, 'class_modified': False,
            'class_does_not_bind_consumers_back': True}


def admit_executed(directory, *, root=ROOT):
    root, directory = Path(root), Path(directory)
    completion = strict_json.load(directory / 'completion.json')
    regional = strict_json.load(directory / 'regional.json')
    approved = compare_receipts(completion, regional, root=root)
    execution = strict_json.load(directory / 'executions.json')
    receipt_hashes = {}
    before = {name: (root / name).read_bytes() for name in (COMPLETION, REGIONAL)}
    for name, value in ((COMPLETION, completion), (REGIONAL, regional)):
        canonical_bytes = (json.dumps(value, indent=2, sort_keys=True) + '\n').encode()
        if hashlib.sha256(canonical_bytes).hexdigest() != execution[name]['admitted_sha256']:
            raise ValueError('saved actual-owner output differs from execution record')
        receipt_hashes[name] = hashlib.sha256(canonical_bytes).hexdigest()
    validate_execution(execution, receipt_hashes, approved, root=root)
    for name, raw in before.items():
        if hashlib.sha256(raw).hexdigest() not in {
                approved['retained_receipts'][name]['retained_sha256'], receipt_hashes[name]}:
            raise ValueError('canonical late receipt changed before reviewed admission')
    authority(root)
    if any((root / name).read_bytes() != raw for name, raw in before.items()):
        raise ValueError('late-owner output changed during admission')
    # All semantic checks precede any receipt mutation. This write follows the
    # actual saved original report/checker executions, not source-pin guessing.
    for name, value in ((COMPLETION, completion), (REGIONAL, regional)):
        (root / name).write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')
    (root / EXECUTIONS).write_text(json.dumps(execution, indent=2, sort_keys=True) + '\n')
    value = record_value(root=root)
    (root / OUTPUT).write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')
    return value


def replay():
    authority()
    from model.particle_identification import completion_certificate
    fresh = completion_certificate.report()
    if live.canonical(fresh) != live.canonical(strict_json.load(ROOT / COMPLETION)):
        raise ValueError('full original completion report differs')
    subprocess.run([sys.executable, '-B', '-W', 'error', '-m', 'model.regional_clock.verify',
                    '--check', str(ROOT / REGIONAL)], cwd=ROOT, check=True)
    return check()


def check():
    value = record_value()
    if strict_json.load(ROOT / OUTPUT) != value:
        raise ValueError('late class-consumer admission record is stale')
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--check', action='store_true')
    mode.add_argument('--replay', action='store_true')
    mode.add_argument('--admit-executed', type=Path)
    args = parser.parse_args()
    value = admit_executed(args.admit_executed) if args.admit_executed else replay() if args.replay else check()
    print(json.dumps(value, indent=2, sort_keys=True))
