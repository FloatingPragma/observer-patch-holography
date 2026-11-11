"""Keep historical execution separate from exact current dependency checks."""
from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import re

from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
AUTHORITY_SHA256 = '5fa922ee926ba79f0d3b99a77111669293e76df02ec1a163cd8d3b3a8b22500b'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def implementation_template_sha(path):
    content = Path(path).read_text()
    content, count = re.subn(r"AUTHORITY_SHA256 = '[0-9a-f]{64}'",
                            "AUTHORITY_SHA256 = '" + '0'*64 + "'", content)
    if count != 1:
        raise ValueError('execution implementation authority slot changed')
    return hashlib.sha256(content.encode()).hexdigest()


def authority(root=ROOT):
    root = Path(root)
    path = root / 'model/retained_identity_execution/authority.json'
    if sha(path) != AUTHORITY_SHA256:
        raise ValueError('retained identity execution authority changed')
    value = strict_json.load(path)
    if implementation_template_sha(root / 'model/retained_identity_execution/live.py') != value['live_template_sha256']:
        raise ValueError('execution implementation template changed')
    for path, expected in value['implementation_sha256'].items():
        if sha(root / path) != expected:
            raise ValueError('execution implementation changed: ' + path)
    for path, expected in value['verification_archive_sha256'].items():
        if sha(root/path) != expected:
            raise ValueError('archived actual execution evidence changed: '+path)
    if (value['test_file_count'] != 81 or len(value['test_files']) != 81
            or value['test_count'] != 1005
            or set(value['test_files']) != set(value['nodeids_by_file'])
            or sum(map(len, value['nodeids_by_file'].values())) != 1005):
        raise ValueError('retained test inventory changed')
    baseline = value['baseline_collection_order']
    if len(baseline) != 1005 or set(baseline) != {n for rows in value['nodeids_by_file'].values() for n in rows}:
        raise ValueError('retained original collection order changed')
    for path, expected in value['test_files'].items():
        if sha(root / path) != expected:
            raise ValueError('original retained test bytes changed: ' + path)
    for path, expected in value['current_validation_sha256'].items():
        if sha(root / path) != expected:
            raise ValueError('current validation input changed: ' + path)
    return value


def compare_transition(row, old, current, *, view, root=ROOT):
    """Only the exact recorded source maps and typed receipt pointers differ."""
    from tools import renew_canonical_live_receipts as live
    root, view = Path(root), Path(view)
    retained = copy.deepcopy(old)
    for pointer in row['typed_nested_pointers']:
        before = live.at_path(retained, pointer['json_path'])
        after = live.at_path(current, pointer['json_path'])
        if (before != pointer['retained_sha256'] or before != sha(view / pointer['target'])
                or after != pointer['current_sha256'] or after != sha(root / pointer['target'])):
            raise ValueError('typed current/historical receipt pointer differs')
        live.set_path(retained, pointer['json_path'], after)
    admitted = live.carry_metadata(row['path'], retained, current,
                                   scope_append=row['source_scope_append'], root=root)
    if live.canonical(admitted) != live.canonical(current):
        raise ValueError('current scientific dependency differs from retained source: ' + row['path'])
    for _, relative, expected, path in live.bindings(root, row['path'], current):
        if sha(path) != expected:
            raise ValueError('current dependency has a stale source: ' + relative)


def validate_view(view, approved, *, root=ROOT):
    """Authenticate every original source input and every original test byte."""
    view, root = Path(view), Path(root)
    for key, count_key in (('calibration', 'calibration_bound_input_count'),
                           ('occupation', 'occupation_bound_input_count')):
        name = approved[key + '_input_receipt']
        expected = approved[key + '_input_receipt_sha256']
        if sha(view / name) != expected or sha(root / name) != expected:
            raise ValueError('immutable input receipt changed: ' + name)
        value = strict_json.load(view / name)
        pins = value['input_sha256' if key == 'calibration' else 'source_sha256']
        if len(pins) != approved[count_key]:
            raise ValueError('immutable dependency inventory changed')
        for path, expected in pins.items():
            if sha(view / path) != expected:
                raise ValueError('historical source dependency changed: ' + path)
    for path, expected in approved['test_files'].items():
        if sha(view / path) != expected:
            raise ValueError('historical test bytes changed: ' + path)
    if sha(view / 'model/completion/program_class_v2.json') != approved['retained_class_sha256']:
        raise ValueError('historical reviewed class changed')


def validate_scientific_custody(root=ROOT):
    """The only alternate domain is the separately guarded documentary view."""
    root = Path(root).resolve()
    from model.completion import field_current_records_custody as custody
    approved = custody.load_authority(root)
    manifest = 'paper/paper_release_manifest.json'
    expected = approved['release_transition_scope']['reviewed_release_source_sha256'][manifest]
    if sha(root/manifest) == expected:
        return {'domain': 'current_source_and_sealed_release_identity',
                'validation': custody.validate(root=root)}
    from model.documentary_release_successor.model import validate_current_science
    return {'domain': 'current_source_with_separately_validated_documentary_successor',
            'validation': validate_current_science(root)}


def validate_current(view, *, root=ROOT):
    """Audit all changed dependencies; historical passes are not this check."""
    root, view = Path(root), Path(view)
    approved = authority(root)
    validate_view(view, approved, root=root)
    from tools import renew_canonical_live_receipts as live
    from model.completion import field_current_records_custody as custody
    from model.completion.program_class_v2 import reviewed_contract_digest
    from model.field_current_late_custody import model as late
    custody_validation = validate_scientific_custody(root)
    late.check()
    base = strict_json.load(root / approved['calibration_input_receipt'])['input_sha256']
    changed = {p for p, expected in base.items() if sha(root / p) != expected}
    transitions = {row['path']: row for row in approved['current_changes']}
    if changed != set(transitions) or len(changed) != 9:
        raise ValueError('unreviewed current scientific/executable dependency change')
    executions = strict_json.load(root / 'model/completion/field_current_records_replays.json')['executed']
    by_receipt = {row.get('receipt'): row for row in executions.values()}
    for path, row in transitions.items():
        if (sha(view / path) != row['retained_sha256'] or sha(root / path) != row['current_sha256']
                or path not in by_receipt):
            raise ValueError('current dependency lacks its executed transition: ' + path)
        executed = by_receipt[path]
        if executed.get('current_sha256', executed.get('fresh_sha256')) != row['current_sha256']:
            raise ValueError('actual current producer does not bind the admitted dependency')
        compare_transition(row, strict_json.load(view / path), strict_json.load(root / path),
                           view=view, root=root)
    fixed = strict_json.load(root / approved['occupation_input_receipt'])['source_sha256']
    if any(sha(root / p) != h for p, h in fixed.items()):
        raise ValueError('current fixed-algebra scientific/executable dependency changed')
    old_class = strict_json.load(view / custody.CLASS)
    current_class = strict_json.load(root / custody.CLASS)
    if reviewed_contract_digest(old_class) != approved['retained_reviewed_class_digest']:
        raise ValueError('retained fixed-class review differs')
    if live.canonical(custody.science_view(custody.CLASS, old_class)) != live.canonical(
            custody.science_view(custody.CLASS, current_class)):
        raise ValueError('current fixed-class scientific clauses changed')
    return {'calibration_bound_paths_checked': len(base), 'changed_receipts': sorted(changed),
            'changed_executable_inputs': [], 'actual_current_producer_records_checked': 9,
            'all_nine_exact_typed_scientific_comparisons_passed': True,
            'fixed_algebra_bound_paths_unchanged': len(fixed),
            'fixed_class_scientific_clauses_equal': True,
            'current_scientific_custody_validation': custody_validation,
            'historical_test_passes_claimed_to_verify_current_physics': False}
