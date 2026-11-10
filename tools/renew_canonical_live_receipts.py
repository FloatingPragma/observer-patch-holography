#!/usr/bin/env python3
"""Replay enumerated canonical live receipts and renew only source metadata.

Every owner has an authenticated archived input and pinned producer module.
The two explicitly recorded parent-custody adapters preserve their original
scientific calculations; their prior producer bytes are archived separately.
Exact payload comparison is the default. The three named binary64 owners use
their existing report/array comparisons, retaining every original number and
NPZ byte. A discrepancy outside those contracts stops before any output write.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import importlib
import json
from pathlib import Path
import sys
import tempfile
import shutil

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from tools import strict_json  # noqa: E402
from tools.replay_model_receipts import bindings, science_payload, source_maps  # noqa: E402

CATALOG = 'tools/canonical_live_receipt_owners_r2063.json'
CATALOG_SHA256 = '5a40afb85402428cf8d168676248d2e39b73cb6175dbf168ba29ae3989638ea8'
CAMPAIGN = 'model/completion/campaign.json'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def verify_executable_source_pins(pins, root=None):
    """Metadata maintenance cannot silently admit changed transitive code."""
    root = ROOT if root is None else Path(root).resolve()
    for name, expected in pins.items():
        path = (root / name).resolve()
        if not path.is_relative_to(root) or not path.is_file() or digest(path) != expected:
            raise ValueError('reviewed transitive executable source changed: ' + name)


def load_catalog():
    if digest(ROOT / CATALOG) != CATALOG_SHA256:
        raise ValueError('live owner authority changed')
    catalog = strict_json.load(ROOT / CATALOG)
    verify_executable_source_pins(catalog['reviewed_executable_source_sha256'])
    return catalog['owners']


def campaign_pointer(value):
    return (isinstance(value, dict) and value.get('path') == CAMPAIGN
            and 'sha256' in value and set(value) <= {'path', 'sha256', 'campaign_id'})


def at_path(value, path):
    for key in path:
        value = value[key]
    return value


def set_path(value, path, replacement):
    at_path(value, path[:-1])[path[-1]] = replacement


def metadata_view(receipt, value, pointers=(), root=None):
    root = ROOT if root is None else Path(root)
    result = science_payload(root, receipt, value)
    def walk(node):
        if campaign_pointer(node):
            node['sha256'] = '<ADMITTED_NATIVE_CONTROL_CUSTODY>'
        if isinstance(node, dict):
            for child in node.values():
                walk(child)
        elif isinstance(node, list):
            for child in node:
                walk(child)
    walk(result)
    for pointer in pointers:
        set_path(result, pointer['json_path'], '<REVIEWED_LIVE_RECEIPT_POINTER>')
    return result


def carry_metadata(receipt, retained, fresh, scope_append=None, pointers=(),
                   allow_retained_pointers=False, root=None):
    """Copy only equally shaped source maps and the admitted campaign pointer."""
    output = copy.deepcopy(retained)
    root = ROOT if root is None else Path(root)
    additions = {} if scope_append is None else scope_append
    for path, expected in additions.items():
        if digest(root / path) != expected:
            raise ValueError('reviewed source-scope input changed')
    left = {(key, where): value for key, value, where in source_maps(output)}
    right = {(key, where): value for key, value, where in source_maps(fresh)}
    if set(left) != set(right):
        raise ValueError('source map inventory changed')
    expanded_locations = set()
    for location, a in left.items():
        b = right[location]
        if type(a) is not dict or type(b) is not dict:
            if canonical(a) != canonical(b):
                raise ValueError('non-map source inventory changed')
            continue
        if set(a) != set(b):
            if not additions or set(b) != set(a) | set(additions):
                raise ValueError('source map scope changed')
            expanded_locations.add(location)
            if any(b[path] != expected for path, expected in additions.items()):
                raise ValueError('source-scope append is not the reviewed input')
        a.update(b)
    def walk(a, b):
        if campaign_pointer(a):
            if not campaign_pointer(b) or b['sha256'] != digest(root / CAMPAIGN):
                raise ValueError('campaign pointer is not the actual admitted source')
            a['sha256'] = b['sha256']
        if isinstance(a, dict):
            for key, child in a.items():
                if key in b:
                    walk(child, b[key])
        elif isinstance(a, list) and isinstance(b, list):
            for x, y in zip(a, b):
                walk(x, y)
    walk(output, fresh)
    for pointer in pointers:
        path = pointer['json_path']
        old_hash = pointer['retained_sha256']
        if (digest(root / pointer['retained_path']) != old_hash
                or at_path(retained, path) != old_hash):
            raise ValueError('typed receipt pointer retained authority changed')
        new_hash = digest(root / pointer['target'])
        allowed = {new_hash, old_hash} if allow_retained_pointers else {new_hash}
        if at_path(fresh, path) not in allowed:
            raise ValueError('typed receipt pointer differs from actual owner')
        set_path(output, path, at_path(fresh, path))
    baseline = copy.deepcopy(retained)
    for key, source_map, where in source_maps(baseline):
        if (key, where) in expanded_locations:
            source_map.update(additions)
    if canonical(metadata_view(receipt, output, pointers, root)) != canonical(metadata_view(receipt, baseline, pointers, root)):
        raise ValueError('custody projection changed scientific content')
    return output


def artifact_name(retained):
    name = retained['artifact'].get('filename', retained['artifact'].get('name'))
    if not name or Path(name).name != name:
        raise ValueError('array artifact must be adjacent')
    return name


def native_mass_owner_check(spec, producer, retained, before, *, renew):
    """Use the complete existing checker, including its four typed summaries."""
    receipt = spec['receipt']
    fresh = copy.deepcopy(retained)
    fresh['source_bindings'] = producer.source_bindings()
    output = carry_metadata(receipt, retained, fresh)
    artifact = (ROOT / receipt).with_suffix('.npz')
    if digest(artifact) != retained['artifact']['sha256']:
        raise ValueError('retained native mass arrays changed')
    with tempfile.TemporaryDirectory(prefix='oph-native-mass-owner-') as scratch:
        candidate = Path(scratch) / Path(receipt).name
        candidate.write_text(json.dumps(output, indent=2, sort_keys=True) + '\n')
        shutil.copyfile(artifact, candidate.with_suffix('.npz'))
        producer.check(candidate)
    if (ROOT / receipt).read_bytes() != before:
        raise ValueError('canonical output changed during replay')
    if renew:
        (ROOT / receipt).write_text(json.dumps(output, indent=2, sort_keys=True) + '\n')
    elif canonical(strict_json.loads(before.decode())) != canonical(output):
        raise ValueError('canonical metadata differs from complete owner replay')
    return {'receipt': receipt, 'producer': spec['module'] + ':check',
            'retained': spec['retained_path'], 'retained_sha256': spec['retained_sha256'],
            'canonical_scientific_payload_exactly_preserved': True,
            'actual_binary64_replay_byte_equality_claimed': False,
            'comparison': 'UNCHANGED_COMPLETE_NATIVE_MASS_OWNER_CHECK',
            'full_array_absolute_admission': producer.ARRAY_ADMISSION,
            'diagnostic_absolute_admission': producer.DIAGNOSTIC_ADMISSION,
            'physical_error_bound_claimed_by_replay_tolerance': False,
            'fresh_sha256': digest(ROOT / receipt)}


def execute(owner, *, renew=False):
    owners = load_catalog()
    if owner not in owners:
        raise ValueError('unknown canonical live owner')
    spec = owners[owner]
    receipt = spec['receipt']
    archive = ROOT / spec['retained_path']
    if digest(archive) != spec['retained_sha256']:
        raise ValueError('retained owner bytes changed')
    if digest(ROOT / spec['producer_source']) != spec['producer_source_sha256']:
        raise ValueError('owner producer implementation changed')
    prior = spec.get('reviewed_custody_interface')
    if prior and digest(ROOT / prior['retained_source']) != prior['retained_source_sha256']:
        raise ValueError('retained owner interface bytes changed')
    old = strict_json.load(archive)
    before = (ROOT / receipt).read_bytes()
    current = strict_json.loads(before.decode())
    pointers = spec.get('reviewed_receipt_pointers', ())
    projected_current = carry_metadata(receipt, old, current,
        spec.get('reviewed_source_scope_append'), pointers, allow_retained_pointers=True)
    if canonical(current) != canonical(projected_current):
        raise ValueError('canonical science differs from its retained authority')
    producer = importlib.import_module(spec['module'])
    if spec['comparison'] == 'native_mass':
        return native_mass_owner_check(spec, producer, old, before, renew=renew)
    fresh = getattr(producer, spec['entrypoint'])()
    arrays = None
    if isinstance(fresh, tuple) and len(fresh) == 2:
        fresh, arrays = fresh
        if 'source_bindings' not in fresh:
            fresh['source_bindings'] = producer.source_bindings()
    if 'artifact' in old:
        fresh['artifact'] = old['artifact']
    admitted = carry_metadata(receipt, old, fresh, spec.get('reviewed_source_scope_append'), pointers)
    left, right = metadata_view(receipt, admitted, pointers), metadata_view(receipt, fresh, pointers)
    exact_report = canonical(left) == canonical(right)
    policy = spec['comparison']
    if policy == 'exact':
        if not exact_report:
            raise ValueError('non-custody scientific change requires review')
    elif policy == 'chiral':
        from model.common_source.replay import compare_report
        compare_report(left, right)
    elif policy in ('hierarchy', 'full_repair'):
        from model.common_source.hierarchy_diagnosis.audit import compare_replay
        compare_replay(left, right)
    elif policy == 'native_mass':
        producer._compare({k: v for k, v in left.items() if k not in ('source_bindings', 'artifact')},
                          {k: v for k, v in right.items() if k not in ('source_bindings', 'artifact')})
    else:
        raise ValueError('unreviewed numerical comparison')
    array_exact, max_delta = True, 0.0
    if arrays is not None:
        import numpy as np
        name = artifact_name(old)
        artifact = (ROOT / receipt).with_name(name)
        if digest(artifact) != old['artifact']['sha256']:
            raise ValueError('retained quantum arrays changed')
        with np.load(artifact, allow_pickle=False) as stored:
            if set(stored.files) != set(arrays):
                raise ValueError('complete array inventory changed')
            for key in stored.files:
                a, b = stored[key], arrays[key]
                if a.shape != b.shape or a.dtype.kind != b.dtype.kind:
                    raise ValueError('quantum array shape/type changed')
                same = np.array_equal(a, b)
                array_exact = array_exact and same
                if a.dtype.kind in 'fc':
                    if not np.isfinite(a).all() or not np.isfinite(b).all():
                        raise ValueError('nonfinite quantum array')
                    delta = float(np.max(np.abs(a-b)))
                    max_delta = max(max_delta, delta)
                    if policy == 'native_mass':
                        array_admitted = delta <= producer.ARRAY_ADMISSION
                    elif policy in ('chiral', 'full_repair'):
                        from model.common_source.hierarchy_diagnosis.audit import REPLAY_ABSOLUTE_TOLERANCE, REPLAY_RELATIVE_TOLERANCE
                        array_admitted = np.allclose(a, b, atol=REPLAY_ABSOLUTE_TOLERANCE,
                                               rtol=REPLAY_RELATIVE_TOLERANCE)
                    else:
                        array_admitted = same
                    if not array_admitted:
                        raise ValueError('complete quantum array replay failed: ' + key)
                elif not same:
                    raise ValueError('discrete quantum array changed')
    output = admitted
    for here, relative, expected, target in bindings(ROOT, receipt, output):
        if digest(target) != expected:
            raise ValueError('producer source custody mismatch: ' + here + ' ' + relative)
    if renew:
        if (ROOT / receipt).read_bytes() != before:
            raise ValueError('canonical output changed during replay')
        (ROOT / receipt).write_text(json.dumps(output, indent=2, allow_nan=False) + '\n')
    elif canonical(current) != canonical(output):
        raise ValueError('canonical metadata differs from actual replay')
    return {'receipt': receipt, 'producer': spec['module'] + ':' + spec['entrypoint'],
            'retained': spec['retained_path'], 'retained_sha256': spec['retained_sha256'],
            'canonical_scientific_payload_exactly_preserved': True,
            'actual_binary64_replay_byte_exact': exact_report and array_exact,
            'comparison': policy, 'reviewed_source_scope_append': spec.get('reviewed_source_scope_append', {}),
            'maximum_absolute_array_replay_difference': max_delta,
            'physical_error_bound_claimed_by_replay_tolerance': False,
            'fresh_sha256': digest(ROOT / receipt)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--renew', action='store_true')
    mode.add_argument('--check', action='store_true')
    parser.add_argument('--owner', required=True, choices=sorted(load_catalog()))
    args = parser.parse_args()
    print(json.dumps(execute(args.owner, renew=args.renew), indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
