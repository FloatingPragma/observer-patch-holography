"""Three executed receipt renewals layered over the immutable atlas authority.

The retained source IDs, original historical replay and physical guards are
unchanged. A private instance of the retained implementation resolves exactly
three additional historical input paths. Their original owners have executed;
this service checks those executions, exact scientific equality, the one
enumerated source append, and every live source binding on each admission.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import importlib.util
import json
from pathlib import Path

from tools import strict_json
from tools import renew_canonical_live_receipts as live
from .retained_source_equivalence import RETAINED_NAMES as _RETAINED_NAMES

ROOT = Path(__file__).resolve().parents[2]
AUTHORITY = 'model/particle_atlas_explainer/current_retained_source_authority.json'
AUTHORITY_SHA256 = '952c01e174966c74c1ef3d0658d2bca03e6fc3eb50fd0c8f38372ceb1d715315'
OUTPUT = 'model/particle_atlas_explainer/current_retained_source_navigation.json'
SOURCE = 'model/particle_atlas_explainer/current_retained_source_equivalence.py'
RETAINED_SOURCE = 'model/particle_atlas_explainer/retained_source_equivalence.py'
RETAINED_NAMES = frozenset(_RETAINED_NAMES)
NAMES = frozenset({
    'model/protected_defect_catalog/catalog.json',
    'model/spin_charge_dispersion/results.json',
    'model/spin_joint_field_limit/results.json',
})


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def authority(root=ROOT):
    root = Path(root)
    if sha(root / AUTHORITY) != AUTHORITY_SHA256:
        raise ValueError('current atlas custody authority changed')
    value = strict_json.load(root / AUTHORITY)
    if set(value['additional_live_receipt_transitions']) != NAMES:
        raise ValueError('current atlas transition inventory changed')
    for name, expected in value['unchanged_source_sha256'].items():
        if sha(root / name) != expected:
            raise ValueError('retained atlas or actual owner evidence changed: ' + name)
    return value


def load_authority(root=ROOT):
    """Expose the original authority to callers reading retained receipt IDs."""
    authority(root)
    from . import retained_source_equivalence as retained
    return retained.load_authority(root)


def compare_transition(name, spec, old, current, *, root=ROOT):
    """Preserve every scientific value and the precisely reviewed source scope."""
    if name not in NAMES or spec['receipt'] != name:
        raise ValueError('unreviewed current atlas receipt')
    additions = spec.get('reviewed_source_scope_append', {})
    admitted = live.carry_metadata(name, old, current, scope_append=additions, root=root)
    if live.canonical(admitted) != live.canonical(current):
        raise ValueError('current atlas source science differs: ' + name)
    for _, pins, _ in live.source_maps(current):
        for path, expected in additions.items():
            if pins.get(path) != expected:
                raise ValueError('reviewed atlas source append changed: ' + path)


def validate_current_transition(name, spec, *, root=ROOT):
    root = Path(root)
    approved = authority(root)
    if approved['additional_live_receipt_transitions'].get(name) != spec:
        raise ValueError('unreviewed current atlas transition')
    archive = root / spec['retained_path']
    if sha(archive) != spec['retained_sha256'] or sha(root / name) != spec['current_sha256']:
        raise ValueError('executed current atlas receipt or retained input changed: ' + name)
    record = strict_json.load(root / approved['actual_owner_replay'])['executed'][spec['owner']]
    expected_record = {k: v for k, v in spec.items()
                       if k not in {'reviewed_nested_receipt_pointers', 'owner'}}
    if 'owner' in record:
        expected_record['owner'] = spec['owner']
    if record != expected_record:
        raise ValueError('current atlas transition lacks its exact owner execution')
    compare_transition(name, spec, strict_json.load(archive), strict_json.load(root / name), root=root)
    for _, relative, expected, path in live.bindings(root, name, strict_json.load(root / name)):
        if sha(path) != expected:
            raise ValueError('current atlas source custody is stale: ' + name + ':' + relative)
    return {'retained_sha256': spec['retained_sha256'],
            'current_sha256': spec['current_sha256'], 'retained_path': spec['retained_path'],
            'non_provenance_payload_equal': True,
            'actual_unchanged_owner_execution_authenticated': True,
            'reviewed_source_scope_append': spec.get('reviewed_source_scope_append', {})}


def service(root=ROOT):
    """Keep the imported retained module and all its callers untouched."""
    root = Path(root)
    authority(root)
    spec = importlib.util.spec_from_file_location(
        'model.particle_atlas_explainer._current_retained_source', root / RETAINED_SOURCE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    original_load = module.load_authority
    original_transition = module.validate_transition
    original_navigation = module.navigation

    def combined_authority(root=root):
        extra = authority(root)
        combined = copy.deepcopy(original_load(root))
        for name, row in extra['additional_live_receipt_transitions'].items():
            if (name in combined['allowed_live_receipt_transitions']
                    or combined['retained_input_sha256'].get(name) != row['retained_sha256']):
                raise ValueError('current atlas transition changes retained input scope')
            combined['allowed_live_receipt_transitions'][name] = row
        return combined

    def transition(name, row, *, root=root):
        if name in NAMES:
            return validate_current_transition(name, row, root=root)
        return original_transition(name, row, root=root)

    def current_navigation(*, root=root):
        value = original_navigation(root=root)
        del value['navigation_digest']
        value.update(schema='oph.current_retained_primary_navigation.v1',
                     current_custody_authority_sha256=AUTHORITY_SHA256,
                     current_custody_source_sha256=sha(Path(root) / SOURCE),
                     additional_executed_receipts=sorted(NAMES),
                     retained_authority_bytes_changed=False)
        value['navigation_digest'] = hashlib.sha256(live.canonical(value).encode()).hexdigest()
        return value

    module.load_authority = combined_authority
    module.validate_transition = transition
    module.navigation = current_navigation
    return module


def navigation(*, root=ROOT):
    return service(root).navigation(root=root)


def historical_input_paths(*, root=ROOT):
    return service(root).historical_input_paths(root=root)


def authenticate_pins(pins, *, root=ROOT):
    module = service(root)
    result = module.authenticate_pins(pins, root=root)
    # The current service remains visible even when every particular input pin
    # already equals its live value; it is never an implicit hash exception.
    result.update({AUTHORITY: AUTHORITY_SHA256, SOURCE: sha(Path(root) / SOURCE)})
    result.update(authority(root)['unchanged_source_sha256'])
    return result


def read_receipt(path, *, root=ROOT):
    return service(root).read_receipt(path, root=root)


def validate_receipt_value(name, value, *, root=ROOT):
    return service(root).validate_receipt_value(name, value, root=root)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--renew', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    value = navigation()
    if args.renew:
        (ROOT / OUTPUT).write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')
    elif strict_json.load(ROOT / OUTPUT) != value:
        raise ValueError('current atlas custody projection is stale')
    print('Verified three executed atlas transitions; retained source IDs unchanged: '
          + value['navigation_digest'])


if __name__ == '__main__':
    main()
