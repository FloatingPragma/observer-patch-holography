"""Replay eight owners for one proof insertion and one enumerated source append.

Two named binary64 producers retain their complete original validators and
canonical scientific values. All other scientific payloads compare exactly.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import importlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from tools import strict_json
from tools import renew_canonical_live_receipts as live

AUTHORITY = 'tools/field_current_documentary_owners.json'
AUTHORITY_SHA256 = '69f75432e4c00cef6fc9771e3fbe14d10ed0b9410a0c22731886424fc91128f6'
OWNERS = ('spin_dynamical_memory', 'charge_channel_lifetime', 'source_energy_effort',
          'spin_charge_dispersion', 'source_charge_scattering', 'spin_boost_stress',
          'spin_energy_polarization', 'spin_multipacket_comparison')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def authority(root=ROOT):
    root = Path(root)
    require(live.digest(root/AUTHORITY) == AUTHORITY_SHA256, 'documentary owner authority changed')
    value = strict_json.load(root/AUTHORITY)
    require(tuple(value['owners']) == OWNERS, 'documentary owner inventory changed')
    delta = value['only_documentary_delta']
    old = (root/delta['retained_path']).read_bytes()
    current = (root/delta['path']).read_bytes()
    require(hashlib.sha256(old).hexdigest() == delta['retained_sha256'], 'retained proof changed')
    require(hashlib.sha256(current).hexdigest() == delta['current_sha256'], 'reviewed proof changed')
    anchor, insertion = delta['unique_anchor'].encode(), delta['inserted_text'].encode()
    require(old.count(anchor) == 1 and insertion and current == old.replace(anchor, insertion+anchor),
            'proof changed outside the one reviewed section insertion')
    live.verify_executable_source_pins(value['unchanged_source_sha256'], root)
    return value


def project(spec, old, fresh, approved, *, root=ROOT):
    delta = approved['only_documentary_delta']
    old_maps = {(key, location): pins for key, pins, location in live.source_maps(old)}
    new_maps = {(key, location): pins for key, pins, location in live.source_maps(fresh)}
    additions = spec.get('reviewed_source_scope_append', {})
    require(old_maps.keys() == new_maps.keys(), 'scientific source map inventory changed')
    for location, before in old_maps.items():
        after = new_maps[location]
        require(type(before) is dict and type(after) is dict
                and set(after) == set(before) | set(additions),
                'scientific source map scope changed')
        for name, expected in additions.items():
            require(after[name] == expected, 'reviewed source-scope addition changed: '+name)
        for name, expected in before.items():
            want = delta['current_sha256'] if name == delta['path'] else expected
            require(after[name] == want, 'unreviewed source digest changed: '+name)
    output = live.carry_metadata(spec['receipt'], old, fresh, scope_append=additions, root=root)
    if spec.get('comparison', 'exact') == 'exact':
        require(live.canonical(output) == live.canonical(fresh), 'scientific owner payload changed')
    else:
        require(spec['comparison'] == 'unchanged_complete_owner_validator'
                and spec['module'] in ('model.source_energy_effort.packet', 'model.spin_boost_stress.packet'),
                'unreviewed numerical owner admission')
    return output


def execute(owner, *, renew=False, root=ROOT):
    root = Path(root)
    approved = authority(root)
    require(owner in OWNERS, 'unknown documentary owner')
    spec = approved['owners'][owner]
    old_path = root/spec['retained_path']
    require(live.digest(old_path) == spec['retained_sha256'], 'retained owner receipt changed')
    old = strict_json.load(old_path)
    output = root/spec['receipt']
    before = output.read_bytes()
    current = strict_json.loads(before.decode())
    # Current canonical output may have either the old or admitted document
    # pin. Every scientific field and all unlisted input pins are immutable.
    candidate = copy.deepcopy(current)
    for _, pins, _ in live.source_maps(candidate):
        if approved['only_documentary_delta']['path'] in pins:
            require(pins[approved['only_documentary_delta']['path']] in {
                approved['only_documentary_delta']['retained_sha256'],
                approved['only_documentary_delta']['current_sha256']}, 'unknown documentary pin')
            pins[approved['only_documentary_delta']['path']] = approved['only_documentary_delta']['current_sha256']
        for name, expected in spec.get('reviewed_source_scope_append', {}).items():
            require(name not in pins or pins[name] == expected, 'existing appended source changed: '+name)
            pins[name] = expected
    canonical_candidate = project(spec, old, candidate, approved, root=root)
    require(live.canonical(canonical_candidate) == live.canonical(candidate),
            'canonical scientific values changed before replay')
    producer = importlib.import_module(spec['module'])
    fresh = getattr(producer, spec['entrypoint'])()
    admitted = project(spec, old, fresh, approved, root=root)
    if spec.get('comparison', 'exact') == 'unchanged_complete_owner_validator':
        producer.validate(admitted)  # Full original comparator and all failure amplitudes.
        if owner == 'source_energy_effort':
            require((output.parent/'ENERGY_READOUT.md').read_text() == fresh['energy_table_markdown']
                    == admitted['energy_table_markdown'], 'original energy rendering changed')
    for _, name, expected, target in live.bindings(root, spec['receipt'], admitted):
        require(live.digest(target) == expected, 'actual owner source changed: '+name)
    authority(root)
    require(output.read_bytes() == before, 'owner output changed during actual replay')
    if renew:
        text = (producer.canonical(admitted) if hasattr(producer, 'canonical') else
                json.dumps(admitted, indent=2, sort_keys=True, allow_nan=False)+'\n')
        output.write_text(text)
    else:
        require(live.canonical(current) == live.canonical(admitted), 'documentary custody stale; run --renew')
    return dict(owner=owner, receipt=spec['receipt'], retained_path=spec['retained_path'],
                retained_sha256=spec['retained_sha256'], current_sha256=live.digest(output),
                actual_unchanged_producer_replayed=spec['module']+':'+spec['entrypoint'],
                scientific_payload_exactly_preserved=True, source_scope_exactly_preserved=not bool(spec.get('reviewed_source_scope_append')),
                reviewed_source_scope_append=spec.get('reviewed_source_scope_append', {}),
                comparison=spec.get('comparison', 'exact'),
                fresh_binary64_science_byte_equality_claimed=spec.get('comparison', 'exact') == 'exact',
                only_admitted_source_change=approved['only_documentary_delta']['path'],
                physical_statuses_changed=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--owner', choices=OWNERS, required=True)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--renew', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    print(json.dumps(execute(args.owner, renew=args.renew), indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
