#!/usr/bin/env python3
"""Replay the particle comparison against its retained c7 field-custody checkpoint.

The original live-owner authority remains unchanged. This successor names the
one nested cross-scale receipt pointer whose previous live value is retained
in c7. Every scientific value, input-map key and executable remains exact.
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
from tools import strict_json  # noqa: E402
from tools import renew_canonical_live_receipts as retained  # noqa: E402

ARCHIVE = 'model/completion/historical/source_records_inputs/'
RECEIPT = 'model/particle_comparison/results.json'
PARENT = 'model/empirical_diagnostics/cross_scale_results.json'
RETAINED_SHA256 = '71d59d98fc8fef3da83a760e78912725357d4c0b44fb8a9d64d3d9eb3c0c1e47'
PARENT_SHA256 = '2fc9458d63f9f45209c4dfb5516cfc52b4560887f06bbeaf41e6b99e98d9bf41'
POINTER = ['actual_data_replay', 'full_replay_sha256']


def policy():
    original = retained.load_catalog()['particle_comparison']
    spec = copy.deepcopy(original)
    spec['retained_path'] = ARCHIVE + RECEIPT
    spec['retained_sha256'] = RETAINED_SHA256
    spec['reviewed_receipt_pointers'] = [{
        'json_path': POINTER, 'target': PARENT,
        'retained_path': ARCHIVE + PARENT, 'retained_sha256': PARENT_SHA256}]
    for name, expected in ((spec['retained_path'], RETAINED_SHA256),
                           (ARCHIVE + PARENT, PARENT_SHA256)):
        if retained.digest(ROOT / name) != expected:
            raise ValueError('retained c7 comparison custody changed')
    return spec


def conserve(old, fresh, *, allow_retained=False, root=ROOT):
    pointers = policy()['reviewed_receipt_pointers']
    projected = retained.carry_metadata(RECEIPT, old, fresh, pointers=pointers,
                                        allow_retained_pointers=allow_retained, root=root)
    if retained.canonical(projected) != retained.canonical(fresh):
        raise ValueError('comparison scientific payload changed')
    return projected


def execute(*, renew=False):
    spec = policy()
    old = strict_json.load(ROOT / spec['retained_path'])
    output = ROOT / RECEIPT
    before = output.read_bytes()
    current = strict_json.loads(before.decode())
    conserve(old, current, allow_retained=True)
    producer = importlib.import_module(spec['module'])
    fresh = getattr(producer, spec['entrypoint'])()
    admitted = conserve(old, fresh)
    for _, name, expected, target in retained.bindings(ROOT, RECEIPT, admitted):
        if retained.digest(target) != expected:
            raise ValueError('actual comparison source changed: ' + name)
    policy()
    if output.read_bytes() != before:
        raise ValueError('comparison changed during producer execution')
    if renew:
        output.write_text(json.dumps(admitted, indent=2, allow_nan=False) + '\n')
    elif retained.canonical(current) != retained.canonical(admitted):
        raise ValueError('comparison live custody stale; run --renew')
    return {'receipt': RECEIPT, 'retained_checkpoint': 'c7b09261853f5a3c5a67e21294c3d9511ed59fc7',
            'retained_sha256': RETAINED_SHA256, 'current_sha256': retained.digest(output),
            'unchanged_producer': spec['module'] + ':' + spec['entrypoint'],
            'actual_producer_replayed': True, 'scientific_payload_exactly_preserved': True,
            'source_scope_exactly_preserved': True, 'physical_statuses_changed': False,
            'typed_nested_receipt_pointer': spec['reviewed_receipt_pointers'][0]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--renew', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    print(json.dumps(execute(renew=args.renew), indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
