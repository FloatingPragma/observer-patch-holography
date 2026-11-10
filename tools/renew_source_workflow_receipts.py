#!/usr/bin/env python3
"""Replay the five named live workflow owners without changing their science.

PYTHONPATH=code python3 tools/renew_source_workflow_receipts.py --renew
PYTHONPATH=code python3 tools/renew_source_workflow_receipts.py --check

The unchanged builders and independent verifiers execute in both modes. The
controls and receipts may change only their top-level custody maps and derived
custody digests. Captured simulator operations and frozen policies are inputs.
No frozen campaign or historical comparison record is written.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'code'), str(ROOT)]
from tools import strict_json  # noqa: E402

OWNERS = ('source_temporal_acceptance', 'source_checkpoint_selection',
          'source_publication_selection', 'source_operation_reads',
          'independent_postdictions')


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()


def conserved_payload(retained, replay, custody_fields):
    """Exact type-sensitive conservation; nested fields are never masked."""
    left = {k: v for k, v in retained.items() if k not in custody_fields}
    right = {k: v for k, v in replay.items() if k not in custody_fields}
    if canonical(left) != canonical(right):
        raise ValueError('scientific payload differs from retained live receipt')
    return hashlib.sha256(canonical(left)).hexdigest()


def replay_owner(owner):
    if owner not in OWNERS:
        raise ValueError('unknown live workflow owner')
    here = ROOT / 'code' / owner
    check = importlib.import_module(owner + '.verify')
    if owner == 'source_operation_reads':
        result = check.verify(check.check.strict_load(here / 'capture.json'))
        return [('receipt.json', result, {'source_pins'}, check.check.canonical)]
    build = importlib.import_module(owner + '.build')
    if owner == 'independent_postdictions':
        result = {'schema': 'oph-independent-postdictions-v1',
                  'sources': check.pins(), 'calculations': build.calculations()}
        check.verify_policy()
        check.verify_calculations(result['calculations'])
        def encode(value):
            return (json.dumps(value, indent=2, sort_keys=True) + '\n').encode('ascii')
        return [('receipt.json', result, {'sources'}, encode)]
    codec = importlib.import_module(owner + '.codec')
    controls = build.packet() if owner == 'source_publication_selection' else build.build()
    receipt = check.verify(controls)
    control_fields = {'source_pins', 'sha256'} if owner == 'source_publication_selection' else {'pins'}
    receipt_fields = {
        'source_temporal_acceptance': {'pins', 'controls_sha256'},
        'source_checkpoint_selection': {'controls_sha256', 'source_pins_sha256'},
        'source_publication_selection': {'verified_controls_sha256'},
    }[owner]
    return [('controls.json', controls, control_fields, codec.canonical),
            ('receipt.json', receipt, receipt_fields, codec.canonical)]


def execute(owner, *, renew=False):
    writes, report = [], []
    for name, fresh, custody_fields, encode in replay_owner(owner):
        path = ROOT / 'code' / owner / name
        original = path.read_bytes()
        retained = strict_json.loads(original.decode())
        science = conserved_payload(retained, fresh, custody_fields)
        encoded = encode(fresh)
        if not renew and original != encoded:
            raise ValueError('live custody differs from producer replay: ' + path.relative_to(ROOT).as_posix())
        writes.append((path, original, encoded))
        report.append({'path': path.relative_to(ROOT).as_posix(),
                       'scientific_payload_sha256': science,
                       'scientific_payload_exactly_preserved': True,
                       'custody_fields': sorted(custody_fields),
                       'sha256': hashlib.sha256(encoded).hexdigest()})
    # The paired controls and receipt both pass conservation before writing.
    if renew:
        for path, original, _ in writes:
            if path.read_bytes() != original:
                raise ValueError('live receipt changed concurrently: ' + str(path))
        for path, _, encoded in writes:
            path.write_bytes(encoded)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--renew', action='store_true')
    mode.add_argument('--check', action='store_true')
    parser.add_argument('--owner', choices=OWNERS, action='append')
    args = parser.parse_args()
    for owner in args.owner or OWNERS:
        print(json.dumps({'owner': owner, 'results': execute(owner, renew=args.renew)},
                         sort_keys=True), flush=True)


if __name__ == '__main__':
    main()
