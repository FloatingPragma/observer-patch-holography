#!/usr/bin/env python3
"""Renew one documentation pin after complete unchanged chiral-source replay.

The canonical numerical report and NPZ stay byte-identical in scientific
content. Binary64 replay is checked with the owning package's existing report
and full-array criteria and reported separately from exact custody conservation.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from tools import strict_json  # noqa: E402

RECEIPT = 'model/common_source/chiral_candidate/results.json'
ARCHIVE = 'model/completion/historical/r2063_live_inputs/' + RECEIPT
ARCHIVE_SHA256 = '930b3caaccf7adfa1d61388911dbb5c417902cf8d609643e24f6ebef087f8d8a'
DOCUMENT = 'docs/research/OPH_PARTICLE_MASS_RESPONSE.md'


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def execute(*, renew):
    import numpy as np
    from model.common_source.chiral_candidate import run
    from model.common_source.replay import compare_report, compare_arrays
    from model.common_source.hierarchy_diagnosis.audit import REPLAY_ABSOLUTE_TOLERANCE, REPLAY_RELATIVE_TOLERANCE
    if digest(ROOT / ARCHIVE) != ARCHIVE_SHA256:
        raise ValueError('retained chiral input changed')
    retained = strict_json.load(ROOT / ARCHIVE)
    current = strict_json.load(ROOT / RECEIPT)
    science = lambda x: canonical({k: v for k, v in x.items() if k != 'source_bindings'})
    if science(current) != science(retained):
        raise ValueError('canonical chiral science differs from retained bytes')
    artifact = ROOT / Path(RECEIPT).parent / retained['artifact']['filename']
    if digest(artifact) != retained['artifact']['sha256']:
        raise ValueError('canonical chiral arrays changed')
    fresh, arrays = run.run()
    fresh['artifact'] = retained['artifact']
    if set(fresh['source_bindings']) != set(retained['source_bindings']):
        raise ValueError('chiral source scope changed')
    changed = {p for p, h in retained['source_bindings'].items() if fresh['source_bindings'][p] != h}
    if not changed <= {DOCUMENT}:
        raise ValueError('chiral executable source changed: ' + repr(sorted(changed - {DOCUMENT})))
    admitted = json.loads(json.dumps(retained))
    admitted['source_bindings'] = fresh['source_bindings']
    compare_report(admitted, fresh)
    with np.load(artifact, allow_pickle=False) as stored:
        compare_arrays(dict(stored), arrays)
        exact = all(np.array_equal(stored[k], arrays[k]) for k in stored.files)
        max_delta = max(float(np.max(np.abs(stored[k] - arrays[k])))
                        for k in stored.files if stored[k].dtype.kind in 'fc')
    if renew:
        # Preserve original formatting and every other value, including the NPZ.
        raw = (ROOT / RECEIPT).read_text()
        before = current['source_bindings'][DOCUMENT]
        after = fresh['source_bindings'][DOCUMENT]
        if raw.count('"' + before + '"') != 1:
            raise ValueError('documentation pin is not a unique serialized value')
        updated = raw.replace('"' + before + '"', '"' + after + '"')
        if strict_json.loads(updated) != admitted:
            raise ValueError('canonical update changed more than its documentary pin')
        (ROOT / RECEIPT).write_text(updated)
    elif canonical(current) != canonical(admitted):
        raise ValueError('chiral documentation custody differs from actual replay')
    run.verify(ROOT / RECEIPT)
    return {'retained_receipt': {'path': ARCHIVE, 'sha256': ARCHIVE_SHA256},
            'current_receipt_sha256': digest(ROOT / RECEIPT),
            'unchanged_npz_sha256': digest(artifact),
            'canonical_scientific_payload_exactly_preserved': True,
            'actual_binary64_replay_byte_exact': exact and science(fresh) == science(retained),
            'existing_complete_report_and_array_comparisons': 'PASS',
            'replay_absolute_tolerance': REPLAY_ABSOLUTE_TOLERANCE,
            'replay_relative_tolerance': REPLAY_RELATIVE_TOLERANCE,
            'maximum_absolute_array_replay_difference': max_delta,
            'physical_error_bound_claimed_by_replay_tolerance': False,
            'scientific_status_changed': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--renew', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    print(json.dumps(execute(renew=args.renew), indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
