#!/usr/bin/env python3
"""Replay the canonical native control and its admitted custody-only parents.

The frozen clock contract explicitly permits the control's one documentation
pin and the completion campaign's corresponding digest to renew. Exact prior
bytes are authenticated under model/completion/historical/r2063_live_inputs.
This operation never changes a producer, a clock campaign or a scientific gate.
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

ARCHIVE = 'model/completion/historical/r2063_live_inputs'
CONTROL = 'model/common_source/native_defect_decision_results.json'
CAMPAIGN = 'model/completion/campaign.json'
QUANTUM = 'model/native_quantum_defects/results.json'
DOCUMENT = 'docs/research/OPH_FIELD_RECONSTRUCTION.md'
EXPECTED = {
    CONTROL: '77065a134b14867e5e5426172249608a17ff19413fa7e8b89043c2377bfcb862',
    CAMPAIGN: '3328d8bc91e58bcb6643a79db451a8cd9217650d860b44aa6edcea9fdc209b75',
    QUANTUM: '3a34a40b3f5857c442d8934f29f7392189b02291f0ecc948c8baa3321a155418',
}


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def retained(name):
    path = ROOT / ARCHIVE / name
    if digest(path) != EXPECTED[name]:
        raise ValueError('retained canonical input changed: ' + name)
    return strict_json.load(path)


def science(value):
    return canonical({k: v for k, v in value.items() if k != 'source_sha256'})


def compare_source_receipt(fresh, name, renewable):
    original = retained(name)
    if science(fresh) != science(original):
        raise ValueError('native scientific payload changed: ' + name)
    before, after = original['source_sha256'], fresh['source_sha256']
    if set(before) != set(after):
        raise ValueError('native source scope changed: ' + name)
    changed = {p for p in before if before[p] != after[p]}
    if not changed <= renewable:
        raise ValueError('native executable source changed: ' + repr(sorted(changed - renewable)))
    for p, expected in after.items():
        source = (ROOT / p).resolve()
        if not source.is_relative_to(ROOT) or not source.is_file() or digest(source) != expected:
            raise ValueError('native source pin mismatch: ' + p)
    return hashlib.sha256(science(fresh).encode()).hexdigest()


def execute(*, renew):
    from model.common_source import native_defect_decision
    from model.native_quantum_defects import certificate
    from model.clock_defect_energy import contract
    from model.completion import family
    previous = {name: (ROOT / name).read_bytes() for name in EXPECTED}
    fresh = native_defect_decision.report()
    control_science = compare_source_receipt(fresh, CONTROL, {DOCUMENT})
    original_campaign = retained(CAMPAIGN)
    current_campaign = strict_json.load(ROOT / CAMPAIGN)
    normalized = json.loads(json.dumps(current_campaign))
    normalized['retained_controls'][0]['sha256'] = original_campaign['retained_controls'][0]['sha256']
    if canonical(normalized) != canonical(original_campaign):
        raise ValueError('completion campaign changed beyond its admitted control digest')
    try:
        if renew:
            (ROOT / CONTROL).write_text(json.dumps(fresh, indent=2) + '\n')
            current_campaign['retained_controls'][0]['sha256'] = digest(ROOT / CONTROL)
            (ROOT / CAMPAIGN).write_text(json.dumps(current_campaign, indent=2) + '\n')
        elif canonical(strict_json.load(ROOT / CONTROL)) != canonical(fresh):
            raise ValueError('canonical native control differs from fresh producer')
        # Both unchanged contracts validate their actual live source views.
        family.load_contract()
        contract.load_contract()
        quantum = certificate.report()
        quantum_science = compare_source_receipt(quantum, QUANTUM, {CONTROL})
        if renew:
            (ROOT / QUANTUM).write_text(json.dumps(quantum, indent=2, sort_keys=True) + '\n')
        elif canonical(strict_json.load(ROOT / QUANTUM)) != canonical(quantum):
            raise ValueError('canonical quantum-defect result differs from fresh producer')
    except Exception:
        if renew:
            for name, data in previous.items():
                (ROOT / name).write_bytes(data)
        raise
    return {'scope': 'canonical live control metadata; unchanged source and complete scientific payloads',
            'scientific_status_changed': False, 'frozen_clock_campaign_changed': False,
            'retained_bytes': {ARCHIVE + '/' + p: h for p, h in EXPECTED.items()},
            'fresh_scientific_payload_sha256': {CONTROL: control_science, QUANTUM: quantum_science},
            'current_outputs': {p: digest(ROOT / p) for p in EXPECTED}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--renew', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    print(json.dumps(execute(renew=args.renew), indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
