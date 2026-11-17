#!/usr/bin/env python3
"""Bind scoped native-record claims to immutable source proofs and receipts."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / 'code/openwave_native_record_rows'
OUTPUT = ROOT / 'claims/openwave_native_record_rows.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def generate():
    spec = importlib.util.spec_from_file_location('native_record_custody', PACKAGE / 'replay.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    leaves = module.verify()
    claims = json.loads((PACKAGE / 'CLAIMS.json').read_text())
    origins = json.loads((PACKAGE / 'ORIGINS.json').read_text())
    assert [r['row'] for r in claims['rows']] == list(range(21, 32))
    assert claims['official_awards_added'] == claims['physical_predictions_added'] == 0
    assert claims['one_interacting_physical_candidate'] is False
    for row in claims['rows']:
        assert row['source_package'] in origins
        assert row['official_status'] == 'untested' and row['G6_awarded'] is False
        row['source_manifest_sha256'] = origins[row['source_package']]['manifest_sha256']
    evidence = PACKAGE / 'evidence'
    unified = json.loads((evidence / 'unified/evidence/unified/RECEIPT.json').read_text())
    triality = json.loads((evidence / 'confinement/RECEIPT.json').read_text())
    decay = json.loads((evidence / 'fixed_charge/VERIFICATION.json').read_text())
    assert unified['all_registered_checks'] == 'pass' and unified['total_native_calls'] == 492837
    assert triality['native_calls'] == 663600 and triality['full_flux_states'] == 535086
    assert decay['pass'] and decay['all_native_calls_accounted'] == 2454480
    assert decay['initial_failed_attempt_retained'] and not decay['physical_spin_or_lepton_assignment']
    replay = json.loads((PACKAGE / 'RESULT.json').read_text())
    assert replay['custody'] == 'pass' and not replay['drafts_executed']
    assert replay['actual_native_calls_in_replay'] == 6582294
    for execution in replay['executions'].values():
        assert all(execution['byte_replay'].values())
    claims.update(evidence_leaves=leaves, packages=origins,
                  representative_replay='python3 code/openwave_native_record_rows/replay.py --mode representative',
                  source_sha256={str(p.relative_to(ROOT)): sha(p) for p in
                     [Path(__file__), PACKAGE/'CLAIMS.json', PACKAGE/'ORIGINS.json',
                      PACKAGE/'SNAPSHOT_MANIFEST.json', PACKAGE/'PROOF.md', PACKAGE/'replay.py',
                      PACKAGE/'RESULT.json']})
    return claims


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    text = json.dumps(generate(), indent=2, sort_keys=True) + '\n'
    if args.write:
        OUTPUT.write_text(text)
    elif OUTPUT.read_text() != text:
        raise SystemExit('Stale native-record row claims')
    print('PASS_NATIVE_RECORD_ROW_CLAIMS')


if __name__ == '__main__':
    main()
