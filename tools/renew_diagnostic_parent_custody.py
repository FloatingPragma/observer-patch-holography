#!/usr/bin/env python3
"""Replay two finite diagnostics after the exact mode-energy custody renewal.

Each original frozen campaign and receipt remains archived. The new live
campaign changes only its mode-energy source hash. Its unchanged producer then
recomputes all numerical diagnostics; only the two resulting source hashes and
the explicit campaign digest may differ. No source path or status is added.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from tools import strict_json  # noqa: E402
from tools import renew_canonical_live_receipts as retained  # noqa: E402

AUTHORITY = 'tools/diagnostic_parent_custody_authority.json'
AUTHORITY_SHA256 = '545f8510a854980d9016344936fb08344d7fee2501fc04c3a87036cf7ff281d0'
OWNERS = ('active_code_energy', 'screened_transport')
PARENT = 'model/clock_defect_energy/mode_energy_results.json'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(name):
    return hashlib.sha256((ROOT / name).read_bytes()).hexdigest()


def exact(left, right):
    """JSON equality retains Boolean/number types as well as every value."""
    return retained.canonical(left) == retained.canonical(right)


def authority():
    require(sha(AUTHORITY) == AUTHORITY_SHA256, 'diagnostic custody authority changed')
    value = strict_json.load(ROOT / AUTHORITY)
    require(set(value['owners']) == set(OWNERS), 'diagnostic owner inventory changed')
    retained.verify_executable_source_pins(value['unchanged_executable_sha256'], ROOT)
    parent = value['parent']
    require(parent['path'] == PARENT and sha(PARENT) == parent['current_sha256'],
            'mode-energy current source differs from reviewed renewal')
    require(sha(parent['original_path']) == parent['original_sha256'], 'retained mode-energy source changed')
    old = strict_json.load(ROOT / parent['original_path'])
    current = strict_json.load(ROOT / PARENT)
    require(parent['only_changed_source_pointer'] == 'Lean/lakefile.lean', 'parent transition scope expanded')
    old['source_sha256']['Lean/lakefile.lean'] = sha('Lean/lakefile.lean')
    require(exact(old, current), 'mode-energy physics or unlisted custody changed')
    for target, expected in current['source_sha256'].items():
        require(sha(target) == expected, 'mode-energy direct source custody changed: ' + target)
    for row in value['owners'].values():
        require(sha(row['original_campaign']) == row['original_campaign_sha256']
                and sha(row['original_receipt']) == row['original_receipt_sha256'],
                'retained diagnostic evidence changed')
    return value


def projected_campaign(original, parent_hash):
    result = copy.deepcopy(original)
    result['source_sha256'][PARENT] = parent_hash
    return result


def conserve_campaign(original, current, parent_hash, *, allow_original=False):
    expected = projected_campaign(original, parent_hash)
    require(exact(current, expected) or (allow_original and exact(current, original)),
            'campaign changed outside its one reviewed mode-energy pin')


def conserve_receipt(spec, original, fresh, parent_hash, campaign_hash):
    require(sorted(fresh['source_sha256']) == spec['exact_source_paths']
            == sorted(original['source_sha256']), 'diagnostic source scope changed')
    expected = copy.deepcopy(original)
    expected['source_sha256'][PARENT] = parent_hash
    expected['source_sha256'][spec['campaign']] = campaign_hash
    expected['contract_sha256'] = campaign_hash
    require(exact(expected, fresh), 'diagnostic science or unlisted metadata changed')
    return expected


def execute(owner, *, renew=False):
    approved = authority()
    require(owner in OWNERS, 'unknown diagnostic owner')
    spec = approved['owners'][owner]
    campaign_path, receipt_path = ROOT / spec['campaign'], ROOT / spec['receipt']
    before_campaign, before_receipt = campaign_path.read_bytes(), receipt_path.read_bytes()
    original_campaign = strict_json.load(ROOT / spec['original_campaign'])
    original_receipt = strict_json.load(ROOT / spec['original_receipt'])
    current_campaign = strict_json.loads(before_campaign.decode())
    current_receipt = strict_json.loads(before_receipt.decode())
    parent_hash = approved['parent']['current_sha256']
    conserve_campaign(original_campaign, current_campaign, parent_hash, allow_original=renew)
    # Current receipt may still bind the original campaign, but its entire
    # scientific content and all other pins must already be retained values.
    if not exact(current_receipt, original_receipt):
        conserve_receipt(spec, original_receipt, current_receipt, parent_hash, sha(spec['campaign']))
    expected_campaign = projected_campaign(original_campaign, parent_hash)
    changed_campaign = False
    try:
        if renew and not exact(current_campaign, expected_campaign):
            campaign_path.write_text(json.dumps(expected_campaign, indent=2, sort_keys=True) + '\n')
            changed_campaign = True
        require(exact(strict_json.load(campaign_path), expected_campaign), 'current campaign custody is stale')
        with tempfile.TemporaryDirectory(prefix='oph-diagnostic-replay-') as folder:
            output = Path(folder) / 'fresh.json'
            command = [sys.executable, '-m', spec['module'], '--write', str(output)]
            replay = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
            require(replay.returncode == 0, 'unchanged producer replay failed:\n' + replay.stdout + replay.stderr)
            fresh = strict_json.load(output)
            admitted = conserve_receipt(spec, original_receipt, fresh, parent_hash, sha(spec['campaign']))
            check = subprocess.run([sys.executable, '-m', spec['module'], '--check', str(output)],
                                   cwd=ROOT, capture_output=True, text=True)
            require(check.returncode == 0, 'unchanged cold checker failed:\n' + check.stdout + check.stderr)
            for target, expected in fresh['source_sha256'].items():
                require(sha(target) == expected, 'actual producer source changed: ' + target)
            authority()
            require(receipt_path.read_bytes() == before_receipt, 'canonical receipt changed during replay')
            if renew:
                # Keep the retained key order to expose the three metadata
                # changes in review. Its complete JSON value equals the
                # actual unchanged producer output, checked above and cold.
                receipt_path.write_text(json.dumps(admitted, indent=2, ensure_ascii=False) + '\n')
            else:
                require(exact(current_receipt, fresh), 'canonical diagnostic custody stale; run --renew')
        return {'owner': owner, 'unchanged_producer': spec['module'],
                'actual_producer_and_cold_checker_executed': True,
                'scientific_payload_exactly_preserved': True, 'source_scope_exactly_preserved': True,
                'physical_statuses_changed': False,
                'retained_campaign_sha256': spec['original_campaign_sha256'],
                'current_campaign_sha256': sha(spec['campaign']),
                'retained_receipt_sha256': spec['original_receipt_sha256'],
                'current_receipt_sha256': sha(spec['receipt']),
                'retained_campaign_identity_claimed_equal_to_current_digest': False}
    except BaseException:
        if changed_campaign:
            campaign_path.write_bytes(before_campaign)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--owner', choices=OWNERS, action='append')
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--renew', action='store_true')
    mode.add_argument('--check', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    # These two unchanged producers consume the authenticated parent receipt.
    # A new execution of its larger documentary ancestry is a distinct task.
    parent = authority()['parent']
    selected = set(args.owner or OWNERS)
    result = {'authority': {'path': AUTHORITY, 'sha256': AUTHORITY_SHA256},
              'authenticated_parent': parent,
              'parent_scientific_payload_exactly_preserved': True,
              'fresh_parent_producer_execution_claimed_by_this_pass': False,
              'owners': [execute(owner, renew=args.renew) for owner in OWNERS if owner in selected]}
    text = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(text)
    print(text, end='')


if __name__ == '__main__':
    main()
