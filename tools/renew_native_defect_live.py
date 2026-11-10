#!/usr/bin/env python3
"""Produce the live native-defect receipt without rewriting retained evidence.

Regenerate: python3 tools/renew_native_defect_live.py --renew
Verify by actual replay: python3 tools/renew_native_defect_live.py --check
The unchanged producer executes in both modes. Only its source_sha256 map
may differ from the immutable retained result; no scientific status changes.
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

RETAINED = 'model/completion/historical/r2063_live_inputs/model/common_source/native_defect_decision_results.json'
RETAINED_SHA256 = '77065a134b14867e5e5426172249608a17ff19413fa7e8b89043c2377bfcb862'
LIVE = 'model/common_source/native_defect_decision_live_results.json'


def validate_live(packet, root=ROOT):
    root = Path(root).resolve()
    payload = (root / RETAINED).read_bytes()
    if hashlib.sha256(payload).hexdigest() != RETAINED_SHA256:
        raise ValueError('retained native-defect evidence changed')
    retained = strict_json.loads(payload.decode())
    science = lambda row: json.dumps({k: v for k, v in row.items() if k != 'source_sha256'},
                                    sort_keys=True, separators=(',', ':'), allow_nan=False)
    if science(packet) != science(retained):
        raise ValueError('live native-defect scientific payload differs from retained evidence')
    pins = packet.get('source_sha256', {})
    if set(pins) != set(retained['source_sha256']):
        raise ValueError('live native-defect source coverage changed')
    for name, expected in pins.items():
        path = (root / name).resolve()
        if (not path.is_relative_to(root) or not path.is_file()
                or hashlib.sha256(path.read_bytes()).hexdigest() != expected):
            raise ValueError('live native-defect source custody changed: ' + name)
    return {'path': RETAINED, 'sha256': RETAINED_SHA256,
            'role': 'RETAINED_EVIDENCE_WITH_HISTORICAL_DOCUMENTATION_PIN',
            'historical_documentation_pin_claimed_current': False,
            'scientific_payload_exactly_preserved': True}


def replay():
    from model.common_source import native_defect_decision
    result = native_defect_decision.report()
    validate_live(result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--renew', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    fresh = replay()
    if args.renew:
        (ROOT / LIVE).write_text(json.dumps(fresh, indent=2) + '\n')
    elif strict_json.load(ROOT / LIVE) != fresh:
        raise ValueError('live native-defect receipt differs from actual producer replay')
    print('Verified live source custody and exact retained native-defect scientific payload')


if __name__ == '__main__':
    main()
